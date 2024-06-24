#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

#include "../cellmodels/Land_2016.hpp"
#include "../cellmodels/Ohara_Rudy_2011.hpp"
#include "../utils/constants.hpp"
#include "glob_funct.hpp"
#include "glob_type.hpp"
#include "gpu.cuh"
#include "param.hpp"

/**
 * @brief Main kernel function to run drug simulation for all samples in parallel.
 *
 * @param d_ic50 Array of IC50 values.
 * @param d_cvar Array of conductance variability values.
 * @param d_conc Array of drug concentrations.
 * @param d_CONSTANTS Array of constants.
 * @param d_STATES Array of states.
 * @param d_RATES Array of rates.
 * @param d_ALGEBRAIC Array of algebraic values.
 * @param d_STATES_RESULT Array to store the result states.
 * @param sample_size Sample size.
 * @param temp_result Temporary result array.
 * @param cipa_result CIPA result array.
 * @param p_param Parameters.
 */
__global__ void kernel_DrugSimulation(double *d_ic50, double *d_cvar, double *d_conc, double *d_CONSTANTS,
                                      double *d_STATES, double *d_STATES_init, double *d_RATES, double *d_ALGEBRAIC,
                                      double *d_mec_CONSTANTS, double *d_mec_RATES, double *d_mec_STATES,
                                      double *d_mec_ALGEBRAIC, double *d_STATES_RESULT, double *time, double *states,
                                      double *out_dt, double *cai_result, double *ina, double *inal, double *ical,
                                      double *ito, double *ikr, double *iks, double *ik1, unsigned int sample_size,
                                      cipa_t *temp_result, cipa_t *cipa_result, param_t *p_param) {
    unsigned short thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_id >= sample_size) return;

    // Local arrays for each sample
    double time_for_each_sample[10000];
    double dt_for_each_sample[10000];

    // Run the drug simulation for each sample
    kernel_DoDrugSim_init(d_ic50, d_cvar, d_conc[thread_id], d_CONSTANTS, d_STATES, d_RATES, d_ALGEBRAIC,
                          d_mec_CONSTANTS, d_mec_RATES, d_mec_STATES, d_mec_ALGEBRAIC, d_STATES_RESULT,
                          time_for_each_sample, dt_for_each_sample, thread_id, sample_size, temp_result, cipa_result,
                          p_param);
}

/**
 * @brief Runs a single drug simulation on the GPU for a given sample.
 *
 * @param d_ic50 Array of IC50 values.
 * @param d_cvar Array of conductance variability values.
 * @param d_conc Drug concentration.
 * @param d_CONSTANTS Array of constants.
 * @param d_STATES Array of states.
 * @param d_RATES Array of rates.
 * @param d_ALGEBRAIC Array of algebraic values.
 * @param d_STATES_RESULT Array to store the result states.
 * @param tcurr Current time array.
 * @param dt Time step array.
 * @param sample_id Sample ID.
 * @param sample_size Sample size.
 * @param temp_result Temporary result array.
 * @param cipa_result CIPA result array.
 * @param p_param Parameters.
 */
__device__ void kernel_DoDrugSim_init(double *d_ic50, double *d_cvar, double d_conc, double *d_CONSTANTS,
                                      double *d_STATES, double *d_RATES, double *d_ALGEBRAIC, double *d_STATES_RESULT,
                                      double *d_mec_CONSTANTS, double *d_mec_RATES, double *d_mec_STATES,
                                      double *d_mec_ALGEBRAIC, double *tcurr, double *dt, unsigned short sample_id,
                                      unsigned int sample_size, cipa_t *temp_result, cipa_t *cipa_result,
                                      param_t *p_param) {
    unsigned int input_counter = 0;

    // Initialize temporary result and CiPA result structures
    auto init_result = [](cipa_t &result, const double *STATES, unsigned int sample_id) {
        result.qnet = 0.;
        result.inal_auc = 0.;
        result.ical_auc = 0.;
        result.dvmdt_repol = -999;
        result.dvmdt_max = -999;
        result.vm_peak = -999;
        result.vm_valley = STATES[(sample_id * ORd_num_of_states) + V];
        result.vm_dia = -999;
        result.apd90 = 0.;
        result.apd50 = 0.;
        result.ca_peak = -999;
        result.ca_valley = STATES[(sample_id * ORd_num_of_states) + cai];
        result.ca_dia = -999;
        result.cad90 = 0.;
        result.cad50 = 0.;
    };

    // Initialize results for this sample
    init_result(temp_result[sample_id], d_STATES, sample_id);
    init_result(cipa_result[sample_id], d_STATES, sample_id);

    // Simulation variables
    bool is_peak = false;
    tcurr[sample_id] = 0.000001;
    dt[sample_id] = p_param->dt;
    double max_time_step = 0.1, time_point = 25.0;
    double dt_set;
    int cipa_datapoint = 0;
    unsigned short pace_count = 0;
    double t_peak_capture = 0.0;
    unsigned short pace_steepest = 0;
    bool init_states_captured = false;
    bool is_eligible_AP;
    const double bcl = p_param->bcl;
    const unsigned short pace_max = p_param->pace_max;
    const unsigned short last_drug_check_pace = p_param->find_steepest_start;
    double tmax = pace_max * bcl;
    double conc = d_conc;
    double type = p_param->celltype;
    double y[7] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
    double epsilon = 10E-14;
    double vm_repol30, vm_repol90;

    // Initialize constants and apply drug effects
    initConsts(d_CONSTANTS, d_STATES, type, conc, d_ic50, d_cvar, p_param->is_dutta, p_param->is_cvar, bcl, sample_id);
    applyDrugEffect(d_CONSTANTS, conc, d_ic50, epsilon, sample_id);
    land_initConsts(false, false, y, d_mec_CONSTANTS, d_mec_RATES, d_mec_STATES, d_mec_ALGEBRAIC, sample_id);

    d_CONSTANTS[BCL + (sample_id * ORd_num_of_constants)] = bcl;

    // Main simulation loop
    // dt_set = 0.001;
    while (tcurr[sample_id] < tmax) {
        // Compute rates
        coupledComputeRates(tcurr[sample_id], d_CONSTANTS, d_RATES, d_STATES, d_ALGEBRAIC, sample_id,
                     d_mec_RATES[TRPN + (sample_id * Land_num_of_rates)]);
        land_computeRates(tcurr[sample_id], d_mec_CONSTANTS, d_mec_RATES, d_mec_STATES, d_mec_ALGEBRAIC, y, sample_id);
        // Set time step (adaptive dt)
        //NOTE: Disabled in Margara
        dt_set = set_time_step(tcurr[sample_id], time_point, max_time_step, d_CONSTANTS, d_RATES, d_STATES, d_ALGEBRAIC,
                              sample_id);
        // dt_set = 0.005;
        // Check if within the same cycle
        if (floor((tcurr[sample_id] + dt_set) / bcl) == floor(tcurr[sample_id] / bcl)) {
            dt[sample_id] = dt_set;
        } else {
            // Handle end of pacing cycle
            dt[sample_id] = (floor(tcurr[sample_id] / bcl) + 1) * bcl - tcurr[sample_id];

            // Update temporary results if this is the steepest pace
            if (temp_result[sample_id].dvmdt_repol > cipa_result[sample_id].dvmdt_repol) {
                pace_steepest = pace_count;
                cipa_result[sample_id] = temp_result[sample_id];
                cipa_result[sample_id].ca_valley = d_STATES[(sample_id * ORd_num_of_states) + cai];
                cipa_result[sample_id].vm_valley = d_STATES[(sample_id * ORd_num_of_states) + V];
                is_peak = true;
                init_states_captured = false;
            } else {
                is_peak = false;
            }

            // Reset variables for next pacing cycle
            t_peak_capture = 0.0;
            init_result(temp_result[sample_id], d_STATES, sample_id);
            pace_count++;
            input_counter = 0;
            cipa_datapoint = 0;
            is_eligible_AP = false;

            // Debug output
            if (sample_id == 0) {
                printf("core: %d pace count: %d t: %lf, steepest: %d, dvmdt_repol: %lf, conc: %lf\n", sample_id,
                       pace_count, tcurr[sample_id], pace_steepest, cipa_result[sample_id].dvmdt_repol, conc);
            }
        }

        // Solve ODEs analytically
        solveAnalytical(d_CONSTANTS, d_STATES, d_ALGEBRAIC, d_RATES, dt[sample_id], sample_id);
        land_solveEuler(dt[sample_id], tcurr[sample_id], d_STATES[cai + (sample_id * ORd_num_of_states)] * 1000.,
                        d_mec_CONSTANTS, d_mec_RATES, d_mec_STATES, sample_id);

        // Perform checks in the last few pacing cycles
        if (pace_count >= pace_max - last_drug_check_pace) {
            if (tcurr[sample_id] > ((d_CONSTANTS[(sample_id * ORd_num_of_constants) + BCL] * pace_count) +
                                    (d_CONSTANTS[(sample_id * ORd_num_of_constants) + stim_start] + 2)) &&
                tcurr[sample_id] < ((d_CONSTANTS[(sample_id * ORd_num_of_constants) + BCL] * pace_count) +
                                    (d_CONSTANTS[(sample_id * ORd_num_of_constants) + stim_start] + 10)) &&
                abs(d_ALGEBRAIC[(sample_id * ORd_num_of_algebraic) + INa]) < 1) {
                if (d_STATES[(sample_id * ORd_num_of_states) + V] > temp_result[sample_id].vm_peak) {
                    temp_result[sample_id].vm_peak = d_STATES[(sample_id * ORd_num_of_states) + V];
                    if (temp_result[sample_id].vm_peak > 0) {
                        vm_repol30 = temp_result[sample_id].vm_peak -
                                     (0.3 * (temp_result[sample_id].vm_peak - temp_result[sample_id].vm_valley));
                        vm_repol90 = temp_result[sample_id].vm_peak -
                                     (0.9 * (temp_result[sample_id].vm_peak - temp_result[sample_id].vm_valley));
                        is_eligible_AP = true;
                        t_peak_capture = tcurr[sample_id];
                    } else {
                        is_eligible_AP = false;
                    }
                }
            } else if (tcurr[sample_id] > ((d_CONSTANTS[(sample_id * ORd_num_of_constants) + BCL] * pace_count) +
                                           (d_CONSTANTS[(sample_id * ORd_num_of_constants) + stim_start] + 10)) &&
                       is_eligible_AP) {
                if (d_RATES[(sample_id * ORd_num_of_rates) + V] > temp_result[sample_id].dvmdt_repol &&
                    d_STATES[(sample_id * ORd_num_of_states) + V] <= vm_repol30 &&
                    d_STATES[(sample_id * ORd_num_of_states) + V] >= vm_repol90) {
                    temp_result[sample_id].dvmdt_repol = d_RATES[(sample_id * ORd_num_of_rates) + V];
                }
            }

            // Capture initial states and data points if in the last few paces
            if ((pace_count >= pace_max - last_drug_check_pace) && (is_peak == true) && (pace_count < pace_max)) {
                if (!init_states_captured) {
                    for (int counter = 0; counter < ORd_num_of_states; counter++) {
                        d_STATES_RESULT[(sample_id * ORd_num_of_states) + counter] =
                            d_STATES[(sample_id * ORd_num_of_states) + counter];
                    }
                    init_states_captured = true;
                }

                input_counter += sample_size;
                cipa_datapoint++;
            }
        }
        tcurr[sample_id] += dt[sample_id];
    }
}
