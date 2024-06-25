#ifndef GPU_CUH
#define GPU_CUH

#include <cuda.h>
#include <cuda_runtime.h>

#include "param.hpp"

/**
 * @brief Kernel function for drug simulation.
 *
 * @param d_ic50 Pointer to IC50 data on the device.
 * @param d_cvar Pointer to conductance variability data on the device.
 * @param d_conc Pointer to drug concentration data on the device.
 * @param d_CONSTANTS Pointer to constants data on the device.
 * @param d_STATES Pointer to states data on the device.
 * @param d_STATES_init Pointer to initial states data on the device.
 * @param d_RATES Pointer to rates data on the device.
 * @param d_ALGEBRAIC Pointer to algebraic data on the device.
 * @param d_STATES_RESULT Pointer to states result data on the device.
 * @param d_all_states Pointer to all states data on the device.
 * @param time Pointer to time data on the device.
 * @param states Pointer to states data on the device.
 * @param out_dt Pointer to output delta time data on the device.
 * @param cai_result Pointer to calcium result data on the device.
 * @param ina Pointer to INa data on the device.
 * @param inal Pointer to INaL data on the device.
 * @param ical Pointer to ICaL data on the device.
 * @param ito Pointer to Ito data on the device.
 * @param ikr Pointer to IKr data on the device.
 * @param iks Pointer to IKs data on the device.
 * @param ik1 Pointer to IK1 data on the device.
 * @param sample_size Size of the sample.
 * @param temp_result Pointer to temporary results on the device.
 * @param cipa_result Pointer to CiPA results on the device.
 * @param p_param Pointer to parameters on the device.
 */
__global__ void kernel_DrugSimulation(double *d_ic50, double *d_cvar, double *d_conc, double *d_CONSTANTS,
                                      double *d_STATES, double *d_STATES_init, double *d_RATES, double *d_ALGEBRAIC,
                                      double *d_mec_CONSTANTS, double *d_mec_RATES, double *d_mec_STATES,
                                      double *d_mec_ALGEBRAIC, double *d_STATES_RESULT, double *time, double *states,
                                      double *out_dt, double *cai_result, double *ina, double *inal, double *ical,
                                      double *ito, double *ikr, double *iks, double *ik1, unsigned int sample_size,
                                      cipa_t *temp_result, cipa_t *cipa_result, param_t *p_param);

/**
 * @brief Device function for initializing drug simulation.
 *
 * @param d_ic50 Pointer to IC50 data on the device.
 * @param d_cvar Pointer to conductance variability data on the device.
 * @param d_conc Drug concentration.
 * @param d_CONSTANTS Pointer to constants data on the device.
 * @param d_STATES Pointer to states data on the device.
 * @param d_RATES Pointer to rates data on the device.
 * @param d_ALGEBRAIC Pointer to algebraic data on the device.
 * @param d_STATES_RESULT Pointer to states result data on the device.
 * @param d_all_states Pointer to all states data on the device.
 * @param tcurr Pointer to current time data on the device.
 * @param dt Pointer to delta time data on the device.
 * @param sample_id Sample ID.
 * @param sample_size Size of the sample.
 * @param temp_result Pointer to temporary results on the device.
 * @param cipa_result Pointer to CiPA results on the device.
 * @param p_param Pointer to parameters on the device.
 */
__device__ void kernel_DoDrugSim_init(double *d_ic50, double *d_cvar, double d_conc, double *d_CONSTANTS,
                                      double *d_STATES, double *d_RATES, double *d_ALGEBRAIC, double *d_mec_CONSTANTS,
                                      double *d_mec_RATES, double *d_mec_STATES, double *d_mec_ALGEBRAIC,
                                      double *d_STATES_RESULT, double *tcurr, double *dt, unsigned short sample_id,
                                      unsigned int sample_size, cipa_t *temp_result, cipa_t *cipa_result,
                                      param_t *p_param);

/**
 * @brief Device function for post-processing drug simulation.
 *
 * @param d_ic50 Pointer to IC50 data on the device.
 * @param d_cvar Pointer to conductance variability data on the device.
 * @param d_conc Drug concentration.
 * @param d_CONSTANTS Pointer to constants data on the device.
 * @param d_STATES Pointer to states data on the device.
 * @param d_STATES_init Pointer to initial states data on the device.
 * @param d_RATES Pointer to rates data on the device.
 * @param d_ALGEBRAIC Pointer to algebraic data on the device.
 * @param time Pointer to time data on the device.
 * @param states Pointer to states data on the device.
 * @param out_dt Pointer to output delta time data on the device.
 * @param cai_result Pointer to calcium result data on the device.
 * @param ina Pointer to INa data on the device.
 * @param inal Pointer to INaL data on the device.
 * @param ical Pointer to ICaL data on the device.
 * @param ito Pointer to Ito data on the device.
 * @param ikr Pointer to IKr data on the device.
 * @param iks Pointer to IKs data on the device.
 * @param ik1 Pointer to IK1 data on the device.
 * @param tcurr Pointer to current time data on the device.
 * @param dt Pointer to delta time data on the device.
 * @param sample_id Sample ID.
 * @param sample_size Size of the sample.
 * @param temp_result Pointer to temporary results on the device.
 * @param cipa_result Pointer to CiPA results on the device.
 * @param p_param Pointer to parameters on the device.
 */
__device__ void kernel_DoDrugSim_post(double *d_ic50, double *d_cvar, double d_conc, double *d_CONSTANTS,
                                      double *d_STATES, double *d_STATES_init, double *d_RATES, double *d_ALGEBRAIC,
                                      double *d_mec_CONSTANTS, double *d_mec_RATES, double *d_mec_STATES,
                                      double *d_mec_ALGEBRAIC, double *time, double *states, double *out_dt,
                                      double *cai_result, double *ina, double *inal, double *ical, double *ito,
                                      double *ikr, double *iks, double *ik1, double *tcurr, double *dt,
                                      unsigned short sample_id, unsigned int sample_size, cipa_t *temp_result,
                                      cipa_t *cipa_result, param_t *p_param);



// __device__ void kernel_DoDrugSim_single(double *d_ic50, double *d_cvar, double d_conc, double *d_CONSTANTS,
//                                         double *d_STATES, double *d_STATES_cache, double *d_RATES, double *d_ALGEBRAIC,
//                                         double *d_mec_CONSTANTS, double *d_mec_STATES, double *d_mec_RATES,
//                                         double *d_mec_ALGEBRAIC, double *time, double *states, double *out_dt,
//                                         double *cai_result, double *ina, double *inal, double *ical, double *ito,
//                                         double *ikr, double *iks, double *ik1, double *tension, double *tcurr, double *dt,
//                                         unsigned short sample_id, unsigned int sample_size, cipa_t *temp_result,
//                                         cipa_t *cipa_result, param_t *p_param);

#endif  // GPU_CUH
