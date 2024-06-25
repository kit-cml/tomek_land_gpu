/**
 * @file main.cu
 * @brief Main entry point for the Drug Simulation using CUDA
 */

#include <sys/stat.h>

#include <cuda.h>
#include <cuda_runtime.h>

#include "modules/glob_funct.hpp"
#include "modules/glob_type.hpp"
#include "modules/gpu.cuh"
#include "utils/constants.hpp"
#include "utils/file_operations.hpp"
#include "utils/gpu_operations.hpp"
#include "utils/timing.hpp"

/**
 * @brief Main function for running the drug simulation
 *
 * @param argc Number of command-line arguments
 * @param argv Array of command-line arguments
 * @return int Exit status of the program
 */
int main(int argc, char **argv) {
    param_t *p_param = new param_t();  // input data for CPU
    param_t *d_p_param;                // input data for GPU parsing
    p_param->init();
    edison_assign_params(argc, argv, p_param);

    double *ic50 = (double *)malloc(14 * sample_limit * sizeof(double));
    double *conc = (double *)malloc(sample_limit * sizeof(double));
    double *cvar = (double *)malloc(18 * sample_limit * sizeof(double));
    char *drug_name = get_drug_name(p_param->hill_file);
    double *d_ic50, *d_conc, *d_cvar, *d_ALGEBRAIC, *d_CONSTANTS, *d_RATES, *d_STATES, *d_STATES_RESULT, *d_STATES_init;
    double *d_mec_ALGEBRAIC, *d_mec_CONSTANTS, *d_mec_RATES, *d_mec_STATES;
    double *time, *dt, *states, *ical, *inal, *cai_result, *ina, *ito, *ikr, *iks, *ik1;
    cipa_t *temp_result, *cipa_result;

    int sample_size = get_IC50_data_from_file(p_param->hill_file, ic50, conc);
    int blocksPerGrid = (sample_size + threadsPerBlock - 1) / threadsPerBlock;
    printf("Sample size: %d\nSet GPU Number: %d\n", sample_size, p_param->gpu_index);

    cudaSetDevice(p_param->gpu_index);

    if (p_param->is_cvar) {
        int cvar_sample = get_cvar_data_from_file(p_param->cvar_file, sample_size, cvar);
        printf("Reading: %d Conductance Variability samples\n", cvar_sample);
    }

    prepingGPUMemory(sample_size, d_ALGEBRAIC, d_CONSTANTS, d_RATES, d_STATES, d_mec_ALGEBRAIC, d_mec_CONSTANTS,
                     d_mec_RATES, d_mec_STATES, d_p_param, temp_result, cipa_result, d_STATES_RESULT, d_ic50, ic50,
                     d_conc, conc, p_param);

    tic();

    if (p_param->is_time_series == 1 ) {

        std::regex pattern("/([a-zA-Z0-9_\.]+)\.csv");
        std::smatch match;
        std::string fname = p_param->hill_file;
        std::regex_search(fname, match, pattern);
        
        printf("%s\n", match[1].str().c_str());

        printf("Using cached initial state from previous result!!!! \n\n");

        const unsigned int datapoint_size = p_param->sampling_limit;
        double *cache;
        cache = (double *)malloc((num_of_states + 2) * sample_limit * sizeof(double));

        double *d_ic50;
        double *d_conc;
        double *d_cvar;
        double *d_ALGEBRAIC;
        double *d_CONSTANTS;
        double *d_RATES;
        double *d_STATES;
        double *d_STATES_cache;
        double *d_mec_CONSTANTS, *d_mec_STATES, *d_mec_RATES, *d_mec_ALGEBRAIC;
        // actually not used but for now, this is only for satisfiying the GPU regulator parameters
        double *d_STATES_RESULT;
        double *d_all_states;

        double *time;
        double *dt;
        double *states;
        double *ical;
        double *inal;
        double *cai_result;
        double *ina;
        double *ito;
        double *ikr;
        double *iks;
        double *ik1;
        double *tension;
        cipa_t *temp_result, *cipa_result;

        static const int CALCIUM_SCALING = 1000000;
        static const int CURRENT_SCALING = 1000;

        int num_of_constants = 146;
        int num_of_states = 42;
        int num_of_algebraic = 199;
        int num_of_rates = 42;

        // snprintf(buffer, sizeof(buffer),
        //   "./drugs/bepridil/IC50_samples.csv"
        //   // "./drugs/bepridil/IC50_optimal.csv"
        //   // "./IC50_samples.csv"
        //   );

        int sample_size = get_IC50_data_from_file(p_param->hill_file, ic50, conc, drug_name);
        if (sample_size == 0)
            printf("Something problem with the IC50 file!\n");
        // else if(sample_size > 2000)
        //     printf("Too much input! Maximum sample data is 2000!\n");
        printf("Sample size: %d\n", sample_size);
        printf("Set GPU Number: %d\n", p_param->gpu_index);

        cudaSetDevice(p_param->gpu_index);

        if (p_param->is_cvar == true) {
            int cvar_sample = get_cvar_data_from_file(p_param->cvar_file, sample_size, cvar);
            printf("Reading: %d Conductance Variability samples\n", cvar_sample);
        }

        printf("preparing GPU memory space \n");

        // char buffer_cvar[255];
        // snprintf(buffer_cvar, sizeof(buffer_cvar),
        // "./result/66_00.csv"
        // // "./drugs/optimized_pop_10k.csv"
        // );
        int cache_num = get_init_data_from_file(p_param->cache_file, cache);
        printf("Found cache for %d samples\n", cache_num);
        // note to self:
        // num of states+2 gave you at the very end of the file (pace number)
        // the very beginning -> the core number
        //   for (int z = 0; z <  num_of_states; z++) {printf("%lf\n", cache[z+1]);}
        //   printf("\n");
        //   for (int z = 0; z <  num_of_states; z++) {printf("%lf\n", cache[ 1*(num_of_states+2) + (z+2)]);}
        //   printf("\n");
        //   for (int z = 0; z <  num_of_states; z++) {printf("%lf\n", cache[ 2*(num_of_states+2) + (z+3)]);}
        // return 0 ;

        cudaMalloc(&d_ALGEBRAIC, num_of_algebraic * sample_size * sizeof(double));
        cudaMalloc(&d_CONSTANTS, num_of_constants * sample_size * sizeof(double));
        cudaMalloc(&d_RATES, num_of_rates * sample_size * sizeof(double));
        cudaMalloc(&d_STATES, num_of_states * sample_size * sizeof(double));
        cudaMalloc(&d_STATES_cache, (num_of_states + 2) * sample_size * sizeof(double));
        cudaMalloc(&d_mec_ALGEBRAIC, 24 * sample_size * sizeof(double));
        cudaMalloc(&d_mec_CONSTANTS, 29 * sample_size * sizeof(double));
        cudaMalloc(&d_mec_RATES, 7 * sample_size * sizeof(double));
        cudaMalloc(&d_mec_STATES, 7 * sample_size * sizeof(double));

        cudaMalloc(&d_p_param, sizeof(param_t));

        // prep for 1 cycle plus a bit (7000 * sample_size)
        cudaMalloc(&temp_result, sample_size * sizeof(cipa_t));
        cudaMalloc(&cipa_result, sample_size * sizeof(cipa_t));

        cudaMalloc(&time, sample_size * datapoint_size * sizeof(double));
        cudaMalloc(&dt, sample_size * datapoint_size * sizeof(double));
        cudaMalloc(&states, sample_size * datapoint_size * sizeof(double));
        cudaMalloc(&ical, sample_size * datapoint_size * sizeof(double));
        cudaMalloc(&inal, sample_size * datapoint_size * sizeof(double));
        cudaMalloc(&cai_result, sample_size * datapoint_size * sizeof(double));
        cudaMalloc(&ina, sample_size * datapoint_size * sizeof(double));
        cudaMalloc(&ito, sample_size * datapoint_size * sizeof(double));
        cudaMalloc(&ikr, sample_size * datapoint_size * sizeof(double));
        cudaMalloc(&iks, sample_size * datapoint_size * sizeof(double));
        cudaMalloc(&ik1, sample_size * datapoint_size * sizeof(double));
        cudaMalloc(&tension, sample_size * datapoint_size * sizeof(double));
        // cudaMalloc(&d_STATES_RESULT, (num_of_states+1) * sample_size * sizeof(double));
        // cudaMalloc(&d_all_states, num_of_states * sample_size * p_param->find_steepest_start * sizeof(double));

        printf("Copying sample files to GPU memory space \n");
        cudaMalloc(&d_ic50, sample_size * 14 * sizeof(double));
        cudaMalloc(&d_cvar, sample_size * 18 * sizeof(double));
        cudaMalloc(&d_conc, sample_size * sizeof(double));
        cudaMemcpy(d_STATES_cache, cache, (num_of_states + 2) * sample_size * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_ic50, ic50, sample_size * 14 * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_cvar, cvar, sample_size * 18 * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_conc, conc, sample_size * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(d_p_param, p_param, sizeof(param_t), cudaMemcpyHostToDevice);

        // // Get the maximum number of active blocks per multiprocessor
        // cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocks, do_drug_sim_analytical, threadsPerBlock);

        // // Calculate the total number of blocks
        // int numTotalBlocks = numBlocks * cudaDeviceGetMultiprocessorCount();

        tic();
        printf("Timer started, doing simulation.... \n\n\nGPU Usage at this moment: \n");
        int thread = 32;
        int block = (sample_size + thread - 1) / thread;
        // int block = (sample_size + thread - 1) / thread;
        if (gpu_check(15 * sample_size * sizeof(double) + sizeof(param_t)) == 1) {
            printf("GPU memory insufficient!\n");
            return 0;
        }
        printf("Sample size: %d\n", sample_size);
        cudaSetDevice(p_param->gpu_index);
        printf("\n   Configuration: \n\n\tblock\t||\tthread\n---------------------------------------\n  \t%d\t||\t%d\n\n\n", block, thread);
        // initscr();
        // printf("[____________________________________________________________________________________________________]  0.00 %% \n");

        kernel_DrugSimulation_postpro<<<block, thread>>>(d_ic50, d_cvar, d_conc, d_CONSTANTS, d_STATES, d_STATES_cache, d_RATES, d_ALGEBRAIC,
                                                 d_mec_CONSTANTS, d_mec_STATES, d_mec_RATES, d_mec_ALGEBRAIC,
                                                 d_STATES_RESULT, d_all_states,
                                                 time, states, dt, cai_result,
                                                 ina, inal,
                                                 ical, ito,
                                                 ikr, iks,
                                                 ik1, tension,
                                                 sample_size,
                                                 temp_result, cipa_result,
                                                 d_p_param);
        // block per grid, threads per block
        // endwin();

        cudaDeviceSynchronize();

        printf("allocating memory for computation result in the CPU, malloc style \n");
        double *h_states, *h_time, *h_dt, *h_ical, *h_inal, *h_cai_result, *h_ina, *h_ito, *h_ikr, *h_iks, *h_ik1, *h_tension;
        cipa_t *h_cipa_result;

        h_states = (double *)malloc(datapoint_size * sample_size * sizeof(double));
        printf("...allocated for STATES, \n");
        h_time = (double *)malloc(datapoint_size * sample_size * sizeof(double));
        printf("...allocated for time, \n");
        h_dt = (double *)malloc(datapoint_size * sample_size * sizeof(double));
        printf("...allocated for dt, \n");
        h_cai_result = (double *)malloc(datapoint_size * sample_size * sizeof(double));
        printf("...allocated for Cai, \n");
        h_ina = (double *)malloc(datapoint_size * sample_size * sizeof(double));
        printf("...allocated for iNa, \n");
        h_ito = (double *)malloc(datapoint_size * sample_size * sizeof(double));
        printf("...allocated for ito, \n");
        h_ikr = (double *)malloc(datapoint_size * sample_size * sizeof(double));
        printf("...allocated for ikr, \n");
        h_iks = (double *)malloc(datapoint_size * sample_size * sizeof(double));
        printf("...allocated for iks, \n");
        h_ik1 = (double *)malloc(datapoint_size * sample_size * sizeof(double));
        printf("...allocated for ik1, \n");
        h_ical = (double *)malloc(datapoint_size * sample_size * sizeof(double));
        printf("...allocated for ICaL, \n");
        h_inal = (double *)malloc(datapoint_size * sample_size * sizeof(double));
        h_tension = (double *)malloc(datapoint_size * sample_size * sizeof(double));
        h_cipa_result = (cipa_t *)malloc(sample_size * sizeof(cipa_t));
        printf("...allocating for INaL and postprocessing, all set!\n");

        ////// copy the data back to CPU, and write them into file ////////
        printf("copying the data back to the CPU \n");

        cudaMemcpy(h_states, states, sample_size * datapoint_size * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_time, time, sample_size * datapoint_size * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_dt, dt, sample_size * datapoint_size * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_ical, ical, sample_size * datapoint_size * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_inal, inal, sample_size * datapoint_size * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_cai_result, cai_result, sample_size * datapoint_size * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_ina, ina, sample_size * datapoint_size * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_ito, ito, sample_size * datapoint_size * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_ikr, ikr, sample_size * datapoint_size * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_iks, iks, sample_size * datapoint_size * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_ik1, ik1, sample_size * datapoint_size * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_tension, tension, sample_size * datapoint_size * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_cipa_result, cipa_result, sample_size * sizeof(cipa_t), cudaMemcpyDeviceToHost);

        cudaFree(d_ALGEBRAIC);
        cudaFree(d_CONSTANTS);
        cudaFree(d_RATES);
        cudaFree(d_STATES);
        cudaFree(d_mec_ALGEBRAIC);
        cudaFree(d_mec_CONSTANTS);
        cudaFree(d_mec_RATES);
        cudaFree(d_mec_STATES);
        cudaFree(d_p_param);
        cudaFree(temp_result);
        cudaFree(cipa_result);
        cudaFree(d_STATES_RESULT);
        cudaFree(d_ic50);
        cudaFree(d_cvar);
        cudaFree(d_conc);
        cudaFree(time);
        cudaFree(dt);
        cudaFree(states);
        cudaFree(ical);
        cudaFree(inal);
        cudaFree(cai_result);
        cudaFree(ina);
        cudaFree(ito);
        cudaFree(ikr);
        cudaFree(iks);
        cudaFree(ik1);
        cudaFree(tension);
    
        FILE *writer;
        int check;
        bool folder_created = false;

        printf("writing to file... \n");
        // sample loop
        for (int sample_id = 0; sample_id < sample_size; sample_id++) {
            // printf("writing sample %d... \n",sample_id);
            char sample_str[ENOUGH];
            char conc_str[ENOUGH];
            char filename[500] = "./result/post_";
            sprintf(sample_str, "%d", sample_id);
            //sprintf(conc_str, "%.2f", conc[sample_id]);
            strcat(filename, match[1].str().c_str());
            strcat(filename, "/");
            if (folder_created == false) {
                check = mkdir(filename, 0777);
                // check if directory is created or not
                if (!check) {
                    printf("Directory created\n");
                } else {
                    printf("Unable to create directory, or the folder is already created, relax mate...\n");
                }
                folder_created = true;
            }

            strcat(filename, sample_str);
            strcat(filename, "_pace.csv");

            writer = fopen(filename, "w");
            fprintf(writer, "Time,Vm,dVm/dt,Cai,INa,INaL,ICaL,IKs,IKr,IK1,Ito,Tension\n");
            for (int datapoint = 1; datapoint < datapoint_size; datapoint++) {
                if (h_time[ sample_id + (datapoint * sample_size)] == 0.0) {break;}
                fprintf(writer, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", // change this into string, or limit the decimal accuracy, so we can decrease filesize
                        h_time[sample_id + (datapoint * sample_size)],
                        h_states[sample_id + (datapoint * sample_size)],
                        h_dt[sample_id + (datapoint * sample_size)],
                        h_cai_result[sample_id + (datapoint * sample_size)],

                        h_ina[sample_id + (datapoint * sample_size)],
                        h_inal[sample_id + (datapoint * sample_size)],

                        h_ical[sample_id + (datapoint * sample_size)],
                        h_iks[sample_id + (datapoint * sample_size)],

                        h_ikr[sample_id + (datapoint * sample_size)],
                        h_ik1[sample_id + (datapoint * sample_size)],

                        h_ito[sample_id + (datapoint * sample_size)],
                        h_tension[sample_id + (datapoint * sample_size)]);
            }
            fclose(writer);
        }

        printf("writing each biomarkers value... \n");
        // sample loop
        // char conc_str[ENOUGH];
        char filename[500] = "./result/post_";
        // sprintf(sample_str, "%d", sample_id);
        // sprintf(conc_str, "%.2f", conc[sample_id]);
        strcat(filename, match[1].str().c_str());
        strcat(filename, "/");
        // printf("creating %s... \n", filename);
        if (folder_created == false) {
            check = mkdir(filename, 0777);
            // check if directory is created or not
            if (!check) {
                printf("Directory created\n");
            } else {
                printf("Unable to create directory, or the folder is already created, relax mate...\n");
            }
            folder_created = true;
        }

        // strcat(filename,sample_str);
        strcat(filename, "_biomarkers.csv");

        writer = fopen(filename, "a");

        fprintf(writer, "sample,qnet,inal_auc,ical_auc,apd90,apd50,apd_tri,cad90,cad50,cad_tri,dvmdt_repol,vm_peak,vm_valley,vm_dia,ca_peak,ca_valley,ca_dia\n");
        for (int sample_id = 0; sample_id < sample_size; sample_id++) {
            // printf("writing sample %d... \n",sample_id);

            fprintf(writer, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", // change this into string, or limit the decimal accuracy, so we can decrease filesize
                    sample_id,
                    h_cipa_result[sample_id].qnet,
                    h_cipa_result[sample_id].inal_auc,
                    h_cipa_result[sample_id].ical_auc,

                    h_cipa_result[sample_id].apd90,
                    h_cipa_result[sample_id].apd50,
                    h_cipa_result[sample_id].apd90 - h_cipa_result[sample_id].apd50,

                    h_cipa_result[sample_id].cad90,
                    h_cipa_result[sample_id].cad50,
                    h_cipa_result[sample_id].cad90 - h_cipa_result[sample_id].cad50,

                    h_cipa_result[sample_id].dvmdt_repol,
                    h_cipa_result[sample_id].vm_peak,
                    h_cipa_result[sample_id].vm_valley,
                    h_cipa_result[sample_id].vm_dia,

                    h_cipa_result[sample_id].ca_peak,
                    h_cipa_result[sample_id].ca_valley,
                    h_cipa_result[sample_id].ca_dia

                    //      temp_result[sample_id].qnet = 0.;
                    // temp_result[sample_id].inal_auc = 0.;
                    // temp_result[sample_id].ical_auc = 0.;

                    // temp_result[sample_id].dvmdt_repol = -999;
                    // temp_result[sample_id].dvmdt_max = -999;
                    // temp_result[sample_id].vm_peak = -999;
                    // temp_result[sample_id].vm_valley = d_STATES[(sample_id * num_of_states) +V];
                    // temp_result[sample_id].vm_dia = -999;

                    // temp_result[sample_id].apd90 = 0.;
                    // temp_result[sample_id].apd50 = 0.;
                    // temp_result[sample_id].ca_peak = -999;
                    // temp_result[sample_id].ca_valley = d_STATES[(sample_id * num_of_states) +cai];
                    // temp_result[sample_id].ca_dia = -999;
                    // temp_result[sample_id].cad90 = 0.;
                    // temp_result[sample_id].cad50 = 0.;
            );
        }
        fclose(writer);

        toc();

        return 0; // k bye
    }

    printf("Timer started, doing simulation.... \n\n\nGPU Usage at this moment: \n");
    if (gpu_check(15 * sample_size * datapoint_size * sizeof(double) + sizeof(param_t)) == 1) {
        printf("GPU memory insufficient!\n");
        return 1;
    }
    printf("\n   Configuration: \n\n\tblock\t||\tthread\n---------------------------------------\n  \t%d\t||\t%d\n\n\n",
           blocksPerGrid, threadsPerBlock);

    kernel_DrugSimulation<<<blocksPerGrid, threadsPerBlock>>>(
        d_ic50, d_cvar, d_conc, d_CONSTANTS, d_STATES, d_STATES_init, d_RATES, d_ALGEBRAIC, d_STATES_RESULT,
        d_mec_CONSTANTS, d_mec_STATES, d_mec_RATES, d_mec_ALGEBRAIC, time, states, dt, cai_result, ina, inal, ical, ito,
        ikr, iks, ik1, sample_size, temp_result, cipa_result, d_p_param);
    cudaDeviceSynchronize();

    double *h_states = (double *)malloc(ORd_num_of_states * sample_size * sizeof(double));
    cipa_t *h_cipa_result = (cipa_t *)malloc(sample_size * sizeof(cipa_t));

    printf("copying the data back to the CPU \n");
    cudaMemcpy(h_states, d_STATES_RESULT, sample_size * ORd_num_of_states * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_cipa_result, cipa_result, sample_size * sizeof(cipa_t), cudaMemcpyDeviceToHost);

    write_results_to_file("./result", drug_name, h_states, h_cipa_result, sample_size, ORd_num_of_states);

    freeingMemory(d_ALGEBRAIC, d_CONSTANTS, d_RATES, d_STATES, d_mec_ALGEBRAIC, d_mec_CONSTANTS, d_mec_RATES,
                  d_mec_STATES, d_p_param, temp_result, cipa_result, d_STATES_RESULT, d_ic50, ic50, conc, h_states,
                  h_cipa_result, p_param);

    toc();

    return 0;
}
