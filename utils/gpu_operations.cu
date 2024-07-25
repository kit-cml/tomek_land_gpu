#include <cuda_runtime.h>

#include <cstdio>

#include "gpu_operations.hpp"
#include "constants.hpp"

/**
 * @brief Prepares GPU memory space and copies initial data from host to device.
 *
 * @param d_ALGEBRAIC Pointer to device memory for algebraic variables.
 * @param Tomek_num_of_algebraic Number of algebraic variables.
 * @param sample_size Number of samples.
 * @param d_CONSTANTS Pointer to device memory for constants.
 * @param Tomek_num_of_constants Number of constants.
 * @param d_RATES Pointer to device memory for rates.
 * @param Tomek_num_of_rates Number of rates.
 * @param d_STATES Pointer to device memory for states.
 * @param Tomek_num_of_states Number of states.
 * @param d_p_param Pointer to device memory for parameters.
 * @param temp_result Pointer to device memory for temporary results.
 * @param cipa_result Pointer to device memory for CIPA results.
 * @param d_STATES_RESULT Pointer to device memory for states results.
 * @param d_ic50 Pointer to device memory for IC50 data.
 * @param ic50 Pointer to host memory for IC50 data.
 * @param d_conc Pointer to device memory for concentration data.
 * @param conc Pointer to host memory for concentration data.
 * @param p_param Pointer to host memory for parameters.
 */
void prepingGPUMemory(int sample_size, double *&d_ALGEBRAIC, double *&d_CONSTANTS, double *&d_RATES, double *&d_STATES,
                      double *&d_mec_ALGEBRAIC, double *&d_mec_CONSTANTS, double *&d_mec_RATES, double *&d_mec_STATES,
                      param_t *&d_p_param, cipa_t *&temp_result, cipa_t *&cipa_result, double *&d_STATES_RESULT, double *&d_ic50, 
                      
                      double *ic50, double *&d_conc, double *conc, param_t *p_param) {
    printf("preparing GPU memory space \n");

    // Allocate memory on the device
    cudaMalloc(&d_ALGEBRAIC, Tomek_num_of_algebraic * sample_size * sizeof(double));
    cudaMalloc(&d_CONSTANTS, Tomek_num_of_constants * sample_size * sizeof(double));
    cudaMalloc(&d_RATES, Tomek_num_of_rates * sample_size * sizeof(double));
    cudaMalloc(&d_STATES, Tomek_num_of_states * sample_size * sizeof(double));
    cudaMalloc(&d_mec_ALGEBRAIC, Land_num_of_algebraic * sample_size * sizeof(double));
    cudaMalloc(&d_mec_CONSTANTS, Land_num_of_constants * sample_size * sizeof(double));
    cudaMalloc(&d_mec_RATES, Land_num_of_rates * sample_size * sizeof(double));
    cudaMalloc(&d_mec_STATES, Land_num_of_states * sample_size * sizeof(double));
    cudaMalloc(&d_p_param, sizeof(param_t));
    cudaMalloc(&temp_result, sample_size * sizeof(cipa_t));
    cudaMalloc(&cipa_result, sample_size * sizeof(cipa_t));
    cudaMalloc(&d_STATES_RESULT, Tomek_num_of_states * sample_size * sizeof(double));

    // Allocate memory for IC50 and concentration data
    cudaMalloc(&d_ic50, sample_size * 14 * sizeof(double));
    cudaMalloc(&d_conc, sample_size * sizeof(double));

    // Copy data from host to device
    printf("Copying sample files to GPU memory space \n");
    cudaMemcpy(d_ic50, ic50, sample_size * 14 * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_conc, conc, sample_size * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_p_param, p_param, sizeof(param_t), cudaMemcpyHostToDevice);
}

void prepingGPUMemoryPostpro(int sample_size, double *&d_ALGEBRAIC, double *&d_CONSTANTS, double *&d_RATES, double *&d_STATES, double *d_STATES_cache,
                      double *&d_mec_ALGEBRAIC, double *&d_mec_CONSTANTS, double *&d_mec_RATES, double *&d_mec_STATES,
                      param_t *&d_p_param, cipa_t *&temp_result, cipa_t *&cipa_result, double *&d_STATES_RESULT, double *&d_ic50, 
                     
                      double *ic50, double *&d_conc, double *conc, param_t *p_param, double *cache,
                      double *time, double *dt, double *states, double *ical, double *inal, double *cai_result, double *ina, double *ito, double *ikr, double *iks, double *ik1, double *tension) {
    printf("preparing GPU memory space \n");

    // Allocate memory on the device
    cudaMalloc(&d_ALGEBRAIC, Tomek_num_of_algebraic * sample_size * sizeof(double));
    cudaMalloc(&d_CONSTANTS, Tomek_num_of_constants * sample_size * sizeof(double));
    cudaMalloc(&d_RATES, Tomek_num_of_rates * sample_size * sizeof(double));
    cudaMalloc(&d_STATES, Tomek_num_of_states * sample_size * sizeof(double));
    cudaMalloc(&d_STATES_cache, (Tomek_num_of_states + 2) * sample_size * sizeof(double));

    cudaMalloc(&d_mec_ALGEBRAIC, Land_num_of_algebraic * sample_size * sizeof(double));
    cudaMalloc(&d_mec_CONSTANTS, Land_num_of_constants * sample_size * sizeof(double));
    cudaMalloc(&d_mec_RATES, Land_num_of_rates * sample_size * sizeof(double));
    cudaMalloc(&d_mec_STATES, Land_num_of_states * sample_size * sizeof(double));

    cudaMalloc(&d_p_param, sizeof(param_t));
    cudaMalloc(&temp_result, sample_size * sizeof(cipa_t));
    cudaMalloc(&cipa_result, sample_size * sizeof(cipa_t));
    cudaMalloc(&d_STATES_RESULT, Tomek_num_of_states * sample_size * sizeof(double)); // check for wat later

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

   printf("Copying sample files to GPU memory space \n");
   cudaMalloc(&d_ic50, sample_size * 14 * sizeof(double));
//    cudaMalloc(&d_cvar, sample_size * 18 * sizeof(double));
   cudaMalloc(&d_conc, sample_size * sizeof(double));

   cudaMemcpy(d_STATES_cache, cache, (Tomek_num_of_states + 2) * sample_size * sizeof(double), cudaMemcpyHostToDevice);
   cudaMemcpy(d_ic50, ic50, sample_size * 14 * sizeof(double), cudaMemcpyHostToDevice);
        
//    cudaMemcpy(d_cvar, cvar, sample_size * 18 * sizeof(double), cudaMemcpyHostToDevice);
   cudaMemcpy(d_conc, conc, sample_size * sizeof(double), cudaMemcpyHostToDevice);
   cudaMemcpy(d_p_param, p_param, sizeof(param_t), cudaMemcpyHostToDevice);
}

/**
 * @brief Frees allocated memory on both the host and device.
 *
 * @param d_ALGEBRAIC Pointer to device memory for algebraic variables.
 * @param d_CONSTANTS Pointer to device memory for constants.
 * @param d_RATES Pointer to device memory for rates.
 * @param d_STATES Pointer to device memory for states.
 * @param d_p_param Pointer to device memory for parameters.
 * @param temp_result Pointer to device memory for temporary results.
 * @param cipa_result Pointer to device memory for CIPA results.
 * @param d_STATES_RESULT Pointer to device memory for states results.
 * @param d_ic50 Pointer to device memory for IC50 data.
 * @param ic50 Pointer to host memory for IC50 data.
 * @param conc Pointer to host memory for concentration data.
 * @param h_states Pointer to host memory for states.
 * @param h_cipa_result Pointer to host memory for CIPA results.
 * @param p_param Pointer to host memory for parameters.
 */
void freeingMemory(double *d_ALGEBRAIC, double *d_CONSTANTS, double *d_RATES, double *d_STATES, double *d_mec_ALGEBRAIC,
                   double *d_mec_CONSTANTS, double *d_mec_RATES, double *d_mec_STATES, param_t *d_p_param,
                   cipa_t *temp_result, cipa_t *cipa_result, double *d_STATES_RESULT, double *d_ic50, double *ic50,
                   double *conc, double *h_states, cipa_t *h_cipa_result, param_t *p_param) {
    // Free GPU memory
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

    // Free CPU memory
    free(ic50);
    free(conc);
    free(h_states);
    free(h_cipa_result);
    delete p_param;
}

/**
 * @brief Checks the available GPU memory.
 *
 * @param datasize Size of the data to be checked against available GPU memory.
 * @return int 0 if successful, 1 if insufficient memory.
 */
int gpu_check(unsigned int datasize) {
    int num_gpus;
    float percent;
    int id;
    size_t free, total;
    cudaGetDeviceCount(&num_gpus);
    for (int gpu_id = 0; gpu_id < num_gpus; gpu_id++) {
        cudaSetDevice(gpu_id);
        cudaGetDevice(&id);
        cudaMemGetInfo(&free, &total);
        percent = (free / (float)total);
        printf("GPU No %d\nFree Memory: %ld, Total Memory: %ld (%f percent free)\n", id, free, total, percent * 100.0);
    }
    percent = 1.0 - (datasize / (float)total);

    return (percent >= 0) ? 0 : 1;
}
