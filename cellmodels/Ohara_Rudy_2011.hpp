/**
 * @file ohara_rudy_2011.hpp
 * @brief Header file for the O'Hara-Rudy 2011 cardiac cell model functions.
 *
 * This file contains the declarations of the functions used in the O'Hara-Rudy 2011 cardiac cell model.
 * The functions are implemented in CUDA for parallel computation on GPU.
 */

#ifndef OHARA_RUDY_2011_HPP
#define OHARA_RUDY_2011_HPP

#include <cuda_runtime.h>
#include "enums/enum_Ohara_Rudy_2011.hpp"

/**
 * @brief Initialize the constants and state variables for the O'Hara-Rudy 2011 model.
 *
 * @param CONSTANTS Pointer to the array of constants.
 * @param STATES Pointer to the array of state variables.
 * @param type Type of the cell.
 * @param conc Drug concentration.
 * @param ic50 Pointer to the array of IC50 values for drug effects.
 * @param cvar Pointer to the array of CVAR values.
 * @param is_dutta Flag indicating if Dutta modifications are applied.
 * @param is_cvar Flag indicating if CVAR modifications are applied.
 * @param offset Offset for accessing specific elements in the arrays.
 */
__device__ void initConsts(double *CONSTANTS, double *STATES, double type, double conc, double *ic50, double *cvar, bool is_dutta, bool is_cvar, double bcl, int offset);

/**
 * @brief Compute the rates of change for the state variables.
 *
 * @param TIME The current time in the simulation.
 * @param CONSTANTS Pointer to the array of constants.
 * @param RATES Pointer to the array where the computed rates of change are stored.
 * @param STATES Pointer to the array of current state variable values.
 * @param ALGEBRAIC Pointer to the array where computed algebraic variables are stored.
 * @param offset Offset for accessing specific elements in the arrays.
 */
__device__ void coupledComputeRates(double TIME, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int offset, double land_trpn);
/**
 * @brief Solve the state variables analytically for the next time step.
 *
 * @param CONSTANTS Pointer to the array of constants.
 * @param STATES Pointer to the array of current state variable values.
 * @param ALGEBRAIC Pointer to the array where computed algebraic variables are stored.
 * @param RATES Pointer to the array where the computed rates of change are stored.
 * @param dt The time step size.
 * @param offset Offset for accessing specific elements in the arrays.
 */
__device__ void solveAnalytical(double *CONSTANTS, double *STATES, double *ALGEBRAIC, double *RATES, double dt, int offset);

/**
 * @brief Determine the appropriate time step for the integration based on the rate of change of the membrane potential.
 *
 * @param TIME The current time in the simulation.
 * @param time_point A specific time point used to decide the time step.
 * @param max_time_step The maximum allowed time step.
 * @param CONSTANTS Pointer to the array of constants.
 * @param RATES Pointer to the array where the computed rates of change are stored.
 * @param STATES Pointer to the array of current state variable values.
 * @param ALGEBRAIC Pointer to the array where computed algebraic variables are stored.
 * @param offset Offset for accessing specific elements in the arrays.
 * @return The calculated time step for the next integration step.
 */
__device__ double set_time_step(double TIME, double time_point, double max_time_step, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC, int offset);

/**
 * @brief Apply drug effects on the constants of the model.
 *
 * @param CONSTANTS Pointer to the array of constants.
 * @param conc Drug concentration.
 * @param ic50 Pointer to the array of IC50 values for drug effects.
 * @param epsilon Small value to prevent division by zero.
 * @param offset Offset for accessing specific elements in the arrays.
 */
__device__ void applyDrugEffect(double *CONSTANTS, double conc, double *ic50, double epsilon, int offset);

/**
 * @brief Apply Dutta modifications to the constants of the model.
 *
 * @param CONSTANTS Pointer to the array of constants.
 * @param offset Offset for accessing specific elements in the arrays.
 */
__device__ void ___applyDutta(double *CONSTANTS, int offset);

/**
 * @brief Apply CVAR modifications to the constants of the model.
 *
 * @param CONSTANTS Pointer to the array of constants.
 * @param cvar Pointer to the array of CVAR values.
 * @param offset Offset for accessing specific elements in the arrays.
 */
__device__ void ___applyCvar(double *CONSTANTS, double *cvar, int offset);

#endif // OHARA_RUDY_2011_HPP
