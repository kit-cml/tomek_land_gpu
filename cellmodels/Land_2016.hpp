#ifndef LAND_2016_HPP
#define LAND_2016_HPP

#include "cellmodel.hpp"
#include "enums/enum_Land_2016.hpp"
#include <cuda_runtime.h>

	// void initConsts(double type);
	// void initConsts(bool is_dutta);
	// void initConsts(double type, double conc, double *ic50, bool is_dutta );
	__device__ double check_max(double a, double b);
	__device__ double check_min(double a, double b);
    __device__ void land_initConsts( bool is_skinned, bool BETA, double* y, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC, int offset);
	__device__ void land_computeRates(double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC, double *y, int offset);
	__device__ void land_solveEuler( double dt, double t, double Cai_input,  double *CONSTANTS, double *RATES, double *STATES, int offset);
	//static double set_time_step(double TIME,double time_point,double max_time_step,
//   double* CONSTANTS,
//   double* RATES,
//   double* STATES,
//   double* ALGEBRAIC);
// private:
	// void ___applyDrugEffect(double conc, double *ic50, double epsilon);
	// void ___applyDutta();
	// void ___initConsts(double type);

#endif