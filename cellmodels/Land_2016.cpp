#include "Land_2016.hpp"

#include <cuda.h>
#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "../modules/glob_funct.hpp"
#include "../utils/constants.hpp"
#include "cellmodel.hpp"

/**
 * @brief Returns the maximum of two values.
 *
 * @param a First value.
 * @param b Second value.
 * @return double Maximum of a and b.
 */
__device__ inline double check_max(double a, double b) { return fmax(a, b); }

/**
 * @brief Returns the minimum of two values.
 *
 * @param a First value.
 * @param b Second value.
 * @return double Minimum of a and b.
 */
__device__ inline double check_min(double a, double b) { return fmin(a, b); }

/**
 * @brief Initialize constants and state variables for the Land model.
 *
 * @param is_skinned Boolean indicating if the model is skinned.
 * @param BETA Boolean indicating if beta parameters are used.
 * @param y Array of initial state values.
 * @param CONSTANTS Array of constant parameters.
 * @param RATES Array of rate values.
 * @param STATES Array of state variables.
 * @param ALGEBRAIC Array of algebraic variables.
 * @param sample_id sample_id for multi-threaded computation.
 */
__device__ void land_initConsts(bool is_skinned, bool BETA, double *y, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC, int sample_id) {
    // printf("Successfully go to inside land_initConsts\n");
    CONSTANTS[sample_id * Land_num_of_constants + dlambda_dt] = 0;
    CONSTANTS[sample_id * Land_num_of_constants + lambda] = 1.0;
    CONSTANTS[sample_id * Land_num_of_constants + Cai] = 0.0;

    RATES[sample_id * Land_num_of_rates + XS] = 0;
    RATES[sample_id * Land_num_of_rates + XW] = 0;
    RATES[sample_id * Land_num_of_rates + TRPN] = 0;
    RATES[sample_id * Land_num_of_rates + TmBlocked] = 1;
    RATES[sample_id * Land_num_of_rates + ZETAS] = 0;
    RATES[sample_id * Land_num_of_rates + ZETAW] = 0;
    RATES[sample_id * Land_num_of_rates + dCd_dt] = 0;

    CONSTANTS[sample_id * Land_num_of_constants + lambda] =
        check_min(1.2, CONSTANTS[sample_id * Land_num_of_constants + lambda]);

    CONSTANTS[sample_id * Land_num_of_constants + perm50] = 0.35;
    CONSTANTS[sample_id * Land_num_of_constants + TRPN_n] = 2;
    CONSTANTS[sample_id * Land_num_of_constants + koff] = 0.1;
    CONSTANTS[sample_id * Land_num_of_constants + dr] = 0.25;
    CONSTANTS[sample_id * Land_num_of_constants + wfrac] = 0.5;
    CONSTANTS[sample_id * Land_num_of_constants + TOT_A] = 25;
    CONSTANTS[sample_id * Land_num_of_constants + ktm_unblock] = 0.04;
    CONSTANTS[sample_id * Land_num_of_constants + beta_1] = -2.4;
    CONSTANTS[sample_id * Land_num_of_constants + beta_0] = 2.3;
    CONSTANTS[sample_id * Land_num_of_constants + gamma_idx] = 0.0085;
    CONSTANTS[sample_id * Land_num_of_constants + gamma_wu] = 0.615;
    CONSTANTS[sample_id * Land_num_of_constants + phi] = 2.23;

    if (is_skinned) {
        CONSTANTS[sample_id * Land_num_of_constants + nperm] = 2.2;
        ALGEBRAIC[sample_id * Land_num_of_algebraic + ca50] = 2.5;
        CONSTANTS[sample_id * Land_num_of_constants + Tref] = 40.5;
        CONSTANTS[sample_id * Land_num_of_constants + nu] = 1;
        CONSTANTS[sample_id * Land_num_of_constants + mu] = 1;
    } else {
        CONSTANTS[sample_id * Land_num_of_constants + nperm] = 2.4;
        ALGEBRAIC[sample_id * Land_num_of_algebraic + ca50] = 0.805;
        CONSTANTS[sample_id * Land_num_of_constants + Tref] = 120.0;
        CONSTANTS[sample_id * Land_num_of_constants + nu] = 7;
        CONSTANTS[sample_id * Land_num_of_constants + mu] = 3;
    }

    if (BETA) {
        // If BETA is true, use beta values (currently commented out)
        // CONSTANTS[sample_id * Land_num_of_constants + beta_1] = beta[1];
        // CONSTANTS[sample_id * Land_num_of_constants + beta_0] = beta[0];
    }

    CONSTANTS[sample_id * Land_num_of_constants + k_ws] = 0.004 * CONSTANTS[sample_id * Land_num_of_constants + mu];
    CONSTANTS[sample_id * Land_num_of_constants + k_uw] = 0.026 * CONSTANTS[sample_id * Land_num_of_constants + nu];

    // Initialize STATES with checks
    STATES[sample_id * Land_num_of_states + XS] = check_max(0, y[0]);
    STATES[sample_id * Land_num_of_states + XW] = check_max(0, y[1]);
    STATES[sample_id * Land_num_of_states + TRPN] = check_max(0, y[2]);
    STATES[sample_id * Land_num_of_states + TmBlocked] = y[3];
    STATES[sample_id * Land_num_of_states + ZETAS] = y[4];
    STATES[sample_id * Land_num_of_states + ZETAW] = y[5];
    STATES[sample_id * Land_num_of_states + dCd_dt] = y[6];

    CONSTANTS[sample_id * Land_num_of_constants + par_k] = 7;
    CONSTANTS[sample_id * Land_num_of_constants + b] = 9.1;
    CONSTANTS[sample_id * Land_num_of_constants + eta_l] = 200;
    CONSTANTS[sample_id * Land_num_of_constants + eta_s] = 20;
    CONSTANTS[sample_id * Land_num_of_constants + land_a] = 2.1;
    // printf("Successfully finish land_initConsts\n");
}

/**
 * @brief Compute the rates for the Land model.
 *
 * @param TIME Current time.
 * @param CONSTANTS Array of constant parameters.
 * @param RATES Array of rate values.
 * @param STATES Array of state variables.
 * @param ALGEBRAIC Array of algebraic variables.
 * @param y Array of state values.
 * @param sample_id sample_id for multi-threaded computation.
 */
__device__ void land_computeRates(double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC,
                                  double *y, int sample_id) {
    CONSTANTS[sample_id * Land_num_of_constants + lambda] =
        check_min(1.2, CONSTANTS[sample_id * Land_num_of_constants + lambda]);
    ALGEBRAIC[sample_id * Land_num_of_algebraic + Lfac] =
        check_max(0, 1 + CONSTANTS[sample_id * Land_num_of_constants + beta_0] *
                             (CONSTANTS[sample_id * Land_num_of_constants + lambda] +
                              check_min(0.87, CONSTANTS[sample_id * Land_num_of_constants + lambda]) - 1.87));

    ALGEBRAIC[sample_id * Land_num_of_algebraic + cdw] =
        CONSTANTS[sample_id * Land_num_of_constants + phi] * CONSTANTS[sample_id * Land_num_of_constants + k_uw] *
        (1 - CONSTANTS[sample_id * Land_num_of_constants + dr]) * (1 - CONSTANTS[sample_id * Land_num_of_constants + wfrac]) /
        ((1 - CONSTANTS[sample_id * Land_num_of_constants + dr]) * CONSTANTS[sample_id * Land_num_of_constants + wfrac]);
    ALGEBRAIC[sample_id * Land_num_of_algebraic + cds] =
        CONSTANTS[sample_id * Land_num_of_constants + phi] * CONSTANTS[sample_id * Land_num_of_constants + k_ws] *
        (1 - CONSTANTS[sample_id * Land_num_of_constants + dr]) * CONSTANTS[sample_id * Land_num_of_constants + wfrac] /
        CONSTANTS[sample_id * Land_num_of_constants + dr];

    ALGEBRAIC[sample_id * Land_num_of_algebraic + k_wu] =
        CONSTANTS[sample_id * Land_num_of_constants + k_uw] * (1 / CONSTANTS[sample_id * Land_num_of_constants + wfrac] - 1) -
        CONSTANTS[sample_id * Land_num_of_constants + k_ws];
    ALGEBRAIC[sample_id * Land_num_of_algebraic + k_su] = CONSTANTS[sample_id * Land_num_of_constants + k_ws] *
                                                       (1 / CONSTANTS[sample_id * Land_num_of_constants + dr] - 1) *
                                                       CONSTANTS[sample_id * Land_num_of_constants + wfrac];
    ALGEBRAIC[sample_id * Land_num_of_algebraic + A] =
        (0.25 * CONSTANTS[sample_id * Land_num_of_constants + TOT_A]) /
        ((1 - CONSTANTS[sample_id * Land_num_of_constants + dr]) * CONSTANTS[sample_id * Land_num_of_constants + wfrac] +
         CONSTANTS[sample_id * Land_num_of_constants + dr]) *
        (CONSTANTS[sample_id * Land_num_of_constants + dr] / 0.25);

    ALGEBRAIC[sample_id * Land_num_of_algebraic + XU] = (1 - STATES[sample_id * Land_num_of_states + TmBlocked]) -
                                                     STATES[sample_id * Land_num_of_states + XW] -
                                                     STATES[sample_id * Land_num_of_states + XS];

    ALGEBRAIC[sample_id * Land_num_of_algebraic + xb_ws] =
        CONSTANTS[sample_id * Land_num_of_constants + k_ws] * STATES[sample_id * Land_num_of_states + XW];
    ALGEBRAIC[sample_id * Land_num_of_algebraic + xb_uw] =
        CONSTANTS[sample_id * Land_num_of_constants + k_uw] * ALGEBRAIC[sample_id * Land_num_of_algebraic + XU];
    ALGEBRAIC[sample_id * Land_num_of_algebraic + xb_wu] =
        ALGEBRAIC[sample_id * Land_num_of_algebraic + k_wu] * STATES[sample_id * Land_num_of_states + XW];
    ALGEBRAIC[sample_id * Land_num_of_algebraic + xb_su] =
        ALGEBRAIC[sample_id * Land_num_of_algebraic + k_su] * STATES[sample_id * Land_num_of_states + XS];

    double temp_zetas1 =
        STATES[sample_id * Land_num_of_states + ZETAS] > 0 ? STATES[sample_id * Land_num_of_states + ZETAS] : 0;
    double temp_zetas2 =
        STATES[sample_id * Land_num_of_states + ZETAS] < -1 ? -STATES[sample_id * Land_num_of_states + ZETAS] - 1 : 0;
    ALGEBRAIC[sample_id * Land_num_of_algebraic + gamma_rate] =
        CONSTANTS[sample_id * Land_num_of_constants + gamma_idx] * check_max(temp_zetas1, temp_zetas2);

    ALGEBRAIC[sample_id * Land_num_of_algebraic + xb_su_gamma] =
        ALGEBRAIC[sample_id * Land_num_of_algebraic + gamma_rate] * STATES[sample_id * Land_num_of_states + XS];
    ALGEBRAIC[sample_id * Land_num_of_algebraic + gamma_rate_w] =
        CONSTANTS[sample_id * Land_num_of_constants + gamma_wu] * fabs(STATES[sample_id * Land_num_of_states + ZETAW]);
    ALGEBRAIC[sample_id * Land_num_of_algebraic + xb_wu_gamma] =
        ALGEBRAIC[sample_id * Land_num_of_algebraic + gamma_rate_w] * STATES[sample_id * Land_num_of_states + XW];

    RATES[sample_id * Land_num_of_rates + XS] = ALGEBRAIC[sample_id * Land_num_of_algebraic + xb_ws] -
                                             ALGEBRAIC[sample_id * Land_num_of_algebraic + xb_su] -
                                             ALGEBRAIC[sample_id * Land_num_of_algebraic + xb_wu_gamma];
    RATES[sample_id * Land_num_of_rates + XW] =
        ALGEBRAIC[sample_id * Land_num_of_algebraic + xb_uw] - ALGEBRAIC[sample_id * Land_num_of_algebraic + xb_wu] -
        ALGEBRAIC[sample_id * Land_num_of_algebraic + xb_ws] - ALGEBRAIC[sample_id * Land_num_of_algebraic + xb_wu_gamma];

    ALGEBRAIC[sample_id * Land_num_of_algebraic + ca50] +=
        CONSTANTS[sample_id * Land_num_of_constants + beta_1] *
        check_min(0.2, CONSTANTS[sample_id * Land_num_of_constants + lambda] - 1);
    RATES[sample_id * Land_num_of_rates + TRPN] =
        CONSTANTS[sample_id * Land_num_of_constants + koff] *
        (pow((CONSTANTS[sample_id * Land_num_of_constants + Cai] / ALGEBRAIC[sample_id * Land_num_of_algebraic + ca50]),
             CONSTANTS[sample_id * Land_num_of_constants + TRPN_n]) *
             (1 - STATES[sample_id * Land_num_of_states + TRPN]) -
         STATES[sample_id * Land_num_of_states + TRPN]);

    ALGEBRAIC[sample_id * Land_num_of_algebraic + XSSS] = CONSTANTS[sample_id * Land_num_of_constants + dr] * 0.5;
    ALGEBRAIC[sample_id * Land_num_of_algebraic + XWSS] =
        (1 - CONSTANTS[sample_id * Land_num_of_constants + dr]) * CONSTANTS[sample_id * Land_num_of_constants + wfrac] * 0.5;
    ALGEBRAIC[sample_id * Land_num_of_algebraic + ktm_block] =
        CONSTANTS[sample_id * Land_num_of_constants + ktm_unblock] *
        (pow(CONSTANTS[sample_id * Land_num_of_constants + perm50], CONSTANTS[sample_id * Land_num_of_constants + nperm]) *
         0.5) /
        (0.5 - ALGEBRAIC[sample_id * Land_num_of_algebraic + XSSS] - ALGEBRAIC[sample_id * Land_num_of_algebraic + XWSS]);

    RATES[sample_id * Land_num_of_rates + TmBlocked] =
        CONSTANTS[sample_id * Land_num_of_constants + ktm_block] *
            check_min(100, pow(STATES[sample_id * Land_num_of_states + TRPN],
                               -(CONSTANTS[sample_id * Land_num_of_constants + nperm] / 2))) *
            ALGEBRAIC[sample_id * Land_num_of_algebraic + XU] -
        CONSTANTS[sample_id * Land_num_of_constants + ktm_unblock] *
            pow(STATES[sample_id * Land_num_of_states + TRPN], (CONSTANTS[sample_id * Land_num_of_constants + nperm] / 2)) *
            STATES[sample_id * Land_num_of_states + TmBlocked];

    RATES[sample_id * Land_num_of_rates + ZETAS] =
        CONSTANTS[sample_id * Land_num_of_constants + A] * CONSTANTS[sample_id * Land_num_of_constants + dlambda_dt] -
        ALGEBRAIC[sample_id * Land_num_of_algebraic + cds] * STATES[sample_id * Land_num_of_states + ZETAS];
    RATES[sample_id * Land_num_of_rates + ZETAW] =
        CONSTANTS[sample_id * Land_num_of_constants + A] * CONSTANTS[sample_id * Land_num_of_constants + dlambda_dt] -
        ALGEBRAIC[sample_id * Land_num_of_algebraic + cdw] * STATES[sample_id * Land_num_of_states + ZETAW];

    CONSTANTS[sample_id * Land_num_of_constants + Cd] = y[6];
    CONSTANTS[sample_id * Land_num_of_constants + C] = CONSTANTS[sample_id * Land_num_of_constants + lambda] - 1;

    CONSTANTS[sample_id * Land_num_of_constants + eta] =
        CONSTANTS[sample_id * Land_num_of_constants + C] - CONSTANTS[sample_id * Land_num_of_constants + Cd] < 0
            ? CONSTANTS[sample_id * Land_num_of_constants + eta_s]
            : CONSTANTS[sample_id * Land_num_of_constants + eta_l];

    STATES[sample_id * Land_num_of_states + dCd_dt] =
        CONSTANTS[sample_id * Land_num_of_constants + par_k] *
        (CONSTANTS[sample_id * Land_num_of_constants + C] - CONSTANTS[sample_id * Land_num_of_constants + Cd]) /
        CONSTANTS[sample_id * Land_num_of_constants + eta];
    RATES[sample_id * Land_num_of_rates + dCd_dt] = STATES[sample_id * Land_num_of_states + dCd_dt];

    ALGEBRAIC[sample_id * Land_num_of_algebraic + Fd] =
        CONSTANTS[sample_id * Land_num_of_constants + eta] * STATES[sample_id * Land_num_of_states + dCd_dt];
    ALGEBRAIC[sample_id * Land_num_of_algebraic + F1] =
        (exp(CONSTANTS[sample_id * Land_num_of_constants + b] * CONSTANTS[sample_id * Land_num_of_constants + C]) - 1);
    ALGEBRAIC[sample_id * Land_num_of_algebraic + Tp] =
        CONSTANTS[sample_id * Land_num_of_constants + land_a] *
        (ALGEBRAIC[sample_id * Land_num_of_algebraic + F1] + ALGEBRAIC[sample_id * Land_num_of_algebraic + Fd]);

    ALGEBRAIC[sample_id * Land_num_of_algebraic + Ta] =
        ALGEBRAIC[sample_id * Land_num_of_algebraic + Lfac] *
        (CONSTANTS[sample_id * Land_num_of_constants + Tref] / CONSTANTS[sample_id * Land_num_of_constants + dr]) *
        ((STATES[sample_id * Land_num_of_states + ZETAS] + 1) * STATES[sample_id * Land_num_of_states + XS] +
         STATES[sample_id * Land_num_of_states + ZETAW] * STATES[sample_id * Land_num_of_states + XW]);
    ALGEBRAIC[sample_id * Land_num_of_algebraic + land_T] =
        ALGEBRAIC[sample_id * Land_num_of_algebraic + Ta] + ALGEBRAIC[sample_id * Land_num_of_algebraic + Tp];
}

/**
 * @brief Solve the Euler method for the Land model.
 *
 * @param dt Time step.
 * @param t Current time.
 * @param Cai_input Calcium input.
 * @param CONSTANTS Array of constant parameters.
 * @param RATES Array of rate values.
 * @param STATES Array of state variables.
 * @param sample_id sample_id for multi-threaded computation.
 */
__device__ void land_solveEuler(double dt, double t, double Cai_input, double *CONSTANTS, double *RATES, double *STATES,
                                int sample_id) {
    CONSTANTS[sample_id * Land_num_of_constants + Cai] = Cai_input;

    STATES[sample_id * Land_num_of_states + XS] += RATES[sample_id * Land_num_of_rates + XS] * dt;
    STATES[sample_id * Land_num_of_states + XW] += RATES[sample_id * Land_num_of_rates + XW] * dt;
    STATES[sample_id * Land_num_of_states + TRPN] += RATES[sample_id * Land_num_of_rates + TRPN] * dt;
    STATES[sample_id * Land_num_of_states + TmBlocked] += RATES[sample_id * Land_num_of_rates + TmBlocked] * dt;
    STATES[sample_id * Land_num_of_states + ZETAS] += RATES[sample_id * Land_num_of_rates + ZETAS] * dt;
    STATES[sample_id * Land_num_of_states + ZETAW] += RATES[sample_id * Land_num_of_rates + ZETAW] * dt;
    STATES[sample_id * Land_num_of_states + dCd_dt] += RATES[sample_id * Land_num_of_rates + dCd_dt] * dt;
}
