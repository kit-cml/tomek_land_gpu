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
 * @param offset Offset for multi-threaded computation.
 */
__device__ void land_initConsts(bool is_skinned, bool BETA, double *y, double *CONSTANTS, double *RATES, double *STATES,
                                double *ALGEBRAIC, int offset) {
    printf("Successfully go to inside land_initConsts\n");
    CONSTANTS[offset * Land_num_of_constants + dlambda_dt] = 0;
    CONSTANTS[offset * Land_num_of_constants + lambda] = 1.0;
    CONSTANTS[offset * Land_num_of_constants + Cai] = 0.0;

    RATES[offset * Land_num_of_rates + XS] = 0;
    RATES[offset * Land_num_of_rates + XW] = 0;
    RATES[offset * Land_num_of_rates + TRPN] = 0;
    RATES[offset * Land_num_of_rates + TmBlocked] = 1;
    RATES[offset * Land_num_of_rates + ZETAS] = 0;
    RATES[offset * Land_num_of_rates + ZETAW] = 0;
    RATES[offset * Land_num_of_rates + dCd_dt] = 0;

    CONSTANTS[offset * Land_num_of_constants + lambda] =
        check_min(1.2, CONSTANTS[offset * Land_num_of_constants + lambda]);

    CONSTANTS[offset * Land_num_of_constants + perm50] = 0.35;
    CONSTANTS[offset * Land_num_of_constants + TRPN_n] = 2;
    CONSTANTS[offset * Land_num_of_constants + koff] = 0.1;
    CONSTANTS[offset * Land_num_of_constants + dr] = 0.25;
    CONSTANTS[offset * Land_num_of_constants + wfrac] = 0.5;
    CONSTANTS[offset * Land_num_of_constants + TOT_A] = 25;
    CONSTANTS[offset * Land_num_of_constants + ktm_unblock] = 0.04;
    CONSTANTS[offset * Land_num_of_constants + beta_1] = -2.4;
    CONSTANTS[offset * Land_num_of_constants + beta_0] = 2.3;
    CONSTANTS[offset * Land_num_of_constants + gamma_idx] = 0.0085;
    CONSTANTS[offset * Land_num_of_constants + gamma_wu] = 0.615;
    CONSTANTS[offset * Land_num_of_constants + phi] = 2.23;

    if (is_skinned) {
        CONSTANTS[offset * Land_num_of_constants + nperm] = 2.2;
        ALGEBRAIC[offset * Land_num_of_algebraic + ca50] = 2.5;
        CONSTANTS[offset * Land_num_of_constants + Tref] = 40.5;
        CONSTANTS[offset * Land_num_of_constants + nu] = 1;
        CONSTANTS[offset * Land_num_of_constants + mu] = 1;
    } else {
        CONSTANTS[offset * Land_num_of_constants + nperm] = 2.4;
        ALGEBRAIC[offset * Land_num_of_algebraic + ca50] = 0.805;
        CONSTANTS[offset * Land_num_of_constants + Tref] = 120.0;
        CONSTANTS[offset * Land_num_of_constants + nu] = 7;
        CONSTANTS[offset * Land_num_of_constants + mu] = 3;
    }

    if (BETA) {
        // If BETA is true, use beta values (currently commented out)
        // CONSTANTS[offset * Land_num_of_constants + beta_1] = beta[1];
        // CONSTANTS[offset * Land_num_of_constants + beta_0] = beta[0];
    }

    CONSTANTS[offset * Land_num_of_constants + k_ws] = 0.004 * CONSTANTS[offset * Land_num_of_constants + mu];
    CONSTANTS[offset * Land_num_of_constants + k_uw] = 0.026 * CONSTANTS[offset * Land_num_of_constants + nu];

    // Initialize STATES with checks
    STATES[offset * Land_num_of_states + XS] = check_max(0, y[0]);
    STATES[offset * Land_num_of_states + XW] = check_max(0, y[1]);
    STATES[offset * Land_num_of_states + TRPN] = check_max(0, y[2]);
    STATES[offset * Land_num_of_states + TmBlocked] = y[3];
    STATES[offset * Land_num_of_states + ZETAS] = y[4];
    STATES[offset * Land_num_of_states + ZETAW] = y[5];
    STATES[offset * Land_num_of_states + dCd_dt] = y[6];

    CONSTANTS[offset * Land_num_of_constants + par_k] = 7;
    CONSTANTS[offset * Land_num_of_constants + b] = 9.1;
    CONSTANTS[offset * Land_num_of_constants + eta_l] = 200;
    CONSTANTS[offset * Land_num_of_constants + eta_s] = 20;
    CONSTANTS[offset * Land_num_of_constants + land_a] = 2.1;
    printf("Successfully finish land_initConsts\n");
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
 * @param offset Offset for multi-threaded computation.
 */
__device__ void land_computeRates(double TIME, double *CONSTANTS, double *RATES, double *STATES, double *ALGEBRAIC,
                                  double *y, int offset) {
    CONSTANTS[offset * Land_num_of_constants + lambda] =
        check_min(1.2, CONSTANTS[offset * Land_num_of_constants + lambda]);
    ALGEBRAIC[offset * Land_num_of_algebraic + Lfac] =
        check_max(0, 1 + CONSTANTS[offset * Land_num_of_constants + beta_0] *
                             (CONSTANTS[offset * Land_num_of_constants + lambda] +
                              check_min(0.87, CONSTANTS[offset * Land_num_of_constants + lambda]) - 1.87));

    ALGEBRAIC[offset * Land_num_of_algebraic + cdw] =
        CONSTANTS[offset * Land_num_of_constants + phi] * CONSTANTS[offset * Land_num_of_constants + k_uw] *
        (1 - CONSTANTS[offset * Land_num_of_constants + dr]) * (1 - CONSTANTS[offset * Land_num_of_constants + wfrac]) /
        ((1 - CONSTANTS[offset * Land_num_of_constants + dr]) * CONSTANTS[offset * Land_num_of_constants + wfrac]);
    ALGEBRAIC[offset * Land_num_of_algebraic + cds] =
        CONSTANTS[offset * Land_num_of_constants + phi] * CONSTANTS[offset * Land_num_of_constants + k_ws] *
        (1 - CONSTANTS[offset * Land_num_of_constants + dr]) * CONSTANTS[offset * Land_num_of_constants + wfrac] /
        CONSTANTS[offset * Land_num_of_constants + dr];

    ALGEBRAIC[offset * Land_num_of_algebraic + k_wu] =
        CONSTANTS[offset * Land_num_of_constants + k_uw] * (1 / CONSTANTS[offset * Land_num_of_constants + wfrac] - 1) -
        CONSTANTS[offset * Land_num_of_constants + k_ws];
    ALGEBRAIC[offset * Land_num_of_algebraic + k_su] = CONSTANTS[offset * Land_num_of_constants + k_ws] *
                                                       (1 / CONSTANTS[offset * Land_num_of_constants + dr] - 1) *
                                                       CONSTANTS[offset * Land_num_of_constants + wfrac];
    ALGEBRAIC[offset * Land_num_of_algebraic + A] =
        (0.25 * CONSTANTS[offset * Land_num_of_constants + TOT_A]) /
        ((1 - CONSTANTS[offset * Land_num_of_constants + dr]) * CONSTANTS[offset * Land_num_of_constants + wfrac] +
         CONSTANTS[offset * Land_num_of_constants + dr]) *
        (CONSTANTS[offset * Land_num_of_constants + dr] / 0.25);

    ALGEBRAIC[offset * Land_num_of_algebraic + XU] = (1 - STATES[offset * Land_num_of_states + TmBlocked]) -
                                                     STATES[offset * Land_num_of_states + XW] -
                                                     STATES[offset * Land_num_of_states + XS];

    ALGEBRAIC[offset * Land_num_of_algebraic + xb_ws] =
        CONSTANTS[offset * Land_num_of_constants + k_ws] * STATES[offset * Land_num_of_states + XW];
    ALGEBRAIC[offset * Land_num_of_algebraic + xb_uw] =
        CONSTANTS[offset * Land_num_of_constants + k_uw] * ALGEBRAIC[offset * Land_num_of_algebraic + XU];
    ALGEBRAIC[offset * Land_num_of_algebraic + xb_wu] =
        ALGEBRAIC[offset * Land_num_of_algebraic + k_wu] * STATES[offset * Land_num_of_states + XW];
    ALGEBRAIC[offset * Land_num_of_algebraic + xb_su] =
        ALGEBRAIC[offset * Land_num_of_algebraic + k_su] * STATES[offset * Land_num_of_states + XS];

    double temp_zetas1 =
        STATES[offset * Land_num_of_states + ZETAS] > 0 ? STATES[offset * Land_num_of_states + ZETAS] : 0;
    double temp_zetas2 =
        STATES[offset * Land_num_of_states + ZETAS] < -1 ? -STATES[offset * Land_num_of_states + ZETAS] - 1 : 0;
    ALGEBRAIC[offset * Land_num_of_algebraic + gamma_rate] =
        CONSTANTS[offset * Land_num_of_constants + gamma_idx] * check_max(temp_zetas1, temp_zetas2);

    ALGEBRAIC[offset * Land_num_of_algebraic + xb_su_gamma] =
        ALGEBRAIC[offset * Land_num_of_algebraic + gamma_rate] * STATES[offset * Land_num_of_states + XS];
    ALGEBRAIC[offset * Land_num_of_algebraic + gamma_rate_w] =
        CONSTANTS[offset * Land_num_of_constants + gamma_wu] * fabs(STATES[offset * Land_num_of_states + ZETAW]);
    ALGEBRAIC[offset * Land_num_of_algebraic + xb_wu_gamma] =
        ALGEBRAIC[offset * Land_num_of_algebraic + gamma_rate_w] * STATES[offset * Land_num_of_states + XW];

    RATES[offset * Land_num_of_rates + XS] = ALGEBRAIC[offset * Land_num_of_algebraic + xb_ws] -
                                             ALGEBRAIC[offset * Land_num_of_algebraic + xb_su] -
                                             ALGEBRAIC[offset * Land_num_of_algebraic + xb_wu_gamma];
    RATES[offset * Land_num_of_rates + XW] =
        ALGEBRAIC[offset * Land_num_of_algebraic + xb_uw] - ALGEBRAIC[offset * Land_num_of_algebraic + xb_wu] -
        ALGEBRAIC[offset * Land_num_of_algebraic + xb_ws] - ALGEBRAIC[offset * Land_num_of_algebraic + xb_wu_gamma];

    ALGEBRAIC[offset * Land_num_of_algebraic + ca50] +=
        CONSTANTS[offset * Land_num_of_constants + beta_1] *
        check_min(0.2, CONSTANTS[offset * Land_num_of_constants + lambda] - 1);
    RATES[offset * Land_num_of_rates + TRPN] =
        CONSTANTS[offset * Land_num_of_constants + koff] *
        (pow((CONSTANTS[offset * Land_num_of_constants + Cai] / ALGEBRAIC[offset * Land_num_of_algebraic + ca50]),
             CONSTANTS[offset * Land_num_of_constants + TRPN_n]) *
             (1 - STATES[offset * Land_num_of_states + TRPN]) -
         STATES[offset * Land_num_of_states + TRPN]);

    ALGEBRAIC[offset * Land_num_of_algebraic + XSSS] = CONSTANTS[offset * Land_num_of_constants + dr] * 0.5;
    ALGEBRAIC[offset * Land_num_of_algebraic + XWSS] =
        (1 - CONSTANTS[offset * Land_num_of_constants + dr]) * CONSTANTS[offset * Land_num_of_constants + wfrac] * 0.5;
    ALGEBRAIC[offset * Land_num_of_algebraic + ktm_block] =
        CONSTANTS[offset * Land_num_of_constants + ktm_unblock] *
        (pow(CONSTANTS[offset * Land_num_of_constants + perm50], CONSTANTS[offset * Land_num_of_constants + nperm]) *
         0.5) /
        (0.5 - ALGEBRAIC[offset * Land_num_of_algebraic + XSSS] - ALGEBRAIC[offset * Land_num_of_algebraic + XWSS]);

    RATES[offset * Land_num_of_rates + TmBlocked] =
        CONSTANTS[offset * Land_num_of_constants + ktm_block] *
            check_min(100, pow(STATES[offset * Land_num_of_states + TRPN],
                               -(CONSTANTS[offset * Land_num_of_constants + nperm] / 2))) *
            ALGEBRAIC[offset * Land_num_of_algebraic + XU] -
        CONSTANTS[offset * Land_num_of_constants + ktm_unblock] *
            pow(STATES[offset * Land_num_of_states + TRPN], (CONSTANTS[offset * Land_num_of_constants + nperm] / 2)) *
            STATES[offset * Land_num_of_states + TmBlocked];

    RATES[offset * Land_num_of_rates + ZETAS] =
        CONSTANTS[offset * Land_num_of_constants + A] * CONSTANTS[offset * Land_num_of_constants + dlambda_dt] -
        ALGEBRAIC[offset * Land_num_of_algebraic + cds] * STATES[offset * Land_num_of_states + ZETAS];
    RATES[offset * Land_num_of_rates + ZETAW] =
        CONSTANTS[offset * Land_num_of_constants + A] * CONSTANTS[offset * Land_num_of_constants + dlambda_dt] -
        ALGEBRAIC[offset * Land_num_of_algebraic + cdw] * STATES[offset * Land_num_of_states + ZETAW];

    CONSTANTS[offset * Land_num_of_constants + Cd] = y[6];
    CONSTANTS[offset * Land_num_of_constants + C] = CONSTANTS[offset * Land_num_of_constants + lambda] - 1;

    CONSTANTS[offset * Land_num_of_constants + eta] =
        CONSTANTS[offset * Land_num_of_constants + C] - CONSTANTS[offset * Land_num_of_constants + Cd] < 0
            ? CONSTANTS[offset * Land_num_of_constants + eta_s]
            : CONSTANTS[offset * Land_num_of_constants + eta_l];

    STATES[offset * Land_num_of_states + dCd_dt] =
        CONSTANTS[offset * Land_num_of_constants + par_k] *
        (CONSTANTS[offset * Land_num_of_constants + C] - CONSTANTS[offset * Land_num_of_constants + Cd]) /
        CONSTANTS[offset * Land_num_of_constants + eta];
    RATES[offset * Land_num_of_rates + dCd_dt] = STATES[offset * Land_num_of_states + dCd_dt];

    ALGEBRAIC[offset * Land_num_of_algebraic + Fd] =
        CONSTANTS[offset * Land_num_of_constants + eta] * STATES[offset * Land_num_of_states + dCd_dt];
    ALGEBRAIC[offset * Land_num_of_algebraic + F1] =
        (exp(CONSTANTS[offset * Land_num_of_constants + b] * CONSTANTS[offset * Land_num_of_constants + C]) - 1);
    ALGEBRAIC[offset * Land_num_of_algebraic + Tp] =
        CONSTANTS[offset * Land_num_of_constants + land_a] *
        (ALGEBRAIC[offset * Land_num_of_algebraic + F1] + ALGEBRAIC[offset * Land_num_of_algebraic + Fd]);

    ALGEBRAIC[offset * Land_num_of_algebraic + Ta] =
        ALGEBRAIC[offset * Land_num_of_algebraic + Lfac] *
        (CONSTANTS[offset * Land_num_of_constants + Tref] / CONSTANTS[offset * Land_num_of_constants + dr]) *
        ((STATES[offset * Land_num_of_states + ZETAS] + 1) * STATES[offset * Land_num_of_states + XS] +
         STATES[offset * Land_num_of_states + ZETAW] * STATES[offset * Land_num_of_states + XW]);
    ALGEBRAIC[offset * Land_num_of_algebraic + land_T] =
        ALGEBRAIC[offset * Land_num_of_algebraic + Ta] + ALGEBRAIC[offset * Land_num_of_algebraic + Tp];
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
 * @param offset Offset for multi-threaded computation.
 */
__device__ void land_solveEuler(double dt, double t, double Cai_input, double *CONSTANTS, double *RATES, double *STATES,
                                int offset) {
    CONSTANTS[offset * Land_num_of_constants + Cai] = Cai_input;

    STATES[offset * Land_num_of_states + XS] += RATES[offset * Land_num_of_rates + XS] * dt;
    STATES[offset * Land_num_of_states + XW] += RATES[offset * Land_num_of_rates + XW] * dt;
    STATES[offset * Land_num_of_states + TRPN] += RATES[offset * Land_num_of_rates + TRPN] * dt;
    STATES[offset * Land_num_of_states + TmBlocked] += RATES[offset * Land_num_of_rates + TmBlocked] * dt;
    STATES[offset * Land_num_of_states + ZETAS] += RATES[offset * Land_num_of_rates + ZETAS] * dt;
    STATES[offset * Land_num_of_states + ZETAW] += RATES[offset * Land_num_of_rates + ZETAW] * dt;
    STATES[offset * Land_num_of_states + dCd_dt] += RATES[offset * Land_num_of_rates + dCd_dt] * dt;
}
