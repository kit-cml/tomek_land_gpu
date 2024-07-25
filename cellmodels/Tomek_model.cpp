/*
   There are a total of 223 entries in the algebraic variable array.
   There are a total of 43 entries in each of the rate and state variable arrays.
   There are a total of 163 entries in the constant variable array.
 */

#include "Tomek_model.hpp"
#include <cmath>
#include <cstdlib>
#include "../enums/enum_Tomek_model.hpp"
#include <iostream>
#include "../utils/constants.hpp"
// #include "../../functions/inputoutput.hpp"

/*
 * TIME is time in component environment (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype] is celltype in component environment (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + nao] is nao in component extracellular (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + cao] is cao in component extracellular (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + ko] is ko in component extracellular (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + clo] is clo in component extracellular (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + R] is R in component physical_constants (joule_per_kilomole_kelvin).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + T] is T in component physical_constants (kelvin).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + F] is F in component physical_constants (coulomb_per_mole).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + zna] is zna in component physical_constants (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + zca] is zca in component physical_constants (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + zk] is zk in component physical_constants (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + zcl] is zcl in component physical_constants (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + L] is L in component cell_geometry (centimeter).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + rad] is rad in component cell_geometry (centimeter).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + vcell] is vcell in component cell_geometry (microliter).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Ageo] is Ageo in component cell_geometry (centimeter_squared).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap] is Acap in component cell_geometry (centimeter_squared).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo] is vmyo in component cell_geometry (microliter).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + vnsr] is vnsr in component cell_geometry (microliter).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + vjsr] is vjsr in component cell_geometry (microliter).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + vss] is vss in component cell_geometry (microliter).
 * STATES[V] is v in component membrane (millivolt).
 * ALGEBRAIC[vffrt] is vffrt in component membrane (coulomb_per_mole).
 * ALGEBRAIC[vfrt] is vfrt in component membrane (dimensionless).
 * ALGEBRAIC[INa] is INa in component INa (microA_per_microF).
 * ALGEBRAIC[INaL] is INaL in component INaL (microA_per_microF).
 * ALGEBRAIC[Ito] is Ito in component Ito (microA_per_microF).
 * ALGEBRAIC[ICaL] is ICaL in component ICaL (microA_per_microF).
 * ALGEBRAIC[ICaNa] is ICaNa in component ICaL (microA_per_microF).
 * ALGEBRAIC[ICaK] is ICaK in component ICaL (microA_per_microF).
 * ALGEBRAIC[IKr] is IKr in component IKr (microA_per_microF).
 * ALGEBRAIC[IKs] is IKs in component IKs (microA_per_microF).
 * ALGEBRAIC[IK1] is IK1 in component IK1 (microA_per_microF).
 * ALGEBRAIC[INaCa_i] is INaCa_i in component INaCa (microA_per_microF).
 * ALGEBRAIC[INaCa_ss] is INaCa_ss in component INaCa (microA_per_microF).
 * ALGEBRAIC[INaK] is INaK in component INaK (microA_per_microF).
 * ALGEBRAIC[INab] is INab in component INab (microA_per_microF).
 * ALGEBRAIC[IKb] is IKb in component IKb (microA_per_microF).
 * ALGEBRAIC[IpCa] is IpCa in component IpCa (microA_per_microF).
 * ALGEBRAIC[ICab] is ICab in component ICab (microA_per_microF).
 * ALGEBRAIC[IClCa] is IClCa in component ICl (microA_per_microF).
 * ALGEBRAIC[IClb] is IClb in component ICl (microA_per_microF).
 * ALGEBRAIC[I_katp] is I_katp in component I_katp (microA_per_microF).
 * ALGEBRAIC[Istim] is Istim in component membrane (microA_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + stim_start] is stim_start in component membrane (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + i_Stim_End] is i_Stim_End in component membrane (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + i_Stim_Amplitude] is i_Stim_Amplitude in component membrane (microA_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + BCL] is BCL in component membrane (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + i_Stim_PulseDuration] is i_Stim_PulseDuration in component membrane (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaMK] is KmCaMK in component CaMK (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + aCaMK] is aCaMK in component CaMK (per_millimolar_per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + bCaMK] is bCaMK in component CaMK (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + CaMKo] is CaMKo in component CaMK (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaM] is KmCaM in component CaMK (millimolar).
 * ALGEBRAIC[CaMKb] is CaMKb in component CaMK (millimolar).
 * ALGEBRAIC[CaMKa] is CaMKa in component CaMK (millimolar).
 * STATES[CaMKt] is CaMKt in component CaMK (millimolar).
 * STATES[cass] is cass in component intracellular_ions (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + cmdnmax_b] is cmdnmax_b in component intracellular_ions (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + cmdnmax] is cmdnmax in component intracellular_ions (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcmdn] is kmcmdn in component intracellular_ions (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + trpnmax] is trpnmax in component intracellular_ions (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + kmtrpn] is kmtrpn in component intracellular_ions (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + BSRmax] is BSRmax in component intracellular_ions (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + KmBSR] is KmBSR in component intracellular_ions (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + BSLmax] is BSLmax in component intracellular_ions (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + KmBSL] is KmBSL in component intracellular_ions (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + csqnmax] is csqnmax in component intracellular_ions (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcsqn] is kmcsqn in component intracellular_ions (millimolar).
 * STATES[nai] is nai in component intracellular_ions (millimolar).
 * STATES[nass] is nass in component intracellular_ions (millimolar).
 * STATES[ki] is ki in component intracellular_ions (millimolar).
 * STATES[kss] is kss in component intracellular_ions (millimolar).
 * STATES[cansr] is cansr in component intracellular_ions (millimolar).
 * STATES[cajsr] is cajsr in component intracellular_ions (millimolar).
 * STATES[cai] is cai in component intracellular_ions (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + cli] is cli in component intracellular_ions (millimolar).
 * ALGEBRAIC[ICaL_ss] is ICaL_ss in component ICaL (microA_per_microF).
 * ALGEBRAIC[ICaNa_ss] is ICaNa_ss in component ICaL (microA_per_microF).
 * ALGEBRAIC[ICaK_ss] is ICaK_ss in component ICaL (microA_per_microF).
 * ALGEBRAIC[ICaL_i] is ICaL_i in component ICaL (microA_per_microF).
 * ALGEBRAIC[ICaNa_i] is ICaNa_i in component ICaL (microA_per_microF).
 * ALGEBRAIC[ICaK_i] is ICaK_i in component ICaL (microA_per_microF).
 * ALGEBRAIC[JdiffNa] is JdiffNa in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[Jdiff] is Jdiff in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[Jup] is Jup in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[JdiffK] is JdiffK in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[Jrel] is Jrel in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[Jtr] is Jtr in component trans_flux (millimolar_per_millisecond).
 * ALGEBRAIC[Bcai] is Bcai in component intracellular_ions (dimensionless).
 * ALGEBRAIC[Bcajsr] is Bcajsr in component intracellular_ions (dimensionless).
 * ALGEBRAIC[Bcass] is Bcass in component intracellular_ions (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PKNa] is PKNa in component reversal_potentials (dimensionless).
 * ALGEBRAIC[ENa] is ENa in component reversal_potentials (millivolt).
 * ALGEBRAIC[EK] is EK in component reversal_potentials (millivolt).
 * ALGEBRAIC[EKs] is EKs in component reversal_potentials (millivolt).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + ECl] is ECl in component reversal_potentials (millivolt).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + gkatp] is gkatp in component I_katp (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + fkatp] is fkatp in component I_katp (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + K_o_n] is K_o_n in component I_katp (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + A_atp] is A_atp in component I_katp (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + K_atp] is K_atp in component I_katp (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + akik] is akik in component I_katp (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + bkik] is bkik in component I_katp (dimensionless).
 * ALGEBRAIC[mss] is mss in component INa (dimensionless).
 * ALGEBRAIC[tm] is tm in component INa (millisecond).
 * STATES[m] is m in component INa (dimensionless).
 * ALGEBRAIC[hss] is hss in component INa (dimensionless).
 * ALGEBRAIC[ah] is ah in component INa (dimensionless).
 * ALGEBRAIC[bh] is bh in component INa (dimensionless).
 * ALGEBRAIC[th] is th in component INa (millisecond).
 * STATES[h] is h in component INa (dimensionless).
 * ALGEBRAIC[jss] is jss in component INa (dimensionless).
 * ALGEBRAIC[aj] is aj in component INa (dimensionless).
 * ALGEBRAIC[bj] is bj in component INa (dimensionless).
 * ALGEBRAIC[tj] is tj in component INa (millisecond).
 * STATES[j] is j in component INa (dimensionless).
 * ALGEBRAIC[hssp] is hssp in component INa (dimensionless).
 * STATES[hp] is hp in component INa (dimensionless).
 * ALGEBRAIC[tjp] is tjp in component INa (millisecond).
 * STATES[jp] is jp in component INa (dimensionless).
 * ALGEBRAIC[fINap] is fINap in component INa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GNa] is GNa in component INa (milliS_per_microF).
 * ALGEBRAIC[mLss] is mLss in component INaL (dimensionless).
 * ALGEBRAIC[tmL] is tmL in component INaL (millisecond).
 * STATES[mL] is mL in component INaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + thL] is thL in component INaL (millisecond).
 * ALGEBRAIC[hLss] is hLss in component INaL (dimensionless).
 * STATES[hL] is hL in component INaL (dimensionless).
 * ALGEBRAIC[hLssp] is hLssp in component INaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + thLp] is thLp in component INaL (millisecond).
 * STATES[hLp] is hLp in component INaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL_b] is GNaL_b in component INaL (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL] is GNaL in component INaL (milliS_per_microF).
 * ALGEBRAIC[fINaLp] is fINaLp in component INaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Gto_b] is Gto_b in component Ito (milliS_per_microF).
 * ALGEBRAIC[ass] is ass in component Ito (dimensionless).
 * ALGEBRAIC[ta] is ta in component Ito (millisecond).
 * STATES[a] is a in component Ito (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift] is EKshift in component Ito (millivolt).
 * ALGEBRAIC[iss] is iss in component Ito (dimensionless).
 * ALGEBRAIC[delta_epi] is delta_epi in component Ito (dimensionless).
 * ALGEBRAIC[tiF_b] is tiF_b in component Ito (millisecond).
 * ALGEBRAIC[tiS_b] is tiS_b in component Ito (millisecond).
 * ALGEBRAIC[tiF] is tiF in component Ito (millisecond).
 * ALGEBRAIC[tiS] is tiS in component Ito (millisecond).
 * ALGEBRAIC[AiF] is AiF in component Ito (dimensionless).
 * ALGEBRAIC[AiS] is AiS in component Ito (dimensionless).
 * STATES[iF] is iF in component Ito (dimensionless).
 * STATES[iS] is iS in component Ito (dimensionless).
 * ALGEBRAIC[i] is i in component Ito (dimensionless).
 * ALGEBRAIC[assp] is assp in component Ito (dimensionless).
 * STATES[ap] is ap in component Ito (dimensionless).
 * ALGEBRAIC[dti_develop] is dti_develop in component Ito (dimensionless).
 * ALGEBRAIC[dti_recover] is dti_recover in component Ito (dimensionless).
 * ALGEBRAIC[tiFp] is tiFp in component Ito (millisecond).
 * ALGEBRAIC[tiSp] is tiSp in component Ito (millisecond).
 * STATES[iFp] is iFp in component Ito (dimensionless).
 * STATES[iSp] is iSp in component Ito (dimensionless).
 * ALGEBRAIC[ip] is ip in component Ito (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Gto] is Gto in component Ito (milliS_per_microF).
 * ALGEBRAIC[fItop] is fItop in component Ito (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Kmn] is Kmn in component ICaL (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] is k2n in component ICaL (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa_b] is PCa_b in component ICaL (dimensionless).
 * ALGEBRAIC[dss] is dss in component ICaL (dimensionless).
 * STATES[d] is d in component ICaL (dimensionless).
 * ALGEBRAIC[fss] is fss in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Aff] is Aff in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Afs] is Afs in component ICaL (dimensionless).
 * STATES[ff] is ff in component ICaL (dimensionless).
 * STATES[fs] is fs in component ICaL (dimensionless).
 * ALGEBRAIC[f] is f in component ICaL (dimensionless).
 * ALGEBRAIC[fcass] is fcass in component ICaL (dimensionless).
 * ALGEBRAIC[jcass] is jcass in component ICaL (dimensionless).
 * ALGEBRAIC[Afcaf] is Afcaf in component ICaL (dimensionless).
 * ALGEBRAIC[Afcas] is Afcas in component ICaL (dimensionless).
 * STATES[fcaf] is fcaf in component ICaL (dimensionless).
 * STATES[fcas] is fcas in component ICaL (dimensionless).
 * ALGEBRAIC[fca] is fca in component ICaL (dimensionless).
 * STATES[jca] is jca in component ICaL (dimensionless).
 * STATES[ffp] is ffp in component ICaL (dimensionless).
 * ALGEBRAIC[fp] is fp in component ICaL (dimensionless).
 * STATES[fcafp] is fcafp in component ICaL (dimensionless).
 * ALGEBRAIC[fcap] is fcap in component ICaL (dimensionless).
 * ALGEBRAIC[km2n] is km2n in component ICaL (per_millisecond).
 * ALGEBRAIC[anca_ss] is anca_ss in component ICaL (dimensionless).
 * STATES[nca_ss] is nca_ss in component ICaL (dimensionless).
 * ALGEBRAIC[anca_i] is anca_i in component ICaL (dimensionless).
 * STATES[nca_i] is nca_i in component ICaL (dimensionless).
 * ALGEBRAIC[PhiCaL_ss] is PhiCaL_ss in component ICaL (dimensionless).
 * ALGEBRAIC[PhiCaNa_ss] is PhiCaNa_ss in component ICaL (dimensionless).
 * ALGEBRAIC[PhiCaK_ss] is PhiCaK_ss in component ICaL (dimensionless).
 * ALGEBRAIC[PhiCaL_i] is PhiCaL_i in component ICaL (dimensionless).
 * ALGEBRAIC[PhiCaNa_i] is PhiCaNa_i in component ICaL (dimensionless).
 * ALGEBRAIC[PhiCaK_i] is PhiCaK_i in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa] is PCa in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCap] is PCap in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaNa] is PCaNa in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaK] is PCaK in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaNap] is PCaNap in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaKp] is PCaKp in component ICaL (dimensionless).
 * ALGEBRAIC[fICaLp] is fICaLp in component ICaL (dimensionless).
 * ALGEBRAIC[td] is td in component ICaL (millisecond).
 * ALGEBRAIC[tff] is tff in component ICaL (millisecond).
 * ALGEBRAIC[tfs] is tfs in component ICaL (millisecond).
 * ALGEBRAIC[tfcaf] is tfcaf in component ICaL (millisecond).
 * ALGEBRAIC[tfcas] is tfcas in component ICaL (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + tjca] is tjca in component ICaL (millisecond).
 * ALGEBRAIC[tffp] is tffp in component ICaL (millisecond).
 * ALGEBRAIC[tfcafp] is tfcafp in component ICaL (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + vShift] is vShift in component ICaL (millivolt).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + sample_id] is sample_id in component ICaL (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Io] is Io in component ICaL (dimensionless).
 * ALGEBRAIC[Iss] is Iss in component ICaL (dimensionless).
 * ALGEBRAIC[Ii] is Ii in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + dielConstant] is dielConstant in component ICaL (per_kelvin).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + constA] is constA in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_cao] is gamma_cao in component ICaL (dimensionless).
 * ALGEBRAIC[gamma_cass] is gamma_cass in component ICaL (dimensionless).
 * ALGEBRAIC[gamma_cai] is gamma_cai in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_nao] is gamma_nao in component ICaL (dimensionless).
 * ALGEBRAIC[gamma_nass] is gamma_nass in component ICaL (dimensionless).
 * ALGEBRAIC[gamma_nai] is gamma_nai in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_ko] is gamma_ko in component ICaL (dimensionless).
 * ALGEBRAIC[gamma_kss] is gamma_kss in component ICaL (dimensionless).
 * ALGEBRAIC[gamma_ki] is gamma_ki in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS] is ICaL_fractionSS in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr_b] is GKr_b in component IKr (milliS_per_microF).
 * STATES[C1] is C1 in component IKr (dimensionless).
 * STATES[C2] is C2 in component IKr (dimensionless).
 * STATES[C3] is C3 in component IKr (dimensionless).
 * STATES[I] is I in component IKr (dimensionless).
 * STATES[O] is O in component IKr (dimensionless).
 * ALGEBRAIC[alpha] is alpha in component IKr (per_millisecond).
 * ALGEBRAIC[beta] is beta in component IKr (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + alpha_1] is alpha_1 in component IKr (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + beta_1] is beta_1 in component IKr (per_millisecond).
 * ALGEBRAIC[alpha_2] is alpha_2 in component IKr (per_millisecond).
 * ALGEBRAIC[beta_2] is beta_2 in component IKr (per_millisecond).
 * ALGEBRAIC[alpha_i] is alpha_i in component IKr (per_millisecond).
 * ALGEBRAIC[beta_i] is beta_i in component IKr (per_millisecond).
 * ALGEBRAIC[alpha_C2ToI] is alpha_C2ToI in component IKr (per_millisecond).
 * ALGEBRAIC[beta_ItoC2] is beta_ItoC2 in component IKr (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr] is GKr in component IKr (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs_b] is GKs_b in component IKs (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs] is GKs in component IKs (milliS_per_microF).
 * ALGEBRAIC[xs1ss] is xs1ss in component IKs (dimensionless).
 * ALGEBRAIC[xs2ss] is xs2ss in component IKs (dimensionless).
 * ALGEBRAIC[txs1] is txs1 in component IKs (millisecond).
 * STATES[xs1] is xs1 in component IKs (dimensionless).
 * STATES[xs2] is xs2 in component IKs (dimensionless).
 * ALGEBRAIC[KsCa] is KsCa in component IKs (dimensionless).
 * ALGEBRAIC[txs2] is txs2 in component IKs (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1] is GK1 in component IK1 (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1_b] is GK1_b in component IK1 (milliS_per_microF).
 * ALGEBRAIC[aK1] is aK1 in component IK1 (dimensionless).
 * ALGEBRAIC[bK1] is bK1 in component IK1 (dimensionless).
 * ALGEBRAIC[K1ss] is K1ss in component IK1 (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + INaCa_fractionSS] is INaCa_fractionSS in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + kna1] is kna1 in component INaCa (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + kna2] is kna2 in component INaCa (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3] is kna3 in component INaCa (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + kasymm] is kasymm in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + wna] is wna in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + wca] is wca in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + wnaca] is wnaca in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + kcaon] is kcaon in component INaCa (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + kcaoff] is kcaoff in component INaCa (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + qna] is qna in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + qca] is qca in component INaCa (dimensionless).
 * ALGEBRAIC[hna] is hna in component INaCa (dimensionless).
 * ALGEBRAIC[hca] is hca in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaAct] is KmCaAct in component INaCa (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Gncx_b] is Gncx_b in component INaCa (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Gncx] is Gncx in component INaCa (milliS_per_microF).
 * ALGEBRAIC[h1_i] is h1_i in component INaCa (dimensionless).
 * ALGEBRAIC[h2_i] is h2_i in component INaCa (dimensionless).
 * ALGEBRAIC[h3_i] is h3_i in component INaCa (dimensionless).
 * ALGEBRAIC[h4_i] is h4_i in component INaCa (dimensionless).
 * ALGEBRAIC[h5_i] is h5_i in component INaCa (dimensionless).
 * ALGEBRAIC[h6_i] is h6_i in component INaCa (dimensionless).
 * ALGEBRAIC[h7_i] is h7_i in component INaCa (dimensionless).
 * ALGEBRAIC[h8_i] is h8_i in component INaCa (dimensionless).
 * ALGEBRAIC[h9_i] is h9_i in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + h10_i] is h10_i in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + h11_i] is h11_i in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + h12_i] is h12_i in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_i] is k1_i in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_i] is k2_i in component INaCa (dimensionless).
 * ALGEBRAIC[k3p_i] is k3p_i in component INaCa (dimensionless).
 * ALGEBRAIC[k3pp_i] is k3pp_i in component INaCa (dimensionless).
 * ALGEBRAIC[k3_i] is k3_i in component INaCa (dimensionless).
 * ALGEBRAIC[k4_i] is k4_i in component INaCa (dimensionless).
 * ALGEBRAIC[k4p_i] is k4p_i in component INaCa (dimensionless).
 * ALGEBRAIC[k4pp_i] is k4pp_i in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_i] is k5_i in component INaCa (dimensionless).
 * ALGEBRAIC[k6_i] is k6_i in component INaCa (dimensionless).
 * ALGEBRAIC[k7_i] is k7_i in component INaCa (dimensionless).
 * ALGEBRAIC[k8_i] is k8_i in component INaCa (dimensionless).
 * ALGEBRAIC[x1_i] is x1_i in component INaCa (dimensionless).
 * ALGEBRAIC[x2_i] is x2_i in component INaCa (dimensionless).
 * ALGEBRAIC[x3_i] is x3_i in component INaCa (dimensionless).
 * ALGEBRAIC[x4_i] is x4_i in component INaCa (dimensionless).
 * ALGEBRAIC[E1_i] is E1_i in component INaCa (dimensionless).
 * ALGEBRAIC[E2_i] is E2_i in component INaCa (dimensionless).
 * ALGEBRAIC[E3_i] is E3_i in component INaCa (dimensionless).
 * ALGEBRAIC[E4_i] is E4_i in component INaCa (dimensionless).
 * ALGEBRAIC[allo_i] is allo_i in component INaCa (dimensionless).
 * ALGEBRAIC[JncxNa_i] is JncxNa_i in component INaCa (millimolar_per_millisecond).
 * ALGEBRAIC[JncxCa_i] is JncxCa_i in component INaCa (millimolar_per_millisecond).
 * ALGEBRAIC[h1_ss] is h1_ss in component INaCa (dimensionless).
 * ALGEBRAIC[h2_ss] is h2_ss in component INaCa (dimensionless).
 * ALGEBRAIC[h3_ss] is h3_ss in component INaCa (dimensionless).
 * ALGEBRAIC[h4_ss] is h4_ss in component INaCa (dimensionless).
 * ALGEBRAIC[h5_ss] is h5_ss in component INaCa (dimensionless).
 * ALGEBRAIC[h6_ss] is h6_ss in component INaCa (dimensionless).
 * ALGEBRAIC[h7_ss] is h7_ss in component INaCa (dimensionless).
 * ALGEBRAIC[h8_ss] is h8_ss in component INaCa (dimensionless).
 * ALGEBRAIC[h9_ss] is h9_ss in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + h10_ss] is h10_ss in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + h11_ss] is h11_ss in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + h12_ss] is h12_ss in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_ss] is k1_ss in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_ss] is k2_ss in component INaCa (dimensionless).
 * ALGEBRAIC[k3p_ss] is k3p_ss in component INaCa (dimensionless).
 * ALGEBRAIC[k3pp_ss] is k3pp_ss in component INaCa (dimensionless).
 * ALGEBRAIC[k3_ss] is k3_ss in component INaCa (dimensionless).
 * ALGEBRAIC[k4_ss] is k4_ss in component INaCa (dimensionless).
 * ALGEBRAIC[k4p_ss] is k4p_ss in component INaCa (dimensionless).
 * ALGEBRAIC[k4pp_ss] is k4pp_ss in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_ss] is k5_ss in component INaCa (dimensionless).
 * ALGEBRAIC[k6_ss] is k6_ss in component INaCa (dimensionless).
 * ALGEBRAIC[k7_ss] is k7_ss in component INaCa (dimensionless).
 * ALGEBRAIC[k8_ss] is k8_ss in component INaCa (dimensionless).
 * ALGEBRAIC[x1_ss] is x1_ss in component INaCa (dimensionless).
 * ALGEBRAIC[x2_ss] is x2_ss in component INaCa (dimensionless).
 * ALGEBRAIC[x3_ss] is x3_ss in component INaCa (dimensionless).
 * ALGEBRAIC[x4_ss] is x4_ss in component INaCa (dimensionless).
 * ALGEBRAIC[E1_ss] is E1_ss in component INaCa (dimensionless).
 * ALGEBRAIC[E2_ss] is E2_ss in component INaCa (dimensionless).
 * ALGEBRAIC[E3_ss] is E3_ss in component INaCa (dimensionless).
 * ALGEBRAIC[E4_ss] is E4_ss in component INaCa (dimensionless).
 * ALGEBRAIC[allo_ss] is allo_ss in component INaCa (dimensionless).
 * ALGEBRAIC[JncxNa_ss] is JncxNa_ss in component INaCa (millimolar_per_millisecond).
 * ALGEBRAIC[JncxCa_ss] is JncxCa_ss in component INaCa (millimolar_per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k1p] is k1p in component INaK (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k1m] is k1m in component INaK (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2p] is k2p in component INaK (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2m] is k2m in component INaK (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k3p] is k3p in component INaK (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k3m] is k3m in component INaK (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k4p] is k4p in component INaK (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k4m] is k4m in component INaK (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Knai0] is Knai0 in component INaK (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Knao0] is Knao0 in component INaK (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + delta] is delta in component INaK (millivolt).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Kki] is Kki in component INaK (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Kko] is Kko in component INaK (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + MgADP] is MgADP in component INaK (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + MgATP] is MgATP in component INaK (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Kmgatp] is Kmgatp in component INaK (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + H] is H in component INaK (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + eP] is eP in component INaK (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Khp] is Khp in component INaK (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Knap] is Knap in component INaK (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Kxkur] is Kxkur in component INaK (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Pnak_b] is Pnak_b in component INaK (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Pnak] is Pnak in component INaK (milliS_per_microF).
 * ALGEBRAIC[Knai] is Knai in component INaK (millimolar).
 * ALGEBRAIC[Knao] is Knao in component INaK (millimolar).
 * ALGEBRAIC[P] is P in component INaK (dimensionless).
 * ALGEBRAIC[a1] is a1 in component INaK (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + b1] is b1 in component INaK (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + a2] is a2 in component INaK (dimensionless).
 * ALGEBRAIC[b2] is b2 in component INaK (dimensionless).
 * ALGEBRAIC[a3] is a3 in component INaK (dimensionless).
 * ALGEBRAIC[b3] is b3 in component INaK (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + a4] is a4 in component INaK (dimensionless).
 * ALGEBRAIC[b4] is b4 in component INaK (dimensionless).
 * ALGEBRAIC[x1] is x1 in component INaK (dimensionless).
 * ALGEBRAIC[x2] is x2 in component INaK (dimensionless).
 * ALGEBRAIC[x3] is x3 in component INaK (dimensionless).
 * ALGEBRAIC[x4] is x4 in component INaK (dimensionless).
 * ALGEBRAIC[E1] is E1 in component INaK (dimensionless).
 * ALGEBRAIC[E2] is E2 in component INaK (dimensionless).
 * ALGEBRAIC[E3] is E3 in component INaK (dimensionless).
 * ALGEBRAIC[E4] is E4 in component INaK (dimensionless).
 * ALGEBRAIC[JnakNa] is JnakNa in component INaK (millimolar_per_millisecond).
 * ALGEBRAIC[JnakK] is JnakK in component INaK (millimolar_per_millisecond).
 * ALGEBRAIC[xkb] is xkb in component IKb (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GKb_b] is GKb_b in component IKb (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GKb] is GKb in component IKb (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PNab] is PNab in component INab (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCab] is PCab in component ICab (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GpCa] is GpCa in component IpCa (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCap] is KmCap in component IpCa (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GClCa] is GClCa in component ICl (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GClb] is GClb in component ICl (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + KdClCa] is KdClCa in component ICl (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Fjunc] is Fjunc in component ICl (dimensionless).
 * ALGEBRAIC[IClCa_junc] is IClCa_junc in component ICl (microA_per_microF).
 * ALGEBRAIC[IClCa_sl] is IClCa_sl in component ICl (microA_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + tauNa] is tauNa in component diff (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + tauK] is tauK in component diff (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + tauCa] is tauCa in component diff (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + bt] is bt in component ryr (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + a_rel] is a_rel in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[Jrel_inf_b] is Jrel_inf_b in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[Jrel_inf] is Jrel_inf in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[tau_rel_b] is tau_rel_b in component ryr (millisecond).
 * ALGEBRAIC[tau_rel] is tau_rel in component ryr (millisecond).
 * STATES[Jrel_np] is Jrel_np in component ryr (millimolar_per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + btp] is btp in component ryr (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + a_relp] is a_relp in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[Jrel_infp_b] is Jrel_infp_b in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[Jrel_infp] is Jrel_infp in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[tau_relp_b] is tau_relp_b in component ryr (millisecond).
 * ALGEBRAIC[tau_relp] is tau_relp in component ryr (millisecond).
 * STATES[Jrel_p] is Jrel_p in component ryr (millimolar_per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + cajsr_half] is cajsr_half in component ryr (millimolar).
 * ALGEBRAIC[fJrelp] is fJrelp in component ryr (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Jrel_b] is Jrel_b in component ryr (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + upScale] is upScale in component SERCA (dimensionless).
 * ALGEBRAIC[Jupnp] is Jupnp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[Jupp] is Jupp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[fJupp] is fJupp in component SERCA (dimensionless).
 * ALGEBRAIC[Jleak] is Jleak in component SERCA (millimolar_per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Jup_b] is Jup_b in component SERCA (dimensionless).
 * RATES[V] is d/dt v in component membrane (millivolt).
 * RATES[CaMKt] is d/dt CaMKt in component CaMK (millimolar).
 * RATES[nai] is d/dt nai in component intracellular_ions (millimolar).
 * RATES[nass] is d/dt nass in component intracellular_ions (millimolar).
 * RATES[ki] is d/dt ki in component intracellular_ions (millimolar).
 * RATES[kss] is d/dt kss in component intracellular_ions (millimolar).
 * RATES[cai] is d/dt cai in component intracellular_ions (millimolar).
 * RATES[cass] is d/dt cass in component intracellular_ions (millimolar).
 * RATES[cansr] is d/dt cansr in component intracellular_ions (millimolar).
 * RATES[cajsr] is d/dt cajsr in component intracellular_ions (millimolar).
 * RATES[m] is d/dt m in component INa (dimensionless).
 * RATES[h] is d/dt h in component INa (dimensionless).
 * RATES[j] is d/dt j in component INa (dimensionless).
 * RATES[hp] is d/dt hp in component INa (dimensionless).
 * RATES[jp] is d/dt jp in component INa (dimensionless).
 * RATES[mL] is d/dt mL in component INaL (dimensionless).
 * RATES[hL] is d/dt hL in component INaL (dimensionless).
 * RATES[hLp] is d/dt hLp in component INaL (dimensionless).
 * RATES[a] is d/dt a in component Ito (dimensionless).
 * RATES[iF] is d/dt iF in component Ito (dimensionless).
 * RATES[iS] is d/dt iS in component Ito (dimensionless).
 * RATES[ap] is d/dt ap in component Ito (dimensionless).
 * RATES[iFp] is d/dt iFp in component Ito (dimensionless).
 * RATES[iSp] is d/dt iSp in component Ito (dimensionless).
 * RATES[d] is d/dt d in component ICaL (dimensionless).
 * RATES[ff] is d/dt ff in component ICaL (dimensionless).
 * RATES[fs] is d/dt fs in component ICaL (dimensionless).
 * RATES[fcaf] is d/dt fcaf in component ICaL (dimensionless).
 * RATES[fcas] is d/dt fcas in component ICaL (dimensionless).
 * RATES[jca] is d/dt jca in component ICaL (dimensionless).
 * RATES[ffp] is d/dt ffp in component ICaL (dimensionless).
 * RATES[fcafp] is d/dt fcafp in component ICaL (dimensionless).
 * RATES[nca_ss] is d/dt nca_ss in component ICaL (dimensionless).
 * RATES[nca_i] is d/dt nca_i in component ICaL (dimensionless).
 * RATES[C3] is d/dt C3 in component IKr (dimensionless).
 * RATES[C2] is d/dt C2 in component IKr (dimensionless).
 * RATES[C1] is d/dt C1 in component IKr (dimensionless).
 * RATES[O] is d/dt O in component IKr (dimensionless).
 * RATES[I] is d/dt I in component IKr (dimensionless).
 * RATES[xs1] is d/dt xs1 in component IKs (dimensionless).
 * RATES[xs2] is d/dt xs2 in component IKs (dimensionless).
 * RATES[Jrel_np] is d/dt Jrel_np in component ryr (millimolar_per_millisecond).
 * RATES[Jrel_p] is d/dt Jrel_p in component ryr (millimolar_per_millisecond).
 */


// Tomek_model::Tomek_model()
// {

// }

// Tomek_model::~Tomek_model()
// {

}

__device__ void ___initConsts(double *CONSTANTS, double *STATES, double type, double bcl, int sample_id)
    {
    CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype] = type;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + nao] = 140.0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + cao] = 1.8;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + ko] = 5.0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + clo] = 150.0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + R] = 8314;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + T] = 310;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + F] = 96485;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + zna] = 1;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + zca] = 2;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + zk] = 1;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + zcl] = -1;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + L] = 0.01;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + rad] = 0.0011;
    STATES[V] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ?  -89.14 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ?  -89.1704 : -88.7638);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + stim_start] = 10;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + i_Stim_End] = 100000000000000000;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + i_Stim_Amplitude] = -53;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + BCL] = bcl;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + i_Stim_PulseDuration] = 1.0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaMK] = 0.15;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + aCaMK] = 0.05;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + bCaMK] = 0.00068;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + CaMKo] = 0.05;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaM] = 0.0015;
    STATES[CaMKt] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.0129 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.0192 : 0.0111);
    STATES[cass] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 5.77E-05 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 6.58E-05 : 7.0305e-5);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + cmdnmax_b] = 0.05;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcmdn] = 0.00238;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + trpnmax] = 0.07;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + kmtrpn] = 0.0005;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + BSRmax] = 0.047;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + KmBSR] = 0.00087;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + BSLmax] = 1.124;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + KmBSL] = 0.0087;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + csqnmax] = 10;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcsqn] = 0.8;
    STATES[nai] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 12.1025 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 15.0038 : 12.1025);
    STATES[nass] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 12.8366 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 15.0043 : 12.1029);
    STATES[ki] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 142.6951 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 143.0403 : 142.3002);
    STATES[kss] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 142.6951 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 143.0402 : 142.3002);
    STATES[cansr] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 1.8119 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 1.9557 : 1.5211);
    STATES[cajsr] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 1.8102 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 1.9593 : 1.5214);
    STATES[cai] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 6.63E-05 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 8.17E-05 : 8.1583e-05);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + cli] = 24.0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + PKNa] = 0.01833;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + gkatp] = 4.3195;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + fkatp] = 0.0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + K_o_n] = 5;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + A_atp] = 2;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + K_atp] = 0.25;
    STATES[m] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 7.43E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 7.38E-04 : 8.0572e-4);
    STATES[h] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.836 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.8365 : 0.8286);
    STATES[j] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.8359 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.8363 : 0.8284);
    STATES[hp] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.6828 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.6838 : 0.6707);
    STATES[jp] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.8357 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.8358 : 0.8281);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GNa] = 11.7802;
    STATES[mL] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 1.52E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 1.51E-04 : 1.629e-4);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + thL] = 200;
    STATES[hL] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.5401 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.5327 : 0.5255);
    STATES[hLp] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.3034 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.2834 : 0.2872);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL_b] = 0.0279;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Gto_b] = 0.16;
    STATES[a] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 9.27E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 9.25E-04 : 9.5098e-4);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift] = 0;
    STATES[iF] = 0.9996;
    STATES[iS] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.9996 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.5671 : 0.5936);
    STATES[ap] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 4.72E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 4.71E-04 : 4.8454e-4);
    STATES[iFp] = 0.9996;
    STATES[iSp] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.9996 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.6261 :0.6538);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Kmn] = 0.002;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] = 500;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa_b] = 8.3757e-05;
    STATES[d] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.0 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.0 : 8.1084e-9);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Aff] = 0.6;
    STATES[ff] = 1.0;
    STATES[fs] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.9485 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.92 : 0.939);
    STATES[fcaf] = 1.0;
    STATES[fcas] = 0.9999;
    STATES[jca] = 1.0;
    STATES[ffp] = 1.0;
    STATES[fcafp] = 1.0;
    STATES[nca_ss] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 3.09E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 5.14E-04 : 6.6462e-4);
    STATES[nca_i] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 5.30E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.0012 : 0.0012);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + tjca] = 75;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + vShift] = 0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + sample_id] = 0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + dielConstant] = 74;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS] = 0.8;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr_b] = 0.0321;
    STATES[C1] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 6.79E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 6.96E-04 : 7.0344e-4);
    STATES[C2] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 8.29E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 8.27E-04 : 8.5109e-4);
    STATES[C3] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.9982 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.9979 : 0.9981);
    STATES[I] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 9.54E-06 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 1.88E-05 : 1.3289e-5);
    STATES[O] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 2.76E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 5.42E-04 : 3.7585e-4);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + alpha_1] = 0.154375;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + beta_1] = 0.1911;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs_b] = 0.0011;
    STATES[xs1] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.2309 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.2653 : 0.248);
    STATES[xs2] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 1.70E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 1.69E-04 : 1.7707e-4);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1_b] = 0.6992;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + INaCa_fractionSS] = 0.35;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + kna1] = 15;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + kna2] = 5;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3] = 88.12;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + kasymm] = 12.5;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + wna] = 6e4;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + wca] = 6e4;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + wnaca] = 5e3;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + kcaon] = 1.5e6;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + kcaoff] = 5e3;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + qna] = 0.5224;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + qca] = 0.167;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaAct] = 150e-6;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Gncx_b] = 0.0034;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k1p] = 949.5;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k1m] = 182.4;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k2p] = 687.2;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k2m] = 39.4;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k3p] = 1899;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k3m] = 79300;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k4p] = 639;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k4m] = 40;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Knai0] = 9.073;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Knao0] = 27.78;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + delta] = -0.155;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Kki] = 0.5;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Kko] = 0.3582;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + MgADP] = 0.05;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + MgATP] = 9.8;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Kmgatp] = 1.698e-7;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + H] = 1e-7;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + eP] = 4.2;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Khp] = 1.698e-7;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Knap] = 224;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Kxkur] = 292;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Pnak_b] = 15.4509;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GKb_b] = 0.0189;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + PNab] = 1.9239e-09;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + PCab] = 5.9194e-08;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GpCa] = 5e-04;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCap] = 0.0005;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GClCa] = 0.2843;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GClb] = 1.98e-3;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + KdClCa] = 0.1;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Fjunc] = 1;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + tauNa] = 2.0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + tauK] = 2.0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + tauCa] = 0.2;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + bt] = 4.75;
    STATES[Jrel_np] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 2.82E-24 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0. : 1.6129e-22);
    STATES[Jrel_p] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0. : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0. : 1.2475e-20);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + cajsr_half] = 1.7;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Jrel_b] = 1.5378;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Jup_b] = 1.0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + vcell] =  1000.00*3.14000*CONSTANTS[(sample_id * Tomek_num_of_constants) + rad]*CONSTANTS[(sample_id * Tomek_num_of_constants) + rad]*CONSTANTS[(sample_id * Tomek_num_of_constants) + L];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + cmdnmax] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + cmdnmax_b]*1.30000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + cmdnmax_b]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + ECl] =  (( CONSTANTS[(sample_id * Tomek_num_of_constants) + R]*CONSTANTS[(sample_id * Tomek_num_of_constants) + T])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + zcl]*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]))*log(CONSTANTS[(sample_id * Tomek_num_of_constants) + clo]/CONSTANTS[(sample_id * Tomek_num_of_constants) + cli]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + akik] = pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]/CONSTANTS[(sample_id * Tomek_num_of_constants) + K_o_n], 0.240000);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + bkik] = 1.00000/(1.00000+pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + A_atp]/CONSTANTS[(sample_id * Tomek_num_of_constants) + K_atp], 2.00000));
    CONSTANTS[(sample_id * Tomek_num_of_constants) + thLp] =  3.00000*CONSTANTS[(sample_id * Tomek_num_of_constants) + thL];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL_b]*0.600000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL_b]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Gto] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + Gto_b]*2.00000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + Gto_b]*2.00000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + Gto_b]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Afs] = 1.00000 - CONSTANTS[(sample_id * Tomek_num_of_constants) + Aff];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa_b]*1.20000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa_b]*2.00000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa_b]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Io] = ( 0.500000*(CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]+CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]+CONSTANTS[(sample_id * Tomek_num_of_constants) + clo]+ 4.00000*CONSTANTS[(sample_id * Tomek_num_of_constants) + cao]))/1000.00;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr_b]*1.30000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr_b]*0.800000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr_b]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs_b]*1.40000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs_b]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1_b]*1.20000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1_b]*1.30000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1_b]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GKb] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + GKb_b]*0.600000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + GKb_b]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + a_rel] = ( 0.500000*CONSTANTS[(sample_id * Tomek_num_of_constants) + bt])/1.00000;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + btp] =  1.25000*CONSTANTS[(sample_id * Tomek_num_of_constants) + bt];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + upScale] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 1.30000 : 1.00000);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Ageo] =  2.00000*3.14000*CONSTANTS[(sample_id * Tomek_num_of_constants) + rad]*CONSTANTS[(sample_id * Tomek_num_of_constants) + rad]+ 2.00000*3.14000*CONSTANTS[(sample_id * Tomek_num_of_constants) + rad]*CONSTANTS[(sample_id * Tomek_num_of_constants) + L];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + PCap] =  1.10000*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaNa] =  0.00125000*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaK] =  0.000357400*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + constA] =  1.82000e+06*pow( CONSTANTS[(sample_id * Tomek_num_of_constants) + dielConstant]*CONSTANTS[(sample_id * Tomek_num_of_constants) + T], - 1.50000);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + a_relp] = ( 0.500000*CONSTANTS[(sample_id * Tomek_num_of_constants) + btp])/1.00000;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap] =  2.00000*CONSTANTS[(sample_id * Tomek_num_of_constants) + Ageo];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaNap] =  0.00125000*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCap];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaKp] =  0.000357400*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCap];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_cao] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*4.00000*( pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + Io], 1.0 / 2)/(1.00000+ pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + Io], 1.0 / 2)) -  0.300000*CONSTANTS[(sample_id * Tomek_num_of_constants) + Io]));
    CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_nao] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*1.00000*( pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + Io], 1.0 / 2)/(1.00000+ pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + Io], 1.0 / 2)) -  0.300000*CONSTANTS[(sample_id * Tomek_num_of_constants) + Io]));
    CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_ko] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*1.00000*( pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + Io], 1.0 / 2)/(1.00000+ pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + Io], 1.0 / 2)) -  0.300000*CONSTANTS[(sample_id * Tomek_num_of_constants) + Io]));
    CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo] =  0.680000*CONSTANTS[(sample_id * Tomek_num_of_constants) + vcell];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + vnsr] =  0.0552000*CONSTANTS[(sample_id * Tomek_num_of_constants) + vcell];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + vjsr] =  0.00480000*CONSTANTS[(sample_id * Tomek_num_of_constants) + vcell];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + vss] =  0.0200000*CONSTANTS[(sample_id * Tomek_num_of_constants) + vcell];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + h10_i] = CONSTANTS[(sample_id * Tomek_num_of_constants) + kasymm]+1.00000+ (CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna1])*(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna2]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + h11_i] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]*CONSTANTS[(sample_id * Tomek_num_of_constants) + nao])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + h10_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kna1]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kna2]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + h12_i] = 1.00000/CONSTANTS[(sample_id * Tomek_num_of_constants) + h10_i];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_i] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + h12_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + cao]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kcaon];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_i] = CONSTANTS[(sample_id * Tomek_num_of_constants) + kcaoff];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_i] = CONSTANTS[(sample_id * Tomek_num_of_constants) + kcaoff];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Gncx] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + Gncx_b]*1.10000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + Gncx_b]*1.40000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + Gncx_b]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + h10_ss] = CONSTANTS[(sample_id * Tomek_num_of_constants) + kasymm]+1.00000+ (CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna1])*(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna2]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + h11_ss] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]*CONSTANTS[(sample_id * Tomek_num_of_constants) + nao])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + h10_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kna1]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kna2]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + h12_ss] = 1.00000/CONSTANTS[(sample_id * Tomek_num_of_constants) + h10_ss];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + h12_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + cao]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kcaon];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_ss] = CONSTANTS[(sample_id * Tomek_num_of_constants) + kcaoff];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_ss] = CONSTANTS[(sample_id * Tomek_num_of_constants) + kcaoff];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + b1] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k1m]*CONSTANTS[(sample_id * Tomek_num_of_constants) + MgADP];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + a2] = CONSTANTS[(sample_id * Tomek_num_of_constants) + k2p];
    CONSTANTS[(sample_id * Tomek_num_of_constants) + a4] = (( CONSTANTS[(sample_id * Tomek_num_of_constants) + k4p]*CONSTANTS[(sample_id * Tomek_num_of_constants) + MgATP])/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kmgatp])/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + MgATP]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kmgatp]);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Pnak] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + Pnak_b]*0.900000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ?  CONSTANTS[(sample_id * Tomek_num_of_constants) + Pnak_b]*0.700000 : CONSTANTS[(sample_id * Tomek_num_of_constants) + Pnak_b]);
}

__device__ void ___applyDrugEffect(double *CONSTANTS, double conc, double *hill, int sample_id)
    {
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1] = CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1] * ((hill[(sample_id * 14) +2] > 10E-14 && hill[(sample_id * 14) +3] > 10E-14) ? 1./(1.+pow(conc/hill[(sample_id * 14) +2],hill[(sample_id * 14) +3])) : 1.);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr] = CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr] * ((hill[(sample_id * 14) +12] > 10E-14 && hill[(sample_id * 14) +13] > 10E-14) ? 1./(1.+pow(conc/hill[(sample_id * 14) +12],hill[(sample_id * 14) +13])) : 1.);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs] = CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs] * ((hill[(sample_id * 14) +4] > 10E-14 && hill[(sample_id * 14) +5] > 10E-14) ? 1./(1.+pow(conc/hill[(sample_id * 14) +4],hill[(sample_id * 14) +5])) : 1.);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL] = CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL] * ((hill[(sample_id * 14) +8] > 10E-14 && hill[(sample_id * 14) +9] > 10E-14) ? 1./(1.+pow(conc/hill[(sample_id * 14) +8],hill[(sample_id * 14) +9])) : 1.);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GNa] = CONSTANTS[(sample_id * Tomek_num_of_constants) + GNa] * ((hill[(sample_id * 14) +6] > 10E-14 && hill[(sample_id * 14) +7] > 10E-14) ? 1./(1.+pow(conc/hill[(sample_id * 14) +6],hill[(sample_id * 14) +7])) : 1.);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Gto] = CONSTANTS[(sample_id * Tomek_num_of_constants) + Gto] * ((hill[(sample_id * 14) +10] > 10E-14 && hill[(sample_id * 14) +11] > 10E-14) ? 1./(1.+pow(conc/hill[(sample_id * 14) +10],hill[(sample_id * 14) +11])) : 1.);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa] = CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa] * ( (hill[(sample_id * 14) +0] > 10E-14 && hill[(sample_id * 14) +1] > 10E-14) ? 1./(1.+pow(conc/hill[(sample_id * 14) +0],hill[(sample_id * 14) +1])) : 1.);
    }

// __device__ void initConsts()
// {
// 	___initConsts(0.);
// }

// __device__ void initConsts(double type)
// {
// 	___initConsts(type);
// }

__device__ void initConsts(double *CONSTANTS, double *STATES, double type, double conc, double *ic50, double *cvar,  bool is_cvar, double bcl, double epsilon, int sample_id)
    {
	___initConsts(type);
	printf("Celltype: %lf\n", CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]);
	#ifndef COMPONENT_PATCH
	printf("Control %lf %lf %lf %lf %lf\n", CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa], CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1], CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs], CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL], CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr]);
	#endif
	___applyDrugEffect(conc, hill);
	#ifndef COMPONENT_PATCH
	printf("After drug %lf %lf %lf %lf %lf\n", CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa], CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1], CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs], CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL], CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr]);
	#endif
    }

__device__ void computeRates(double TIME, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int sample_id,  double land_trpn)
    {
    ALGEBRAIC[hLss] = 1.00000/(1.00000+exp((STATES[V]+87.6100)/7.48800));
    ALGEBRAIC[hLssp] = 1.00000/(1.00000+exp((STATES[V]+93.8100)/7.48800));
    ALGEBRAIC[jcass] = 1.00000/(1.00000+exp((STATES[V]+18.0800)/2.79160));
    ALGEBRAIC[mss] = 1.00000/pow(1.00000+exp(- (STATES[V]+56.8600)/9.03000), 2.00000);
    ALGEBRAIC[tm] =  0.129200*exp(- pow((STATES[V]+45.7900)/15.5400, 2.00000))+ 0.0648700*exp(- pow((STATES[V] - 4.82300)/51.1200, 2.00000));
    ALGEBRAIC[mLss] = 1.00000/(1.00000+exp(- (STATES[V]+42.8500)/5.26400));
    ALGEBRAIC[tmL] =  0.129200*exp(- pow((STATES[V]+45.7900)/15.5400, 2.00000))+ 0.0648700*exp(- pow((STATES[V] - 4.82300)/51.1200, 2.00000));
    ALGEBRAIC[ass] = 1.00000/(1.00000+exp(- ((STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]) - 14.3400)/14.8200));
    ALGEBRAIC[ta] = 1.05150/(1.00000/( 1.20890*(1.00000+exp(- ((STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]) - 18.4099)/29.3814)))+3.50000/(1.00000+exp((STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+100.000)/29.3814)));
    ALGEBRAIC[dss] = (STATES[V]>=31.4978 ? 1.00000 :  1.07630*exp( - 1.00700*exp( - 0.0829000*STATES[V])));
    ALGEBRAIC[td] = CONSTANTS[(sample_id * Tomek_num_of_constants) + sample_id]+0.600000+1.00000/(exp( - 0.0500000*(STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + vShift]+6.00000))+exp( 0.0900000*(STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + vShift]+14.0000)));
    ALGEBRAIC[fss] = 1.00000/(1.00000+exp((STATES[V]+19.5800)/3.69600));
    ALGEBRAIC[tff] = 7.00000+1.00000/( 0.00450000*exp(- (STATES[V]+20.0000)/10.0000)+ 0.00450000*exp((STATES[V]+20.0000)/10.0000));
    ALGEBRAIC[tfs] = 1000.00+1.00000/( 3.50000e-05*exp(- (STATES[V]+5.00000)/4.00000)+ 3.50000e-05*exp((STATES[V]+5.00000)/6.00000));
    ALGEBRAIC[km2n] =  STATES[jca]*1.00000;
    ALGEBRAIC[anca_ss] = 1.00000/(CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n]/ALGEBRAIC[km2n]+pow(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + Kmn]/STATES[cass], 4.00000));
    ALGEBRAIC[anca_i] = 1.00000/(CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n]/ALGEBRAIC[km2n]+pow(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + Kmn]/STATES[cai], 4.00000));
    ALGEBRAIC[xs1ss] = 1.00000/(1.00000+exp(- (STATES[V]+11.6000)/8.93200));
    ALGEBRAIC[txs1] = 817.300+1.00000/( 0.000232600*exp((STATES[V]+48.2800)/17.8000)+ 0.00129200*exp(- (STATES[V]+210.000)/230.000));
    ALGEBRAIC[assp] = 1.00000/(1.00000+exp(- ((STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]) - 24.3400)/14.8200));
    ALGEBRAIC[fcass] = ALGEBRAIC[fss];
    ALGEBRAIC[tfcaf] = 7.00000+1.00000/( 0.0400000*exp(- (STATES[V] - 4.00000)/7.00000)+ 0.0400000*exp((STATES[V] - 4.00000)/7.00000));
    ALGEBRAIC[tfcas] = 100.000+1.00000/( 0.000120000*exp(- STATES[V]/3.00000)+ 0.000120000*exp(STATES[V]/7.00000));
    ALGEBRAIC[tffp] =  2.50000*ALGEBRAIC[tff];
    ALGEBRAIC[xs2ss] = ALGEBRAIC[xs1ss];
    ALGEBRAIC[txs2] = 1.00000/( 0.0100000*exp((STATES[V] - 50.0000)/20.0000)+ 0.0193000*exp(- (STATES[V]+66.5400)/31.0000));
    ALGEBRAIC[CaMKb] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + CaMKo]*(1.00000 - STATES[CaMKt]))/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaM]/STATES[cass]);
    ALGEBRAIC[hss] = 1.00000/pow(1.00000+exp((STATES[V]+71.5500)/7.43000), 2.00000);
    ALGEBRAIC[ah] = (STATES[V]>=- 40.0000 ? 0.00000 :  0.0570000*exp(- (STATES[V]+80.0000)/6.80000));
    ALGEBRAIC[bh] = (STATES[V]>=- 40.0000 ? 0.770000/( 0.130000*(1.00000+exp(- (STATES[V]+10.6600)/11.1000))) :  2.70000*exp( 0.0790000*STATES[V])+ 310000.*exp( 0.348500*STATES[V]));
    ALGEBRAIC[th] = 1.00000/(ALGEBRAIC[ah]+ALGEBRAIC[bh]);
    ALGEBRAIC[tfcafp] =  2.50000*ALGEBRAIC[tfcaf];
    ALGEBRAIC[jss] = ALGEBRAIC[hss];
    ALGEBRAIC[aj] = (STATES[V]>=- 40.0000 ? 0.00000 : ( ( - 25428.0*exp( 0.244400*STATES[V]) -  6.94800e-06*exp( - 0.0439100*STATES[V]))*(STATES[V]+37.7800))/(1.00000+exp( 0.311000*(STATES[V]+79.2300))));
    ALGEBRAIC[bj] = (STATES[V]>=- 40.0000 ? ( 0.600000*exp( 0.0570000*STATES[V]))/(1.00000+exp( - 0.100000*(STATES[V]+32.0000))) : ( 0.0242400*exp( - 0.0105200*STATES[V]))/(1.00000+exp( - 0.137800*(STATES[V]+40.1400))));
    ALGEBRAIC[tj] = 1.00000/(ALGEBRAIC[aj]+ALGEBRAIC[bj]);
    ALGEBRAIC[hssp] = 1.00000/pow(1.00000+exp((STATES[V]+77.5500)/7.43000), 2.00000);
    ALGEBRAIC[iss] = 1.00000/(1.00000+exp((STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+43.9400)/5.71100));
    ALGEBRAIC[delta_epi] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 1.00000 - 0.950000/(1.00000+exp((STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+70.0000)/5.00000)) : 1.00000);
    ALGEBRAIC[tiF_b] = 4.56200+1.00000/( 0.393300*exp(- (STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+100.000)/100.000)+ 0.0800400*exp((STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+50.0000)/16.5900));
    ALGEBRAIC[tiF] =  ALGEBRAIC[tiF_b]*ALGEBRAIC[delta_epi];
    ALGEBRAIC[vfrt] = ( STATES[V]*CONSTANTS[(sample_id * Tomek_num_of_constants) + F])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + R]*CONSTANTS[(sample_id * Tomek_num_of_constants) + T]);
    ALGEBRAIC[alpha] =  0.116100*exp( 0.299000*ALGEBRAIC[vfrt]);
    ALGEBRAIC[beta] =  0.244200*exp( - 1.60400*ALGEBRAIC[vfrt]);
    ALGEBRAIC[tjp] =  1.46000*ALGEBRAIC[tj];
    ALGEBRAIC[tiS_b] = 23.6200+1.00000/( 0.00141600*exp(- (STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+96.5200)/59.0500)+ 1.78000e-08*exp((STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+114.100)/8.07900));
    ALGEBRAIC[tiS] =  ALGEBRAIC[tiS_b]*ALGEBRAIC[delta_epi];
    ALGEBRAIC[alpha_2] =  0.0578000*exp( 0.971000*ALGEBRAIC[vfrt]);
    ALGEBRAIC[beta_2] =  0.000349000*exp( - 1.06200*ALGEBRAIC[vfrt]);
    ALGEBRAIC[alpha_i] =  0.253300*exp( 0.595300*ALGEBRAIC[vfrt]);
    ALGEBRAIC[beta_i] =  0.0652500*exp( - 0.820900*ALGEBRAIC[vfrt]);
    ALGEBRAIC[dti_develop] = 1.35400+0.000100000/(exp(((STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]) - 167.400)/15.8900)+exp(- ((STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]) - 12.2300)/0.215400));
    ALGEBRAIC[dti_recover] = 1.00000 - 0.500000/(1.00000+exp((STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+70.0000)/20.0000));
    ALGEBRAIC[tiFp] =  ALGEBRAIC[dti_develop]*ALGEBRAIC[dti_recover]*ALGEBRAIC[tiF];
    ALGEBRAIC[tiSp] =  ALGEBRAIC[dti_develop]*ALGEBRAIC[dti_recover]*ALGEBRAIC[tiS];
    ALGEBRAIC[alpha_C2ToI] =  5.20000e-05*exp( 1.52500*ALGEBRAIC[vfrt]);
    ALGEBRAIC[beta_ItoC2] = ( ALGEBRAIC[beta_2]*ALGEBRAIC[beta_i]*ALGEBRAIC[alpha_C2ToI])/( ALGEBRAIC[alpha_2]*ALGEBRAIC[alpha_i]);
    ALGEBRAIC[f] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Aff]*STATES[ff]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + Afs]*STATES[fs];
    ALGEBRAIC[Afcaf] = 0.300000+0.600000/(1.00000+exp((STATES[V] - 10.0000)/10.0000));
    ALGEBRAIC[Afcas] = 1.00000 - ALGEBRAIC[Afcaf];
    ALGEBRAIC[fca] =  ALGEBRAIC[Afcaf]*STATES[fcaf]+ ALGEBRAIC[Afcas]*STATES[fcas];
    ALGEBRAIC[fp] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Aff]*STATES[ffp]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + Afs]*STATES[fs];
    ALGEBRAIC[fcap] =  ALGEBRAIC[Afcaf]*STATES[fcafp]+ ALGEBRAIC[Afcas]*STATES[fcas];
    ALGEBRAIC[vffrt] = ( STATES[V]*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + F])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + R]*CONSTANTS[(sample_id * Tomek_num_of_constants) + T]);
    ALGEBRAIC[Iss] = ( 0.500000*(STATES[nass]+STATES[kss]+CONSTANTS[(sample_id * Tomek_num_of_constants) + cli]+ 4.00000*STATES[cass]))/1000.00;
    ALGEBRAIC[gamma_cass] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*4.00000*( pow(ALGEBRAIC[Iss], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[Iss], 1.0 / 2)) -  0.300000*ALGEBRAIC[Iss]));
    ALGEBRAIC[PhiCaL_ss] = ( 4.00000*ALGEBRAIC[vffrt]*( ALGEBRAIC[gamma_cass]*STATES[cass]*exp( 2.00000*ALGEBRAIC[vfrt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_cao]*CONSTANTS[(sample_id * Tomek_num_of_constants) + cao]))/(exp( 2.00000*ALGEBRAIC[vfrt]) - 1.00000);
    ALGEBRAIC[CaMKa] = ALGEBRAIC[CaMKb]+STATES[CaMKt];
    ALGEBRAIC[fICaLp] = 1.00000/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaMK]/ALGEBRAIC[CaMKa]);
    ALGEBRAIC[ICaL_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS]*( (1.00000 - ALGEBRAIC[fICaLp])*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa]*ALGEBRAIC[PhiCaL_ss]*STATES[d]*( ALGEBRAIC[f]*(1.00000 - STATES[nca_ss])+ STATES[jca]*ALGEBRAIC[fca]*STATES[nca_ss])+ ALGEBRAIC[fICaLp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCap]*ALGEBRAIC[PhiCaL_ss]*STATES[d]*( ALGEBRAIC[fp]*(1.00000 - STATES[nca_ss])+ STATES[jca]*ALGEBRAIC[fcap]*STATES[nca_ss]));
    ALGEBRAIC[Jrel_inf_b] = (( - CONSTANTS[(sample_id * Tomek_num_of_constants) + a_rel]*ALGEBRAIC[ICaL_ss])/1.00000)/(1.00000+pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + cajsr_half]/STATES[cajsr], 8.00000));
    ALGEBRAIC[Jrel_inf] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ?  ALGEBRAIC[Jrel_inf_b]*1.70000 : ALGEBRAIC[Jrel_inf_b]);
    ALGEBRAIC[tau_rel_b] = CONSTANTS[(sample_id * Tomek_num_of_constants) + bt]/(1.00000+0.0123000/STATES[cajsr]);
    ALGEBRAIC[tau_rel] = (ALGEBRAIC[tau_rel_b]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_rel_b]);
    ALGEBRAIC[Jrel_infp_b] = (( - CONSTANTS[(sample_id * Tomek_num_of_constants) + a_relp]*ALGEBRAIC[ICaL_ss])/1.00000)/(1.00000+pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + cajsr_half]/STATES[cajsr], 8.00000));
    ALGEBRAIC[Jrel_infp] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ?  ALGEBRAIC[Jrel_infp_b]*1.70000 : ALGEBRAIC[Jrel_infp_b]);
    ALGEBRAIC[tau_relp_b] = CONSTANTS[(sample_id * Tomek_num_of_constants) + btp]/(1.00000+0.0123000/STATES[cajsr]);
    ALGEBRAIC[tau_relp] = (ALGEBRAIC[tau_relp_b]<0.00100000 ? 0.00100000 : ALGEBRAIC[tau_relp_b]);
    ALGEBRAIC[EK] =  (( CONSTANTS[(sample_id * Tomek_num_of_constants) + R]*CONSTANTS[(sample_id * Tomek_num_of_constants) + T])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + zk]*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]))*log(CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]/STATES[ki]);
    ALGEBRAIC[AiF] = 1.00000/(1.00000+exp(((STATES[V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]) - 213.600)/151.200));
    ALGEBRAIC[AiS] = 1.00000 - ALGEBRAIC[AiF];
    ALGEBRAIC[i] =  ALGEBRAIC[AiF]*STATES[iF]+ ALGEBRAIC[AiS]*STATES[iS];
    ALGEBRAIC[ip] =  ALGEBRAIC[AiF]*STATES[iFp]+ ALGEBRAIC[AiS]*STATES[iSp];
    ALGEBRAIC[fItop] = 1.00000/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaMK]/ALGEBRAIC[CaMKa]);
    ALGEBRAIC[Ito] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Gto]*(STATES[V] - ALGEBRAIC[EK])*( (1.00000 - ALGEBRAIC[fItop])*STATES[a]*ALGEBRAIC[i]+ ALGEBRAIC[fItop]*STATES[ap]*ALGEBRAIC[ip]);
    ALGEBRAIC[IKr] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr]* pow((CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]/5.00000), 1.0 / 2)*STATES[O]*(STATES[V] - ALGEBRAIC[EK]);
    ALGEBRAIC[EKs] =  (( CONSTANTS[(sample_id * Tomek_num_of_constants) + R]*CONSTANTS[(sample_id * Tomek_num_of_constants) + T])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + zk]*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]))*log((CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + PKNa]*CONSTANTS[(sample_id * Tomek_num_of_constants) + nao])/(STATES[ki]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + PKNa]*STATES[nai]));
    ALGEBRAIC[KsCa] = 1.00000+0.600000/(1.00000+pow(3.80000e-05/STATES[cai], 1.40000));
    ALGEBRAIC[IKs] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs]*ALGEBRAIC[KsCa]*STATES[xs1]*STATES[xs2]*(STATES[V] - ALGEBRAIC[EKs]);
    ALGEBRAIC[aK1] = 4.09400/(1.00000+exp( 0.121700*((STATES[V] - ALGEBRAIC[EK]) - 49.9340)));
    ALGEBRAIC[bK1] = ( 15.7200*exp( 0.0674000*((STATES[V] - ALGEBRAIC[EK]) - 3.25700))+exp( 0.0618000*((STATES[V] - ALGEBRAIC[EK]) - 594.310)))/(1.00000+exp( - 0.162900*((STATES[V] - ALGEBRAIC[EK])+14.2070)));
    ALGEBRAIC[K1ss] = ALGEBRAIC[aK1]/(ALGEBRAIC[aK1]+ALGEBRAIC[bK1]);
    ALGEBRAIC[IK1] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1]* pow((CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]/5.00000), 1.0 / 2)*ALGEBRAIC[K1ss]*(STATES[V] - ALGEBRAIC[EK]);
    ALGEBRAIC[Knao] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Knao0]*exp(( (1.00000 - CONSTANTS[(sample_id * Tomek_num_of_constants) + delta])*ALGEBRAIC[vfrt])/3.00000);
    ALGEBRAIC[a3] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + k3p]*pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kko], 2.00000))/((pow(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/ALGEBRAIC[Knao], 3.00000)+pow(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kko], 2.00000)) - 1.00000);
    ALGEBRAIC[P] = CONSTANTS[(sample_id * Tomek_num_of_constants) + eP]/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + H]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Khp]+STATES[nai]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Knap]+STATES[ki]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kxkur]);
    ALGEBRAIC[b3] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + k3m]*ALGEBRAIC[P]*CONSTANTS[(sample_id * Tomek_num_of_constants) + H])/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + MgATP]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kmgatp]);
    ALGEBRAIC[Knai] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Knai0]*exp(( CONSTANTS[(sample_id * Tomek_num_of_constants) + delta]*ALGEBRAIC[vfrt])/3.00000);
    ALGEBRAIC[a1] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + k1p]*pow(STATES[nai]/ALGEBRAIC[Knai], 3.00000))/((pow(1.00000+STATES[nai]/ALGEBRAIC[Knai], 3.00000)+pow(1.00000+STATES[ki]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kki], 2.00000)) - 1.00000);
    ALGEBRAIC[b2] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + k2m]*pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/ALGEBRAIC[Knao], 3.00000))/((pow(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/ALGEBRAIC[Knao], 3.00000)+pow(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kko], 2.00000)) - 1.00000);
    ALGEBRAIC[b4] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + k4m]*pow(STATES[ki]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kki], 2.00000))/((pow(1.00000+STATES[nai]/ALGEBRAIC[Knai], 3.00000)+pow(1.00000+STATES[ki]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kki], 2.00000)) - 1.00000);
    ALGEBRAIC[x1] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + a4]*ALGEBRAIC[a1]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a2]+ ALGEBRAIC[b2]*ALGEBRAIC[b4]*ALGEBRAIC[b3]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + a2]*ALGEBRAIC[b4]*ALGEBRAIC[b3]+ ALGEBRAIC[b3]*ALGEBRAIC[a1]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a2];
    ALGEBRAIC[x2] =  ALGEBRAIC[b2]*CONSTANTS[(sample_id * Tomek_num_of_constants) + b1]*ALGEBRAIC[b4]+ ALGEBRAIC[a1]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a2]*ALGEBRAIC[a3]+ ALGEBRAIC[a3]*CONSTANTS[(sample_id * Tomek_num_of_constants) + b1]*ALGEBRAIC[b4]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + a2]*ALGEBRAIC[a3]*ALGEBRAIC[b4];
    ALGEBRAIC[x3] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + a2]*ALGEBRAIC[a3]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a4]+ ALGEBRAIC[b3]*ALGEBRAIC[b2]*CONSTANTS[(sample_id * Tomek_num_of_constants) + b1]+ ALGEBRAIC[b2]*CONSTANTS[(sample_id * Tomek_num_of_constants) + b1]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a4]+ ALGEBRAIC[a3]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a4]*CONSTANTS[(sample_id * Tomek_num_of_constants) + b1];
    ALGEBRAIC[x4] =  ALGEBRAIC[b4]*ALGEBRAIC[b3]*ALGEBRAIC[b2]+ ALGEBRAIC[a3]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a4]*ALGEBRAIC[a1]+ ALGEBRAIC[b2]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a4]*ALGEBRAIC[a1]+ ALGEBRAIC[b3]*ALGEBRAIC[b2]*ALGEBRAIC[a1];
    ALGEBRAIC[E1] = ALGEBRAIC[x1]/(ALGEBRAIC[x1]+ALGEBRAIC[x2]+ALGEBRAIC[x3]+ALGEBRAIC[x4]);
    ALGEBRAIC[E2] = ALGEBRAIC[x2]/(ALGEBRAIC[x1]+ALGEBRAIC[x2]+ALGEBRAIC[x3]+ALGEBRAIC[x4]);
    ALGEBRAIC[JnakNa] =  3.00000*( ALGEBRAIC[E1]*ALGEBRAIC[a3] -  ALGEBRAIC[E2]*ALGEBRAIC[b3]);
    ALGEBRAIC[E3] = ALGEBRAIC[x3]/(ALGEBRAIC[x1]+ALGEBRAIC[x2]+ALGEBRAIC[x3]+ALGEBRAIC[x4]);
    ALGEBRAIC[E4] = ALGEBRAIC[x4]/(ALGEBRAIC[x1]+ALGEBRAIC[x2]+ALGEBRAIC[x3]+ALGEBRAIC[x4]);
    ALGEBRAIC[JnakK] =  2.00000*( ALGEBRAIC[E4]*CONSTANTS[(sample_id * Tomek_num_of_constants) + b1] -  ALGEBRAIC[E3]*ALGEBRAIC[a1]);
    ALGEBRAIC[INaK] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Pnak]*( CONSTANTS[(sample_id * Tomek_num_of_constants) + zna]*ALGEBRAIC[JnakNa]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + zk]*ALGEBRAIC[JnakK]);
    ALGEBRAIC[xkb] = 1.00000/(1.00000+exp(- (STATES[V] - 10.8968)/23.9871));
    ALGEBRAIC[IKb] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + GKb]*ALGEBRAIC[xkb]*(STATES[V] - ALGEBRAIC[EK]);
    ALGEBRAIC[I_katp] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + fkatp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + gkatp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + akik]*CONSTANTS[(sample_id * Tomek_num_of_constants) + bkik]*(STATES[V] - ALGEBRAIC[EK]);
    ALGEBRAIC[Istim] = (TIME>=CONSTANTS[(sample_id * Tomek_num_of_constants) + stim_start]&&TIME<=CONSTANTS[(sample_id * Tomek_num_of_constants) + i_Stim_End]&&(TIME - CONSTANTS[(sample_id * Tomek_num_of_constants) + stim_start]) -  floor((TIME - CONSTANTS[(sample_id * Tomek_num_of_constants) + stim_start])/CONSTANTS[(sample_id * Tomek_num_of_constants) + BCL])*CONSTANTS[(sample_id * Tomek_num_of_constants) + BCL]<=CONSTANTS[(sample_id * Tomek_num_of_constants) + i_Stim_PulseDuration] ? CONSTANTS[(sample_id * Tomek_num_of_constants) + i_Stim_Amplitude] : 0.00000);
    ALGEBRAIC[Ii] = ( 0.500000*(STATES[nai]+STATES[ki]+CONSTANTS[(sample_id * Tomek_num_of_constants) + cli]+ 4.00000*STATES[cai]))/1000.00;
    ALGEBRAIC[gamma_ki] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*1.00000*( pow(ALGEBRAIC[Ii], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[Ii], 1.0 / 2)) -  0.300000*ALGEBRAIC[Ii]));
    ALGEBRAIC[PhiCaK_i] = ( 1.00000*ALGEBRAIC[vffrt]*( ALGEBRAIC[gamma_ki]*STATES[ki]*exp( 1.00000*ALGEBRAIC[vfrt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_ko]*CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]))/(exp( 1.00000*ALGEBRAIC[vfrt]) - 1.00000);
    ALGEBRAIC[ICaK_i] =  (1.00000 - CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS])*( (1.00000 - ALGEBRAIC[fICaLp])*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaK]*ALGEBRAIC[PhiCaK_i]*STATES[d]*( ALGEBRAIC[f]*(1.00000 - STATES[nca_i])+ STATES[jca]*ALGEBRAIC[fca]*STATES[nca_i])+ ALGEBRAIC[fICaLp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaKp]*ALGEBRAIC[PhiCaK_i]*STATES[d]*( ALGEBRAIC[fp]*(1.00000 - STATES[nca_i])+ STATES[jca]*ALGEBRAIC[fcap]*STATES[nca_i]));
    ALGEBRAIC[JdiffK] = (STATES[kss] - STATES[ki])/CONSTANTS[(sample_id * Tomek_num_of_constants) + tauK];
    ALGEBRAIC[gamma_kss] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*1.00000*( pow(ALGEBRAIC[Iss], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[Iss], 1.0 / 2)) -  0.300000*ALGEBRAIC[Iss]));
    ALGEBRAIC[PhiCaK_ss] = ( 1.00000*ALGEBRAIC[vffrt]*( ALGEBRAIC[gamma_kss]*STATES[kss]*exp( 1.00000*ALGEBRAIC[vfrt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_ko]*CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]))/(exp( 1.00000*ALGEBRAIC[vfrt]) - 1.00000);
    ALGEBRAIC[ICaK_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS]*( (1.00000 - ALGEBRAIC[fICaLp])*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaK]*ALGEBRAIC[PhiCaK_ss]*STATES[d]*( ALGEBRAIC[f]*(1.00000 - STATES[nca_ss])+ STATES[jca]*ALGEBRAIC[fca]*STATES[nca_ss])+ ALGEBRAIC[fICaLp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaKp]*ALGEBRAIC[PhiCaK_ss]*STATES[d]*( ALGEBRAIC[fp]*(1.00000 - STATES[nca_ss])+ STATES[jca]*ALGEBRAIC[fcap]*STATES[nca_ss]));
    ALGEBRAIC[ENa] =  (( CONSTANTS[(sample_id * Tomek_num_of_constants) + R]*CONSTANTS[(sample_id * Tomek_num_of_constants) + T])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + zna]*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]))*log(CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/STATES[nai]);
    ALGEBRAIC[fINap] = 1.00000/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaMK]/ALGEBRAIC[CaMKa]);
    ALGEBRAIC[INa] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + GNa]*(STATES[V] - ALGEBRAIC[ENa])*pow(STATES[m], 3.00000)*( (1.00000 - ALGEBRAIC[fINap])*STATES[h]*STATES[j]+ ALGEBRAIC[fINap]*STATES[hp]*STATES[jp]);
    ALGEBRAIC[fINaLp] = 1.00000/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaMK]/ALGEBRAIC[CaMKa]);
    ALGEBRAIC[INaL] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL]*(STATES[V] - ALGEBRAIC[ENa])*STATES[mL]*( (1.00000 - ALGEBRAIC[fINaLp])*STATES[hL]+ ALGEBRAIC[fINaLp]*STATES[hLp]);
    ALGEBRAIC[allo_i] = 1.00000/(1.00000+pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaAct]/STATES[cai], 2.00000));
    ALGEBRAIC[hna] = exp( CONSTANTS[(sample_id * Tomek_num_of_constants) + qna]*ALGEBRAIC[vfrt]);
    ALGEBRAIC[h7_i] = 1.00000+ (CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3])*(1.00000+1.00000/ALGEBRAIC[hna]);
    ALGEBRAIC[h8_i] = CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/( CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3]*ALGEBRAIC[hna]*ALGEBRAIC[h7_i]);
    ALGEBRAIC[k3pp_i] =  ALGEBRAIC[h8_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wnaca];
    ALGEBRAIC[h1_i] = 1.00000+ (STATES[nai]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3])*(1.00000+ALGEBRAIC[hna]);
    ALGEBRAIC[h2_i] = ( STATES[nai]*ALGEBRAIC[hna])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3]*ALGEBRAIC[h1_i]);
    ALGEBRAIC[k4pp_i] =  ALGEBRAIC[h2_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wnaca];
    ALGEBRAIC[h4_i] = 1.00000+ (STATES[nai]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna1])*(1.00000+STATES[nai]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna2]);
    ALGEBRAIC[h5_i] = ( STATES[nai]*STATES[nai])/( ALGEBRAIC[h4_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kna1]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kna2]);
    ALGEBRAIC[k7_i] =  ALGEBRAIC[h5_i]*ALGEBRAIC[h2_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wna];
    ALGEBRAIC[k8_i] =  ALGEBRAIC[h8_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + h11_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wna];
    ALGEBRAIC[h9_i] = 1.00000/ALGEBRAIC[h7_i];
    ALGEBRAIC[k3p_i] =  ALGEBRAIC[h9_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wca];
    ALGEBRAIC[k3_i] = ALGEBRAIC[k3p_i]+ALGEBRAIC[k3pp_i];
    ALGEBRAIC[hca] = exp( CONSTANTS[(sample_id * Tomek_num_of_constants) + qca]*ALGEBRAIC[vfrt]);
    ALGEBRAIC[h3_i] = 1.00000/ALGEBRAIC[h1_i];
    ALGEBRAIC[k4p_i] = ( ALGEBRAIC[h3_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wca])/ALGEBRAIC[hca];
    ALGEBRAIC[k4_i] = ALGEBRAIC[k4p_i]+ALGEBRAIC[k4pp_i];
    ALGEBRAIC[h6_i] = 1.00000/ALGEBRAIC[h4_i];
    ALGEBRAIC[k6_i] =  ALGEBRAIC[h6_i]*STATES[cai]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kcaon];
    ALGEBRAIC[x1_i] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_i]*ALGEBRAIC[k4_i]*(ALGEBRAIC[k7_i]+ALGEBRAIC[k6_i])+ CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_i]*ALGEBRAIC[k7_i]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_i]+ALGEBRAIC[k3_i]);
    ALGEBRAIC[x2_i] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_i]*ALGEBRAIC[k7_i]*(ALGEBRAIC[k4_i]+CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_i])+ ALGEBRAIC[k4_i]*ALGEBRAIC[k6_i]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_i]+ALGEBRAIC[k8_i]);
    ALGEBRAIC[x3_i] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_i]*ALGEBRAIC[k3_i]*(ALGEBRAIC[k7_i]+ALGEBRAIC[k6_i])+ ALGEBRAIC[k8_i]*ALGEBRAIC[k6_i]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_i]+ALGEBRAIC[k3_i]);
    ALGEBRAIC[x4_i] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_i]*ALGEBRAIC[k8_i]*(ALGEBRAIC[k4_i]+CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_i])+ ALGEBRAIC[k3_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_i]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_i]+ALGEBRAIC[k8_i]);
    ALGEBRAIC[E1_i] = ALGEBRAIC[x1_i]/(ALGEBRAIC[x1_i]+ALGEBRAIC[x2_i]+ALGEBRAIC[x3_i]+ALGEBRAIC[x4_i]);
    ALGEBRAIC[E2_i] = ALGEBRAIC[x2_i]/(ALGEBRAIC[x1_i]+ALGEBRAIC[x2_i]+ALGEBRAIC[x3_i]+ALGEBRAIC[x4_i]);
    ALGEBRAIC[E3_i] = ALGEBRAIC[x3_i]/(ALGEBRAIC[x1_i]+ALGEBRAIC[x2_i]+ALGEBRAIC[x3_i]+ALGEBRAIC[x4_i]);
    ALGEBRAIC[E4_i] = ALGEBRAIC[x4_i]/(ALGEBRAIC[x1_i]+ALGEBRAIC[x2_i]+ALGEBRAIC[x3_i]+ALGEBRAIC[x4_i]);
    ALGEBRAIC[JncxNa_i] = ( 3.00000*( ALGEBRAIC[E4_i]*ALGEBRAIC[k7_i] -  ALGEBRAIC[E1_i]*ALGEBRAIC[k8_i])+ ALGEBRAIC[E3_i]*ALGEBRAIC[k4pp_i]) -  ALGEBRAIC[E2_i]*ALGEBRAIC[k3pp_i];
    ALGEBRAIC[JncxCa_i] =  ALGEBRAIC[E2_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_i] -  ALGEBRAIC[E1_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_i];
    ALGEBRAIC[INaCa_i] =  (1.00000 - CONSTANTS[(sample_id * Tomek_num_of_constants) + INaCa_fractionSS])*CONSTANTS[(sample_id * Tomek_num_of_constants) + Gncx]*ALGEBRAIC[allo_i]*( CONSTANTS[(sample_id * Tomek_num_of_constants) + zna]*ALGEBRAIC[JncxNa_i]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + zca]*ALGEBRAIC[JncxCa_i]);
    ALGEBRAIC[INab] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + PNab]*ALGEBRAIC[vffrt]*( STATES[nai]*exp(ALGEBRAIC[vfrt]) - CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]))/(exp(ALGEBRAIC[vfrt]) - 1.00000);
    ALGEBRAIC[gamma_nai] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*1.00000*( pow(ALGEBRAIC[Ii], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[Ii], 1.0 / 2)) -  0.300000*ALGEBRAIC[Ii]));
    ALGEBRAIC[PhiCaNa_i] = ( 1.00000*ALGEBRAIC[vffrt]*( ALGEBRAIC[gamma_nai]*STATES[nai]*exp( 1.00000*ALGEBRAIC[vfrt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_nao]*CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]))/(exp( 1.00000*ALGEBRAIC[vfrt]) - 1.00000);
    ALGEBRAIC[ICaNa_i] =  (1.00000 - CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS])*( (1.00000 - ALGEBRAIC[fICaLp])*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaNa]*ALGEBRAIC[PhiCaNa_i]*STATES[d]*( ALGEBRAIC[f]*(1.00000 - STATES[nca_i])+ STATES[jca]*ALGEBRAIC[fca]*STATES[nca_i])+ ALGEBRAIC[fICaLp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaNap]*ALGEBRAIC[PhiCaNa_i]*STATES[d]*( ALGEBRAIC[fp]*(1.00000 - STATES[nca_i])+ STATES[jca]*ALGEBRAIC[fcap]*STATES[nca_i]));
    ALGEBRAIC[JdiffNa] = (STATES[nass] - STATES[nai])/CONSTANTS[(sample_id * Tomek_num_of_constants) + tauNa];
    ALGEBRAIC[allo_ss] = 1.00000/(1.00000+pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaAct]/STATES[cass], 2.00000));
    ALGEBRAIC[h7_ss] = 1.00000+ (CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3])*(1.00000+1.00000/ALGEBRAIC[hna]);
    ALGEBRAIC[h8_ss] = CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/( CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3]*ALGEBRAIC[hna]*ALGEBRAIC[h7_ss]);
    ALGEBRAIC[k3pp_ss] =  ALGEBRAIC[h8_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wnaca];
    ALGEBRAIC[h1_ss] = 1.00000+ (STATES[nass]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3])*(1.00000+ALGEBRAIC[hna]);
    ALGEBRAIC[h2_ss] = ( STATES[nass]*ALGEBRAIC[hna])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3]*ALGEBRAIC[h1_ss]);
    ALGEBRAIC[k4pp_ss] =  ALGEBRAIC[h2_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wnaca];
    ALGEBRAIC[h4_ss] = 1.00000+ (STATES[nass]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna1])*(1.00000+STATES[nass]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna2]);
    ALGEBRAIC[h5_ss] = ( STATES[nass]*STATES[nass])/( ALGEBRAIC[h4_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kna1]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kna2]);
    ALGEBRAIC[k7_ss] =  ALGEBRAIC[h5_ss]*ALGEBRAIC[h2_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wna];
    ALGEBRAIC[k8_ss] =  ALGEBRAIC[h8_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + h11_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wna];
    ALGEBRAIC[h9_ss] = 1.00000/ALGEBRAIC[h7_ss];
    ALGEBRAIC[k3p_ss] =  ALGEBRAIC[h9_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wca];
    ALGEBRAIC[k3_ss] = ALGEBRAIC[k3p_ss]+ALGEBRAIC[k3pp_ss];
    ALGEBRAIC[h3_ss] = 1.00000/ALGEBRAIC[h1_ss];
    ALGEBRAIC[k4p_ss] = ( ALGEBRAIC[h3_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wca])/ALGEBRAIC[hca];
    ALGEBRAIC[k4_ss] = ALGEBRAIC[k4p_ss]+ALGEBRAIC[k4pp_ss];
    ALGEBRAIC[h6_ss] = 1.00000/ALGEBRAIC[h4_ss];
    ALGEBRAIC[k6_ss] =  ALGEBRAIC[h6_ss]*STATES[cass]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kcaon];
    ALGEBRAIC[x1_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_ss]*ALGEBRAIC[k4_ss]*(ALGEBRAIC[k7_ss]+ALGEBRAIC[k6_ss])+ CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_ss]*ALGEBRAIC[k7_ss]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_ss]+ALGEBRAIC[k3_ss]);
    ALGEBRAIC[x2_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_ss]*ALGEBRAIC[k7_ss]*(ALGEBRAIC[k4_ss]+CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_ss])+ ALGEBRAIC[k4_ss]*ALGEBRAIC[k6_ss]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_ss]+ALGEBRAIC[k8_ss]);
    ALGEBRAIC[x3_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_ss]*ALGEBRAIC[k3_ss]*(ALGEBRAIC[k7_ss]+ALGEBRAIC[k6_ss])+ ALGEBRAIC[k8_ss]*ALGEBRAIC[k6_ss]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_ss]+ALGEBRAIC[k3_ss]);
    ALGEBRAIC[x4_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_ss]*ALGEBRAIC[k8_ss]*(ALGEBRAIC[k4_ss]+CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_ss])+ ALGEBRAIC[k3_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_ss]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_ss]+ALGEBRAIC[k8_ss]);
    ALGEBRAIC[E1_ss] = ALGEBRAIC[x1_ss]/(ALGEBRAIC[x1_ss]+ALGEBRAIC[x2_ss]+ALGEBRAIC[x3_ss]+ALGEBRAIC[x4_ss]);
    ALGEBRAIC[E2_ss] = ALGEBRAIC[x2_ss]/(ALGEBRAIC[x1_ss]+ALGEBRAIC[x2_ss]+ALGEBRAIC[x3_ss]+ALGEBRAIC[x4_ss]);
    ALGEBRAIC[E3_ss] = ALGEBRAIC[x3_ss]/(ALGEBRAIC[x1_ss]+ALGEBRAIC[x2_ss]+ALGEBRAIC[x3_ss]+ALGEBRAIC[x4_ss]);
    ALGEBRAIC[E4_ss] = ALGEBRAIC[x4_ss]/(ALGEBRAIC[x1_ss]+ALGEBRAIC[x2_ss]+ALGEBRAIC[x3_ss]+ALGEBRAIC[x4_ss]);
    ALGEBRAIC[JncxNa_ss] = ( 3.00000*( ALGEBRAIC[E4_ss]*ALGEBRAIC[k7_ss] -  ALGEBRAIC[E1_ss]*ALGEBRAIC[k8_ss])+ ALGEBRAIC[E3_ss]*ALGEBRAIC[k4pp_ss]) -  ALGEBRAIC[E2_ss]*ALGEBRAIC[k3pp_ss];
    ALGEBRAIC[JncxCa_ss] =  ALGEBRAIC[E2_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_ss] -  ALGEBRAIC[E1_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_ss];
    ALGEBRAIC[INaCa_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + INaCa_fractionSS]*CONSTANTS[(sample_id * Tomek_num_of_constants) + Gncx]*ALGEBRAIC[allo_ss]*( CONSTANTS[(sample_id * Tomek_num_of_constants) + zna]*ALGEBRAIC[JncxNa_ss]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + zca]*ALGEBRAIC[JncxCa_ss]);
    ALGEBRAIC[gamma_nass] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*1.00000*( pow(ALGEBRAIC[Iss], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[Iss], 1.0 / 2)) -  0.300000*ALGEBRAIC[Iss]));
    ALGEBRAIC[PhiCaNa_ss] = ( 1.00000*ALGEBRAIC[vffrt]*( ALGEBRAIC[gamma_nass]*STATES[nass]*exp( 1.00000*ALGEBRAIC[vfrt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_nao]*CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]))/(exp( 1.00000*ALGEBRAIC[vfrt]) - 1.00000);
    ALGEBRAIC[ICaNa_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS]*( (1.00000 - ALGEBRAIC[fICaLp])*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaNa]*ALGEBRAIC[PhiCaNa_ss]*STATES[d]*( ALGEBRAIC[f]*(1.00000 - STATES[nca_ss])+ STATES[jca]*ALGEBRAIC[fca]*STATES[nca_ss])+ ALGEBRAIC[fICaLp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaNap]*ALGEBRAIC[PhiCaNa_ss]*STATES[d]*( ALGEBRAIC[fp]*(1.00000 - STATES[nca_ss])+ STATES[jca]*ALGEBRAIC[fcap]*STATES[nca_ss]));
    ALGEBRAIC[Jdiff] = (STATES[cass] - STATES[cai])/CONSTANTS[(sample_id * Tomek_num_of_constants) + tauCa];
    ALGEBRAIC[fJrelp] = 1.00000/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaMK]/ALGEBRAIC[CaMKa]);
    ALGEBRAIC[Jrel] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Jrel_b]*( (1.00000 - ALGEBRAIC[fJrelp])*STATES[Jrel_np]+ ALGEBRAIC[fJrelp]*STATES[Jrel_p]);
    ALGEBRAIC[Bcass] = 1.00000/(1.00000+( CONSTANTS[(sample_id * Tomek_num_of_constants) + BSRmax]*CONSTANTS[(sample_id * Tomek_num_of_constants) + KmBSR])/pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + KmBSR]+STATES[cass], 2.00000)+( CONSTANTS[(sample_id * Tomek_num_of_constants) + BSLmax]*CONSTANTS[(sample_id * Tomek_num_of_constants) + KmBSL])/pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + KmBSL]+STATES[cass], 2.00000));
    ALGEBRAIC[gamma_cai] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*4.00000*( pow(ALGEBRAIC[Ii], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[Ii], 1.0 / 2)) -  0.300000*ALGEBRAIC[Ii]));
    ALGEBRAIC[PhiCaL_i] = ( 4.00000*ALGEBRAIC[vffrt]*( ALGEBRAIC[gamma_cai]*STATES[cai]*exp( 2.00000*ALGEBRAIC[vfrt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_cao]*CONSTANTS[(sample_id * Tomek_num_of_constants) + cao]))/(exp( 2.00000*ALGEBRAIC[vfrt]) - 1.00000);
    ALGEBRAIC[ICaL_i] =  (1.00000 - CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS])*( (1.00000 - ALGEBRAIC[fICaLp])*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa]*ALGEBRAIC[PhiCaL_i]*STATES[d]*( ALGEBRAIC[f]*(1.00000 - STATES[nca_i])+ STATES[jca]*ALGEBRAIC[fca]*STATES[nca_i])+ ALGEBRAIC[fICaLp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCap]*ALGEBRAIC[PhiCaL_i]*STATES[d]*( ALGEBRAIC[fp]*(1.00000 - STATES[nca_i])+ STATES[jca]*ALGEBRAIC[fcap]*STATES[nca_i]));
    ALGEBRAIC[ICaL] = ALGEBRAIC[ICaL_ss]+ALGEBRAIC[ICaL_i];
    ALGEBRAIC[ICaNa] = ALGEBRAIC[ICaNa_ss]+ALGEBRAIC[ICaNa_i];
    ALGEBRAIC[ICaK] = ALGEBRAIC[ICaK_ss]+ALGEBRAIC[ICaK_i];
    ALGEBRAIC[IpCa] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + GpCa]*STATES[cai])/(CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCap]+STATES[cai]);
    ALGEBRAIC[ICab] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + PCab]*4.00000*ALGEBRAIC[vffrt]*( ALGEBRAIC[gamma_cai]*STATES[cai]*exp( 2.00000*ALGEBRAIC[vfrt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_cao]*CONSTANTS[(sample_id * Tomek_num_of_constants) + cao]))/(exp( 2.00000*ALGEBRAIC[vfrt]) - 1.00000);
    ALGEBRAIC[IClCa_junc] =  (( CONSTANTS[(sample_id * Tomek_num_of_constants) + Fjunc]*CONSTANTS[(sample_id * Tomek_num_of_constants) + GClCa])/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KdClCa]/STATES[cass]))*(STATES[V] - CONSTANTS[(sample_id * Tomek_num_of_constants) + ECl]);
    ALGEBRAIC[IClCa_sl] =  (( (1.00000 - CONSTANTS[(sample_id * Tomek_num_of_constants) + Fjunc])*CONSTANTS[(sample_id * Tomek_num_of_constants) + GClCa])/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KdClCa]/STATES[cai]))*(STATES[V] - CONSTANTS[(sample_id * Tomek_num_of_constants) + ECl]);
    ALGEBRAIC[IClCa] = ALGEBRAIC[IClCa_junc]+ALGEBRAIC[IClCa_sl];
    ALGEBRAIC[IClb] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + GClb]*(STATES[V] - CONSTANTS[(sample_id * Tomek_num_of_constants) + ECl]);
    ALGEBRAIC[Jupnp] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + upScale]*0.00542500*STATES[cai])/(STATES[cai]+0.000920000);
    ALGEBRAIC[Jupp] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + upScale]*2.75000*0.00542500*STATES[cai])/((STATES[cai]+0.000920000) - 0.000170000);
    ALGEBRAIC[fJupp] = 1.00000/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaMK]/ALGEBRAIC[CaMKa]);
    ALGEBRAIC[Jleak] = ( 0.00488250*STATES[cansr])/15.0000;
    ALGEBRAIC[Jup] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Jup_b]*(( (1.00000 - ALGEBRAIC[fJupp])*ALGEBRAIC[Jupnp]+ ALGEBRAIC[fJupp]*ALGEBRAIC[Jupp]) - ALGEBRAIC[Jleak]);
    // ALGEBRAIC[Bcai] = 1.00000/(1.00000+( CONSTANTS[(sample_id * Tomek_num_of_constants) + cmdnmax]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcmdn])/pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcmdn]+STATES[cai], 2.00000)+( CONSTANTS[(sample_id * Tomek_num_of_constants) + trpnmax]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kmtrpn])/pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + kmtrpn]+STATES[cai], 2.00000));
    ALGEBRAIC[Bcai] = 1.00000/(1.00000+( CONSTANTS[(sample_id * Tomek_num_of_constants) + cmdnmax]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcmdn])/pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcmdn]+STATES[cai], 2.00000)); //modified
    ALGEBRAIC[Jtr] = (STATES[cansr] - STATES[cajsr])/60.0000;
    ALGEBRAIC[Bcajsr] = 1.00000/(1.00000+( CONSTANTS[(sample_id * Tomek_num_of_constants) + csqnmax]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcsqn])/pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcsqn]+STATES[cajsr], 2.00000));

    RATES[hL] = (ALGEBRAIC[hLss] - STATES[hL])/CONSTANTS[(sample_id * Tomek_num_of_constants) + thL];
    RATES[hLp] = (ALGEBRAIC[hLssp] - STATES[hLp])/CONSTANTS[(sample_id * Tomek_num_of_constants) + thLp];
    RATES[jca] = (ALGEBRAIC[jcass] - STATES[jca])/CONSTANTS[(sample_id * Tomek_num_of_constants) + tjca];
    RATES[m] = (ALGEBRAIC[mss] - STATES[m])/ALGEBRAIC[tm];
    RATES[mL] = (ALGEBRAIC[mLss] - STATES[mL])/ALGEBRAIC[tmL];
    RATES[a] = (ALGEBRAIC[ass] - STATES[a])/ALGEBRAIC[ta];
    RATES[d] = (ALGEBRAIC[dss] - STATES[d])/ALGEBRAIC[td];
    RATES[ff] = (ALGEBRAIC[fss] - STATES[ff])/ALGEBRAIC[tff];
    RATES[fs] = (ALGEBRAIC[fss] - STATES[fs])/ALGEBRAIC[tfs];
    RATES[nca_ss] =  ALGEBRAIC[anca_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] -  STATES[nca_ss]*ALGEBRAIC[km2n];
    RATES[nca_i] =  ALGEBRAIC[anca_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] -  STATES[nca_i]*ALGEBRAIC[km2n];
    RATES[xs1] = (ALGEBRAIC[xs1ss] - STATES[xs1])/ALGEBRAIC[txs1];
    RATES[ap] = (ALGEBRAIC[assp] - STATES[ap])/ALGEBRAIC[ta];
    RATES[fcaf] = (ALGEBRAIC[fcass] - STATES[fcaf])/ALGEBRAIC[tfcaf];
    RATES[fcas] = (ALGEBRAIC[fcass] - STATES[fcas])/ALGEBRAIC[tfcas];
    RATES[ffp] = (ALGEBRAIC[fss] - STATES[ffp])/ALGEBRAIC[tffp];
    RATES[xs2] = (ALGEBRAIC[xs2ss] - STATES[xs2])/ALGEBRAIC[txs2];
    RATES[CaMKt] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + aCaMK]*ALGEBRAIC[CaMKb]*(ALGEBRAIC[CaMKb]+STATES[CaMKt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + bCaMK]*STATES[CaMKt];
    RATES[h] = (ALGEBRAIC[hss] - STATES[h])/ALGEBRAIC[th];
    RATES[fcafp] = (ALGEBRAIC[fcass] - STATES[fcafp])/ALGEBRAIC[tfcafp];
    RATES[j] = (ALGEBRAIC[jss] - STATES[j])/ALGEBRAIC[tj];
    RATES[hp] = (ALGEBRAIC[hssp] - STATES[hp])/ALGEBRAIC[th];
    RATES[iF] = (ALGEBRAIC[iss] - STATES[iF])/ALGEBRAIC[tiF];
    RATES[C3] =  ALGEBRAIC[beta]*STATES[C2] -  ALGEBRAIC[alpha]*STATES[C3];
    RATES[C2] = ( ALGEBRAIC[alpha]*STATES[C3]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + beta_1]*STATES[C1]) -  (ALGEBRAIC[beta]+CONSTANTS[(sample_id * Tomek_num_of_constants) + alpha_1])*STATES[C2];
    RATES[jp] = (ALGEBRAIC[jss] - STATES[jp])/ALGEBRAIC[tjp];
    RATES[iS] = (ALGEBRAIC[iss] - STATES[iS])/ALGEBRAIC[tiS];
    RATES[O] = ( ALGEBRAIC[alpha_2]*STATES[C1]+ ALGEBRAIC[beta_i]*STATES[I]) -  (ALGEBRAIC[beta_2]+ALGEBRAIC[alpha_i])*STATES[O];
    RATES[iFp] = (ALGEBRAIC[iss] - STATES[iFp])/ALGEBRAIC[tiFp];
    RATES[iSp] = (ALGEBRAIC[iss] - STATES[iSp])/ALGEBRAIC[tiSp];
    RATES[C1] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + alpha_1]*STATES[C2]+ ALGEBRAIC[beta_2]*STATES[O]+ ALGEBRAIC[beta_ItoC2]*STATES[I]) -  (CONSTANTS[(sample_id * Tomek_num_of_constants) + beta_1]+ALGEBRAIC[alpha_2]+ALGEBRAIC[alpha_C2ToI])*STATES[C1];
    RATES[I] = ( ALGEBRAIC[alpha_C2ToI]*STATES[C1]+ ALGEBRAIC[alpha_i]*STATES[O]) -  (ALGEBRAIC[beta_ItoC2]+ALGEBRAIC[beta_i])*STATES[I];
    RATES[Jrel_np] = (ALGEBRAIC[Jrel_inf] - STATES[Jrel_np])/ALGEBRAIC[tau_rel];
    RATES[Jrel_p] = (ALGEBRAIC[Jrel_infp] - STATES[Jrel_p])/ALGEBRAIC[tau_relp];
    RATES[ki] = ( - (((ALGEBRAIC[Ito]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[IKb]+ALGEBRAIC[I_katp]+ALGEBRAIC[Istim]) -  2.00000*ALGEBRAIC[INaK])+ALGEBRAIC[ICaK_i])*CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo])+( ALGEBRAIC[JdiffK]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vss])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo];
    RATES[kss] = ( - ALGEBRAIC[ICaK_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vss]) - ALGEBRAIC[JdiffK];
    RATES[nai] = ( - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ 3.00000*ALGEBRAIC[INaCa_i]+ALGEBRAIC[ICaNa_i]+ 3.00000*ALGEBRAIC[INaK]+ALGEBRAIC[INab])*CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo])+( ALGEBRAIC[JdiffNa]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vss])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo];
    RATES[nass] = ( - (ALGEBRAIC[ICaNa_ss]+ 3.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vss]) - ALGEBRAIC[JdiffNa];
    RATES[cass] =  ALGEBRAIC[Bcass]*((( - (ALGEBRAIC[ICaL_ss] -  2.00000*ALGEBRAIC[INaCa_ss])*CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap])/( 2.00000*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vss])+( ALGEBRAIC[Jrel]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vjsr])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vss]) - ALGEBRAIC[Jdiff]);
    RATES[V] = - (ALGEBRAIC[INa]+ALGEBRAIC[INaL]+ALGEBRAIC[Ito]+ALGEBRAIC[ICaL]+ALGEBRAIC[ICaNa]+ALGEBRAIC[ICaK]+ALGEBRAIC[IKr]+ALGEBRAIC[IKs]+ALGEBRAIC[IK1]+ALGEBRAIC[INaCa_i]+ALGEBRAIC[INaCa_ss]+ALGEBRAIC[INaK]+ALGEBRAIC[INab]+ALGEBRAIC[IKb]+ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]+ALGEBRAIC[IClCa]+ALGEBRAIC[IClb]+ALGEBRAIC[I_katp]+ALGEBRAIC[Istim]);
    // new for coupling
    RATES[ca_trpn] = CONSTANTS[(sample_id * Tomek_num_of_constants) + trpnmax] * land_trpn;
    // RATES[cai] =  ALGEBRAIC[Bcai]*((( - ((ALGEBRAIC[ICaL_i]+ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]) -  2.00000*ALGEBRAIC[INaCa_i])*CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap])/( 2.00000*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo]) - ( ALGEBRAIC[Jup]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vnsr])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo])+( ALGEBRAIC[Jdiff]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vss])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo]);
    RATES[cai] =  ALGEBRAIC[Bcai]*((( - ((ALGEBRAIC[IpCa]+ALGEBRAIC[ICab]) -  2.00000*ALGEBRAIC[INaCa_i])* /*CONSTANTS[(sample_id * Tomek_num_of_constants) + cm]* */CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap])/( 2.00000*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo]) - ( ALGEBRAIC[Jup]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vnsr])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo])+( ALGEBRAIC[Jdiff]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vss])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo] - RATES[ca_trpn]); // modified -> cm is unknown here
    RATES[cansr] = ALGEBRAIC[Jup] - ( ALGEBRAIC[Jtr]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vjsr])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vnsr];
    RATES[cajsr] =  ALGEBRAIC[Bcajsr]*(ALGEBRAIC[Jtr] - ALGEBRAIC[Jrel]);

    }

//// freeze rk4 since unused for now
// __device__ void solveRK4(double TIME, double dt)
// {
// 	double k1[43],k23[43];
// 	double yk123[43];
// 	int idx;


// 	// assuming first computeRates() have been executed
// 	computeRates( TIME, CONSTANTS, RATES, STATES, ALGEBRAIC );
// 	for( idx = 0; idx < states_size; idx++ ) {
// 		k1[idx] = RATES[idx];
// 		yk123[idx] = STATES[idx] + (k1[idx]*dt*0.5);
// 	}
// 	computeRates( TIME+(dt*0.5), CONSTANTS, RATES, yk123, ALGEBRAIC );
// 	for( idx = 0; idx < states_size; idx++ ) {
// 		k23[idx] = RATES[idx];
// 		yk123[idx] = STATES[idx] + (k23[idx]*dt*0.5);
// 	}
// 	computeRates( TIME+(dt*0.5), CONSTANTS, RATES, yk123, ALGEBRAIC );
//   for( idx = 0; idx < states_size; idx++ ) {
//     k23[idx] += RATES[idx];
//     yk123[idx] = STATES[idx] + (k23[idx]*dt);
//   }
//   computeRates( TIME+dt, CONSTANTS, RATES, yk123, ALGEBRAIC );
// 	for( idx = 0; idx < states_size; idx++ ) {
// 		STATES[idx] += (k1[idx]+(2*k23[idx])+RATES[idx])/6. * dt;
//   }


// }

__device__ void solveAnalytical(double *CONSTANTS, double *STATES, double *ALGEBRAIC, double *RATES, double dt, int sample_id)
    {
    #ifdef EULER
    STATES[V] = STATES[V] + RATES[V] * dt;
    STATES[CaMKt] = STATES[CaMKt] + RATES[CaMKt] * dt;
    STATES[cass] = STATES[cass] + RATES[cass] * dt;
    STATES[nai] = STATES[nai] + RATES[nai] * dt;
    STATES[nass] = STATES[nass] + RATES[nass] * dt;
    STATES[ki] = STATES[ki] + RATES[ki] * dt;
    STATES[kss] = STATES[kss] + RATES[kss] * dt;
    STATES[cansr] = STATES[cansr] + RATES[cansr] * dt;
    STATES[cajsr] = STATES[cajsr] + RATES[cajsr] * dt;
    STATES[cai] = STATES[cai] + RATES[cai] * dt;
    STATES[m] = STATES[m] + RATES[m] * dt;
    STATES[h] = STATES[h] + RATES[h] * dt;
    STATES[j] = STATES[j] + RATES[j] * dt;
    STATES[hp] = STATES[hp] + RATES[hp] * dt;
    STATES[jp] = STATES[jp] + RATES[jp] * dt;
    STATES[mL] = STATES[mL] + RATES[mL] * dt;
    STATES[hL] = STATES[hL] + RATES[hL] * dt;
    STATES[hLp] = STATES[hLp] + RATES[hLp] * dt;
    STATES[a] = STATES[a] + RATES[a] * dt;
    STATES[iF] = STATES[iF] + RATES[iF] * dt;
    STATES[iS] = STATES[iS] + RATES[iS] * dt;
    STATES[ap] = STATES[ap] + RATES[ap] * dt;
    STATES[iFp] = STATES[iFp] + RATES[iFp] * dt;
    STATES[iSp] = STATES[iSp] + RATES[iSp] * dt;
    STATES[d] = STATES[d] + RATES[d] * dt;
    STATES[ff] = STATES[ff] + RATES[ff] * dt;
    STATES[fs] = STATES[fs] + RATES[fs] * dt;
    STATES[fcaf] = STATES[fcaf] + RATES[fcaf] * dt;
    STATES[fcas] = STATES[fcas] + RATES[fcas] * dt;
    STATES[jca] = STATES[jca] + RATES[jca] * dt;
    STATES[ffp] = STATES[ffp] + RATES[ffp] * dt;
    STATES[fcafp] = STATES[fcafp] + RATES[fcafp] * dt;
    STATES[nca_ss] = STATES[nca_ss] + RATES[nca_ss] * dt;
    STATES[nca_i] = STATES[nca_i] + RATES[nca_i] * dt;
    STATES[O] = STATES[O] + RATES[O] * dt;
    STATES[I] = STATES[I] + RATES[I] * dt;
        STATES[C3] = STATES[C3] + RATES[C3] * dt;
        STATES[C2] = STATES[C2] + RATES[C2] * dt;
        STATES[C1] = STATES[C1] + RATES[C1] * dt;
    STATES[xs1] = STATES[xs1] + RATES[xs1] * dt;
    STATES[xs2] = STATES[xs2] + RATES[xs2] * dt;
    STATES[Jrel_np] = STATES[Jrel_np] + RATES[Jrel_np] * dt;
    STATES[Jrel_p] = STATES[Jrel_p] + RATES[Jrel_p] * dt;
    #else
    ////==============
    ////Exact solution
    ////==============
    ////INa
    STATES[m] = ALGEBRAIC[mss] - (ALGEBRAIC[mss] - STATES[m]) * exp(-dt / ALGEBRAIC[tm]);
    STATES[h] = ALGEBRAIC[hss] - (ALGEBRAIC[hss] - STATES[h]) * exp(-dt / ALGEBRAIC[th]);
    STATES[j] = ALGEBRAIC[jss] - (ALGEBRAIC[jss] - STATES[j]) * exp(-dt / ALGEBRAIC[tj]);
    STATES[hp] = ALGEBRAIC[hssp] - (ALGEBRAIC[hssp] - STATES[hp]) * exp(-dt / ALGEBRAIC[th]);
    STATES[jp] = ALGEBRAIC[jss] - (ALGEBRAIC[jss] - STATES[jp]) * exp(-dt / ALGEBRAIC[tjp]);
    STATES[mL] = ALGEBRAIC[mLss] - (ALGEBRAIC[mLss] - STATES[mL]) * exp(-dt / ALGEBRAIC[tmL]);
    STATES[hL] = ALGEBRAIC[hLss] - (ALGEBRAIC[hLss] - STATES[hL]) * exp(-dt / CONSTANTS[(sample_id * Tomek_num_of_constants) + thL]);
    STATES[hLp] = ALGEBRAIC[hLssp] - (ALGEBRAIC[hLssp] - STATES[hLp]) * exp(-dt / CONSTANTS[(sample_id * Tomek_num_of_constants) + thLp]);
    ////Ito
    STATES[a] = ALGEBRAIC[ass] - (ALGEBRAIC[ass] - STATES[a]) * exp(-dt / ALGEBRAIC[ta]);
    STATES[iF] = ALGEBRAIC[iss] - (ALGEBRAIC[iss] - STATES[iF]) * exp(-dt / ALGEBRAIC[tiF]);
    STATES[iS] = ALGEBRAIC[iss] - (ALGEBRAIC[iss] - STATES[iS]) * exp(-dt / ALGEBRAIC[tiS]);
    STATES[ap] = ALGEBRAIC[assp] - (ALGEBRAIC[assp] - STATES[ap]) * exp(-dt / ALGEBRAIC[ta]);
    STATES[iFp] = ALGEBRAIC[iss] - (ALGEBRAIC[iss] - STATES[iFp]) * exp(-dt / ALGEBRAIC[tiFp]);
    STATES[iSp] = ALGEBRAIC[iss] - (ALGEBRAIC[iss] - STATES[iSp]) * exp(-dt / ALGEBRAIC[tiSp]);
    ////ICaL
    STATES[d] = ALGEBRAIC[dss] - (ALGEBRAIC[dss] - STATES[d]) * exp(-dt / ALGEBRAIC[td]);
    STATES[ff] = ALGEBRAIC[fss] - (ALGEBRAIC[fss] - STATES[ff]) * exp(-dt / ALGEBRAIC[tff]);
    STATES[fs] = ALGEBRAIC[fss] - (ALGEBRAIC[fss] - STATES[fs]) * exp(-dt / ALGEBRAIC[tfs]);
    STATES[fcaf] = ALGEBRAIC[fcass] - (ALGEBRAIC[fcass] - STATES[fcaf]) * exp(-dt / ALGEBRAIC[tfcaf]);
    STATES[fcas] = ALGEBRAIC[fcass] - (ALGEBRAIC[fcass] - STATES[fcas]) * exp(-dt / ALGEBRAIC[tfcas]);
    STATES[jca] = ALGEBRAIC[jcass] - (ALGEBRAIC[jcass] - STATES[jca]) * exp(- dt / CONSTANTS[(sample_id * Tomek_num_of_constants) + tjca]);
    STATES[ffp] = ALGEBRAIC[fss] - (ALGEBRAIC[fss] - STATES[ffp]) * exp(-dt / ALGEBRAIC[tffp]);
    STATES[fcafp] = ALGEBRAIC[fcass] - (ALGEBRAIC[fcass] - STATES[fcafp]) * exp(-d / ALGEBRAIC[tfcafp]);
        STATES[nca_i] = STATES[nca_i] + RATES[nca_i]*dt;
        STATES[nca_ss] = STATES[nca_ss] + RATES[nca_ss]*dt;
    //  STATES[nca_i] = ALGEBRAIC[anca_i] * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] / ALGEBRAIC[km2n] -
    //      (ALGEBRAIC[anca_i] * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] / ALGEBRAIC[km2n] - STATES[nca_i]) * exp(-ALGEBRAIC[km2n] * dt);
    //  STATES[nca_ss] = ALGEBRAIC[anca_ss] * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] / ALGEBRAIC[km2n] -
    //      (ALGEBRAIC[anca_ss] * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] / ALGEBRAIC[km2n] - STATES[nca_ss]) * exp(-ALGEBRAIC[km2n] * dt);
    ////IKr
    //STATES[O] = STATES[O] + RATES[O] * dt;
    //STATES[I] = STATES[I] + RATES[I] * dt;
    //STATES[C3] = STATES[C3] + RATES[C3] * dt;
    //STATES[C2] = STATES[C2] + RATES[C2] * dt;
    //STATES[C1] = STATES[C1] + RATES[C1] * dt;
    double* coeffs = new double[15];
    coeffs[0] = -  (ALGEBRAIC[beta_2]+ALGEBRAIC[alpha_i]);
    coeffs[1] = ALGEBRAIC[beta_i];
    coeffs[2] = ALGEBRAIC[alpha_2];
    coeffs[3] = ALGEBRAIC[alpha_i];
    coeffs[4] = -  (ALGEBRAIC[beta_ItoC2]+ALGEBRAIC[beta_i]);
    coeffs[5] = ALGEBRAIC[alpha_C2ToI];
    coeffs[6] = ALGEBRAIC[beta_2];
    coeffs[7] = ALGEBRAIC[beta_ItoC2];
    coeffs[8] = -  (CONSTANTS[(sample_id * Tomek_num_of_constants) + beta_1]+ALGEBRAIC[alpha_2]+ALGEBRAIC[alpha_C2ToI]);
    coeffs[9] = CONSTANTS[(sample_id * Tomek_num_of_constants) + alpha_1];
    coeffs[10] = CONSTANTS[(sample_id * Tomek_num_of_constants) + beta_1];
    coeffs[11] = -  (ALGEBRAIC[beta]+CONSTANTS[(sample_id * Tomek_num_of_constants) + alpha_1]);
    coeffs[12] = ALGEBRAIC[alpha];
    coeffs[13] = ALGEBRAIC[beta];
    coeffs[14] = -  ALGEBRAIC[alpha];
    int m = 5;
    double* a = new double[m*m]; // Flattened a
    a[0 * m + 0] = 1.0 - dt * coeffs[0];   a[0 * m + 1] = - dt * coeffs[1];     a[0 * m + 2] = - dt * coeffs[2];     a[0 * m + 3] = 0.0;                      a[0 * m + 4] = 0.0;
    a[1 * m + 0] = - dt * coeffs[3];       a[1 * m + 1] = 1.0 - dt * coeffs[4]; a[1 * m + 2] = - dt * coeffs[5];     a[1 * m + 3] = 0.0;                      a[1 * m + 4] = 0.0;
    a[2 * m + 0] = - dt * coeffs[6];       a[2 * m + 1] = - dt * coeffs[7];     a[2 * m + 2] = 1.0 - dt * coeffs[8]; a[2 * m + 3] = - dt * coeffs[9];         a[2 * m + 4] = 0.0;
    a[3 * m + 0] = 0.0;                    a[3 * m + 1] = 0.0;                  a[3 * m + 2] = - dt * coeffs[10];    a[3 * m + 3] = 1.0 - dt * coeffs[11];    a[3 * m + 4] = - dt * coeffs[12];
    a[4 * m + 0] = 0.0;                    a[4 * m + 1] = 0.0;                  a[4 * m + 2] = 0.0;                  a[4 * m + 3] = - dt * coeffs[13];;       a[4 * m + 4] = 1.0 - dt * coeffs[14];
    double* b = new double[m];
    b[0] = STATES[O];
    b[1] = STATES[I];
    b[2] = STATES[C1];
    b[3] = STATES[C2];
    b[4] = STATES[C3];
    double* x = new double[m];
    for(int i = 0; i < m; i++){
        x[i] = 0.0;
    }
    ___gaussElimination(a,b,x,m);
    STATES[O] = x[0];
    STATES[I] = x[1];
    STATES[C1] = x[2];
    STATES[C2] = x[3];
    STATES[C3] = x[4];
    delete[] coeffs;
    delete[] a;
    delete[] b;
    delete[] x;
    
    ////IKs
    STATES[xs1] = ALGEBRAIC[xs1ss] - (ALGEBRAIC[xs1ss] - STATES[xs1]) * exp(-dt / ALGEBRAIC[txs1]);
    STATES[xs2] = ALGEBRAIC[xs2ss] - (ALGEBRAIC[xs2ss] - STATES[xs2]) * exp(-dt / ALGEBRAIC[txs2]);
    ////IK1
    ////RyR receptors
    STATES[Jrel_np] = ALGEBRAIC[Jrel_inf] - (ALGEBRAIC[Jrel_inf] - STATES[Jrel_np]) * exp(-dt / ALGEBRAIC[tau_rel]);
    STATES[Jrel_p] = ALGEBRAIC[Jrel_infp] - (ALGEBRAIC[Jrel_infp] - STATES[Jrel_p]) * exp(-dt / ALGEBRAIC[tau_relp]);
    ////=============================
    ////Approximated solution (Euler)
    ////=============================
    ////CaMK
    STATES[CaMKt] = STATES[CaMKt] + RATES[CaMKt] * dt;
    ////Membrane potential
    STATES[V] = STATES[V] + RATES[V] * dt;
    ////Ion Concentrations and Buffers
    STATES[nai] = STATES[nai] + RATES[nai] * dt;
    STATES[nass] = STATES[nass] + RATES[nass] * dt;
    STATES[ki] = STATES[ki] + RATES[ki] * dt;
    STATES[kss] = STATES[kss] + RATES[kss] * dt;
    STATES[cai] = STATES[cai] + RATES[cai] * dt;
    STATES[cass] = STATES[cass] + RATES[cass] * dt;
    STATES[cansr] = STATES[cansr] + RATES[cansr] * dt;
    STATES[cajsr] = STATES[cajsr] + RATES[cajsr] * dt;
    #endif

    }

__device__ void ___gaussElimination(double *A, double *b, double *x, int N) {
        // Using A as a flat array to represent an N x N matrix
    for (int i = 0; i < N; i++) {
        // Search for maximum in this column
        double maxEl = fabs(A[i*N + i]);
        int maxRow = i;
        for (int k = i + 1; k < N; k++) {
            if (fabs(A[k*N + i]) > maxEl) {
                maxEl = fabs(A[k*N + i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k = i; k < N; k++) {
            double tmp = A[maxRow*N + k];
            A[maxRow*N + k] = A[i*N + k];
            A[i*N + k] = tmp;
        }
        double tmp = b[maxRow];
        b[maxRow] = b[i];
        b[i] = tmp;

        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < N; k++) {
            double c = -A[k*N + i] / A[i*N + i];
            for (int j = i; j < N; j++) {
                if (i == j) {
                    A[k*N + j] = 0;
                } else {
                    A[k*N + j] += c * A[i*N + j];
                }
            }
            b[k] += c * b[i];
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    for (int i = N - 1; i >= 0; i--) {
        x[i] = b[i] / A[i*N + i];
        for (int k = i - 1; k >= 0; k--) {
            b[k] -= A[k*N + i] * x[i];
        }
    }
}

__device__ double set_time_step (double TIME, double time_point, double max_time_step, double *CONSTANTS, double *RATES, int sample_id) {
 double min_time_step = 0.005;
 double time_step = min_time_step;
 double min_dV = 0.2;
 double max_dV = 0.8;

 
 if (TIME <= time_point || (TIME - floor(TIME / CONSTANTS[(sample_id * Tomek_num_of_constants) + BCL]) * CONSTANTS[(sample_id * Tomek_num_of_constants) + BCL]) <= time_point) {
    //printf("TIME <= time_point ms\n");
    return time_step;
    //printf("TIME = %E, dV = %E, time_step = %E\n",TIME, RATES[V] * time_step, time_step);
  }
  else {
    //printf("TIME > time_point ms\n");
    if (std::abs(RATES[(rates_size * sample_id) + V] * time_step) <= min_dV) {//Slow changes in V
        //printf("dV/dt <= 0.2\n");
        time_step = std::abs(max_dV / RATES[(rates_size * sample_id) + V]);
        //Make sure time_step is between min time step and max_time_step
        if (time_step < min_time_step) {
            time_step = min_time_step;
        }
        else if (time_step > max_time_step) {
            time_step = max_time_step;
        }
        //printf("TIME = %E, dV = %E, time_step = %E\n",TIME, RATES[V] * time_step, time_step);
    }
    else if (std::abs(RATES[V] * time_step) >= max_dV) {//Fast changes in V
        //printf("dV/dt >= 0.8\n");
        time_step = std::abs(min_dV / RATES[V]);
        //Make sure time_step is not less than 0.005
        if (time_step < min_time_step) {
            time_step = min_time_step;
        }
        //printf("TIME = %E, dV = %E, time_step = %E\n",TIME, RATES[V] * time_step, time_step);
    } else {
        time_step = min_time_step;
    }
    return time_step;
  }
}