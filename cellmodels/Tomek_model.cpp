/*
   There are a total of 223 entries in the algebraic variable array.
   There are a total of 43 entries in each of the rate and state variable arrays.
   There are a total of 163 entries in the constant variable array.
 */

#include "Tomek_model.hpp"
#include <cmath>
#include <cstdlib>
#include "./enums/enum_Tomek_model.hpp"
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
 * STATES[(sample_id * Tomek_num_of_states) + V] is v in component membrane (millivolt).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vffrt] is vffrt in component membrane (coulomb_per_mole).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt] is vfrt in component membrane (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INa] is INa in component INa (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaL] is INaL in component INaL (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ito] is Ito in component Ito (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaL] is ICaL in component ICaL (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaNa] is ICaNa in component ICaL (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaK] is ICaK in component ICaL (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IKr] is IKr in component IKr (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IKs] is IKs in component IKs (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IK1] is IK1 in component IK1 (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaCa_i] is INaCa_i in component INaCa (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaCa_ss] is INaCa_ss in component INaCa (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaK] is INaK in component INaK (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INab] is INab in component INab (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IKb] is IKb in component IKb (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IpCa] is IpCa in component IpCa (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICab] is ICab in component ICab (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IClCa] is IClCa in component ICl (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IClb] is IClb in component ICl (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + I_katp] is I_katp in component I_katp (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Istim] is Istim in component membrane (microA_per_microF).
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
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + CaMKb] is CaMKb in component CaMK (millimolar).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + CaMKa] is CaMKa in component CaMK (millimolar).
 * STATES[(sample_id * Tomek_num_of_states) + CaMKt] is CaMKt in component CaMK (millimolar).
 * STATES[(sample_id * Tomek_num_of_states) + cass] is cass in component intracellular_ions (millimolar).
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
 * STATES[(sample_id * Tomek_num_of_states) + nai] is nai in component intracellular_ions (millimolar).
 * STATES[(sample_id * Tomek_num_of_states) + nass] is nass in component intracellular_ions (millimolar).
 * STATES[(sample_id * Tomek_num_of_states) + ki] is ki in component intracellular_ions (millimolar).
 * STATES[(sample_id * Tomek_num_of_states) + kss] is kss in component intracellular_ions (millimolar).
 * STATES[(sample_id * Tomek_num_of_states) + cansr] is cansr in component intracellular_ions (millimolar).
 * STATES[(sample_id * Tomek_num_of_states) + cajsr] is cajsr in component intracellular_ions (millimolar).
 * STATES[(sample_id * Tomek_num_of_states) + cai] is cai in component intracellular_ions (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + cli] is cli in component intracellular_ions (millimolar).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaL_ss] is ICaL_ss in component ICaL (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaNa_ss] is ICaNa_ss in component ICaL (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaK_ss] is ICaK_ss in component ICaL (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaL_i] is ICaL_i in component ICaL (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaNa_i] is ICaNa_i in component ICaL (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaK_i] is ICaK_i in component ICaL (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JdiffNa] is JdiffNa in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jdiff] is Jdiff in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jup] is Jup in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JdiffK] is JdiffK in component diff (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel] is Jrel in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jtr] is Jtr in component trans_flux (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Bcai] is Bcai in component intracellular_ions (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Bcajsr] is Bcajsr in component intracellular_ions (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Bcass] is Bcass in component intracellular_ions (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PKNa] is PKNa in component reversal_potentials (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ENa] is ENa in component reversal_potentials (millivolt).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + EK] is EK in component reversal_potentials (millivolt).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + EKs] is EKs in component reversal_potentials (millivolt).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + ECl] is ECl in component reversal_potentials (millivolt).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + gkatp] is gkatp in component I_katp (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + fkatp] is fkatp in component I_katp (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + K_o_n] is K_o_n in component I_katp (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + A_atp] is A_atp in component I_katp (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + K_atp] is K_atp in component I_katp (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + akik] is akik in component I_katp (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + bkik] is bkik in component I_katp (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + mss] is mss in component INa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tm] is tm in component INa (millisecond).
 * STATES[(sample_id * Tomek_num_of_states) + m] is m in component INa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hss] is hss in component INa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ah] is ah in component INa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + bh] is bh in component INa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + th] is th in component INa (millisecond).
 * STATES[(sample_id * Tomek_num_of_states) + h] is h in component INa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + jss] is jss in component INa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + aj] is aj in component INa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + bj] is bj in component INa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tj] is tj in component INa (millisecond).
 * STATES[(sample_id * Tomek_num_of_states) + j] is j in component INa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hssp] is hssp in component INa (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + hp] is hp in component INa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tjp] is tjp in component INa (millisecond).
 * STATES[(sample_id * Tomek_num_of_states) + jp] is jp in component INa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fINap] is fINap in component INa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GNa] is GNa in component INa (milliS_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + mLss] is mLss in component INaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tmL] is tmL in component INaL (millisecond).
 * STATES[(sample_id * Tomek_num_of_states) + mL] is mL in component INaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + thL] is thL in component INaL (millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hLss] is hLss in component INaL (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + hL] is hL in component INaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hLssp] is hLssp in component INaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + thLp] is thLp in component INaL (millisecond).
 * STATES[(sample_id * Tomek_num_of_states) + hLp] is hLp in component INaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL_b] is GNaL_b in component INaL (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL] is GNaL in component INaL (milliS_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fINaLp] is fINaLp in component INaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Gto_b] is Gto_b in component Ito (milliS_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ass] is ass in component Ito (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ta] is ta in component Ito (millisecond).
 * STATES[(sample_id * Tomek_num_of_states) + a] is a in component Ito (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift] is EKshift in component Ito (millivolt).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + iss] is iss in component Ito (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + delta_epi] is delta_epi in component Ito (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiF_b] is tiF_b in component Ito (millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiS_b] is tiS_b in component Ito (millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiF] is tiF in component Ito (millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiS] is tiS in component Ito (millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + AiF] is AiF in component Ito (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + AiS] is AiS in component Ito (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + iF] is iF in component Ito (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + iS] is iS in component Ito (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + i] is i in component Ito (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + assp] is assp in component Ito (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + ap] is ap in component Ito (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + dti_develop] is dti_develop in component Ito (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + dti_recover] is dti_recover in component Ito (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiFp] is tiFp in component Ito (millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiSp] is tiSp in component Ito (millisecond).
 * STATES[(sample_id * Tomek_num_of_states) + iFp] is iFp in component Ito (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + iSp] is iSp in component Ito (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ip] is ip in component Ito (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Gto] is Gto in component Ito (milliS_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fItop] is fItop in component Ito (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Kmn] is Kmn in component ICaL (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] is k2n in component ICaL (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa_b] is PCa_b in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + dss] is dss in component ICaL (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + d] is d in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fss] is fss in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Aff] is Aff in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Afs] is Afs in component ICaL (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + ff] is ff in component ICaL (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + fs] is fs in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + f] is f in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcass] is fcass in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + jcass] is jcass in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Afcaf] is Afcaf in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Afcas] is Afcas in component ICaL (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + fcaf] is fcaf in component ICaL (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + fcas] is fcas in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fca] is fca in component ICaL (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + jca] is jca in component ICaL (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + ffp] is ffp in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fp] is fp in component ICaL (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + fcafp] is fcafp in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcap] is fcap in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + km2n] is km2n in component ICaL (per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + anca_ss] is anca_ss in component ICaL (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + nca_ss] is nca_ss in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + anca_i] is anca_i in component ICaL (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + nca_i] is nca_i in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaL_ss] is PhiCaL_ss in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaNa_ss] is PhiCaNa_ss in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaK_ss] is PhiCaK_ss in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaL_i] is PhiCaL_i in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaNa_i] is PhiCaNa_i in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaK_i] is PhiCaK_i in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa] is PCa in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCap] is PCap in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaNa] is PCaNa in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaK] is PCaK in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaNap] is PCaNap in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaKp] is PCaKp in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fICaLp] is fICaLp in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + td] is td in component ICaL (millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tff] is tff in component ICaL (millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfs] is tfs in component ICaL (millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfcaf] is tfcaf in component ICaL (millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfcas] is tfcas in component ICaL (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + tjca] is tjca in component ICaL (millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tffp] is tffp in component ICaL (millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfcafp] is tfcafp in component ICaL (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + vShift] is vShift in component ICaL (millivolt).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + sample_id] is sample_id in component ICaL (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Io] is Io in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Iss] is Iss in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ii] is Ii in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + dielConstant] is dielConstant in component ICaL (per_kelvin).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + constA] is constA in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_cao] is gamma_cao in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_cass] is gamma_cass in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_cai] is gamma_cai in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_nao] is gamma_nao in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_nass] is gamma_nass in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_nai] is gamma_nai in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_ko] is gamma_ko in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_kss] is gamma_kss in component ICaL (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_ki] is gamma_ki in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS] is ICaL_fractionSS in component ICaL (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr_b] is GKr_b in component IKr (milliS_per_microF).
 * STATES[(sample_id * Tomek_num_of_states) + C1] is C1 in component IKr (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + C2] is C2 in component IKr (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + C3] is C3 in component IKr (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + I] is I in component IKr (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + O] is O in component IKr (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha] is alpha in component IKr (per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta] is beta in component IKr (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + alpha_1] is alpha_1 in component IKr (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + beta_1] is beta_1 in component IKr (per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_2] is alpha_2 in component IKr (per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_2] is beta_2 in component IKr (per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_i] is alpha_i in component IKr (per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_i] is beta_i in component IKr (per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_C2ToI] is alpha_C2ToI in component IKr (per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_ItoC2] is beta_ItoC2 in component IKr (per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr] is GKr in component IKr (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs_b] is GKs_b in component IKs (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs] is GKs in component IKs (milliS_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + xs1ss] is xs1ss in component IKs (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + xs2ss] is xs2ss in component IKs (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + txs1] is txs1 in component IKs (millisecond).
 * STATES[(sample_id * Tomek_num_of_states) + xs1] is xs1 in component IKs (dimensionless).
 * STATES[(sample_id * Tomek_num_of_states) + xs2] is xs2 in component IKs (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + KsCa] is KsCa in component IKs (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + txs2] is txs2 in component IKs (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1] is GK1 in component IK1 (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1_b] is GK1_b in component IK1 (milliS_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + aK1] is aK1 in component IK1 (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + bK1] is bK1 in component IK1 (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + K1ss] is K1ss in component IK1 (dimensionless).
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
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hna] is hna in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hca] is hca in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaAct] is KmCaAct in component INaCa (millimolar).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Gncx_b] is Gncx_b in component INaCa (milliS_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Gncx] is Gncx in component INaCa (milliS_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h1_i] is h1_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h2_i] is h2_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h3_i] is h3_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h4_i] is h4_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h5_i] is h5_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h6_i] is h6_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h7_i] is h7_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h8_i] is h8_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h9_i] is h9_i in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + h10_i] is h10_i in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + h11_i] is h11_i in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + h12_i] is h12_i in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_i] is k1_i in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_i] is k2_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3p_i] is k3p_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3pp_i] is k3pp_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3_i] is k3_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4_i] is k4_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4p_i] is k4p_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4pp_i] is k4pp_i in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_i] is k5_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k6_i] is k6_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k7_i] is k7_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k8_i] is k8_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1_i] is x1_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2_i] is x2_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3_i] is x3_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4_i] is x4_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E1_i] is E1_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E2_i] is E2_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E3_i] is E3_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E4_i] is E4_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + allo_i] is allo_i in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JncxNa_i] is JncxNa_i in component INaCa (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JncxCa_i] is JncxCa_i in component INaCa (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h1_ss] is h1_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h2_ss] is h2_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h3_ss] is h3_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h4_ss] is h4_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h5_ss] is h5_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h6_ss] is h6_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h7_ss] is h7_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h8_ss] is h8_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h9_ss] is h9_ss in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + h10_ss] is h10_ss in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + h11_ss] is h11_ss in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + h12_ss] is h12_ss in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_ss] is k1_ss in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_ss] is k2_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3p_ss] is k3p_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3pp_ss] is k3pp_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3_ss] is k3_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4_ss] is k4_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4p_ss] is k4p_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4pp_ss] is k4pp_ss in component INaCa (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_ss] is k5_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k6_ss] is k6_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k7_ss] is k7_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k8_ss] is k8_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1_ss] is x1_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2_ss] is x2_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3_ss] is x3_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4_ss] is x4_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E1_ss] is E1_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E2_ss] is E2_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E3_ss] is E3_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E4_ss] is E4_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + allo_ss] is allo_ss in component INaCa (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JncxNa_ss] is JncxNa_ss in component INaCa (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JncxCa_ss] is JncxCa_ss in component INaCa (millimolar_per_millisecond).
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
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Knai] is Knai in component INaK (millimolar).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Knao] is Knao in component INaK (millimolar).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + P] is P in component INaK (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a1] is a1 in component INaK (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + b1] is b1 in component INaK (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + a2] is a2 in component INaK (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b2] is b2 in component INaK (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a3] is a3 in component INaK (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b3] is b3 in component INaK (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + a4] is a4 in component INaK (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b4] is b4 in component INaK (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1] is x1 in component INaK (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2] is x2 in component INaK (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3] is x3 in component INaK (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4] is x4 in component INaK (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E1] is E1 in component INaK (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E2] is E2 in component INaK (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E3] is E3 in component INaK (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E4] is E4 in component INaK (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JnakNa] is JnakNa in component INaK (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JnakK] is JnakK in component INaK (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + xkb] is xkb in component IKb (dimensionless).
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
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IClCa_junc] is IClCa_junc in component ICl (microA_per_microF).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IClCa_sl] is IClCa_sl in component ICl (microA_per_microF).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + tauNa] is tauNa in component diff (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + tauK] is tauK in component diff (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + tauCa] is tauCa in component diff (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + bt] is bt in component ryr (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + a_rel] is a_rel in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_inf_b] is Jrel_inf_b in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_inf] is Jrel_inf in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_rel_b] is tau_rel_b in component ryr (millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_rel] is tau_rel in component ryr (millisecond).
 * STATES[(sample_id * Tomek_num_of_states) + Jrel_np] is Jrel_np in component ryr (millimolar_per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + btp] is btp in component ryr (millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + a_relp] is a_relp in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_infp_b] is Jrel_infp_b in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_infp] is Jrel_infp in component ryr (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_relp_b] is tau_relp_b in component ryr (millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_relp] is tau_relp in component ryr (millisecond).
 * STATES[(sample_id * Tomek_num_of_states) + Jrel_p] is Jrel_p in component ryr (millimolar_per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + cajsr_half] is cajsr_half in component ryr (millimolar).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fJrelp] is fJrelp in component ryr (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Jrel_b] is Jrel_b in component ryr (dimensionless).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + upScale] is upScale in component SERCA (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jupnp] is Jupnp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jupp] is Jupp in component SERCA (millimolar_per_millisecond).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fJupp] is fJupp in component SERCA (dimensionless).
 * ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jleak] is Jleak in component SERCA (millimolar_per_millisecond).
 * CONSTANTS[(sample_id * Tomek_num_of_constants) + Jup_b] is Jup_b in component SERCA (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + V] is d/dt v in component membrane (millivolt).
 * RATES[(sample_id * Tomek_num_of_rates) + CaMKt] is d/dt CaMKt in component CaMK (millimolar).
 * RATES[(sample_id * Tomek_num_of_rates) + nai] is d/dt nai in component intracellular_ions (millimolar).
 * RATES[(sample_id * Tomek_num_of_rates) + nass] is d/dt nass in component intracellular_ions (millimolar).
 * RATES[(sample_id * Tomek_num_of_rates) + ki] is d/dt ki in component intracellular_ions (millimolar).
 * RATES[(sample_id * Tomek_num_of_rates) + kss] is d/dt kss in component intracellular_ions (millimolar).
 * RATES[(sample_id * Tomek_num_of_rates) + cai] is d/dt cai in component intracellular_ions (millimolar).
 * RATES[(sample_id * Tomek_num_of_rates) + cass] is d/dt cass in component intracellular_ions (millimolar).
 * RATES[(sample_id * Tomek_num_of_rates) + cansr] is d/dt cansr in component intracellular_ions (millimolar).
 * RATES[(sample_id * Tomek_num_of_rates) + cajsr] is d/dt cajsr in component intracellular_ions (millimolar).
 * RATES[(sample_id * Tomek_num_of_rates) + m] is d/dt m in component INa (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + h] is d/dt h in component INa (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + j] is d/dt j in component INa (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + hp] is d/dt hp in component INa (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + jp] is d/dt jp in component INa (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + mL] is d/dt mL in component INaL (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + hL] is d/dt hL in component INaL (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + hLp] is d/dt hLp in component INaL (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + a] is d/dt a in component Ito (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + iF] is d/dt iF in component Ito (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + iS] is d/dt iS in component Ito (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + ap] is d/dt ap in component Ito (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + iFp] is d/dt iFp in component Ito (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + iSp] is d/dt iSp in component Ito (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + d] is d/dt d in component ICaL (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + ff] is d/dt ff in component ICaL (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + fs] is d/dt fs in component ICaL (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + fcaf] is d/dt fcaf in component ICaL (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + fcas] is d/dt fcas in component ICaL (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + jca] is d/dt jca in component ICaL (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + ffp] is d/dt ffp in component ICaL (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + fcafp] is d/dt fcafp in component ICaL (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + nca_ss] is d/dt nca_ss in component ICaL (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + nca_i] is d/dt nca_i in component ICaL (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + C3] is d/dt C3 in component IKr (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + C2] is d/dt C2 in component IKr (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + C1] is d/dt C1 in component IKr (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + O] is d/dt O in component IKr (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + I] is d/dt I in component IKr (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + xs1] is d/dt xs1 in component IKs (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + xs2] is d/dt xs2 in component IKs (dimensionless).
 * RATES[(sample_id * Tomek_num_of_rates) + Jrel_np] is d/dt Jrel_np in component ryr (millimolar_per_millisecond).
 * RATES[(sample_id * Tomek_num_of_rates) + Jrel_p] is d/dt Jrel_p in component ryr (millimolar_per_millisecond).
 */


// Tomek_model::Tomek_model()
// {

// }

// Tomek_model::~Tomek_model()
// {

// }

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
    STATES[(sample_id * Tomek_num_of_states) + V] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ?  -89.14 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ?  -89.1704 : -88.7638);
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
    STATES[(sample_id * Tomek_num_of_states) + CaMKt] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.0129 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.0192 : 0.0111);
    STATES[(sample_id * Tomek_num_of_states) + cass] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 5.77E-05 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 6.58E-05 : 7.0305e-5);
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
    STATES[(sample_id * Tomek_num_of_states) + nai] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 12.1025 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 15.0038 : 12.1025);
    STATES[(sample_id * Tomek_num_of_states) + nass] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 12.8366 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 15.0043 : 12.1029);
    STATES[(sample_id * Tomek_num_of_states) + ki] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 142.6951 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 143.0403 : 142.3002);
    STATES[(sample_id * Tomek_num_of_states) + kss] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 142.6951 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 143.0402 : 142.3002);
    STATES[(sample_id * Tomek_num_of_states) + cansr] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 1.8119 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 1.9557 : 1.5211);
    STATES[(sample_id * Tomek_num_of_states) + cajsr] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 1.8102 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 1.9593 : 1.5214);
    STATES[(sample_id * Tomek_num_of_states) + cai] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 6.63E-05 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 8.17E-05 : 8.1583e-05);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + cli] = 24.0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + PKNa] = 0.01833;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + gkatp] = 4.3195;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + fkatp] = 0.0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + K_o_n] = 5;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + A_atp] = 2;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + K_atp] = 0.25;
    STATES[(sample_id * Tomek_num_of_states) + m] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 7.43E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 7.38E-04 : 8.0572e-4);
    STATES[(sample_id * Tomek_num_of_states) + h] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.836 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.8365 : 0.8286);
    STATES[(sample_id * Tomek_num_of_states) + j] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.8359 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.8363 : 0.8284);
    STATES[(sample_id * Tomek_num_of_states) + hp] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.6828 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.6838 : 0.6707);
    STATES[(sample_id * Tomek_num_of_states) + jp] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.8357 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.8358 : 0.8281);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GNa] = 11.7802;
    STATES[(sample_id * Tomek_num_of_states) + mL] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 1.52E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 1.51E-04 : 1.629e-4);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + thL] = 200;
    STATES[(sample_id * Tomek_num_of_states) + hL] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.5401 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.5327 : 0.5255);
    STATES[(sample_id * Tomek_num_of_states) + hLp] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.3034 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.2834 : 0.2872);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL_b] = 0.0279;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Gto_b] = 0.16;
    STATES[(sample_id * Tomek_num_of_states) + a] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 9.27E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 9.25E-04 : 9.5098e-4);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift] = 0;
    STATES[(sample_id * Tomek_num_of_states) + iF] = 0.9996;
    STATES[(sample_id * Tomek_num_of_states) + iS] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.9996 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.5671 : 0.5936);
    STATES[(sample_id * Tomek_num_of_states) + ap] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 4.72E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 4.71E-04 : 4.8454e-4);
    STATES[(sample_id * Tomek_num_of_states) + iFp] = 0.9996;
    STATES[(sample_id * Tomek_num_of_states) + iSp] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.9996 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.6261 :0.6538);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Kmn] = 0.002;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] = 500;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa_b] = 8.3757e-05;
    STATES[(sample_id * Tomek_num_of_states) + d] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.0 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.0 : 8.1084e-9);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + Aff] = 0.6;
    STATES[(sample_id * Tomek_num_of_states) + ff] = 1.0;
    STATES[(sample_id * Tomek_num_of_states) + fs] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.9485 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.92 : 0.939);
    STATES[(sample_id * Tomek_num_of_states) + fcaf] = 1.0;
    STATES[(sample_id * Tomek_num_of_states) + fcas] = 0.9999;
    STATES[(sample_id * Tomek_num_of_states) + jca] = 1.0;
    STATES[(sample_id * Tomek_num_of_states) + ffp] = 1.0;
    STATES[(sample_id * Tomek_num_of_states) + fcafp] = 1.0;
    STATES[(sample_id * Tomek_num_of_states) + nca_ss] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 3.09E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 5.14E-04 : 6.6462e-4);
    STATES[(sample_id * Tomek_num_of_states) + nca_i] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 5.30E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.0012 : 0.0012);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + tjca] = 75;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + vShift] = 0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + sample_id] = 0;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + dielConstant] = 74;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS] = 0.8;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr_b] = 0.0321;
    STATES[(sample_id * Tomek_num_of_states) + C1] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 6.79E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 6.96E-04 : 7.0344e-4);
    STATES[(sample_id * Tomek_num_of_states) + C2] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 8.29E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 8.27E-04 : 8.5109e-4);
    STATES[(sample_id * Tomek_num_of_states) + C3] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.9982 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.9979 : 0.9981);
    STATES[(sample_id * Tomek_num_of_states) + I] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 9.54E-06 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 1.88E-05 : 1.3289e-5);
    STATES[(sample_id * Tomek_num_of_states) + O] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 2.76E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 5.42E-04 : 3.7585e-4);
    CONSTANTS[(sample_id * Tomek_num_of_constants) + alpha_1] = 0.154375;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + beta_1] = 0.1911;
    CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs_b] = 0.0011;
    STATES[(sample_id * Tomek_num_of_states) + xs1] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0.2309 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0.2653 : 0.248);
    STATES[(sample_id * Tomek_num_of_states) + xs2] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 1.70E-04 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 1.69E-04 : 1.7707e-4);
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
    STATES[(sample_id * Tomek_num_of_states) + Jrel_np] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 2.82E-24 : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0. : 1.6129e-22);
    STATES[(sample_id * Tomek_num_of_states) + Jrel_p] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 0. : CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ? 0. : 1.2475e-20);
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
	___initConsts(CONSTANTS, STATES, type, bcl, sample_id);
	printf("Celltype: %lf\n", CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]);
	#ifndef COMPONENT_PATCH
	printf("Control %lf %lf %lf %lf %lf\n", CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa], CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1], CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs], CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL], CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr]);
	#endif
	___applyDrugEffect(CONSTANTS, conc, ic50, sample_id);
	#ifndef COMPONENT_PATCH
	printf("After drug %lf %lf %lf %lf %lf\n", CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa], CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1], CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs], CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL], CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr]);
	#endif
    }

__device__ void computeRates(double TIME, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int sample_id,  double land_trpn)
    {
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hLss] = 1.00000/(1.00000+exp((STATES[(sample_id * Tomek_num_of_states) + V]+87.6100)/7.48800));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hLssp] = 1.00000/(1.00000+exp((STATES[(sample_id * Tomek_num_of_states) + V]+93.8100)/7.48800));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + jcass] = 1.00000/(1.00000+exp((STATES[(sample_id * Tomek_num_of_states) + V]+18.0800)/2.79160));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + mss] = 1.00000/pow(1.00000+exp(- (STATES[(sample_id * Tomek_num_of_states) + V]+56.8600)/9.03000), 2.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tm] =  0.129200*exp(- pow((STATES[(sample_id * Tomek_num_of_states) + V]+45.7900)/15.5400, 2.00000))+ 0.0648700*exp(- pow((STATES[(sample_id * Tomek_num_of_states) + V] - 4.82300)/51.1200, 2.00000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + mLss] = 1.00000/(1.00000+exp(- (STATES[(sample_id * Tomek_num_of_states) + V]+42.8500)/5.26400));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tmL] =  0.129200*exp(- pow((STATES[(sample_id * Tomek_num_of_states) + V]+45.7900)/15.5400, 2.00000))+ 0.0648700*exp(- pow((STATES[(sample_id * Tomek_num_of_states) + V] - 4.82300)/51.1200, 2.00000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ass] = 1.00000/(1.00000+exp(- ((STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]) - 14.3400)/14.8200));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ta] = 1.05150/(1.00000/( 1.20890*(1.00000+exp(- ((STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]) - 18.4099)/29.3814)))+3.50000/(1.00000+exp((STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+100.000)/29.3814)));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + dss] = (STATES[(sample_id * Tomek_num_of_states) + V]>=31.4978 ? 1.00000 :  1.07630*exp( - 1.00700*exp( - 0.0829000*STATES[(sample_id * Tomek_num_of_states) + V])));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + td] = CONSTANTS[(sample_id * Tomek_num_of_constants) + sample_id]+0.600000+1.00000/(exp( - 0.0500000*(STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + vShift]+6.00000))+exp( 0.0900000*(STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + vShift]+14.0000)));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fss] = 1.00000/(1.00000+exp((STATES[(sample_id * Tomek_num_of_states) + V]+19.5800)/3.69600));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tff] = 7.00000+1.00000/( 0.00450000*exp(- (STATES[(sample_id * Tomek_num_of_states) + V]+20.0000)/10.0000)+ 0.00450000*exp((STATES[(sample_id * Tomek_num_of_states) + V]+20.0000)/10.0000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfs] = 1000.00+1.00000/( 3.50000e-05*exp(- (STATES[(sample_id * Tomek_num_of_states) + V]+5.00000)/4.00000)+ 3.50000e-05*exp((STATES[(sample_id * Tomek_num_of_states) + V]+5.00000)/6.00000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + km2n] =  STATES[(sample_id * Tomek_num_of_states) + jca]*1.00000;
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + anca_ss] = 1.00000/(CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n]/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + km2n]+pow(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + Kmn]/STATES[(sample_id * Tomek_num_of_states) + cass], 4.00000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + anca_i] = 1.00000/(CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n]/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + km2n]+pow(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + Kmn]/STATES[(sample_id * Tomek_num_of_states) + cai], 4.00000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + xs1ss] = 1.00000/(1.00000+exp(- (STATES[(sample_id * Tomek_num_of_states) + V]+11.6000)/8.93200));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + txs1] = 817.300+1.00000/( 0.000232600*exp((STATES[(sample_id * Tomek_num_of_states) + V]+48.2800)/17.8000)+ 0.00129200*exp(- (STATES[(sample_id * Tomek_num_of_states) + V]+210.000)/230.000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + assp] = 1.00000/(1.00000+exp(- ((STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]) - 24.3400)/14.8200));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcass] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fss];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfcaf] = 7.00000+1.00000/( 0.0400000*exp(- (STATES[(sample_id * Tomek_num_of_states) + V] - 4.00000)/7.00000)+ 0.0400000*exp((STATES[(sample_id * Tomek_num_of_states) + V] - 4.00000)/7.00000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfcas] = 100.000+1.00000/( 0.000120000*exp(- STATES[(sample_id * Tomek_num_of_states) + V]/3.00000)+ 0.000120000*exp(STATES[(sample_id * Tomek_num_of_states) + V]/7.00000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tffp] =  2.50000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tff];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + xs2ss] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + xs1ss];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + txs2] = 1.00000/( 0.0100000*exp((STATES[(sample_id * Tomek_num_of_states) + V] - 50.0000)/20.0000)+ 0.0193000*exp(- (STATES[(sample_id * Tomek_num_of_states) + V]+66.5400)/31.0000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + CaMKb] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + CaMKo]*(1.00000 - STATES[(sample_id * Tomek_num_of_states) + CaMKt]))/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaM]/STATES[(sample_id * Tomek_num_of_states) + cass]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hss] = 1.00000/pow(1.00000+exp((STATES[(sample_id * Tomek_num_of_states) + V]+71.5500)/7.43000), 2.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ah] = (STATES[(sample_id * Tomek_num_of_states) + V]>=- 40.0000 ? 0.00000 :  0.0570000*exp(- (STATES[(sample_id * Tomek_num_of_states) + V]+80.0000)/6.80000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + bh] = (STATES[(sample_id * Tomek_num_of_states) + V]>=- 40.0000 ? 0.770000/( 0.130000*(1.00000+exp(- (STATES[(sample_id * Tomek_num_of_states) + V]+10.6600)/11.1000))) :  2.70000*exp( 0.0790000*STATES[(sample_id * Tomek_num_of_states) + V])+ 310000.*exp( 0.348500*STATES[(sample_id * Tomek_num_of_states) + V]));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + th] = 1.00000/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ah]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + bh]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfcafp] =  2.50000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfcaf];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + jss] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hss];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + aj] = (STATES[(sample_id * Tomek_num_of_states) + V]>=- 40.0000 ? 0.00000 : ( ( - 25428.0*exp( 0.244400*STATES[(sample_id * Tomek_num_of_states) + V]) -  6.94800e-06*exp( - 0.0439100*STATES[(sample_id * Tomek_num_of_states) + V]))*(STATES[(sample_id * Tomek_num_of_states) + V]+37.7800))/(1.00000+exp( 0.311000*(STATES[(sample_id * Tomek_num_of_states) + V]+79.2300))));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + bj] = (STATES[(sample_id * Tomek_num_of_states) + V]>=- 40.0000 ? ( 0.600000*exp( 0.0570000*STATES[(sample_id * Tomek_num_of_states) + V]))/(1.00000+exp( - 0.100000*(STATES[(sample_id * Tomek_num_of_states) + V]+32.0000))) : ( 0.0242400*exp( - 0.0105200*STATES[(sample_id * Tomek_num_of_states) + V]))/(1.00000+exp( - 0.137800*(STATES[(sample_id * Tomek_num_of_states) + V]+40.1400))));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tj] = 1.00000/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + aj]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + bj]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hssp] = 1.00000/pow(1.00000+exp((STATES[(sample_id * Tomek_num_of_states) + V]+77.5500)/7.43000), 2.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + iss] = 1.00000/(1.00000+exp((STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+43.9400)/5.71100));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + delta_epi] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==1.00000 ? 1.00000 - 0.950000/(1.00000+exp((STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+70.0000)/5.00000)) : 1.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiF_b] = 4.56200+1.00000/( 0.393300*exp(- (STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+100.000)/100.000)+ 0.0800400*exp((STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+50.0000)/16.5900));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiF] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiF_b]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + delta_epi];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt] = ( STATES[(sample_id * Tomek_num_of_states) + V]*CONSTANTS[(sample_id * Tomek_num_of_constants) + F])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + R]*CONSTANTS[(sample_id * Tomek_num_of_constants) + T]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha] =  0.116100*exp( 0.299000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta] =  0.244200*exp( - 1.60400*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tjp] =  1.46000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tj];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiS_b] = 23.6200+1.00000/( 0.00141600*exp(- (STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+96.5200)/59.0500)+ 1.78000e-08*exp((STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+114.100)/8.07900));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiS] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiS_b]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + delta_epi];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_2] =  0.0578000*exp( 0.971000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_2] =  0.000349000*exp( - 1.06200*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_i] =  0.253300*exp( 0.595300*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_i] =  0.0652500*exp( - 0.820900*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + dti_develop] = 1.35400+0.000100000/(exp(((STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]) - 167.400)/15.8900)+exp(- ((STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]) - 12.2300)/0.215400));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + dti_recover] = 1.00000 - 0.500000/(1.00000+exp((STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]+70.0000)/20.0000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiFp] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + dti_develop]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + dti_recover]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiF];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiSp] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + dti_develop]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + dti_recover]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiS];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_C2ToI] =  5.20000e-05*exp( 1.52500*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_ItoC2] = ( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_2]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_i]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_C2ToI])/( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_2]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_i]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + f] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Aff]*STATES[(sample_id * Tomek_num_of_states) + ff]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + Afs]*STATES[(sample_id * Tomek_num_of_states) + fs];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Afcaf] = 0.300000+0.600000/(1.00000+exp((STATES[(sample_id * Tomek_num_of_states) + V] - 10.0000)/10.0000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Afcas] = 1.00000 - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Afcaf];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fca] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Afcaf]*STATES[(sample_id * Tomek_num_of_states) + fcaf]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Afcas]*STATES[(sample_id * Tomek_num_of_states) + fcas];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fp] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Aff]*STATES[(sample_id * Tomek_num_of_states) + ffp]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + Afs]*STATES[(sample_id * Tomek_num_of_states) + fs];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcap] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Afcaf]*STATES[(sample_id * Tomek_num_of_states) + fcafp]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Afcas]*STATES[(sample_id * Tomek_num_of_states) + fcas];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vffrt] = ( STATES[(sample_id * Tomek_num_of_states) + V]*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + F])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + R]*CONSTANTS[(sample_id * Tomek_num_of_constants) + T]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Iss] = ( 0.500000*(STATES[(sample_id * Tomek_num_of_states) + nass]+STATES[(sample_id * Tomek_num_of_states) + kss]+CONSTANTS[(sample_id * Tomek_num_of_constants) + cli]+ 4.00000*STATES[(sample_id * Tomek_num_of_states) + cass]))/1000.00;
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_cass] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*4.00000*( pow(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Iss], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Iss], 1.0 / 2)) -  0.300000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Iss]));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaL_ss] = ( 4.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vffrt]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_cass]*STATES[(sample_id * Tomek_num_of_states) + cass]*exp( 2.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_cao]*CONSTANTS[(sample_id * Tomek_num_of_constants) + cao]))/(exp( 2.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) - 1.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + CaMKa] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + CaMKb]+STATES[(sample_id * Tomek_num_of_states) + CaMKt];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fICaLp] = 1.00000/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaMK]/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + CaMKa]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaL_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS]*( (1.00000 - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fICaLp])*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaL_ss]*STATES[(sample_id * Tomek_num_of_states) + d]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + f]*(1.00000 - STATES[(sample_id * Tomek_num_of_states) + nca_ss])+ STATES[(sample_id * Tomek_num_of_states) + jca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fca]*STATES[(sample_id * Tomek_num_of_states) + nca_ss])+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fICaLp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCap]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaL_ss]*STATES[(sample_id * Tomek_num_of_states) + d]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fp]*(1.00000 - STATES[(sample_id * Tomek_num_of_states) + nca_ss])+ STATES[(sample_id * Tomek_num_of_states) + jca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcap]*STATES[(sample_id * Tomek_num_of_states) + nca_ss]));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_inf_b] = (( - CONSTANTS[(sample_id * Tomek_num_of_constants) + a_rel]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaL_ss])/1.00000)/(1.00000+pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + cajsr_half]/STATES[(sample_id * Tomek_num_of_states) + cajsr], 8.00000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_inf] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ?  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_inf_b]*1.70000 : ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_inf_b]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_rel_b] = CONSTANTS[(sample_id * Tomek_num_of_constants) + bt]/(1.00000+0.0123000/STATES[(sample_id * Tomek_num_of_states) + cajsr]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_rel] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_rel_b]<0.00100000 ? 0.00100000 : ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_rel_b]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_infp_b] = (( - CONSTANTS[(sample_id * Tomek_num_of_constants) + a_relp]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaL_ss])/1.00000)/(1.00000+pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + cajsr_half]/STATES[(sample_id * Tomek_num_of_states) + cajsr], 8.00000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_infp] = (CONSTANTS[(sample_id * Tomek_num_of_constants) + celltype]==2.00000 ?  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_infp_b]*1.70000 : ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_infp_b]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_relp_b] = CONSTANTS[(sample_id * Tomek_num_of_constants) + btp]/(1.00000+0.0123000/STATES[(sample_id * Tomek_num_of_states) + cajsr]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_relp] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_relp_b]<0.00100000 ? 0.00100000 : ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_relp_b]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + EK] =  (( CONSTANTS[(sample_id * Tomek_num_of_constants) + R]*CONSTANTS[(sample_id * Tomek_num_of_constants) + T])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + zk]*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]))*log(CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]/STATES[(sample_id * Tomek_num_of_states) + ki]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + AiF] = 1.00000/(1.00000+exp(((STATES[(sample_id * Tomek_num_of_states) + V]+CONSTANTS[(sample_id * Tomek_num_of_constants) + EKshift]) - 213.600)/151.200));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + AiS] = 1.00000 - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + AiF];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + i] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + AiF]*STATES[(sample_id * Tomek_num_of_states) + iF]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + AiS]*STATES[(sample_id * Tomek_num_of_states) + iS];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ip] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + AiF]*STATES[(sample_id * Tomek_num_of_states) + iFp]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + AiS]*STATES[(sample_id * Tomek_num_of_states) + iSp];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fItop] = 1.00000/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaMK]/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + CaMKa]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ito] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Gto]*(STATES[(sample_id * Tomek_num_of_states) + V] - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + EK])*( (1.00000 - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fItop])*STATES[(sample_id * Tomek_num_of_states) + a]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + i]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fItop]*STATES[(sample_id * Tomek_num_of_states) + ap]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ip]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IKr] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + GKr]* pow((CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]/5.00000), 1.0 / 2)*STATES[(sample_id * Tomek_num_of_states) + O]*(STATES[(sample_id * Tomek_num_of_states) + V] - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + EK]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + EKs] =  (( CONSTANTS[(sample_id * Tomek_num_of_constants) + R]*CONSTANTS[(sample_id * Tomek_num_of_constants) + T])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + zk]*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]))*log((CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + PKNa]*CONSTANTS[(sample_id * Tomek_num_of_constants) + nao])/(STATES[(sample_id * Tomek_num_of_states) + ki]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + PKNa]*STATES[(sample_id * Tomek_num_of_states) + nai]));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + KsCa] = 1.00000+0.600000/(1.00000+pow(3.80000e-05/STATES[(sample_id * Tomek_num_of_states) + cai], 1.40000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IKs] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + GKs]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + KsCa]*STATES[(sample_id * Tomek_num_of_states) + xs1]*STATES[(sample_id * Tomek_num_of_states) + xs2]*(STATES[(sample_id * Tomek_num_of_states) + V] - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + EKs]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + aK1] = 4.09400/(1.00000+exp( 0.121700*((STATES[(sample_id * Tomek_num_of_states) + V] - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + EK]) - 49.9340)));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + bK1] = ( 15.7200*exp( 0.0674000*((STATES[(sample_id * Tomek_num_of_states) + V] - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + EK]) - 3.25700))+exp( 0.0618000*((STATES[(sample_id * Tomek_num_of_states) + V] - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + EK]) - 594.310)))/(1.00000+exp( - 0.162900*((STATES[(sample_id * Tomek_num_of_states) + V] - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + EK])+14.2070)));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + K1ss] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + aK1]/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + aK1]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + bK1]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IK1] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + GK1]* pow((CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]/5.00000), 1.0 / 2)*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + K1ss]*(STATES[(sample_id * Tomek_num_of_states) + V] - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + EK]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Knao] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Knao0]*exp(( (1.00000 - CONSTANTS[(sample_id * Tomek_num_of_constants) + delta])*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt])/3.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a3] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + k3p]*pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kko], 2.00000))/((pow(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Knao], 3.00000)+pow(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kko], 2.00000)) - 1.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + P] = CONSTANTS[(sample_id * Tomek_num_of_constants) + eP]/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + H]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Khp]+STATES[(sample_id * Tomek_num_of_states) + nai]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Knap]+STATES[(sample_id * Tomek_num_of_states) + ki]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kxkur]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b3] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + k3m]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + P]*CONSTANTS[(sample_id * Tomek_num_of_constants) + H])/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + MgATP]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kmgatp]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Knai] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Knai0]*exp(( CONSTANTS[(sample_id * Tomek_num_of_constants) + delta]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt])/3.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a1] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + k1p]*pow(STATES[(sample_id * Tomek_num_of_states) + nai]/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Knai], 3.00000))/((pow(1.00000+STATES[(sample_id * Tomek_num_of_states) + nai]/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Knai], 3.00000)+pow(1.00000+STATES[(sample_id * Tomek_num_of_states) + ki]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kki], 2.00000)) - 1.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b2] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + k2m]*pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Knao], 3.00000))/((pow(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Knao], 3.00000)+pow(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kko], 2.00000)) - 1.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b4] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + k4m]*pow(STATES[(sample_id * Tomek_num_of_states) + ki]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kki], 2.00000))/((pow(1.00000+STATES[(sample_id * Tomek_num_of_states) + nai]/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Knai], 3.00000)+pow(1.00000+STATES[(sample_id * Tomek_num_of_states) + ki]/CONSTANTS[(sample_id * Tomek_num_of_constants) + Kki], 2.00000)) - 1.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + a4]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a1]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a2]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b2]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b4]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b3]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + a2]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b4]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b3]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b3]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a1]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a2];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b2]*CONSTANTS[(sample_id * Tomek_num_of_constants) + b1]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b4]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a1]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a2]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a3]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a3]*CONSTANTS[(sample_id * Tomek_num_of_constants) + b1]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b4]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + a2]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a3]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b4];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + a2]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a3]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a4]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b3]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b2]*CONSTANTS[(sample_id * Tomek_num_of_constants) + b1]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b2]*CONSTANTS[(sample_id * Tomek_num_of_constants) + b1]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a4]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a3]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a4]*CONSTANTS[(sample_id * Tomek_num_of_constants) + b1];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b4]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b3]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b2]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a3]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a4]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a1]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b2]*CONSTANTS[(sample_id * Tomek_num_of_constants) + a4]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a1]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b3]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b2]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a1];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E1] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1]/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E2] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2]/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JnakNa] =  3.00000*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E1]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a3] -  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E2]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + b3]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E3] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3]/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E4] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4]/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JnakK] =  2.00000*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E4]*CONSTANTS[(sample_id * Tomek_num_of_constants) + b1] -  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E3]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + a1]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaK] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Pnak]*( CONSTANTS[(sample_id * Tomek_num_of_constants) + zna]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JnakNa]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + zk]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JnakK]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + xkb] = 1.00000/(1.00000+exp(- (STATES[(sample_id * Tomek_num_of_states) + V] - 10.8968)/23.9871));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IKb] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + GKb]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + xkb]*(STATES[(sample_id * Tomek_num_of_states) + V] - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + EK]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + I_katp] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + fkatp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + gkatp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + akik]*CONSTANTS[(sample_id * Tomek_num_of_constants) + bkik]*(STATES[(sample_id * Tomek_num_of_states) + V] - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + EK]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Istim] = (TIME>=CONSTANTS[(sample_id * Tomek_num_of_constants) + stim_start]&&TIME<=CONSTANTS[(sample_id * Tomek_num_of_constants) + i_Stim_End]&&(TIME - CONSTANTS[(sample_id * Tomek_num_of_constants) + stim_start]) -  floor((TIME - CONSTANTS[(sample_id * Tomek_num_of_constants) + stim_start])/CONSTANTS[(sample_id * Tomek_num_of_constants) + BCL])*CONSTANTS[(sample_id * Tomek_num_of_constants) + BCL]<=CONSTANTS[(sample_id * Tomek_num_of_constants) + i_Stim_PulseDuration] ? CONSTANTS[(sample_id * Tomek_num_of_constants) + i_Stim_Amplitude] : 0.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ii] = ( 0.500000*(STATES[(sample_id * Tomek_num_of_states) + nai]+STATES[(sample_id * Tomek_num_of_states) + ki]+CONSTANTS[(sample_id * Tomek_num_of_constants) + cli]+ 4.00000*STATES[(sample_id * Tomek_num_of_states) + cai]))/1000.00;
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_ki] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*1.00000*( pow(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ii], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ii], 1.0 / 2)) -  0.300000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ii]));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaK_i] = ( 1.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vffrt]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_ki]*STATES[(sample_id * Tomek_num_of_states) + ki]*exp( 1.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_ko]*CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]))/(exp( 1.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) - 1.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaK_i] =  (1.00000 - CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS])*( (1.00000 - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fICaLp])*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaK]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaK_i]*STATES[(sample_id * Tomek_num_of_states) + d]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + f]*(1.00000 - STATES[(sample_id * Tomek_num_of_states) + nca_i])+ STATES[(sample_id * Tomek_num_of_states) + jca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fca]*STATES[(sample_id * Tomek_num_of_states) + nca_i])+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fICaLp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaKp]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaK_i]*STATES[(sample_id * Tomek_num_of_states) + d]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fp]*(1.00000 - STATES[(sample_id * Tomek_num_of_states) + nca_i])+ STATES[(sample_id * Tomek_num_of_states) + jca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcap]*STATES[(sample_id * Tomek_num_of_states) + nca_i]));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JdiffK] = (STATES[(sample_id * Tomek_num_of_states) + kss] - STATES[(sample_id * Tomek_num_of_states) + ki])/CONSTANTS[(sample_id * Tomek_num_of_constants) + tauK];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_kss] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*1.00000*( pow(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Iss], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Iss], 1.0 / 2)) -  0.300000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Iss]));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaK_ss] = ( 1.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vffrt]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_kss]*STATES[(sample_id * Tomek_num_of_states) + kss]*exp( 1.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_ko]*CONSTANTS[(sample_id * Tomek_num_of_constants) + ko]))/(exp( 1.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) - 1.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaK_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS]*( (1.00000 - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fICaLp])*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaK]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaK_ss]*STATES[(sample_id * Tomek_num_of_states) + d]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + f]*(1.00000 - STATES[(sample_id * Tomek_num_of_states) + nca_ss])+ STATES[(sample_id * Tomek_num_of_states) + jca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fca]*STATES[(sample_id * Tomek_num_of_states) + nca_ss])+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fICaLp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaKp]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaK_ss]*STATES[(sample_id * Tomek_num_of_states) + d]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fp]*(1.00000 - STATES[(sample_id * Tomek_num_of_states) + nca_ss])+ STATES[(sample_id * Tomek_num_of_states) + jca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcap]*STATES[(sample_id * Tomek_num_of_states) + nca_ss]));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ENa] =  (( CONSTANTS[(sample_id * Tomek_num_of_constants) + R]*CONSTANTS[(sample_id * Tomek_num_of_constants) + T])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + zna]*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]))*log(CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/STATES[(sample_id * Tomek_num_of_states) + nai]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fINap] = 1.00000/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaMK]/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + CaMKa]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INa] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + GNa]*(STATES[(sample_id * Tomek_num_of_states) + V] - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ENa])*pow(STATES[(sample_id * Tomek_num_of_states) + m], 3.00000)*( (1.00000 - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fINap])*STATES[(sample_id * Tomek_num_of_states) + h]*STATES[(sample_id * Tomek_num_of_states) + j]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fINap]*STATES[(sample_id * Tomek_num_of_states) + hp]*STATES[(sample_id * Tomek_num_of_states) + jp]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fINaLp] = 1.00000/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaMK]/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + CaMKa]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaL] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + GNaL]*(STATES[(sample_id * Tomek_num_of_states) + V] - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ENa])*STATES[(sample_id * Tomek_num_of_states) + mL]*( (1.00000 - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fINaLp])*STATES[(sample_id * Tomek_num_of_states) + hL]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fINaLp]*STATES[(sample_id * Tomek_num_of_states) + hLp]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + allo_i] = 1.00000/(1.00000+pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaAct]/STATES[(sample_id * Tomek_num_of_states) + cai], 2.00000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hna] = exp( CONSTANTS[(sample_id * Tomek_num_of_constants) + qna]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h7_i] = 1.00000+ (CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3])*(1.00000+1.00000/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hna]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h8_i] = CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/( CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hna]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h7_i]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3pp_i] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h8_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wnaca];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h1_i] = 1.00000+ (STATES[(sample_id * Tomek_num_of_states) + nai]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3])*(1.00000+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hna]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h2_i] = ( STATES[(sample_id * Tomek_num_of_states) + nai]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hna])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h1_i]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4pp_i] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h2_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wnaca];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h4_i] = 1.00000+ (STATES[(sample_id * Tomek_num_of_states) + nai]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna1])*(1.00000+STATES[(sample_id * Tomek_num_of_states) + nai]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna2]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h5_i] = ( STATES[(sample_id * Tomek_num_of_states) + nai]*STATES[(sample_id * Tomek_num_of_states) + nai])/( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h4_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kna1]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kna2]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k7_i] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h5_i]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h2_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wna];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k8_i] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h8_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + h11_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wna];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h9_i] = 1.00000/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h7_i];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3p_i] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h9_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wca];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3_i] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3p_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3pp_i];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hca] = exp( CONSTANTS[(sample_id * Tomek_num_of_constants) + qca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h3_i] = 1.00000/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h1_i];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4p_i] = ( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h3_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wca])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hca];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4_i] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4p_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4pp_i];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h6_i] = 1.00000/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h4_i];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k6_i] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h6_i]*STATES[(sample_id * Tomek_num_of_states) + cai]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kcaon];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1_i] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_i]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4_i]*(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k7_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k6_i])+ CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_i]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k7_i]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3_i]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2_i] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_i]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k7_i]*(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4_i]+CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_i])+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4_i]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k6_i]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k8_i]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3_i] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_i]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3_i]*(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k7_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k6_i])+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k8_i]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k6_i]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3_i]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4_i] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_i]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k8_i]*(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4_i]+CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_i])+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_i]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k8_i]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E1_i] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1_i]/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4_i]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E2_i] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2_i]/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4_i]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E3_i] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3_i]/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4_i]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E4_i] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4_i]/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4_i]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JncxNa_i] = ( 3.00000*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E4_i]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k7_i] -  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E1_i]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k8_i])+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E3_i]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4pp_i]) -  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E2_i]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3pp_i];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JncxCa_i] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E2_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_i] -  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E1_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_i];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaCa_i] =  (1.00000 - CONSTANTS[(sample_id * Tomek_num_of_constants) + INaCa_fractionSS])*CONSTANTS[(sample_id * Tomek_num_of_constants) + Gncx]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + allo_i]*( CONSTANTS[(sample_id * Tomek_num_of_constants) + zna]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JncxNa_i]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + zca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JncxCa_i]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INab] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + PNab]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vffrt]*( STATES[(sample_id * Tomek_num_of_states) + nai]*exp(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) - CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]))/(exp(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) - 1.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_nai] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*1.00000*( pow(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ii], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ii], 1.0 / 2)) -  0.300000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ii]));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaNa_i] = ( 1.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vffrt]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_nai]*STATES[(sample_id * Tomek_num_of_states) + nai]*exp( 1.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_nao]*CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]))/(exp( 1.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) - 1.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaNa_i] =  (1.00000 - CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS])*( (1.00000 - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fICaLp])*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaNa]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaNa_i]*STATES[(sample_id * Tomek_num_of_states) + d]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + f]*(1.00000 - STATES[(sample_id * Tomek_num_of_states) + nca_i])+ STATES[(sample_id * Tomek_num_of_states) + jca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fca]*STATES[(sample_id * Tomek_num_of_states) + nca_i])+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fICaLp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaNap]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaNa_i]*STATES[(sample_id * Tomek_num_of_states) + d]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fp]*(1.00000 - STATES[(sample_id * Tomek_num_of_states) + nca_i])+ STATES[(sample_id * Tomek_num_of_states) + jca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcap]*STATES[(sample_id * Tomek_num_of_states) + nca_i]));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JdiffNa] = (STATES[(sample_id * Tomek_num_of_states) + nass] - STATES[(sample_id * Tomek_num_of_states) + nai])/CONSTANTS[(sample_id * Tomek_num_of_constants) + tauNa];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + allo_ss] = 1.00000/(1.00000+pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaAct]/STATES[(sample_id * Tomek_num_of_states) + cass], 2.00000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h7_ss] = 1.00000+ (CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3])*(1.00000+1.00000/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hna]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h8_ss] = CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]/( CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hna]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h7_ss]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3pp_ss] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h8_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wnaca];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h1_ss] = 1.00000+ (STATES[(sample_id * Tomek_num_of_states) + nass]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3])*(1.00000+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hna]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h2_ss] = ( STATES[(sample_id * Tomek_num_of_states) + nass]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hna])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + kna3]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h1_ss]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4pp_ss] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h2_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wnaca];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h4_ss] = 1.00000+ (STATES[(sample_id * Tomek_num_of_states) + nass]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna1])*(1.00000+STATES[(sample_id * Tomek_num_of_states) + nass]/CONSTANTS[(sample_id * Tomek_num_of_constants) + kna2]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h5_ss] = ( STATES[(sample_id * Tomek_num_of_states) + nass]*STATES[(sample_id * Tomek_num_of_states) + nass])/( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h4_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kna1]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kna2]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k7_ss] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h5_ss]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h2_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wna];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k8_ss] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h8_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + h11_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wna];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h9_ss] = 1.00000/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h7_ss];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3p_ss] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h9_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wca];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3_ss] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3p_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3pp_ss];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h3_ss] = 1.00000/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h1_ss];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4p_ss] = ( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h3_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + wca])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hca];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4_ss] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4p_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4pp_ss];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h6_ss] = 1.00000/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h4_ss];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k6_ss] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + h6_ss]*STATES[(sample_id * Tomek_num_of_states) + cass]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kcaon];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_ss]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4_ss]*(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k7_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k6_ss])+ CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_ss]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k7_ss]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3_ss]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_ss]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k7_ss]*(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4_ss]+CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_ss])+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4_ss]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k6_ss]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k8_ss]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_ss]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3_ss]*(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k7_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k6_ss])+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k8_ss]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k6_ss]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3_ss]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_ss]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k8_ss]*(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4_ss]+CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_ss])+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k5_ss]*(CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k8_ss]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E1_ss] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1_ss]/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4_ss]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E2_ss] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2_ss]/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4_ss]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E3_ss] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3_ss]/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4_ss]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E4_ss] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4_ss]/(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x1_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x2_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x3_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + x4_ss]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JncxNa_ss] = ( 3.00000*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E4_ss]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k7_ss] -  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E1_ss]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k8_ss])+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E3_ss]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k4pp_ss]) -  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E2_ss]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + k3pp_ss];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JncxCa_ss] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E2_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k2_ss] -  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + E1_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k1_ss];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaCa_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + INaCa_fractionSS]*CONSTANTS[(sample_id * Tomek_num_of_constants) + Gncx]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + allo_ss]*( CONSTANTS[(sample_id * Tomek_num_of_constants) + zna]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JncxNa_ss]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + zca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JncxCa_ss]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_nass] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*1.00000*( pow(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Iss], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Iss], 1.0 / 2)) -  0.300000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Iss]));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaNa_ss] = ( 1.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vffrt]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_nass]*STATES[(sample_id * Tomek_num_of_states) + nass]*exp( 1.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_nao]*CONSTANTS[(sample_id * Tomek_num_of_constants) + nao]))/(exp( 1.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) - 1.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaNa_ss] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS]*( (1.00000 - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fICaLp])*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaNa]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaNa_ss]*STATES[(sample_id * Tomek_num_of_states) + d]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + f]*(1.00000 - STATES[(sample_id * Tomek_num_of_states) + nca_ss])+ STATES[(sample_id * Tomek_num_of_states) + jca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fca]*STATES[(sample_id * Tomek_num_of_states) + nca_ss])+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fICaLp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCaNap]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaNa_ss]*STATES[(sample_id * Tomek_num_of_states) + d]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fp]*(1.00000 - STATES[(sample_id * Tomek_num_of_states) + nca_ss])+ STATES[(sample_id * Tomek_num_of_states) + jca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcap]*STATES[(sample_id * Tomek_num_of_states) + nca_ss]));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jdiff] = (STATES[(sample_id * Tomek_num_of_states) + cass] - STATES[(sample_id * Tomek_num_of_states) + cai])/CONSTANTS[(sample_id * Tomek_num_of_constants) + tauCa];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fJrelp] = 1.00000/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaMK]/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + CaMKa]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Jrel_b]*( (1.00000 - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fJrelp])*STATES[(sample_id * Tomek_num_of_states) + Jrel_np]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fJrelp]*STATES[(sample_id * Tomek_num_of_states) + Jrel_p]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Bcass] = 1.00000/(1.00000+( CONSTANTS[(sample_id * Tomek_num_of_constants) + BSRmax]*CONSTANTS[(sample_id * Tomek_num_of_constants) + KmBSR])/pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + KmBSR]+STATES[(sample_id * Tomek_num_of_states) + cass], 2.00000)+( CONSTANTS[(sample_id * Tomek_num_of_constants) + BSLmax]*CONSTANTS[(sample_id * Tomek_num_of_constants) + KmBSL])/pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + KmBSL]+STATES[(sample_id * Tomek_num_of_states) + cass], 2.00000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_cai] = exp( - CONSTANTS[(sample_id * Tomek_num_of_constants) + constA]*4.00000*( pow(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ii], 1.0 / 2)/(1.00000+ pow(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ii], 1.0 / 2)) -  0.300000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ii]));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaL_i] = ( 4.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vffrt]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_cai]*STATES[(sample_id * Tomek_num_of_states) + cai]*exp( 2.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_cao]*CONSTANTS[(sample_id * Tomek_num_of_constants) + cao]))/(exp( 2.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) - 1.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaL_i] =  (1.00000 - CONSTANTS[(sample_id * Tomek_num_of_constants) + ICaL_fractionSS])*( (1.00000 - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fICaLp])*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCa]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaL_i]*STATES[(sample_id * Tomek_num_of_states) + d]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + f]*(1.00000 - STATES[(sample_id * Tomek_num_of_states) + nca_i])+ STATES[(sample_id * Tomek_num_of_states) + jca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fca]*STATES[(sample_id * Tomek_num_of_states) + nca_i])+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fICaLp]*CONSTANTS[(sample_id * Tomek_num_of_constants) + PCap]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + PhiCaL_i]*STATES[(sample_id * Tomek_num_of_states) + d]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fp]*(1.00000 - STATES[(sample_id * Tomek_num_of_states) + nca_i])+ STATES[(sample_id * Tomek_num_of_states) + jca]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcap]*STATES[(sample_id * Tomek_num_of_states) + nca_i]));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaL] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaL_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaL_i];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaNa] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaNa_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaNa_i];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaK] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaK_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaK_i];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IpCa] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + GpCa]*STATES[(sample_id * Tomek_num_of_states) + cai])/(CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCap]+STATES[(sample_id * Tomek_num_of_states) + cai]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICab] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + PCab]*4.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vffrt]*( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + gamma_cai]*STATES[(sample_id * Tomek_num_of_states) + cai]*exp( 2.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + gamma_cao]*CONSTANTS[(sample_id * Tomek_num_of_constants) + cao]))/(exp( 2.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + vfrt]) - 1.00000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IClCa_junc] =  (( CONSTANTS[(sample_id * Tomek_num_of_constants) + Fjunc]*CONSTANTS[(sample_id * Tomek_num_of_constants) + GClCa])/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KdClCa]/STATES[(sample_id * Tomek_num_of_states) + cass]))*(STATES[(sample_id * Tomek_num_of_states) + V] - CONSTANTS[(sample_id * Tomek_num_of_constants) + ECl]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IClCa_sl] =  (( (1.00000 - CONSTANTS[(sample_id * Tomek_num_of_constants) + Fjunc])*CONSTANTS[(sample_id * Tomek_num_of_constants) + GClCa])/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KdClCa]/STATES[(sample_id * Tomek_num_of_states) + cai]))*(STATES[(sample_id * Tomek_num_of_states) + V] - CONSTANTS[(sample_id * Tomek_num_of_constants) + ECl]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IClCa] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IClCa_junc]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IClCa_sl];
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IClb] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + GClb]*(STATES[(sample_id * Tomek_num_of_states) + V] - CONSTANTS[(sample_id * Tomek_num_of_constants) + ECl]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jupnp] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + upScale]*0.00542500*STATES[(sample_id * Tomek_num_of_states) + cai])/(STATES[(sample_id * Tomek_num_of_states) + cai]+0.000920000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jupp] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + upScale]*2.75000*0.00542500*STATES[(sample_id * Tomek_num_of_states) + cai])/((STATES[(sample_id * Tomek_num_of_states) + cai]+0.000920000) - 0.000170000);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fJupp] = 1.00000/(1.00000+CONSTANTS[(sample_id * Tomek_num_of_constants) + KmCaMK]/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + CaMKa]);
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jleak] = ( 0.00488250*STATES[(sample_id * Tomek_num_of_states) + cansr])/15.0000;
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jup] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + Jup_b]*(( (1.00000 - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fJupp])*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jupnp]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fJupp]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jupp]) - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jleak]);
    // ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Bcai] = 1.00000/(1.00000+( CONSTANTS[(sample_id * Tomek_num_of_constants) + cmdnmax]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcmdn])/pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcmdn]+STATES[(sample_id * Tomek_num_of_states) + cai], 2.00000)+( CONSTANTS[(sample_id * Tomek_num_of_constants) + trpnmax]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kmtrpn])/pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + kmtrpn]+STATES[(sample_id * Tomek_num_of_states) + cai], 2.00000));
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Bcai] = 1.00000/(1.00000+( CONSTANTS[(sample_id * Tomek_num_of_constants) + cmdnmax]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcmdn])/pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcmdn]+STATES[(sample_id * Tomek_num_of_states) + cai], 2.00000)); //modified
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jtr] = (STATES[(sample_id * Tomek_num_of_states) + cansr] - STATES[(sample_id * Tomek_num_of_states) + cajsr])/60.0000;
    ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Bcajsr] = 1.00000/(1.00000+( CONSTANTS[(sample_id * Tomek_num_of_constants) + csqnmax]*CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcsqn])/pow(CONSTANTS[(sample_id * Tomek_num_of_constants) + kmcsqn]+STATES[(sample_id * Tomek_num_of_states) + cajsr], 2.00000));

    RATES[(sample_id * Tomek_num_of_rates) + hL] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hLss] - STATES[(sample_id * Tomek_num_of_states) + hL])/CONSTANTS[(sample_id * Tomek_num_of_constants) + thL];
    RATES[(sample_id * Tomek_num_of_rates) + hLp] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hLssp] - STATES[(sample_id * Tomek_num_of_states) + hLp])/CONSTANTS[(sample_id * Tomek_num_of_constants) + thLp];
    RATES[(sample_id * Tomek_num_of_rates) + jca] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + jcass] - STATES[(sample_id * Tomek_num_of_states) + jca])/CONSTANTS[(sample_id * Tomek_num_of_constants) + tjca];
    RATES[(sample_id * Tomek_num_of_rates) + m] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + mss] - STATES[(sample_id * Tomek_num_of_states) + m])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tm];
    RATES[(sample_id * Tomek_num_of_rates) + mL] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + mLss] - STATES[(sample_id * Tomek_num_of_states) + mL])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tmL];
    RATES[(sample_id * Tomek_num_of_rates) + a] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ass] - STATES[(sample_id * Tomek_num_of_states) + a])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ta];
    RATES[(sample_id * Tomek_num_of_rates) + d] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + dss] - STATES[(sample_id * Tomek_num_of_states) + d])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + td];
    RATES[(sample_id * Tomek_num_of_rates) + ff] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fss] - STATES[(sample_id * Tomek_num_of_states) + ff])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tff];
    RATES[(sample_id * Tomek_num_of_rates) + fs] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fss] - STATES[(sample_id * Tomek_num_of_states) + fs])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfs];
    RATES[(sample_id * Tomek_num_of_rates) + nca_ss] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + anca_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] -  STATES[(sample_id * Tomek_num_of_states) + nca_ss]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + km2n];
    RATES[(sample_id * Tomek_num_of_rates) + nca_i] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + anca_i]*CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] -  STATES[(sample_id * Tomek_num_of_states) + nca_i]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + km2n];
    RATES[(sample_id * Tomek_num_of_rates) + xs1] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + xs1ss] - STATES[(sample_id * Tomek_num_of_states) + xs1])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + txs1];
    RATES[(sample_id * Tomek_num_of_rates) + ap] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + assp] - STATES[(sample_id * Tomek_num_of_states) + ap])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ta];
    RATES[(sample_id * Tomek_num_of_rates) + fcaf] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcass] - STATES[(sample_id * Tomek_num_of_states) + fcaf])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfcaf];
    RATES[(sample_id * Tomek_num_of_rates) + fcas] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcass] - STATES[(sample_id * Tomek_num_of_states) + fcas])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfcas];
    RATES[(sample_id * Tomek_num_of_rates) + ffp] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fss] - STATES[(sample_id * Tomek_num_of_states) + ffp])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tffp];
    RATES[(sample_id * Tomek_num_of_rates) + xs2] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + xs2ss] - STATES[(sample_id * Tomek_num_of_states) + xs2])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + txs2];
    RATES[(sample_id * Tomek_num_of_rates) + CaMKt] =  CONSTANTS[(sample_id * Tomek_num_of_constants) + aCaMK]*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + CaMKb]*(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + CaMKb]+STATES[(sample_id * Tomek_num_of_states) + CaMKt]) -  CONSTANTS[(sample_id * Tomek_num_of_constants) + bCaMK]*STATES[(sample_id * Tomek_num_of_states) + CaMKt];
    RATES[(sample_id * Tomek_num_of_rates) + h] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hss] - STATES[(sample_id * Tomek_num_of_states) + h])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + th];
    RATES[(sample_id * Tomek_num_of_rates) + fcafp] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcass] - STATES[(sample_id * Tomek_num_of_states) + fcafp])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfcafp];
    RATES[(sample_id * Tomek_num_of_rates) + j] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + jss] - STATES[(sample_id * Tomek_num_of_states) + j])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tj];
    RATES[(sample_id * Tomek_num_of_rates) + hp] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hssp] - STATES[(sample_id * Tomek_num_of_states) + hp])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + th];
    RATES[(sample_id * Tomek_num_of_rates) + iF] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + iss] - STATES[(sample_id * Tomek_num_of_states) + iF])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiF];
    RATES[(sample_id * Tomek_num_of_rates) + C3] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta]*STATES[(sample_id * Tomek_num_of_states) + C2] -  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha]*STATES[(sample_id * Tomek_num_of_states) + C3];
    RATES[(sample_id * Tomek_num_of_rates) + C2] = ( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha]*STATES[(sample_id * Tomek_num_of_states) + C3]+ CONSTANTS[(sample_id * Tomek_num_of_constants) + beta_1]*STATES[(sample_id * Tomek_num_of_states) + C1]) -  (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta]+CONSTANTS[(sample_id * Tomek_num_of_constants) + alpha_1])*STATES[(sample_id * Tomek_num_of_states) + C2];
    RATES[(sample_id * Tomek_num_of_rates) + jp] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + jss] - STATES[(sample_id * Tomek_num_of_states) + jp])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tjp];
    RATES[(sample_id * Tomek_num_of_rates) + iS] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + iss] - STATES[(sample_id * Tomek_num_of_states) + iS])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiS];
    RATES[(sample_id * Tomek_num_of_rates) + O] = ( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_2]*STATES[(sample_id * Tomek_num_of_states) + C1]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_i]*STATES[(sample_id * Tomek_num_of_states) + I]) -  (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_2]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_i])*STATES[(sample_id * Tomek_num_of_states) + O];
    RATES[(sample_id * Tomek_num_of_rates) + iFp] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + iss] - STATES[(sample_id * Tomek_num_of_states) + iFp])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiFp];
    RATES[(sample_id * Tomek_num_of_rates) + iSp] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + iss] - STATES[(sample_id * Tomek_num_of_states) + iSp])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiSp];
    RATES[(sample_id * Tomek_num_of_rates) + C1] = ( CONSTANTS[(sample_id * Tomek_num_of_constants) + alpha_1]*STATES[(sample_id * Tomek_num_of_states) + C2]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_2]*STATES[(sample_id * Tomek_num_of_states) + O]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_ItoC2]*STATES[(sample_id * Tomek_num_of_states) + I]) -  (CONSTANTS[(sample_id * Tomek_num_of_constants) + beta_1]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_2]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_C2ToI])*STATES[(sample_id * Tomek_num_of_states) + C1];
    RATES[(sample_id * Tomek_num_of_rates) + I] = ( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_C2ToI]*STATES[(sample_id * Tomek_num_of_states) + C1]+ ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_i]*STATES[(sample_id * Tomek_num_of_states) + O]) -  (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_ItoC2]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_i])*STATES[(sample_id * Tomek_num_of_states) + I];
    RATES[(sample_id * Tomek_num_of_rates) + Jrel_np] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_inf] - STATES[(sample_id * Tomek_num_of_states) + Jrel_np])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_rel];
    RATES[(sample_id * Tomek_num_of_rates) + Jrel_p] = (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_infp] - STATES[(sample_id * Tomek_num_of_states) + Jrel_p])/ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_relp];
    RATES[(sample_id * Tomek_num_of_rates) + ki] = ( - (((ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ito]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IKr]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IKs]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IK1]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IKb]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + I_katp]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Istim]) -  2.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaK])+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaK_i])*CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo])+( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JdiffK]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vss])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo];
    RATES[(sample_id * Tomek_num_of_rates) + kss] = ( - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaK_ss]*CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vss]) - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JdiffK];
    RATES[(sample_id * Tomek_num_of_rates) + nai] = ( - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INa]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaL]+ 3.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaCa_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaNa_i]+ 3.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaK]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INab])*CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo])+( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JdiffNa]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vss])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo];
    RATES[(sample_id * Tomek_num_of_rates) + nass] = ( - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaNa_ss]+ 3.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaCa_ss])*CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap])/( CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vss]) - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + JdiffNa];
    RATES[(sample_id * Tomek_num_of_rates) + cass] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Bcass]*((( - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaL_ss] -  2.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaCa_ss])*CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap])/( 2.00000*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vss])+( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vjsr])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vss]) - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jdiff]);
    RATES[(sample_id * Tomek_num_of_rates) + V] = - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INa]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaL]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Ito]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaL]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaNa]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaK]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IKr]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IKs]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IK1]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaCa_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaCa_ss]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaK]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INab]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IKb]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IpCa]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICab]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IClCa]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IClb]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + I_katp]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Istim]);
    // new for coupling
    RATES[(sample_id * Tomek_num_of_rates) + ca_trpn] = CONSTANTS[(sample_id * Tomek_num_of_constants) + trpnmax] * land_trpn;
    // RATES[(sample_id * Tomek_num_of_rates) + cai] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Bcai]*((( - ((ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICaL_i]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IpCa]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICab]) -  2.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaCa_i])*CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap])/( 2.00000*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo]) - ( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jup]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vnsr])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo])+( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jdiff]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vss])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo]);
    RATES[(sample_id * Tomek_num_of_rates) + cai] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Bcai]*((( - ((ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + IpCa]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ICab]) -  2.00000*ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + INaCa_i])* /*CONSTANTS[(sample_id * Tomek_num_of_constants) + cm]* */CONSTANTS[(sample_id * Tomek_num_of_constants) + Acap])/( 2.00000*CONSTANTS[(sample_id * Tomek_num_of_constants) + F]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo]) - ( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jup]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vnsr])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo])+( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jdiff]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vss])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vmyo] - RATES[(sample_id * Tomek_num_of_rates) + ca_trpn]); // modified -> cm is unknown here
    RATES[(sample_id * Tomek_num_of_rates) + cansr] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jup] - ( ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jtr]*CONSTANTS[(sample_id * Tomek_num_of_constants) + vjsr])/CONSTANTS[(sample_id * Tomek_num_of_constants) + vnsr];
    RATES[(sample_id * Tomek_num_of_rates) + cajsr] =  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Bcajsr]*(ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jtr] - ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel]);

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
// 		k1[idx] = RATES[(sample_id * Tomek_num_of_rates) + idx];
// 		yk123[idx] = STATES[(sample_id * Tomek_num_of_states) + idx] + (k1[idx]*dt*0.5);
// 	}
// 	computeRates( TIME+(dt*0.5), CONSTANTS, RATES, yk123, ALGEBRAIC );
// 	for( idx = 0; idx < states_size; idx++ ) {
// 		k23[idx] = RATES[(sample_id * Tomek_num_of_rates) + idx];
// 		yk123[idx] = STATES[(sample_id * Tomek_num_of_states) + idx] + (k23[idx]*dt*0.5);
// 	}
// 	computeRates( TIME+(dt*0.5), CONSTANTS, RATES, yk123, ALGEBRAIC );
//   for( idx = 0; idx < states_size; idx++ ) {
//     k23[idx] += RATES[(sample_id * Tomek_num_of_rates) + idx];
//     yk123[idx] = STATES[(sample_id * Tomek_num_of_states) + idx] + (k23[idx]*dt);
//   }
//   computeRates( TIME+dt, CONSTANTS, RATES, yk123, ALGEBRAIC );
// 	for( idx = 0; idx < states_size; idx++ ) {
// 		STATES[(sample_id * Tomek_num_of_states) + idx] += (k1[idx]+(2*k23[idx])+RATES[(sample_id * Tomek_num_of_rates) + idx])/6. * dt;
//   }


// }

__device__ void solveAnalytical(double *CONSTANTS, double *STATES, double *ALGEBRAIC, double *RATES, double dt, int sample_id)
    {
    #ifdef EULER
    STATES[(sample_id * Tomek_num_of_states) + V] = STATES[(sample_id * Tomek_num_of_states) + V] + RATES[(sample_id * Tomek_num_of_rates) + V] * dt;
    STATES[(sample_id * Tomek_num_of_states) + CaMKt] = STATES[(sample_id * Tomek_num_of_states) + CaMKt] + RATES[(sample_id * Tomek_num_of_rates) + CaMKt] * dt;
    STATES[(sample_id * Tomek_num_of_states) + cass] = STATES[(sample_id * Tomek_num_of_states) + cass] + RATES[(sample_id * Tomek_num_of_rates) + cass] * dt;
    STATES[(sample_id * Tomek_num_of_states) + nai] = STATES[(sample_id * Tomek_num_of_states) + nai] + RATES[(sample_id * Tomek_num_of_rates) + nai] * dt;
    STATES[(sample_id * Tomek_num_of_states) + nass] = STATES[(sample_id * Tomek_num_of_states) + nass] + RATES[(sample_id * Tomek_num_of_rates) + nass] * dt;
    STATES[(sample_id * Tomek_num_of_states) + ki] = STATES[(sample_id * Tomek_num_of_states) + ki] + RATES[(sample_id * Tomek_num_of_rates) + ki] * dt;
    STATES[(sample_id * Tomek_num_of_states) + kss] = STATES[(sample_id * Tomek_num_of_states) + kss] + RATES[(sample_id * Tomek_num_of_rates) + kss] * dt;
    STATES[(sample_id * Tomek_num_of_states) + cansr] = STATES[(sample_id * Tomek_num_of_states) + cansr] + RATES[(sample_id * Tomek_num_of_rates) + cansr] * dt;
    STATES[(sample_id * Tomek_num_of_states) + cajsr] = STATES[(sample_id * Tomek_num_of_states) + cajsr] + RATES[(sample_id * Tomek_num_of_rates) + cajsr] * dt;
    STATES[(sample_id * Tomek_num_of_states) + cai] = STATES[(sample_id * Tomek_num_of_states) + cai] + RATES[(sample_id * Tomek_num_of_rates) + cai] * dt;
    STATES[(sample_id * Tomek_num_of_states) + m] = STATES[(sample_id * Tomek_num_of_states) + m] + RATES[(sample_id * Tomek_num_of_rates) + m] * dt;
    STATES[(sample_id * Tomek_num_of_states) + h] = STATES[(sample_id * Tomek_num_of_states) + h] + RATES[(sample_id * Tomek_num_of_rates) + h] * dt;
    STATES[(sample_id * Tomek_num_of_states) + j] = STATES[(sample_id * Tomek_num_of_states) + j] + RATES[(sample_id * Tomek_num_of_rates) + j] * dt;
    STATES[(sample_id * Tomek_num_of_states) + hp] = STATES[(sample_id * Tomek_num_of_states) + hp] + RATES[(sample_id * Tomek_num_of_rates) + hp] * dt;
    STATES[(sample_id * Tomek_num_of_states) + jp] = STATES[(sample_id * Tomek_num_of_states) + jp] + RATES[(sample_id * Tomek_num_of_rates) + jp] * dt;
    STATES[(sample_id * Tomek_num_of_states) + mL] = STATES[(sample_id * Tomek_num_of_states) + mL] + RATES[(sample_id * Tomek_num_of_rates) + mL] * dt;
    STATES[(sample_id * Tomek_num_of_states) + hL] = STATES[(sample_id * Tomek_num_of_states) + hL] + RATES[(sample_id * Tomek_num_of_rates) + hL] * dt;
    STATES[(sample_id * Tomek_num_of_states) + hLp] = STATES[(sample_id * Tomek_num_of_states) + hLp] + RATES[(sample_id * Tomek_num_of_rates) + hLp] * dt;
    STATES[(sample_id * Tomek_num_of_states) + a] = STATES[(sample_id * Tomek_num_of_states) + a] + RATES[(sample_id * Tomek_num_of_rates) + a] * dt;
    STATES[(sample_id * Tomek_num_of_states) + iF] = STATES[(sample_id * Tomek_num_of_states) + iF] + RATES[(sample_id * Tomek_num_of_rates) + iF] * dt;
    STATES[(sample_id * Tomek_num_of_states) + iS] = STATES[(sample_id * Tomek_num_of_states) + iS] + RATES[(sample_id * Tomek_num_of_rates) + iS] * dt;
    STATES[(sample_id * Tomek_num_of_states) + ap] = STATES[(sample_id * Tomek_num_of_states) + ap] + RATES[(sample_id * Tomek_num_of_rates) + ap] * dt;
    STATES[(sample_id * Tomek_num_of_states) + iFp] = STATES[(sample_id * Tomek_num_of_states) + iFp] + RATES[(sample_id * Tomek_num_of_rates) + iFp] * dt;
    STATES[(sample_id * Tomek_num_of_states) + iSp] = STATES[(sample_id * Tomek_num_of_states) + iSp] + RATES[(sample_id * Tomek_num_of_rates) + iSp] * dt;
    STATES[(sample_id * Tomek_num_of_states) + d] = STATES[(sample_id * Tomek_num_of_states) + d] + RATES[(sample_id * Tomek_num_of_rates) + d] * dt;
    STATES[(sample_id * Tomek_num_of_states) + ff] = STATES[(sample_id * Tomek_num_of_states) + ff] + RATES[(sample_id * Tomek_num_of_rates) + ff] * dt;
    STATES[(sample_id * Tomek_num_of_states) + fs] = STATES[(sample_id * Tomek_num_of_states) + fs] + RATES[(sample_id * Tomek_num_of_rates) + fs] * dt;
    STATES[(sample_id * Tomek_num_of_states) + fcaf] = STATES[(sample_id * Tomek_num_of_states) + fcaf] + RATES[(sample_id * Tomek_num_of_rates) + fcaf] * dt;
    STATES[(sample_id * Tomek_num_of_states) + fcas] = STATES[(sample_id * Tomek_num_of_states) + fcas] + RATES[(sample_id * Tomek_num_of_rates) + fcas] * dt;
    STATES[(sample_id * Tomek_num_of_states) + jca] = STATES[(sample_id * Tomek_num_of_states) + jca] + RATES[(sample_id * Tomek_num_of_rates) + jca] * dt;
    STATES[(sample_id * Tomek_num_of_states) + ffp] = STATES[(sample_id * Tomek_num_of_states) + ffp] + RATES[(sample_id * Tomek_num_of_rates) + ffp] * dt;
    STATES[(sample_id * Tomek_num_of_states) + fcafp] = STATES[(sample_id * Tomek_num_of_states) + fcafp] + RATES[(sample_id * Tomek_num_of_rates) + fcafp] * dt;
    STATES[(sample_id * Tomek_num_of_states) + nca_ss] = STATES[(sample_id * Tomek_num_of_states) + nca_ss] + RATES[(sample_id * Tomek_num_of_rates) + nca_ss] * dt;
    STATES[(sample_id * Tomek_num_of_states) + nca_i] = STATES[(sample_id * Tomek_num_of_states) + nca_i] + RATES[(sample_id * Tomek_num_of_rates) + nca_i] * dt;
    STATES[(sample_id * Tomek_num_of_states) + O] = STATES[(sample_id * Tomek_num_of_states) + O] + RATES[(sample_id * Tomek_num_of_rates) + O] * dt;
    STATES[(sample_id * Tomek_num_of_states) + I] = STATES[(sample_id * Tomek_num_of_states) + I] + RATES[(sample_id * Tomek_num_of_rates) + I] * dt;
        STATES[(sample_id * Tomek_num_of_states) + C3] = STATES[(sample_id * Tomek_num_of_states) + C3] + RATES[(sample_id * Tomek_num_of_rates) + C3] * dt;
        STATES[(sample_id * Tomek_num_of_states) + C2] = STATES[(sample_id * Tomek_num_of_states) + C2] + RATES[(sample_id * Tomek_num_of_rates) + C2] * dt;
        STATES[(sample_id * Tomek_num_of_states) + C1] = STATES[(sample_id * Tomek_num_of_states) + C1] + RATES[(sample_id * Tomek_num_of_rates) + C1] * dt;
    STATES[(sample_id * Tomek_num_of_states) + xs1] = STATES[(sample_id * Tomek_num_of_states) + xs1] + RATES[(sample_id * Tomek_num_of_rates) + xs1] * dt;
    STATES[(sample_id * Tomek_num_of_states) + xs2] = STATES[(sample_id * Tomek_num_of_states) + xs2] + RATES[(sample_id * Tomek_num_of_rates) + xs2] * dt;
    STATES[(sample_id * Tomek_num_of_states) + Jrel_np] = STATES[(sample_id * Tomek_num_of_states) + Jrel_np] + RATES[(sample_id * Tomek_num_of_rates) + Jrel_np] * dt;
    STATES[(sample_id * Tomek_num_of_states) + Jrel_p] = STATES[(sample_id * Tomek_num_of_states) + Jrel_p] + RATES[(sample_id * Tomek_num_of_rates) + Jrel_p] * dt;
    #else
    ////==============
    ////Exact solution
    ////==============
    ////INa
    STATES[(sample_id * Tomek_num_of_states) + m] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + mss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + mss] - STATES[(sample_id * Tomek_num_of_states) + m]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tm]);
    STATES[(sample_id * Tomek_num_of_states) + h] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hss] - STATES[(sample_id * Tomek_num_of_states) + h]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + th]);
    STATES[(sample_id * Tomek_num_of_states) + j] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + jss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + jss] - STATES[(sample_id * Tomek_num_of_states) + j]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tj]);
    STATES[(sample_id * Tomek_num_of_states) + hp] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hssp] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hssp] - STATES[(sample_id * Tomek_num_of_states) + hp]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + th]);
    STATES[(sample_id * Tomek_num_of_states) + jp] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + jss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + jss] - STATES[(sample_id * Tomek_num_of_states) + jp]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tjp]);
    STATES[(sample_id * Tomek_num_of_states) + mL] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + mLss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + mLss] - STATES[(sample_id * Tomek_num_of_states) + mL]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tmL]);
    STATES[(sample_id * Tomek_num_of_states) + hL] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hLss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hLss] - STATES[(sample_id * Tomek_num_of_states) + hL]) * exp(-dt / CONSTANTS[(sample_id * Tomek_num_of_constants) + thL]);
    STATES[(sample_id * Tomek_num_of_states) + hLp] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hLssp] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + hLssp] - STATES[(sample_id * Tomek_num_of_states) + hLp]) * exp(-dt / CONSTANTS[(sample_id * Tomek_num_of_constants) + thLp]);
    ////Ito
    STATES[(sample_id * Tomek_num_of_states) + a] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ass] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ass] - STATES[(sample_id * Tomek_num_of_states) + a]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ta]);
    STATES[(sample_id * Tomek_num_of_states) + iF] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + iss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + iss] - STATES[(sample_id * Tomek_num_of_states) + iF]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiF]);
    STATES[(sample_id * Tomek_num_of_states) + iS] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + iss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + iss] - STATES[(sample_id * Tomek_num_of_states) + iS]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiS]);
    STATES[(sample_id * Tomek_num_of_states) + ap] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + assp] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + assp] - STATES[(sample_id * Tomek_num_of_states) + ap]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + ta]);
    STATES[(sample_id * Tomek_num_of_states) + iFp] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + iss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + iss] - STATES[(sample_id * Tomek_num_of_states) + iFp]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiFp]);
    STATES[(sample_id * Tomek_num_of_states) + iSp] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + iss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + iss] - STATES[(sample_id * Tomek_num_of_states) + iSp]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tiSp]);
    ////ICaL
    STATES[(sample_id * Tomek_num_of_states) + d] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + dss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + dss] - STATES[(sample_id * Tomek_num_of_states) + d]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + td]);
    STATES[(sample_id * Tomek_num_of_states) + ff] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fss] - STATES[(sample_id * Tomek_num_of_states) + ff]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tff]);
    STATES[(sample_id * Tomek_num_of_states) + fs] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fss] - STATES[(sample_id * Tomek_num_of_states) + fs]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfs]);
    STATES[(sample_id * Tomek_num_of_states) + fcaf] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcass] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcass] - STATES[(sample_id * Tomek_num_of_states) + fcaf]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfcaf]);
    STATES[(sample_id * Tomek_num_of_states) + fcas] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcass] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcass] - STATES[(sample_id * Tomek_num_of_states) + fcas]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfcas]);
    STATES[(sample_id * Tomek_num_of_states) + jca] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + jcass] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + jcass] - STATES[(sample_id * Tomek_num_of_states) + jca]) * exp(- dt / CONSTANTS[(sample_id * Tomek_num_of_constants) + tjca]);
    STATES[(sample_id * Tomek_num_of_states) + ffp] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fss] - STATES[(sample_id * Tomek_num_of_states) + ffp]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tffp]);
    STATES[(sample_id * Tomek_num_of_states) + fcafp] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcass] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + fcass] - STATES[(sample_id * Tomek_num_of_states) + fcafp]) * exp(-d / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tfcafp]);
        STATES[(sample_id * Tomek_num_of_states) + nca_i] = STATES[(sample_id * Tomek_num_of_states) + nca_i] + RATES[(sample_id * Tomek_num_of_rates) + nca_i]*dt;
        STATES[(sample_id * Tomek_num_of_states) + nca_ss] = STATES[(sample_id * Tomek_num_of_states) + nca_ss] + RATES[(sample_id * Tomek_num_of_rates) + nca_ss]*dt;
    //  STATES[(sample_id * Tomek_num_of_states) + nca_i] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + anca_i] * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + km2n] -
    //      (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + anca_i] * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + km2n] - STATES[(sample_id * Tomek_num_of_states) + nca_i]) * exp(-ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + km2n] * dt);
    //  STATES[(sample_id * Tomek_num_of_states) + nca_ss] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + anca_ss] * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + km2n] -
    //      (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + anca_ss] * CONSTANTS[(sample_id * Tomek_num_of_constants) + k2n] / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + km2n] - STATES[(sample_id * Tomek_num_of_states) + nca_ss]) * exp(-ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + km2n] * dt);
    ////IKr
    //STATES[(sample_id * Tomek_num_of_states) + O] = STATES[(sample_id * Tomek_num_of_states) + O] + RATES[(sample_id * Tomek_num_of_rates) + O] * dt;
    //STATES[(sample_id * Tomek_num_of_states) + I] = STATES[(sample_id * Tomek_num_of_states) + I] + RATES[(sample_id * Tomek_num_of_rates) + I] * dt;
    //STATES[(sample_id * Tomek_num_of_states) + C3] = STATES[(sample_id * Tomek_num_of_states) + C3] + RATES[(sample_id * Tomek_num_of_rates) + C3] * dt;
    //STATES[(sample_id * Tomek_num_of_states) + C2] = STATES[(sample_id * Tomek_num_of_states) + C2] + RATES[(sample_id * Tomek_num_of_rates) + C2] * dt;
    //STATES[(sample_id * Tomek_num_of_states) + C1] = STATES[(sample_id * Tomek_num_of_states) + C1] + RATES[(sample_id * Tomek_num_of_rates) + C1] * dt;
    double* coeffs = new double[15];
    coeffs[0] = -  (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_2]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_i]);
    coeffs[1] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_i];
    coeffs[2] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_2];
    coeffs[3] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_i];
    coeffs[4] = -  (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_ItoC2]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_i]);
    coeffs[5] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_C2ToI];
    coeffs[6] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_2];
    coeffs[7] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta_ItoC2];
    coeffs[8] = -  (CONSTANTS[(sample_id * Tomek_num_of_constants) + beta_1]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_2]+ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha_C2ToI]);
    coeffs[9] = CONSTANTS[(sample_id * Tomek_num_of_constants) + alpha_1];
    coeffs[10] = CONSTANTS[(sample_id * Tomek_num_of_constants) + beta_1];
    coeffs[11] = -  (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta]+CONSTANTS[(sample_id * Tomek_num_of_constants) + alpha_1]);
    coeffs[12] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha];
    coeffs[13] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + beta];
    coeffs[14] = -  ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + alpha];
    int m = 5;
    double* a = new double[m*m]; // Flattened a
    a[0 * m + 0] = 1.0 - dt * coeffs[0];   a[0 * m + 1] = - dt * coeffs[1];     a[0 * m + 2] = - dt * coeffs[2];     a[0 * m + 3] = 0.0;                      a[0 * m + 4] = 0.0;
    a[1 * m + 0] = - dt * coeffs[3];       a[1 * m + 1] = 1.0 - dt * coeffs[4]; a[1 * m + 2] = - dt * coeffs[5];     a[1 * m + 3] = 0.0;                      a[1 * m + 4] = 0.0;
    a[2 * m + 0] = - dt * coeffs[6];       a[2 * m + 1] = - dt * coeffs[7];     a[2 * m + 2] = 1.0 - dt * coeffs[8]; a[2 * m + 3] = - dt * coeffs[9];         a[2 * m + 4] = 0.0;
    a[3 * m + 0] = 0.0;                    a[3 * m + 1] = 0.0;                  a[3 * m + 2] = - dt * coeffs[10];    a[3 * m + 3] = 1.0 - dt * coeffs[11];    a[3 * m + 4] = - dt * coeffs[12];
    a[4 * m + 0] = 0.0;                    a[4 * m + 1] = 0.0;                  a[4 * m + 2] = 0.0;                  a[4 * m + 3] = - dt * coeffs[13];;       a[4 * m + 4] = 1.0 - dt * coeffs[14];
    double* b = new double[m];
    b[0] = STATES[(sample_id * Tomek_num_of_states) + O];
    b[1] = STATES[(sample_id * Tomek_num_of_states) + I];
    b[2] = STATES[(sample_id * Tomek_num_of_states) + C1];
    b[3] = STATES[(sample_id * Tomek_num_of_states) + C2];
    b[4] = STATES[(sample_id * Tomek_num_of_states) + C3];
    double* x = new double[m];
    for(int i = 0; i < m; i++){
        x[i] = 0.0;
    }
    ___gaussElimination(a,b,x,m);
    STATES[(sample_id * Tomek_num_of_states) + O] = x[0];
    STATES[(sample_id * Tomek_num_of_states) + I] = x[1];
    STATES[(sample_id * Tomek_num_of_states) + C1] = x[2];
    STATES[(sample_id * Tomek_num_of_states) + C2] = x[3];
    STATES[(sample_id * Tomek_num_of_states) + C3] = x[4];
    delete[] coeffs;
    delete[] a;
    delete[] b;
    delete[] x;
    
    ////IKs
    STATES[(sample_id * Tomek_num_of_states) + xs1] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + xs1ss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + xs1ss] - STATES[(sample_id * Tomek_num_of_states) + xs1]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + txs1]);
    STATES[(sample_id * Tomek_num_of_states) + xs2] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + xs2ss] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + xs2ss] - STATES[(sample_id * Tomek_num_of_states) + xs2]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + txs2]);
    ////IK1
    ////RyR receptors
    STATES[(sample_id * Tomek_num_of_states) + Jrel_np] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_inf] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_inf] - STATES[(sample_id * Tomek_num_of_states) + Jrel_np]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_rel]);
    STATES[(sample_id * Tomek_num_of_states) + Jrel_p] = ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_infp] - (ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + Jrel_infp] - STATES[(sample_id * Tomek_num_of_states) + Jrel_p]) * exp(-dt / ALGEBRAIC[(sample_id * Tomek_num_of_algebraic) + tau_relp]);
    ////=============================
    ////Approximated solution (Euler)
    ////=============================
    ////CaMK
    STATES[(sample_id * Tomek_num_of_states) + CaMKt] = STATES[(sample_id * Tomek_num_of_states) + CaMKt] + RATES[(sample_id * Tomek_num_of_rates) + CaMKt] * dt;
    ////Membrane potential
    STATES[(sample_id * Tomek_num_of_states) + V] = STATES[(sample_id * Tomek_num_of_states) + V] + RATES[(sample_id * Tomek_num_of_rates) + V] * dt;
    ////Ion Concentrations and Buffers
    STATES[(sample_id * Tomek_num_of_states) + nai] = STATES[(sample_id * Tomek_num_of_states) + nai] + RATES[(sample_id * Tomek_num_of_rates) + nai] * dt;
    STATES[(sample_id * Tomek_num_of_states) + nass] = STATES[(sample_id * Tomek_num_of_states) + nass] + RATES[(sample_id * Tomek_num_of_rates) + nass] * dt;
    STATES[(sample_id * Tomek_num_of_states) + ki] = STATES[(sample_id * Tomek_num_of_states) + ki] + RATES[(sample_id * Tomek_num_of_rates) + ki] * dt;
    STATES[(sample_id * Tomek_num_of_states) + kss] = STATES[(sample_id * Tomek_num_of_states) + kss] + RATES[(sample_id * Tomek_num_of_rates) + kss] * dt;
    STATES[(sample_id * Tomek_num_of_states) + cai] = STATES[(sample_id * Tomek_num_of_states) + cai] + RATES[(sample_id * Tomek_num_of_rates) + cai] * dt;
    STATES[(sample_id * Tomek_num_of_states) + cass] = STATES[(sample_id * Tomek_num_of_states) + cass] + RATES[(sample_id * Tomek_num_of_rates) + cass] * dt;
    STATES[(sample_id * Tomek_num_of_states) + cansr] = STATES[(sample_id * Tomek_num_of_states) + cansr] + RATES[(sample_id * Tomek_num_of_rates) + cansr] * dt;
    STATES[(sample_id * Tomek_num_of_states) + cajsr] = STATES[(sample_id * Tomek_num_of_states) + cajsr] + RATES[(sample_id * Tomek_num_of_rates) + cajsr] * dt;
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
    //printf("TIME = %E, dV = %E, time_step = %E\n",TIME, RATES[(sample_id * Tomek_num_of_rates) + V] * time_step, time_step);
  }
  else {
    //printf("TIME > time_point ms\n");
    if (std::abs(RATES[(sample_id * Tomek_num_of_rates) + V] * time_step) <= min_dV) {//Slow changes in V
        //printf("dV/dt <= 0.2\n");
        time_step = std::abs(max_dV / RATES[(sample_id * Tomek_num_of_rates)  + V]);
        //Make sure time_step is between min time step and max_time_step
        if (time_step < min_time_step) {
            time_step = min_time_step;
        }
        else if (time_step > max_time_step) {
            time_step = max_time_step;
        }
        //printf("TIME = %E, dV = %E, time_step = %E\n",TIME, RATES[(sample_id * Tomek_num_of_rates) + V] * time_step, time_step);
    }
    else if (std::abs(RATES[(sample_id * Tomek_num_of_rates) + V] * time_step) >= max_dV) {//Fast changes in V
        //printf("dV/dt >= 0.8\n");
        time_step = std::abs(min_dV / RATES[(sample_id * Tomek_num_of_rates) + V]);
        //Make sure time_step is not less than 0.005
        if (time_step < min_time_step) {
            time_step = min_time_step;
        }
        //printf("TIME = %E, dV = %E, time_step = %E\n",TIME, RATES[(sample_id * Tomek_num_of_rates) + V] * time_step, time_step);
    } else {
        time_step = min_time_step;
    }
    return time_step;
  }
}