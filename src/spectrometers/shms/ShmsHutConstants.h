// src/spectrometers/shms/ShmsHutConstants.h (private)
// Fortran: shms/hut.inc — PARAMETERs of the SHMS detector hut (included by mc_shms.f and mc_shms_hut.f),
// copied verbatim. real*8 under -fdefault-real-8; parameter expressions folded in double (constexpr here).
// Positions are relative to the focal plane (z = 0), cm.
#ifndef SIMC_SRC_SPECTROMETERS_SHMS_HUTCONSTANTS_H
#define SIMC_SRC_SPECTROMETERS_SHMS_HUTCONSTANTS_H

namespace simc {
namespace shms {
namespace hut {

// spectrometer exit window
constexpr double hfoil_exit_radlen = 8.89;
constexpr double hfoil_exit_thick = 0.020 * 2.54;
// air and helium bag
constexpr double hair_radlen = 30420.;
constexpr double helbag_hel_radlen = 754560.;
constexpr double helbag_mylar_radlen = 28.53;
constexpr double helbag_mylar_thick = .010 * 2.54;
constexpr double helbag_al_radlen = 8.89;
constexpr double helbag_al_thick = 3e-4;
// drift chambers
constexpr double hdc_entr_radlen = 28.7;
constexpr double hdc_entr_thick = 0.001 * 2.54;
constexpr double hdc_radlen = 16700.0;
constexpr double hdc_thick = 0.125 * 2.54;
constexpr double hdc_wire_radlen = 0.35;
constexpr double hdc_wire_thick = 0.0000354;
constexpr double hdc_cath_radlen = 28.6;
constexpr double hdc_cath_thick = 0.001 * 2.54;
constexpr double hdc_exit_radlen = 28.7;
constexpr double hdc_exit_thick = 0.001 * 2.54;
// hodoscopes
constexpr double hscin_radlen = 42.4;
// heavy-gas Cherenkov (1)
constexpr double hcer_entr_radlen = 19.63;
constexpr double hcer_entr_thick = 0.002 * 2.54;
constexpr double hcer_1_radlen = 11700.0;
constexpr double hcer_mirglass_radlen = 12.29;
constexpr double hcer_mirglass_thick = 0.3;
constexpr double hcer_exit_radlen = 19.63;
constexpr double hcer_exit_thick = 0.002 * 2.54;
// noble-gas Cherenkov (2)
constexpr double hcer_2_entr_radlen = 8.90;
constexpr double hcer_2_entr_thick = 0.040 * 2.54;
constexpr double hcer_2_radlen = 1202.5;
constexpr double hcer_mir_radlen = 400.;
constexpr double hcer_mir_thick = 2.00;
constexpr double hcer_2_exit_radlen = 8.90;
constexpr double hcer_2_exit_thick = 0.040 * 2.54;
// Wire chamber resolutions: real*8 hdc_sigma(1:12) / 12*0.020 /
constexpr double hdc_sigma[12] = {0.020, 0.020, 0.020, 0.020, 0.020, 0.020,
                                  0.020, 0.020, 0.020, 0.020, 0.020, 0.020};
constexpr int hdc_nr_cham = 2;
constexpr int hdc_nr_plan = 6;
constexpr double hdc_1_zpos = -40.656;
constexpr double hdc_2_zpos = 39.332;
constexpr double hdc_1_left = 40.0, hdc_1_right = -40.0, hdc_1y_offset = 0.0;
constexpr double hdc_1_top = -40., hdc_1_bot = 40., hdc_1x_offset = 0.0;
constexpr double hdc_2_left = 40.0, hdc_2_right = -40.0, hdc_2y_offset = 0.0;
constexpr double hdc_2_top = -40., hdc_2_bot = 40., hdc_2x_offset = 0.;
// hodoscope positions / sizes
constexpr double hscin_1x_zpos = 52.1;
constexpr double hscin_1y_zpos = 61.7;
constexpr double hscin_2x_zpos = 271.4;
constexpr double hscin_2y_zpos = 282.4;
constexpr double hscin_1x_thick = 1.000 * 1.067;
constexpr double hscin_1y_thick = 1.000 * 1.067;
constexpr double hscin_2x_thick = 1.000 * 1.067;
constexpr double hscin_2y_thick = 1.000 * 1.067;
constexpr double hscin_1x_left = 50., hscin_1x_right = -50., hscin_1x_offset = 0.0;
constexpr double hscin_1y_top = -45., hscin_1y_bot = 45., hscin_1y_offset = 0.0;
constexpr double hscin_2x_left = 55., hscin_2x_right = -55., hscin_2x_offset = 0.;
constexpr double hscin_2y_top = -62.5, hscin_2y_bot = 62.5, hscin_2y_offset = 0;
// Cherenkov positions
constexpr double hcer_1_zentrance = -291.700;
constexpr double hcer_1_zmirror = -84.900;
constexpr double hcer_1_zexit = -61.700;
constexpr double hcer_2_zentrance = 72.600;
constexpr double hcer_2_zmirror = 179.400;
constexpr double hcer_2_zexit = 202.600;
// calorimeter
constexpr double hcal_4ta_zpos = 341.0;
constexpr double hcal_left = 63.00, hcal_right = -63.00, hcal_top = -70.00, hcal_bottom = 70.00;

}  // namespace hut
}  // namespace shms
}  // namespace simc

#endif
