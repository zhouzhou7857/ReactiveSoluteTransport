/*
 * Chemistry.h
 *
 *  Created on: 2025/12/11
 *      Author: delphine and Wenyu
 */

#ifndef CHEMISTRY_H_
#define CHEMISTRY_H_

#include <algorithm>

// Active chemistry core:
// no external chemistry input file is read by the PT executable at runtime.
// The benchmark script converts PHREEQC output into a fixed cumulative
// mineral-volume law F_ref(t). Each particle contribution is scaled as:
// DeltaV_particle = DeltaV_ref * V_particle / Vref.
// For the current gypsum precipitation benchmark, Vref is the total injected
// supersaturated water volume and the active mode is rational_precipitation.
// Legacy concentration-based values are retained only for inactive helpers.
// Legacy default reference law:
// F_ref(t) =
//   A1 * (1 - exp(-k1 * t)) +
//   A2 * (1 - exp(-k2 * t)) +
//   L  * t
// Here A1, A2, and L are stored in m^3 units.
// Converted from PHREEQC-fit coefficients in cm^3 using 1 cm^3 = 1.0E-6 m^3.
// Positive DeltaV denotes dissolution; negative DeltaV denotes precipitation.
// The active precipitation benchmark overrides this legacy law with:
// F_ref(t) = Pinf * [1 - w*(1+a10*t)^(-1/9) - (1-w)*(1+a2*t)^(-1)].
const double DEFAULT_DELTA_V_FAST1_AMPLITUDE = 2.02192E-5;
const double DEFAULT_DELTA_V_FAST1_RATE = 2.14128E-3;
const double DEFAULT_DELTA_V_FAST2_AMPLITUDE = 2.50122E-6;
const double DEFAULT_DELTA_V_FAST2_RATE = 2.72730E-5;
const double DEFAULT_DELTA_V_SLOW_LINEAR_RATE = 2.48263E-16;
const double DEFAULT_DELTA_V_REFERENCE_WATER_VOLUME = 1.0E-3;
// change to 5 independent variables in volume function, with 5 parameters to fit (A1,k1,A2,k2,L) - DR on 2025/12/11
const double DEFAULT_MINIMUM_FRACTURE_APERTURE = 1.0E-12;
const double DEFAULT_INITIAL_REACTIVE_CONCENTRATION = 500.0; // not used in active DeltaV logic
const double DEFAULT_REACTIVE_CONCENTRATION_DECAY = 5.0E-3; // not used in active DeltaV logic
const double DEFAULT_REACTIVE_TO_MINERAL_STOICH = 1.0; // not used in active DeltaV logic
const double DEFAULT_MINERAL_MOLAR_VOLUME = 7.4E-5; // not used in active DeltaV logic
const double DEFAULT_FRACTURE_OUT_OF_PLANE_THICKNESS = 1.0;

// Particle-based cumulative mineral volume law parameters.
// DELTA_V_* store the legacy two-exponential law; RATIONAL_PRECIP_* store the
// active supersaturated gypsum precipitation benchmark law.
extern double DELTA_V_FAST1_AMPLITUDE;
extern double DELTA_V_FAST1_RATE;
extern double DELTA_V_FAST2_AMPLITUDE;
extern double DELTA_V_FAST2_RATE;
extern double DELTA_V_SLOW_LINEAR_RATE;
extern double DELTA_V_REFERENCE_WATER_VOLUME;
extern double RATIONAL_PRECIP_PINF;
extern double RATIONAL_PRECIP_WEIGHT;
extern double RATIONAL_PRECIP_A10_RATE;
extern double RATIONAL_PRECIP_A2_RATE;
extern double EFFECTIVE_DIFFUSION_COEFFICIENT;
extern double EFFECTIVE_DIFFUSION_TIME;
extern double MINIMUM_FRACTURE_APERTURE;
extern bool USE_RATIONAL_PRECIPITATION_LAW;
extern bool USE_VP_WIDTH_CORRECTION;
extern bool USE_EFFECTIVE_DIFFUSION_HEIGHT_FACTOR;
extern bool CHEMISTRY_DEBUG_LOGGING;

// Deprecated / inactive legacy chemistry parameters.
extern double INITIAL_REACTIVE_CONCENTRATION; // not used in active DeltaV logic
extern double REACTIVE_CONCENTRATION_DECAY; // not used in active DeltaV logic
extern double REACTIVE_TO_MINERAL_STOICH; // not used in active DeltaV logic
extern double MINERAL_MOLAR_VOLUME; // not used in active DeltaV logic
extern double FRACTURE_OUT_OF_PLANE_THICKNESS;

// Deprecated / inactive empirical particle-count coupling parameter.
const double K_REACTION = -1.5E-10;

double ComputeAperture(double b, double t, double Nb, double Nt, bool clamp_nonnegative = true);
double ComputeDeltaAperture(double b, double t, double Nb, double Nt, bool clamp_nonnegative = true);
double safe_ratio(double Nb, double Nt);

// Deprecated / inactive result type retained for compatibility with the legacy
// concentration-based coupling helpers.
class ChemistryStepResult{
public:
	double concentration_out;
	double reacted_reactive_moles;
	double mineral_moles_changed;
	double mineral_volume_change;
	double local_aperture_change;
	double segment_aperture_change;
	double traversed_fraction;
	bool geometry_active;
public:
	ChemistryStepResult():
		concentration_out(0.0),
		reacted_reactive_moles(0.0),
		mineral_moles_changed(0.0),
		mineral_volume_change(0.0),
		local_aperture_change(0.0),
		segment_aperture_change(0.0),
		traversed_fraction(0.0),
		geometry_active(false){}
};

double UpdateReactiveConcentration(double concentration_in,double residence_time);
ChemistryStepResult EvaluateReactiveStep(double concentration_in,double residence_time,double particle_representative_volume,double segment_length,double traversed_length);
void ResetChemistryParametersToDefaults();
void ConfigureChemistryParameters(double initial_concentration,double concentration_decay,
		double stoich,double molar_volume,double fracture_thickness,
		double minimum_fracture_aperture,bool chemistry_debug_logging);
// Set the legacy cumulative PHREEQC-derived reference law:
// F_ref(t) = A1 * (1 - exp(-k1 * t))
//          + A2 * (1 - exp(-k2 * t))
//          + L * t
void ConfigureMineralVolumeLaw(double fast1_amplitude,double fast1_rate,
		double fast2_amplitude,double fast2_rate,double slow_linear_rate);
// Set the active cumulative rational precipitation reference law:
// F_ref(t) = Pinf * [1 - w * (1 + a10*t)^(-1/9)
//                     - (1-w) * (1 + a2*t)^(-1)]
// Pinf is signed; use negative Pinf for precipitation.
void ConfigureRationalPrecipitationLaw(double pinf,double weight,
		double a10_rate,double a2_rate);
// Set the PHREEQC reference water volume Vref used in the scaling factor:
// scale = V_particle / Vref
void ConfigureMineralVolumeScaling(double reference_water_volume);
// Set the out-of-plane fracture thickness used in the volume-to-aperture map:
// delta_b = delta_V / (segment_length * thickness * 2)
void ConfigureFractureOutOfPlaneThickness(double fracture_thickness);
void ConfigureVpWidthCorrection(bool enable_width_correction);
void ConfigureEffectiveDiffusionHeightFactor(bool enable_factor,double diffusion_coefficient,double diffusion_time);

// Return the fixed PHREEQC reference cumulative mineral-volume law.
// If USE_RATIONAL_PRECIPITATION_LAW is true, the rational precipitation law is
// used; otherwise the legacy two-exponential-plus-linear law is used:
// F_ref(t) = A1 * (1 - exp(-k1 * t))
//          + A2 * (1 - exp(-k2 * t))
//          + L * t
double GetCumulativeMineralVolumeChange(double t_particle);
// Return the particle-volume scaling factor:
// scale = V_particle / Vref by default. Optional correction flags can multiply
// this default by the effective diffusion height factor and/or width correction.
double ComputeParticleVolumeScalingFactor(double particle_representative_volume,double aperture_height);
// Return the mineral-volume increment over one particle time interval:
// DeltaV_particle[t_start,t_end]
//   = (DeltaV_ref(t_end) - DeltaV_ref(t_start)) * (V_particle / Vref)
double ComputeParticleSegmentMineralVolumeChange(double t_particle_start,double t_particle_end,
		double particle_representative_volume,double aperture_height);
double ComputeApertureChangeFromMineralVolume(double mineral_volume_change,double segment_length);
double ClampApertureToMinimumThreshold(double aperture_value);

#endif /* CHEMISTRY_H_ */
