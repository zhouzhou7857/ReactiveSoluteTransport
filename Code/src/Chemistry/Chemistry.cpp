/*
 * Chemistry.cpp
 *
 *  Created on: 2025/12/11
 *      Author: delphine and Wenyu
 */

#include "Chemistry.h"
#include "../Utilitaries/Constantes.h"
#include <cmath>

using namespace std;

double INITIAL_REACTIVE_CONCENTRATION = DEFAULT_INITIAL_REACTIVE_CONCENTRATION;
double REACTIVE_CONCENTRATION_DECAY = DEFAULT_REACTIVE_CONCENTRATION_DECAY;
double REACTIVE_TO_MINERAL_STOICH = DEFAULT_REACTIVE_TO_MINERAL_STOICH;
double MINERAL_MOLAR_VOLUME = DEFAULT_MINERAL_MOLAR_VOLUME;
double FRACTURE_OUT_OF_PLANE_THICKNESS = DEFAULT_FRACTURE_OUT_OF_PLANE_THICKNESS;
double DELTA_V_FAST1_AMPLITUDE = DEFAULT_DELTA_V_FAST1_AMPLITUDE;
double DELTA_V_FAST1_RATE = DEFAULT_DELTA_V_FAST1_RATE;
double DELTA_V_FAST2_AMPLITUDE = DEFAULT_DELTA_V_FAST2_AMPLITUDE;
double DELTA_V_FAST2_RATE = DEFAULT_DELTA_V_FAST2_RATE;
double DELTA_V_SLOW_LINEAR_RATE = DEFAULT_DELTA_V_SLOW_LINEAR_RATE;
double DELTA_V_REFERENCE_WATER_VOLUME = DEFAULT_DELTA_V_REFERENCE_WATER_VOLUME;
double MINIMUM_FRACTURE_APERTURE = DEFAULT_MINIMUM_FRACTURE_APERTURE;
bool CHEMISTRY_DEBUG_LOGGING = false;

// DFN-PT-V4 chemistry core:
// the active model uses particle age and a cumulative mineral volume change
// closure instead of the older concentration-to-geometry input-file workflow.

// added by DR on 2025/12/11 - function that returns the change in aperture for a given time spent in the fracture mesh
/*double ChemicalFunction(double fract_time,int nb_part){
  double delta_b=0;
  return delta_b;
}*/
// added by Wenyu on 2025/12/30 - function that updates aperture over one reaction time step
#include <algorithm> // std::max
// Avoid negative aperture
//namespace chemistry {

/*static inline double safe_ratio(int Nb, int Nt)
{
    if (Nt <= 0) return 0.0;
    if (Nb <= 0) return 0.0;
    return static_cast<double>(Nb) / static_cast<double>(Nt);
}*/
// Legacy helper for the old empirical particle-count coupling.
// 在 DFN-PT-V3 当前主路径中并未使用。
double safe_ratio(double Nb, double Nt)
{
    if (Nt <= 0) return 0.0;
    if (Nb <= 0) return 0.0;
    return std::min(1.0, (double)Nb / (double)Nt);
}
// Legacy empirical aperture update using the effective segment throughput Nb
// normalized by a reference injected particle count Nt.
// 在 DFN-PT-V3 当前主路径中并未使用。
double ComputeAperture(double b_old, double t, double Nb, double Nt, bool clamp_nonnegative)
{
    if (t < 0.0) t = 0.0;
    if (Nt <= 0) return b_old;
    if (Nb < 0) Nb = 0;

    const double ratio = safe_ratio(Nb, Nt);
    const double delta = K_REACTION * t * ratio;     // k * t * (Nb/Nt)
    double b_new = b_old - delta;

    if (clamp_nonnegative) b_new = std::max(0.0, b_new);
    if (b_new < 1e-8) b_new = 0.0;
    return b_new;
}
// Legacy aperture-delta helper.
// 在 DFN-PT-V3 当前主路径中并未使用。
double ComputeDeltaAperture(double b, double t, double Nb, double Nt, bool clamp_nonnegative)
{
    const double b_new = ComputeAperture(b, t, Nb, Nt, clamp_nonnegative);
    return b_new - b;
}

void ResetChemistryParametersToDefaults()
{
	INITIAL_REACTIVE_CONCENTRATION = DEFAULT_INITIAL_REACTIVE_CONCENTRATION;
	REACTIVE_CONCENTRATION_DECAY = DEFAULT_REACTIVE_CONCENTRATION_DECAY;
	REACTIVE_TO_MINERAL_STOICH = DEFAULT_REACTIVE_TO_MINERAL_STOICH;
	MINERAL_MOLAR_VOLUME = DEFAULT_MINERAL_MOLAR_VOLUME;
	FRACTURE_OUT_OF_PLANE_THICKNESS = DEFAULT_FRACTURE_OUT_OF_PLANE_THICKNESS;
	DELTA_V_FAST1_AMPLITUDE = DEFAULT_DELTA_V_FAST1_AMPLITUDE;
	DELTA_V_FAST1_RATE = DEFAULT_DELTA_V_FAST1_RATE;
	DELTA_V_FAST2_AMPLITUDE = DEFAULT_DELTA_V_FAST2_AMPLITUDE;
	DELTA_V_FAST2_RATE = DEFAULT_DELTA_V_FAST2_RATE;
	DELTA_V_SLOW_LINEAR_RATE = DEFAULT_DELTA_V_SLOW_LINEAR_RATE;
	DELTA_V_REFERENCE_WATER_VOLUME = DEFAULT_DELTA_V_REFERENCE_WATER_VOLUME;
	MINIMUM_FRACTURE_APERTURE = DEFAULT_MINIMUM_FRACTURE_APERTURE;
	CHEMISTRY_DEBUG_LOGGING = false;
}

void ConfigureChemistryParameters(double initial_concentration,double concentration_decay,
		double stoich,double molar_volume,double fracture_thickness,
		double minimum_fracture_aperture,bool chemistry_debug_logging)
{
	INITIAL_REACTIVE_CONCENTRATION = initial_concentration;
	REACTIVE_CONCENTRATION_DECAY = concentration_decay;
	REACTIVE_TO_MINERAL_STOICH = stoich;
	MINERAL_MOLAR_VOLUME = molar_volume;
	FRACTURE_OUT_OF_PLANE_THICKNESS = fracture_thickness;
	MINIMUM_FRACTURE_APERTURE = minimum_fracture_aperture;
	CHEMISTRY_DEBUG_LOGGING = chemistry_debug_logging;
}

void ConfigureMineralVolumeLaw(double fast1_amplitude,double fast1_rate,
		double fast2_amplitude,double fast2_rate,double slow_linear_rate)
{
	DELTA_V_FAST1_AMPLITUDE = fast1_amplitude;
	DELTA_V_FAST1_RATE = fast1_rate;
	DELTA_V_FAST2_AMPLITUDE = fast2_amplitude;
	DELTA_V_FAST2_RATE = fast2_rate;
	DELTA_V_SLOW_LINEAR_RATE = slow_linear_rate;
}

// Set the PHREEQC calibration basis water volume Vref in:
// scale = V_particle / Vref
void ConfigureMineralVolumeScaling(double reference_water_volume)
{
	if (reference_water_volume>0.0){
		DELTA_V_REFERENCE_WATER_VOLUME = reference_water_volume;
	}
}

void ConfigureFractureOutOfPlaneThickness(double fracture_thickness)
{
	if (fracture_thickness>0.0){
		FRACTURE_OUT_OF_PLANE_THICKNESS = fracture_thickness;
	}
}

// DFN-PT-V4 chemistry core functions below. Modified by WZ on 10/04/2026
// The active PHREEQC-derived law is expressed as a cumulative mineral-volume
// change for one particle parcel:
// DeltaV(t) =
//   A1 * (1 - exp(-k1 * t)) +
//   A2 * (1 - exp(-k2 * t)) +
//   L  * t
// where the current default parameters are:
//   A1 = 20.22E-6 m^3, k1 = 2.141E-3
//   A2 = 2.50E-6 m^3,  k2 = 2.728E-5
//   L  = 2.43E-16 m^3/s
// with positive values denoting dissolution and negative values denoting
// precipitation. The linear term includes both dolomite and calcite contributions + residual effects.
double GetCumulativeMineralVolumeChange(double t_particle)
{
	if (t_particle<=0.0){return 0.0;}
	double t = std::max(0.0,t_particle);
	double fast1_amplitude = DELTA_V_FAST1_AMPLITUDE;
	double fast1_rate = std::max(0.0,DELTA_V_FAST1_RATE);
	double fast2_amplitude = DELTA_V_FAST2_AMPLITUDE;
	double fast2_rate = std::max(0.0,DELTA_V_FAST2_RATE);
	double slow_linear_rate = DELTA_V_SLOW_LINEAR_RATE;
	double fast1_term = fast1_amplitude*(1.0-std::exp(-fast1_rate*t));
	double fast2_term = fast2_amplitude*(1.0-std::exp(-fast2_rate*t));
	double slow_term = slow_linear_rate*t;
	return fast1_term + fast2_term + slow_term;
}

// Particle-specific volume scaling relative to the PHREEQC calibration water:
// scale = V_particle / Vref ; Vref is set by phreeqc model, here we used 1000 cm3 as the default value for Vref, which corresponds to 1 liter of water.
double ComputeParticleVolumeScalingFactor(double particle_representative_volume)
{
	if (particle_representative_volume<=0.0){return 0.0;}
	double reference_volume = DELTA_V_REFERENCE_WATER_VOLUME;
	if (reference_volume<=0.0){return 0.0;}
	return particle_representative_volume/reference_volume;
}

double ComputeParticleSegmentMineralVolumeChange(double t_particle_start,double t_particle_end,
		double particle_representative_volume)
{
	// DFN-PT-V4 segment chemistry contribution:
	// DeltaV_particle[t_start,t_end]
	//   = (DeltaV_ref(t_end) - DeltaV_ref(t_start)) * (V_particle / Vref)
	double t_start = std::max(0.0,t_particle_start);
	double t_end = std::max(t_start,t_particle_end);
	double cumulative_delta_v =
		GetCumulativeMineralVolumeChange(t_end)-GetCumulativeMineralVolumeChange(t_start);
	return cumulative_delta_v*ComputeParticleVolumeScalingFactor(particle_representative_volume);
}

double ComputeApertureChangeFromMineralVolume(double mineral_volume_change,double segment_length)
{
	// DFN-PT-V4 geometry mapping for a 2D segment with unit thickness and
	// symmetric wall dissolution/precipitation.
	if (segment_length<=0.0){return 0.0;}
	return mineral_volume_change/(segment_length*FRACTURE_OUT_OF_PLANE_THICKNESS*2.0);
}

double ClampApertureToMinimumThreshold(double aperture_value) // Modified by WZ on 10/04/2026
{
	// Numerical safeguard for aperture updates:
	// 1. protect the model from erroneous negative aperture values caused by
	//    overshoot or accumulated update errors;
	// 2. enforce a minimum positive threshold, so aperture cannot smoothly
	//    approach zero through a sequence of tiny positive values. With
	//    MINIMUM_FRACTURE_APERTURE > 0, the aperture either stays above the
	//    threshold or is clamped directly to zero when the update crosses zero.
	double clamped_aperture = std::max(0.0,aperture_value);
	if (clamped_aperture>0.0 && clamped_aperture<MINIMUM_FRACTURE_APERTURE){
		clamped_aperture = MINIMUM_FRACTURE_APERTURE;
	}
	return clamped_aperture;
}
// Not used in V4
double UpdateReactiveConcentration(double concentration_in,double residence_time)
{
	if (concentration_in<=0.0){return 0.0;}
	if (residence_time<=0.0){return concentration_in;}
	double concentration_out = concentration_in*std::exp(-REACTIVE_CONCENTRATION_DECAY*residence_time);
	if (concentration_out<0.0){concentration_out = 0.0;}
	if (concentration_out<1.0E-20){concentration_out = 0.0;}
	return concentration_out;
}
// Not used in V4
ChemistryStepResult EvaluateReactiveStep(double concentration_in,double residence_time,double particle_representative_volume,double segment_length,double traversed_length)
{
	ChemistryStepResult result;
	if (concentration_in<=0.0){
		// Zero-concentration particles remain transport-active if the caller keeps
		// moving them, but they no longer modify geometry.
		result.concentration_out = 0.0;
		result.geometry_active = false;
		return result;
	}
	if (residence_time<=0.0 || segment_length<=0.0 || traversed_length<=0.0){
		// Assumption: if no finite segment fraction is traversed during the update
		// interval, neither concentration nor geometry is changed.
		result.concentration_out = concentration_in;
		result.geometry_active = true;
		return result;
	}
	result.concentration_out = UpdateReactiveConcentration(concentration_in,residence_time);
	result.geometry_active = (result.concentration_out>0.0);
	// Mass conservation is enforced by converting the loss of particle-carried
	// reactive concentration into reacted moles over the dynamic particle-carried
	// fluid volume assigned at injection from the instantaneous inlet flux. The
	// same reacted moles are then mapped to mineral volume change through the
	// stoichiometric ratio and mineral molar volume.
	double delta_concentration = concentration_in-result.concentration_out;
	if (delta_concentration<0.0){delta_concentration = 0.0;}
	if (particle_representative_volume<0.0){particle_representative_volume = 0.0;}
	result.reacted_reactive_moles = delta_concentration*particle_representative_volume;
	result.mineral_moles_changed = result.reacted_reactive_moles*REACTIVE_TO_MINERAL_STOICH;
	result.mineral_volume_change = result.mineral_moles_changed*MINERAL_MOLAR_VOLUME;
	double capped_traversed_length = std::min(std::max(traversed_length,0.0),segment_length);
	result.traversed_fraction = capped_traversed_length/segment_length;
	double affected_length = std::max(capped_traversed_length,EPSILON);
	// local_aperture_change is the aperture increment over only the traversed
	// fraction of the fracture. segment_aperture_change is the equivalent
	// segment-averaged increment applied by this code, because the current data
	// model stores one aperture value per segment rather than a sub-segment field.
	result.local_aperture_change = result.mineral_volume_change/(affected_length*FRACTURE_OUT_OF_PLANE_THICKNESS);
	result.segment_aperture_change = result.traversed_fraction*result.local_aperture_change;
	return result;
}
