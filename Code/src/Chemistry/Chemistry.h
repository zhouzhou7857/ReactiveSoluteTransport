/*
 * Chemistry.h
 *
 *  Created on: 2025/12/11
 *      Author: delphine and Wenyu
 */

#ifndef CHEMISTRY_H_
#define CHEMISTRY_H_

#include <algorithm>

// Default chemistry parameters used when no external chemistry input file is
// provided. The active values below are runtime-configurable and may be
// overwritten from Input/Chemistry_files/*.txt.
const double DEFAULT_INITIAL_REACTIVE_CONCENTRATION = 500.0;
const double DEFAULT_REACTIVE_CONCENTRATION_DECAY = 5.0E-3;
const double DEFAULT_REACTIVE_TO_MINERAL_STOICH = 1.0;
const double DEFAULT_MINERAL_MOLAR_VOLUME = 7.4E-5;
const double DEFAULT_FRACTURE_OUT_OF_PLANE_THICKNESS = 1.0;

// Assumed inlet reactive concentration carried by each particle [mol/m^3].
// This codebase does not currently read chemistry parameters from input files, so
// the concentration evolution is parameterized here to keep chemistry modular.
// Assumption: the reactive concentration represents a finite reagent inventory
// carried by each particle. As it is consumed, mineral dissolves and aperture
// opens. A precipitation-oriented closure is not parameterized in the current
// input model, so the default active chemistry is reagent-limited dissolution.
extern double INITIAL_REACTIVE_CONCENTRATION;
// First-order depletion rate of the particle-carried reactive concentration [1/s].
extern double REACTIVE_CONCENTRATION_DECAY;
// Stoichiometric conversion from consumed reactive moles to dissolved/precipitated
// mineral moles [mol mineral / mol reactive species].
extern double REACTIVE_TO_MINERAL_STOICH;
// Mineral molar volume used to convert reacted mineral moles to mineral volume
// change [m^3/mol]. The default value is an assumed gypsum-scale molar volume. （1 unit of H2O convert 1 unit of Ca and SO4）
extern double MINERAL_MOLAR_VOLUME;
// Unit out-of-plane thickness for the 2D DFN representation [m] (1 m for simplifying).
extern double FRACTURE_OUT_OF_PLANE_THICKNESS;

// Legacy empirical particle-count coupling parameter.
// 在 DFN-PT-V3 当前主路径中并未使用，仅保留给旧的 aperture 更新接口。
const double K_REACTION = -1.5E-10; // Reaction factor calculated from trend line slope of Delta_V [MODIFIED!!][exp = 1.5e-10 m/s]
// Gypsum dissolving rate is about 6e-9 m/s， dividing by porosity correction would be 6.7E-9 m/s
/**
 * Legacy empirical aperture update interface.
 * 在 DFN-PT-V3 当前主路径中并未使用，仅保留给旧版本/兼容性代码。
 * Compute aperture b from the empirical formula:
 *   b_new = b - K_REACTION * t * (Nb / Nt)
 * with the ratio (Nb / Nt) capped at 1.
 * @param b  current aperture
 * @param t   reaction update time step (s)
 * @param Nb  effective particle count through the segment during the current reaction time step
 * @param Nt  reference injected particle count for the current reaction time step
 * @param clamp_nonnegative if true, b is clamped to >= 0 (Delta_b can be negative)
 * @return b_new  aperture (new aperture)
 */
double ComputeAperture(double b, double t, double Nb, double Nt, bool clamp_nonnegative = true);

/**
 * Legacy empirical aperture-delta interface.
 * 在 DFN-PT-V3 当前主路径中并未使用，仅保留给旧版本/兼容性代码。
 * Compute delta aperture for updating an existing aperture:
 *   delta_b = b_new - b_old
 *
 * @param b current aperture
 * @param t     reaction update time step
 * @param Nb    effective particles through the segment during the current reaction time step
 * @param Nt    reference injected particle count for the current reaction time step
 * @param clamp_nonnegative if true, b_new is clamped to >= 0 before differencing
 * @return delta_b
 */
double ComputeDeltaAperture(double b, double t, double Nb, double Nt, bool clamp_nonnegative = true);
// Legacy helper used only by the old particle-count aperture update path.
// 在 DFN-PT-V3 当前主路径中并未使用。
double safe_ratio(double Nb, double Nt);

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
		double stoich,double molar_volume,double fracture_thickness);

#endif /* CHEMISTRY_H_ */
