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

double UpdateReactiveConcentration(double concentration_in,double residence_time)
{
	if (concentration_in<=0.0){return 0.0;}
	if (residence_time<=0.0){return concentration_in;}
	double concentration_out = concentration_in*std::exp(-REACTIVE_CONCENTRATION_DECAY*residence_time);
	if (concentration_out<0.0){concentration_out = 0.0;}
	if (concentration_out<1.0E-20){concentration_out = 0.0;}
	return concentration_out;
}

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
