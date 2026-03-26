/*
 * Chemistry.h
 *
 *  Created on: 2025/12/11
 *      Author: delphine and Wenyu
 */

#ifndef CHEMISTRY_H_
#define CHEMISTRY_H_

//namespace chemistry {

const double K_REACTION = -1.5E-10; // Reaction factor calculated from trend line slope of Delta_V [MODIFIED!!][exp = 1.5e-10 m/s]
// Gypsum dissolving rate is about 6e-9 m/s， dividing by porosity correction would be 6.7E-9 m/s
/**
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
double safe_ratio(double Nb, double Nt);

//}

#endif /* CHEMISTRY_H_ */
