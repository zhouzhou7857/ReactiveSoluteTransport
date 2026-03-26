/*
 * Chemistry.h
 *
 *  Created on: 2025/12/11
 *      Author: delphine and Wenyu
 */

#ifndef CHEMISTRY_H_
#define CHEMISTRY_H_

//namespace chemistry {

const double K_REACTION = -6E-7; // Reaction factor calculated from trend line slope of Delta_V [MODIFIED!!][exp = 1.5e-10 m/s]
// Gypsum dissolving rate is about 6e-9 m/s， dividing by porosity correction would be 6.7E-9 m/s
/**
 * Compute aperture b from the empirical formula:
 *   b_new = b - K_REACTION * t * (Nb / Nt)
 * @param b  current aperture
 * @param t   reaction time (s)
 * @param Nb  particle number through the segment
 * @param Nt  total particles number
 * @param clamp_nonnegative if true, b is clamped to >= 0 (Delta_b can be negative)
 * @return b_new  aperture (new aperture)
 */
double ComputeAperture(double b, double t, double Nb, double Nt, bool clamp_nonnegative = true);

/**
 * Compute delta aperture for updating an existing aperture:
 *   delta_b = b_new - b_old
 *
 * @param b current aperture
 * @param t     reaction time
 * @param Nb    particles through the segment
 * @param Nt    total particles
 * @param clamp_nonnegative if true, b_new is clamped to >= 0 before differencing
 * @return delta_b
 */
double ComputeDeltaAperture(double b, double t, double Nb, double Nt, bool clamp_nonnegative = true);
double safe_ratio(double Nb, double Nt);

//}

#endif /* CHEMISTRY_H_ */
