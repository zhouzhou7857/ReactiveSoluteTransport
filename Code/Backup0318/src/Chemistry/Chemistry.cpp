/*
 * Chemistry.cpp
 *
 *  Created on: 2025/12/11
 *      Author: delphine and Wenyu
 */

#include "Chemistry.h"

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
double safe_ratio(double Nb, double Nt)
{
    if (Nt <= 0) return 0.0;
    if (Nb <= 0) return 0.0;
    return (double)Nb / (double)Nt;
}
// Aperture update using the effective segment throughput Nb normalized by
// a reference injected particle count Nt for the same reaction time step.
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
// Compute Delta aperture, maybe useful in following calculation
double ComputeDeltaAperture(double b, double t, double Nb, double Nt, bool clamp_nonnegative)
{
    const double b_new = ComputeAperture(b, t, Nb, Nt, clamp_nonnegative);
    return b_new - b;
}
//}
