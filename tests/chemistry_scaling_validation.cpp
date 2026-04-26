#include "../Code/src/Chemistry/Chemistry.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

bool NearlyEqual(double lhs,double rhs,double relative_tolerance = 1.0E-12)
{
	const double scale = std::max(1.0,std::max(std::fabs(lhs),std::fabs(rhs)));
	return std::fabs(lhs-rhs)<=relative_tolerance*scale;
}

void RequireNearlyEqual(const std::string& name,double actual,double expected)
{
	if (!NearlyEqual(actual,expected)){
		std::cerr << name << " failed: actual=" << actual
		          << " expected=" << expected << std::endl;
		std::exit(1);
	}
}

void RequireZero(const std::string& name,double actual)
{
	RequireNearlyEqual(name,actual,0.0);
}

} // namespace

int main()
{
	ResetChemistryParametersToDefaults();
	ConfigureMineralVolumeScaling(2.0);
	ConfigureVpWidthCorrection(false);
	ConfigureEffectiveDiffusionHeightFactor(false,0.0,0.0);

	const double reference_volume = DELTA_V_REFERENCE_WATER_VOLUME;
	const double aperture_height = 1.0E-3;
	const double t_start = -5.0;
	const double t_end = 5000.0;
	const double clamped_t_start = 0.0;
	const double dV_ref = GetCumulativeMineralVolumeChange(t_end)
		- GetCumulativeMineralVolumeChange(clamped_t_start);

	const double single_particle_dV = ComputeParticleSegmentMineralVolumeChange(
		t_start,t_end,reference_volume,aperture_height);
	RequireNearlyEqual("single particle V_particle=Vref matches dV_ref",
		single_particle_dV,dV_ref);

	const int particle_count = 25;
	double summed_particle_dV = 0.0;
	for (int i=0;i<particle_count;i++){
		summed_particle_dV += ComputeParticleSegmentMineralVolumeChange(
			t_start,t_end,reference_volume/(double)particle_count,aperture_height);
	}
	RequireNearlyEqual("N particles with V_particle=Vref/N sum to dV_ref",
		summed_particle_dV,dV_ref);

	const int coarse_particle_count = 5;
	const int fine_particle_count = 125;
	double coarse_total = 0.0;
	for (int i=0;i<coarse_particle_count;i++){
		coarse_total += ComputeParticleSegmentMineralVolumeChange(
			100.0,t_end,reference_volume/(double)coarse_particle_count,aperture_height);
	}
	double fine_total = 0.0;
	for (int i=0;i<fine_particle_count;i++){
		fine_total += ComputeParticleSegmentMineralVolumeChange(
			100.0,t_end,reference_volume/(double)fine_particle_count,aperture_height);
	}
	RequireNearlyEqual("fixed total volume is independent of particle count",
		fine_total,coarse_total);

	RequireZero("zero representative volume gives zero particle dV",
		ComputeParticleSegmentMineralVolumeChange(0.0,t_end,0.0,aperture_height));
	RequireZero("negative representative volume gives zero particle dV",
		ComputeParticleSegmentMineralVolumeChange(0.0,t_end,-reference_volume,aperture_height));

	return 0;
}
