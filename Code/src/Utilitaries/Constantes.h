/*
 * Constantes.h
 *
 *  Created on: Jul 2, 2012
 *      Author: delphine
 */

#ifndef CONSTANTES_H_
#define CONSTANTES_H_

#include "Pointcpp.h"
#include <math.h>

#ifndef EPSILON
#define EPSILON 1E-10		//epsilon value
#endif

#ifndef EPSILON_DIST
#define EPSILON_DIST 1E-2		//epsilon value
#endif

#ifndef MAX_TIME
#define MAX_TIME 50*365*24*60*60	//3000//140000//	//epsilon value
#endif

#ifndef RHO
#define RHO 1E3		//fluid density
#endif

#ifndef G
#define G 9.8		//gravity
#endif

#ifndef MU
#define MU 1E-3		//dynamic fluid viscosity
#endif

static const double NOT_DEFINED=-999;
static const double PI=4.*atan(1.0);

const int SCALE_LIN = 0;
const int SCALE_LOG = 1;

const int INFINITE_MATRIX = 0;
const int FINITE_MATRIX = 1;

// Boundary conditions definition
const std::string DIRICHLET="DIRICHLET";
const std::string NEUMANN="NEUMANN";
const std::string DIR_GRADIENT="DIR_GRADIENT";	// Dirichlet conditions where the value along the border follows a gradient
const std::string MIXED_BC="MIXED_BC";
const std::string SOURCE_TERM="SOURCE_TERM";

// Domain borders
const std::string LEFT_BORDER="Left_Border";
const std::string RIGHT_BORDER="Right_Border";
const std::string BOTTOM_BORDER="Bottom_Border";
const std::string TOP_BORDER="Top_Border";
const std::string NO_BORDER="No_Border";

// Boundary conditions configuration
const std::string CONFIG1="CONFIG1";	// Configuration for ellipse permeability evaluation with impervious conditions for the lateral borders
const std::string CONFIG2="CONFIG2";	// Configuration for ellipse permeability evaluation with gradient potential on the lateral borders
const std::string CONFIG3="CONFIG3";	// Configuration for ellipse permeability evaluation with gradient potential on the transverse borders

const std::string NO_POSITION="No_Position";

#endif /* CONTANTES_H_ */
