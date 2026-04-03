
#ifndef RNGSTREAM_H
#define RNGSTREAM_H

#include <string>
#include <iostream>

/** @addtogroup RNG RANDOM NUMBER GENERATOR
 *  @{
 */

class RngStream
{
public:

RngStream (const char *name = "");

void setName(const char *name_){ name = name_;}


static bool SetPackageSeed (const unsigned long seed[6]);


void ResetStartStream ();


void ResetStartSubstream ();


void ResetNextSubstream ();


void SetAntithetic (bool a);


void IncreasedPrecis (bool incp);


bool SetSeed (const unsigned long seed[6]);


void AdvanceState (long e, long c);


void GetState (unsigned long seed[6]) const;


void WriteState () const;


void WriteStateFull () const;


double RandU01 ();


int RandInt (int i, int j);



private:

double Cg[6], Bg[6], Ig[6];


bool anti, incPrec;


std::string name;


static double nextSeed[6];


double U01 ();


double U01d ();

friend std::ostream & operator << ( std::ostream & sortie, RngStream & rng);
friend std::istream & operator >> ( std::istream & entree, RngStream & rng);

};


/** @} */ // end of RNG

#endif



