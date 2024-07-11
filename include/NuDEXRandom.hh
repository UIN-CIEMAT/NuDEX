#ifndef NUDEXRANDOM_HH
#define NUDEXRANDOM_HH 1

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

//COMPILATIONTYPE==1 compile with ROOT
//COMPILATIONTYPE==2 compile with GEANT4

#define COMPILATIONTYPE 1

#if COMPILATIONTYPE == 1
//------------------------------------------------------------
// ROOT
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "TRandom2.h"
#pragma GCC diagnostic pop
//------------------------------------------------------------
#elif COMPILATIONTYPE == 2
//------------------------------------------------------------
// GEANT4
#include "Randomize.hh"
#include "globals.hh"
#include "G4Exception.hh"
//------------------------------------------------------------
#else
  #error Unsupported COMPILATIONTYPE setting
#endif

void NuDEXException(const char* originOfException,const char* exceptionCode,const char* description);

class NuDEXRandom{

public:
  NuDEXRandom(unsigned int seed);
  ~NuDEXRandom();

public:
  void SetSeed(unsigned int seed);
  unsigned int GetSeed();
  double Uniform(double Xmin=0,double Xmax=1);
  unsigned int Integer(unsigned int IntegerMax);
  double Exp(double tau);
  double Gaus(double mean=0,double sigma=1);
  int Poisson(double mean);

private:

#if COMPILATIONTYPE == 1
  TRandom2* theRandom;
#elif COMPILATIONTYPE == 2
  CLHEP::HepJamesRandom* theEngine;
  CLHEP::RandFlat* theRandFlat;
  CLHEP::RandExponential* theRandExponential;
  CLHEP::RandGauss* theRandGauss;
  CLHEP::RandPoisson* theRandPoisson;
#endif
};


#endif




