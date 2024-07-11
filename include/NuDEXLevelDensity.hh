#ifndef NUDEXLEVELDENSITY_HH
#define NUDEXLEVELDENSITY_HH 1

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

#include "NuDEXRandom.hh"

//Level densities as they are defined in the RIPL-3 manual

//LDTYPE=1,2,3 --> Back-Shifted-Fermi-Gas model, Constant Temperature, Back-shifted: Egidy

#define DEFAULTLDTYPE 1

//using namespace std;

class NuDEXLevelDensity{

public:
  NuDEXLevelDensity(int aZ,int aA,int ldtype=DEFAULTLDTYPE);
  ~NuDEXLevelDensity(){}


  int ReadLDParameters(const char* dirname,const char* inputfname=0,const char* defaultinputfname=0);
  int CalculateLDParameters_BSFG(const char* dirname);
  int SearchLDParametersInInputFile(const char* inputfname);
  void GetSnD0I0Vals(double &aSn,double &aD0,double &aI0){aSn=Sn; aD0=D0; aI0=I0;}

  int GetLDType(){return LDType;}
  double GetNucleusTemperature(double ExcEnergy);
  double GetLevelDensity(double ExcEnergy_MeV,double spin,bool parity,bool TotalLevelDensity=false);
  double EstimateInverse(double LevDen_iMeV,double spin,bool parity); //an approximate value of ExcEnergy(rho), the inverse function of rho(ExcEnergy) - iMeV means 1/MeV
  double Integrate(double Emin,double Emax,double spin,bool parity);

  void PrintParameters(std::ostream &out);
  void PrintParametersInInputFileFormat(std::ostream &out);

private:

  //General info:
  int A_Int,Z_Int;
  int LDType; //=1,2,3 --> Back-Shifted-Fermi-Gas model, Constant Temperature, Back-shifted: Egidy
  double Sn,D0,I0; //I0 es el del nucleo A-1 (el que captura)
  double Ed;

  bool HasData;

  //Level density parameters:
  double A_mass,ainf_ldpar,gamma_ldpar,dW_ldpar,Delta_ldpar,T_ldpar,E0_ldpar,Ex_ldpar;


};



#endif

