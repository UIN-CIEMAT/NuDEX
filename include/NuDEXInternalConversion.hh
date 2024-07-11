#ifndef NUDEXINTERNALCONVERSION_HH
#define NUDEXINTERNALCONVERSION_HH 1


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

#include "NuDEXRandom.hh"

#define ICC_MAXNSHELLS 40
#define ICC_NMULTIP 5
#define MINZINTABLES 10 //below this value, the alpha is always 0

/*
Class to manage the internal conversion factors and the generation of converted-e-
Still not included the fluorescence-auger effects, i.e., what happens with the hole
We read the occ factors from a file, and they are stored in a matrix
The total Icc are in index=0 (data from the libraries) and index=NShells (sum of the partials)
Data are taken from: https://doi.org/10.1006/adnd.2002.0884
*/

class NuDEXInternalConversion{

public:
  NuDEXInternalConversion(int Z);
  ~NuDEXInternalConversion();
  void Init(const char* fname);
  void PrintICC(std::ostream &out);
  double GetICC(double Ene,int multipolarity,int i_shell=-1);
  bool SampleInternalConversion(double Ene,int multipolarity,double alpha=-1,bool CalculateProducts=true);
  void FillElectronHole(int i_shell); //Fluorescence/auger
  void SetRandom4Seed(unsigned int seed){theRandom4->SetSeed(seed);}


private:
  double Interpolate(double val,int npoints,double* x,double* y);
  void MakeTotal();


private:
  int theZ,NShells;
  double BindingEnergy[ICC_MAXNSHELLS];
  double *Eg[ICC_MAXNSHELLS],*Icc_E[ICC_NMULTIP][ICC_MAXNSHELLS],*Icc_M[ICC_NMULTIP][ICC_MAXNSHELLS];
  int np[ICC_MAXNSHELLS];
  std::string OrbitalName[ICC_MAXNSHELLS];
  NuDEXRandom* theRandom4;  

public:
  int Ne,Ng;
  double Eele[100],Egam[100];
};



#endif

