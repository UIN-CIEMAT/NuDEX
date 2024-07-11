#ifndef NUDEXPSF_HH
#define NUDEXPSF_HH 1

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

//using namespace std;

#include "NuDEXLevelDensity.hh"

/*
All energies in MeV
PSF are defined as in RIPL-3: PSF=Eg**(-2L-1) x Gamma width x level density
JL defines PSF x Eg**(2L+1) instead

  PSFType=0 --> SLO
  PSFType=1 --> EGLO, as defined in RIPL-3, but using always Tf in the formula
  PSFType=2 --> SMLO, as defined in RIPL-3
  PSFType=3 --> GLO  (like EGLO, but k1=k2=1)
  PSFType=4 --> MGLO (like EGLO, but k2=1)
  PSFType=5 --> KMF
  PSFType=6 --> GH
  PSFType=7 --> EGLO, but the k parameter is provided (MEGLO)
  PSFType=8 --> EGLO, but the "k1" and "k2" parameters are provided (MEGLO)
  PSFType=9 --> EGLO, but the k parameter and a constant temperature of the nucleus is provided (MEGLO)
  PSFType=10 --> EGLO, but the "k1" and "k2" parameters and a constant temperature of the nucleus are provided (MEGLO)
  PSFType=11 --> SMLO, as defined in Eur. Phys. J. A (2019) 55: 172
  PSFType=20 --> gaussian (to simulate small bumps or resonances)
  PSFType=21 --> expo -->  C*exp(-eta*Eg). It is defined with three entries: C eta dummy
  PSFType=40 --> pointwise function type 1 (only input file)
  PSFType=41 --> pointwise function type 2 (only input file)

Procedure to obtain the PSF, in order of hierarchy:
  - Get the data from inputfname
  - Get the data from PSF_param.dat file
  - Get the data from IAEA-2019 PSF values (if PSFflag==0)
  - Get the data from RIPL-3 experimental MLO values --> gdr-parameters&errors-exp-MLO.dat
  - Get the data from RIPL-3 Theorethical values --> gdr-parameters-theor.dat
  - Use RIPL-3 and RIPL-2 theoretical formulas
*/



class NuDEXPSF{

public:
  NuDEXPSF(int aZ,int aA);
  ~NuDEXPSF();

  int Init(const char* dirname,NuDEXLevelDensity* aLD,const char* inputfname=0,const char* defaultinputfname=0,int PSFflag=0);
  double GetE1(double Eg,double ExcitationEnergy);
  double GetM1(double Eg,double ExcitationEnergy);
  double GetE2(double Eg,double ExcitationEnergy);
  void PrintPSFParameters(std::ostream &out);
  void PrintPSFParametersInInputFileFormat(std::ostream &out);

private:

  bool TakePSFFromInputFile(const char* fname);
  bool TakePSFFromDetailedParFile(const char* fname);
  bool TakePSFFromIAEA01(const char* fname); // IAEA - PSF values 2019
  bool TakePSFFromRIPL01(const char* fname); // RIPL3-MLO values
  bool TakePSFFromRIPL02(const char* fname); // RIPL3-Theorethical values
  void GenerateM1AndE2FromE1(); // From RIPL-3 and RIPL-2 recommendations


  //Shapes:
  //Typical ones:
  double SLO(double Eg,double Er,double Gr,double sr);                          //PSFType=0
  double EGLO(double Eg,double Er,double Gr,double sr,double ExcitationEnergy); //PSFType=1
  double SMLO(double Eg,double Er,double Gr,double sr,double ExcitationEnergy); //PSFType=2
  double GLO(double Eg,double Er,double Gr,double sr,double ExcitationEnergy);  //PSFType=3
  double MGLO(double Eg,double Er,double Gr,double sr,double ExcitationEnergy); //PSFType=4
  double KMF(double Eg,double Er,double Gr,double sr,double ExcitationEnergy);  //PSFType=5
  double GH(double Eg,double Er,double Gr,double sr,double ExcitationEnergy);   //PSFType=6
  double MEGLO(double Eg,double Er,double Gr,double sr,double ExcitationEnergy,double k_param1,double k_param2,double Temp=-1);//PSFType=6,7,8,9,10
  double SMLO_v2(double Eg,double Er,double Gr,double sr,double ExcitationEnergy); //PSFType=11


  double Gauss(double Eg,double Er,double Gr,double sr); //PSFType=20
  double Expo(double Eg,double C,double eta); //PSFType=21

  //PSFType=40, PSFType=41  are pointwise defined functions
  
  //------------------------------
  double EGLO_GLO_MGLO(double Eg,double Er,double Gr,double sr,double ExcitationEnergy,int Opt);
  double FlexibleGLOType(double Eg,double Er,double Gr,double sr,double Temp1,double k_param1,double Temp2,double k_param2);
  double Gamma_k(double Eg,double Er,double Gr,double Temp,double k_param);

private:
  int Z_Int,A_Int;

  int nR_E1,nR_M1,nR_E2;
  int PSFType_E1[10], PSFType_M1[10], PSFType_E2[10];
  double E_E1[10],G_E1[10],s_E1[10],p1_E1[10],p2_E1[10],p3_E1[10]; 
  double E_M1[10],G_M1[10],s_M1[10],p1_M1[10],p2_M1[10],p3_M1[10]; 
  double E_E2[10],G_E2[10],s_E2[10],p1_E2[10],p2_E2[10],p3_E2[10]; 

  //-----------------------------------------------
  //PSF pointwise defined PSF --> PSFType=3,4,6
  int np_E1,np_M1,np_E2;
  double *x_E1,*y_E1;
  double *x_M1,*y_M1;
  double *x_E2,*y_E2;
  double E1_normFac,M1_normFac,E2_normFac;
  double NormEmin,NormEmax;
  //-----------------------------------------------

  double ScaleFactor_E1,ScaleFactor_M1,ScaleFactor_E2;

  double EvaluateFunction(double xval,int np,double* x,double* y);
  void Renormalize();

  NuDEXLevelDensity* theLD;
};




#endif

