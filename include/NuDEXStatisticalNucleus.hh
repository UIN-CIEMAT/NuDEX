#ifndef NUDEXSTATISTICALNUCLEUS_HH
#define NUDEXSTATISTICALNUCLEUS_HH 1

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "NuDEXRandom.hh"
#include "NuDEXLevelDensity.hh"
#include "NuDEXInternalConversion.hh"
#include "NuDEXPSF.hh"



//This define remains:
//#define GENERATEEXPLICITLYALLLEVELSCHEME 1

//Class to obtain the level density for each excitation energy, spin, and parity
//All energies in MeV, all times in s
//Some of the class methods could be functions out of the class

struct Level{
  double Energy;
  int spinx2;
  bool parity; //true/false --> positive,negative
  unsigned int seed;
  int KnownLevelID;
  int NLevels;
  double Width;
};



//multipolarity of a transition is ...,-2,-1,0,1,2,... --> ...,M2,M1,Unk,E1,E2,...

struct KnownLevel{
  int id;
  double Energy;
  int spinx2;
  bool parity; //true/false --> positive,negative
  double T12; //half life - seconds
  int Ndecays;
  double* decayFraction;
  //std::string* decayMode;
  std::vector<std::string> decayMode;
  int NGammas;
  int *FinalLevelID,*multipolarity;
  double *Eg,*cumulPtot,*Pg,*Pe,*Icc;
};



int ComparisonLevels(const void* va, const void* vb);
void CopyLevel(Level* a,Level* b);
void CopyLevel(KnownLevel* a,Level* b);


class NuDEXStatisticalNucleus{

public:
  NuDEXStatisticalNucleus(int Z,int A);
  ~NuDEXStatisticalNucleus();

public:
  //Initialize everything. All the required files should be in dirname.
  //some of the data could also be in inputfname
  int Init(const char* dirname,const char* inputfname=0);

  //If InitialLevel==-1 then we start from the thermal capture level
  //If ExcitationEnergy>0 then is the excitation energy of the nucleus
  //If ExcitationEnergy<0 then is a capture reaction of a neutron with energy -ExcitationEnergy (MeV)
  int GenerateCascade(int InitialLevel,double ExcitationEnergy,std::vector<char>& pType,std::vector<double>& pEnergy,std::vector<double>& pTime);

  int GetClosestLevel(double Energy,int spinx2,bool parity); //if spinx2<0, then retrieves the closest level of any spin and parity
  double GetLevelEnergy(int i_level);
  void GetSnAndI0(double &sn,double &i0){sn=Sn; i0=I0;}
  Level* GetLevel(int i_level);
  void ChangeLevelSpinParityAndBR(int i_level,int newspinx2,bool newParity,int nlevels,double width,unsigned int seed=0); //if nlevels or width are negative they don't change. If seed (to generate the BR) is 0 it does not change.
  void ChangeThermalCaptureLevelBR(double LevelEnergy,double absoluteIntensity);

  void SetSomeInitalParameters(int LDtype=-1,int PSFFlag=-1,double MaxSpin=-1,int minlevelsperband=-1,double BandWidth_MeV=0,double maxExcEnergy=0,int BrOption=-1,int sampleGammaWidths=-1,unsigned int aseed1=0,unsigned int aseed2=0,unsigned int aseed3=0);
  void SetInitialParameters02(int knownLevelsFlag=-1,int electronConversionFlag=-1,double primGamNormFactor=-1,double primGamEcut=-1,double ecrit=-1);
  void SetBandWidth(double bandWidth){ if(bandWidth==0){bandWidth=-1;} BandWidth=bandWidth;} //So it is not re-written with the lib-params.
  void SetBrOption(int BrOption){BROpt=BrOption;}
  void SetRandom1Seed(unsigned int seed){theRandom1->SetSeed(seed); Rand1seedProvided=true;}
  void SetRandom2Seed(unsigned int seed){theRandom2->SetSeed(seed); Rand2seedProvided=true;}
  void SetRandom3Seed(unsigned int seed){theRandom3->SetSeed(seed); Rand3seedProvided=true;}
  
  NuDEXRandom* GetRandom3(){return theRandom3;}
  bool HasBeenInitialized(){return hasBeenInitialized;}


  //-------------------------------------------------------
  //Print:
  void PrintAll(std::ostream &out);
  
  void PrintParameters(std::ostream &out);
  void PrintKnownLevels(std::ostream &out);
  void PrintLevelDensity(std::ostream &out);
  void PrintLevelScheme(std::ostream &out);
  void PrintThermalPrimaryTransitions(std::ostream &out);
  void PrintPSF(std::ostream &out);
  void PrintICC(std::ostream &out);
  void PrintTotalCumulBR(int i_level,std::ostream &out);
  void PrintBR(int i_level,double MaxExcEneToPrint_MeV,std::ostream &out);
  void PrintInput01(std::ostream &out);
  //----------------
  void PrintKnownLevelsInDEGENformat(std::ostream &out);
  void PrintLevelSchemeInDEGENformat(const char* fname,int MaxLevelID=-1);
  //-------------------------------------------------------

  
private:
  //-------------------------------------------------------
  //Used by Init():
  //Read different data from files (do it in this order). If returnval<0 --> error reading file or nucleus not present in the file:
  int ReadSpecialInputFile(const char* fname);
  int ReadGeneralStatNuclParameters(const char* fname);
  double ReadEcrit(const char* fname);
  double ReadKnownLevels(const char* fname);
  void CreateLevelScheme();
  int InsertHighEnergyKnownLevels();
  void ComputeKnownLevelsMissingBR();
  void MakeSomeParameterChecks01();
  //-------------------------------------------------------
  double TakeTargetNucleiI0(const char* fname,int& check);
  void CreateThermalCaptureLevel(unsigned int seed=0); //If seed (to generate the BR) is 0 it does not change.
  void GenerateThermalCaptureLevelBR(const char* dirname);
  //-------------------------------------------------------

  //-------------------------------------------------------
  //cascade generation:
  double ComputeDecayIntensities(int i_level,double* cumulativeBR=0,double randnumber=-1,double TotGR=-1,bool AllowE1=false);
  int SampleFinalLevel(int i_level,int& multipolarity,double &icc_fac,int nTransition);
  int GetMultipolarity(Level* theInitialLevel,Level* theFinalLevel);
  //-------------------------------------------------------


private:
  //-------------------------------------------------------
  //Used to create the unknown Levels:
  int GenerateLevelsInBigRange(double Emin,double Emax,int spinx2,bool parity,Level* someLevels,int MaxNLevelsToFill); //salen sin ordenar
  int GenerateLevelsInSmallRange(double Emin,double Emax,int spinx2,bool parity,Level* someLevels,int MaxNLevelsToFill); //salen sin ordenar
  int GenerateWignerLevels(double Emin,double Emax,int spinx2,bool parity,Level* someLevels,int MaxNLevelsToFill); //salen ordenados
  int GenerateBandLevels(int bandmin,int bandmax,int spinx2,bool parity,Level* someLevels,int MaxNLevelsToFill);
  int GenerateAllUnknownLevels(Level* someLevels,int MaxNLevelsToFill); //salen ordenados
  int CreateBandsFromLevels(int thisNLevels,Level* someLevels,int spinx2,bool parity); 
  int EstimateNumberOfLevelsToFill(); //to estimate the length of "theLevels" vector
  //-------------------------------------------------------


private:

  //General info:
  int A_Int,Z_Int;
  double Sn,D0,I0; //I0 es el del nucleo A-1 (el que captura)
  bool hasBeenInitialized;
  std::string theLibDir;

  NuDEXRandom* theRandom1;  //To generate the unknown level scheme
  NuDEXRandom* theRandom2;  //To calculate the Gamma-rho values (i.e. to generate the branching ratios)
  NuDEXRandom* theRandom3;  //To generate the cascades
  unsigned int seed1,seed2,seed3;
  bool Rand1seedProvided,Rand2seedProvided,Rand3seedProvided;

  //--------------------------------------------------------------------------
  //Parameters which will define how the level scheme will be created:
  double Ecrit; //Energy between the known and unknown levels
  double MaxExcEnergy,BandWidth;
  int maxspinx2,NBands,MinLevelsPerBand; //maximum spin (x2) to consider, number of bands used to "rebin" the stat. part
  int LevelDensityType; //if negative or cero, use the default one.
  int PSFflag; // use IAEA PSF-data (PSFflag==0), use RIPL-3 data (PSFflag==1)
  double E_unk_min,E_unk_max; //min and max energy where the statistical part will be generated
  double Emin_bands,Emax_bands; //limites de energia para calcular las bandas de niveles
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //Level scheme:
  Level* theLevels; //known+unknown levels
  KnownLevel* theKnownLevels; // known levels
  int NKnownLevels,NUnknownLevels,NLevels,KnownLevelsVectorSize;
  Level theThermalCaptureLevel;
  int NLevelsBelowThermalCaptureLevel; //excluding the last one
  int KnownLevelsFlag;
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //Branching ratios:
  int BROpt,SampleGammaWidths;
  double* TotalGammaRho;
  double* theThermalCaptureLevelCumulBR;
  double** TotalCumulBR; //all BR
  double PrimaryGammasIntensityNormFactor;
  double PrimaryGammasEcut; //This variable can be used to avoid generating transitions close to the "Primary Gammas" region
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //LD,ICC, PSF:
  int ElectronConversionFlag;
  NuDEXLevelDensity* theLD;
  NuDEXInternalConversion* theICC;
  NuDEXPSF* thePSF;
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  //for internal use, when generating the cascades:
  int theSampledLevel,theSampledMultipolarity;
  //--------------------------------------------------------------------------
};

//***************************************************************************************************************
//***************************************************************************************************************




#endif




