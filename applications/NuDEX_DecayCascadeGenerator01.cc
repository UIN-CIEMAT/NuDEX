

#include "NuDEXStatisticalNucleus.hh"
#include <cstring>
#include <iomanip>

using namespace std;

#define MAXNSTARTINGLEVELS 10000

/*

Program to compute (n,g) reaction cascades using NuDEX 

*/

int main(int argc,char** argv){


  if(argc<3){
    std::cout<<" #########################################################################  "<<std::endl;
    std::cout<<" This program can be executed as (two possibilities): "<<std::endl;
    std::cout<<"    NuDEX_DecayCascadeGenerator01 1 outputfilebase LIBDIR ZA EXCENE_MEV [keyname1] [val1] [keyname2] [val2] ..."<<std::endl;
    std::cout<<"    NuDEX_DecayCascadeGenerator01 2 outputfilebase inputfile [keyname1] [val1] [keyname2] [val2] ..."<<std::endl;
    std::cout<<" #########################################################################  "<<std::endl;
    return 1;
  }


  //-------------------------------------------------------------------------------------
  //Define variables:
  int LDtype=-1;
  int PSFflag=-1;
  double MaxSpin=-1;
  int minlevelsperband=-1;
  double BandWidth_MeV=0;
  double MaxExcEnergy=0; // MeV 
  int BrOption=-1;
  int sampleGammaWidths=-1;
  unsigned int seed1=0;
  unsigned int seed2=0;
  unsigned int seed3=0;
  int knownLevelsFlag=-1;
  int electronConversionFlag=-1;
  double primGamNormFactor=-1;
  double primGamEcut=-1;
  double ecrit=-1;
  //----------------------------
  char LibDir[200];
  int ZA=0; //ZA of the  nucleus
  int NCascades=100; //number of cascades to be generated
  //----------------------------
  // Where the cascade starts:
  int NStartingLevels=0;
  double ExcitationEnergy[MAXNSTARTINGLEVELS]; //in MeV
  double Spin[MAXNSTARTINGLEVELS];  // spin of the level
  bool parity[MAXNSTARTINGLEVELS];  // parity of the level (1)
  double FeedingProb[MAXNSTARTINGLEVELS];
  double FeedingCumulProb[MAXNSTARTINGLEVELS];
  int StartingLevelID[MAXNSTARTINGLEVELS];
  //-----------------------------
  double TimeWindow=10e20; //ns - gammas and electrons emitted after TimeWindow are not written in the output
  //-----------------------------
  // (1) true --> positive, false --> negative. 
  //-------------------------------------------------------------------------------------

  for(int i=0;i<MAXNSTARTINGLEVELS;i++){
    ExcitationEnergy[i]=-1;
    Spin[i]=-1;
    parity[i]=true;
    FeedingProb[i]=0;
    FeedingCumulProb[i]=0;
    StartingLevelID[i]=0;
  }

  //--------------------------------------------------------
  char* inputfname=0;
  int exemode=std::atoi(argv[1]);
  char* outfname=argv[2];
  if(exemode==1){
    sprintf(LibDir,"%s",argv[3]);
    ZA=std::atoi(argv[4]);
    ExcitationEnergy[0]=std::atof(argv[5]);
    FeedingProb[0]=1;
    FeedingCumulProb[0]=1;
    NStartingLevels=1;
  }
  else if(exemode==2){
    inputfname=argv[3];
  }
  else{
    cout<<" ############ ERROR: Unknown input mode ---> "<<exemode<<"  ############"<<std::endl; return 1;
  }
  char outfname_out[500],outfname_cas[500],outfname_inp[500];
  sprintf(outfname_out,"%s.out",outfname);
  sprintf(outfname_cas,"%s.cas",outfname);
  sprintf(outfname_inp,"%s.inp",outfname);
  //--------------------------------------------------------

  //--------------------------------------------------------
  //read input:
  if(inputfname!=0){
    string word;
    ifstream in(inputfname);
    double par;
    while(in>>word){
      if(word.c_str()[0]=='#'){in.ignore(10000,'\n');}
      if(word==string("END")){break;}
      else if(word==string("LIBDIR")){in>>LibDir;}
      else if(word==string("ZA")){in>>ZA;}
      else if(word==string("NCASCADES")){in>>NCascades;}

      else if(word==string("STARTINGLEVEL")){
	in>>ExcitationEnergy[NStartingLevels]>>Spin[NStartingLevels]>>par>>FeedingProb[NStartingLevels];
	if(par>0){parity[NStartingLevels]=true;}else{parity[NStartingLevels]=false;}
	NStartingLevels++;
      }
      
      else if(word==string("LEVELDENSITYTYPE")){in>>LDtype;}
      else if(word==string("MAXSPIN")){in>>MaxSpin;}
      else if(word==string("MINLEVELSPERBAND")){in>>minlevelsperband;}
      else if(word==string("BANDWIDTH_MEV")){in>>BandWidth_MeV;}
      else if(word==string("MAXEXCENERGY_MEV")){in>>MaxExcEnergy;}
      else if(word==string("ECRIT_MEV")){in>>ecrit;}
      else if(word==string("KNOWNLEVELSFLAG")){in>>knownLevelsFlag;}

      else if(word==string("PSF_FLAG")){in>>PSFflag;}
      else if(word==string("BROPTION")){in>>BrOption;}
      else if(word==string("SAMPLEGAMMAWIDTHS")){in>>sampleGammaWidths;}
      
      else if(word==string("SEED1")){in>>seed1;}
      else if(word==string("SEED2")){in>>seed2;}
      else if(word==string("SEED3")){in>>seed3;}

      else if(word==string("ELECTRONCONVERSIONFLAG")){in>>electronConversionFlag;}
      else if(word==string("TIMEWINDOW_NS")){in>>TimeWindow;}
    }
    in.close();
  }
  //--------------------------------------------------------


  //--------------------------------------------------------
  //input parameters from the input line:
  int nparameters=0;
  int i_firstpar=5;
  if(exemode==1){
    nparameters=(argc-5)/2;
    i_firstpar=5;
  }
  else{
    nparameters=(argc-3)/2;
    i_firstpar=3;
  }
  if(nparameters>0){
    cout<<" List of input parameters from exe line ("<<nparameters<<"):"<<std::endl;
  }
  for(int i=0;i<nparameters;i++){
    char* parname=argv[i_firstpar+2*i];
    if(string(parname)==string("LIBDIR")){sprintf(LibDir,"%s",argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<LibDir<<std::endl;}
    else if(string(parname)==string("ZA")){ZA=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<ZA<<std::endl;}
    else if(string(parname)==string("NCASCADES")){NCascades=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<NCascades<<std::endl;}
    
    else if(string(parname)==string("LEVELDENSITYTYPE")){LDtype=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<LDtype<<std::endl;}
    else if(string(parname)==string("MAXSPIN")){MaxSpin=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<MaxSpin<<std::endl;}
    else if(string(parname)==string("MINLEVELSPERBAND")){minlevelsperband=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<minlevelsperband<<std::endl;}
    else if(string(parname)==string("BANDWIDTH_MEV")){BandWidth_MeV=std::atof(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<BandWidth_MeV<<std::endl;}
    else if(string(parname)==string("MAXEXCENERGY_MEV")){MaxExcEnergy=std::atof(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<MaxExcEnergy<<std::endl;}
    else if(string(parname)==string("ECRIT_MEV")){ecrit=std::atof(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<ecrit<<std::endl;}
    else if(string(parname)==string("KNOWNLEVELSFLAG")){knownLevelsFlag=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<knownLevelsFlag<<std::endl;}
    
    else if(string(parname)==string("PSF_FLAG")){PSFflag=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<PSFflag<<std::endl;}
    else if(string(parname)==string("BROPTION")){BrOption=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<BrOption<<std::endl;}
    else if(string(parname)==string("SAMPLEGAMMAWIDTHS")){sampleGammaWidths=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<sampleGammaWidths<<std::endl;}
    
    else if(string(parname)==string("SEED1")){seed1=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<seed1<<std::endl;}
    else if(string(parname)==string("SEED2")){seed2=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<seed2<<std::endl;}
    else if(string(parname)==string("SEED3")){seed3=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<seed3<<std::endl;}
    
    else if(string(parname)==string("ELECTRONCONVERSIONFLAG")){electronConversionFlag=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<electronConversionFlag<<std::endl;}
    else if(string(parname)==string("TIMEWINDOW_NS")){TimeWindow=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<TimeWindow<<std::endl;}

    else{
      cout<<" ############ ERROR: Unknown input parameter ---> "<<parname<<"  ############"<<std::endl; return 1;
    }
  }
  //--------------------------------------------------------

  
  //--------------------------------------------------------
  int Z=ZA/1000;
  int A=ZA-1000*Z;
  //--------------------------------------------------------
  //--------------------------------------------------------
  if(MaxExcEnergy==0){ //if not initalized
    for(int i=0;i<NStartingLevels;i++){
      if(MaxExcEnergy<ExcitationEnergy[i]){MaxExcEnergy=ExcitationEnergy[i];}
    }
    MaxExcEnergy+=1.0; //One MeV more than the feeding level with the maximum energy
  }
  //--------------------------------------------------------
  //--------------------------------------------------------
  //Create statistical nucleus
  NuDEXStatisticalNucleus* theStatisticalNucleus=new NuDEXStatisticalNucleus(Z,A);
  theStatisticalNucleus->SetSomeInitalParameters(LDtype,PSFflag,MaxSpin,minlevelsperband,BandWidth_MeV,MaxExcEnergy,BrOption,sampleGammaWidths,seed1,seed2,seed3);
  theStatisticalNucleus->SetInitialParameters02(knownLevelsFlag,electronConversionFlag,primGamNormFactor,primGamEcut,ecrit);
  int check=theStatisticalNucleus->Init(LibDir,inputfname);
  if(check<0){
    std::cout<<" Error initializing StatisticalNucleus with Z = "<<Z<<" , A = "<<A<<std::endl;
    return 0;
  }
  //--------------------------------------------------------
  //--------------------------------------------------------
  //Print the output file:
  std::ofstream outp(outfname_out);
  theStatisticalNucleus->PrintParameters(outp);
  theStatisticalNucleus->PrintKnownLevels(outp);
  theStatisticalNucleus->PrintLevelDensity(outp);
  //theStatisticalNucleus->PrintLevelScheme(outp);
  theStatisticalNucleus->PrintThermalPrimaryTransitions(outp);
  theStatisticalNucleus->PrintPSF(outp);
  //theStatisticalNucleus->PrintICC(outp);
  if(!outp.good()){
    std::cout<<" ############# Error opening or writting in "<<outfname_out<<" #############"<<std::endl; exit(1);
  }
  outp.close();
  //--------------------------------------------------------


  //--------------------------------------------------------
  //Starting levels:
  double TotalProb=0;
  for(int i=0;i<NStartingLevels;i++){
    TotalProb+=FeedingProb[i];
  }
  if(TotalProb<=0){
    std::cout<<" ############# Error defining feeding probabilities. Sum of probabilities is "<<TotalProb<<" #############"<<std::endl; exit(1);
  }
  for(int i=0;i<NStartingLevels;i++){
    FeedingProb[i]=FeedingProb[i]/TotalProb;
    if(i==0){
      FeedingCumulProb[i]=FeedingProb[i];
    }
    else{
      FeedingCumulProb[i]=FeedingCumulProb[i-1]+FeedingProb[i];
    }
    if(Spin[i]<0){
      StartingLevelID[i]=theStatisticalNucleus->GetClosestLevel(ExcitationEnergy[i],-1,true);
    }
    else{
      StartingLevelID[i]=theStatisticalNucleus->GetClosestLevel(ExcitationEnergy[i],(int)(Spin[i]*2+0.1),parity[i]);
    }
  }
  //Print:
  std::cout<<" Number of starting levels: "<<NStartingLevels<<std::endl;
  std::cout<<" ================================================================================ "<<std::endl;
    std::cout<<std::setw(15)<<"ENERGY"<<" "<<std::setw(10)<<"SPIN"<<" "<<std::setw(10)<<"PARITY"<<" "<<std::setw(15)<<"FEEDING_PROB"<<" "<<std::setw(15)<<"FEEDING_CUMULPROB"<<std::endl;
  std::cout<<" ================================================================================ "<<std::endl;
  for(int i=0;i<NStartingLevels;i++){
    Level* iLevel=theStatisticalNucleus->GetLevel(StartingLevelID[i]);
    std::cout<<std::setw(15)<<iLevel->Energy<<" "<<std::setw(10)<<iLevel->spinx2/2.<<" "<<std::setw(10)<<iLevel->parity<<" "<<std::setw(15)<<FeedingProb[i]<<" "<<std::setw(15)<<FeedingCumulProb[i]<<std::endl;
  }
  std::cout<<" ================================================================================ "<<std::endl;
  //--------------------------------------------------------
  


  //--------------------------------------------------------
  //Write input:
  std::ofstream outi(outfname_inp);
  outi<<std::endl;
  outi<<"LIBDIR "<<LibDir<<std::endl;
  outi<<"ZA "<<ZA<<std::endl;
  outi<<"NCASCADES "<<NCascades<<std::endl;
  outi<<"TIMEWINDOW_NS "<<TimeWindow<<std::endl;
  outi<<std::endl;
  for(int i=0;i<NStartingLevels;i++){
    double par=1;
    if(!parity[i]){par=-1;}
    outi<<"STARTINGLEVEL "<<ExcitationEnergy[i]<<"  "<<Spin[i]<<"  "<<par<<"  "<<FeedingProb[i]<<std::endl;
  }
  outi<<std::endl;
  theStatisticalNucleus->PrintInput01(outi);
  if(!outi.good()){
    std::cout<<" ############# Error opening or writting in "<<outfname_inp<<" #############"<<std::endl; exit(1);
  }
  outi.close();
  //--------------------------------------------------------

     
  //--------------------------------------------------------
  //Generate the cascades:
  int Npar,Npar2,i_lev;
  std::vector<char> pType;
  std::vector<double> pEnergy,pTime;
  double rand;
  NuDEXRandom* theRand=theStatisticalNucleus->GetRandom3();
  std::ofstream out(outfname_cas);
  if(!out.good()){
    std::cout<<" ######## Error opening "<<outfname_cas<<" ########"<<std::endl; exit(1);
  }
  for(int i=0;i<NCascades;i++){
    rand=theRand->Uniform();
    i_lev=NStartingLevels-1;
    for(int j=0;j<NStartingLevels;j++){
      if(rand<FeedingCumulProb[j]){
	i_lev=j; break;
      }
    }
    Level* iLevel=theStatisticalNucleus->GetLevel(StartingLevelID[i_lev]);
    Npar=theStatisticalNucleus->GenerateCascade(StartingLevelID[i_lev],iLevel->Energy,pType,pEnergy,pTime);
    Npar2=0; 
    for(int j=0;j<Npar;j++){
      if(pTime.at(j)<TimeWindow*1.e-9){Npar2++;}
    }
    out<<Npar2;
    for(int j=0;j<Npar;j++){
      if(pTime.at(j)<TimeWindow*1.e-9){out<<"  "<<pType.at(j)<<"  "<<pEnergy.at(j)<<"  "<<pTime.at(j);}
    }
    out<<std::endl;
    if(((i+1)*10)%NCascades==0){
      std::cout<<(i+1.)/NCascades*100.<<" % done"<<std::endl;
    }
  }
  out.close();
  //--------------------------------------------------------

  delete theStatisticalNucleus;

  return 0;

}


















