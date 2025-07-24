

#include "NuDEXStatisticalNucleus.hh"
#include <cstring>

using namespace std;

/*

Program to compute (n,g) reaction cascades using NuDEX 

*/

int main(int argc,char** argv){


  if(argc<3){
    std::cout<<" #########################################################################  "<<std::endl;
    std::cout<<" This program can be executed as (two possibilities): "<<std::endl;
    std::cout<<"    NuDEX_NCaptureCascadeGenerator01 outputfilebase LIBDIR ZA [keyname1] [val1] [keyname2] [val2] ..."<<std::endl;
    std::cout<<"    NuDEX_NCaptureCascadeGenerator01 outputfilebase inputfile [keyname1] [val1] [keyname2] [val2] ..."<<std::endl;
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
  double MaxExcEnergy=0; // MeV --> this value should be larger than "Sn+(A-1)/A*NeutronEnergy" (2)
  int BrOption=-1;
  int sampleGammaWidths=-1;
  unsigned int seed1=0;
  unsigned int seed2=0;
  unsigned int seed3=0;
  unsigned int seed4=1234567; // to get the BR in the capture level
  int knownLevelsFlag=-1;
  int electronConversionFlag=-1;
  double primGamNormFactor=-1;
  double primGamEcut=-1;
  double ecrit=-1;
  //----------------------------
  char LibDir[200];
  int ZA=0; //ZA of the target nucleus
  int NCascades=100; //number of cascades to be generated
  bool DoThermal=true; //If true, create thermal cascades.
  //----------------------------
  // Where the capture reaction is produced:
  double NeutronEnergy=1.e-6; //in MeV
  double JSpin=-0.5;  // spin of the resonance
  bool parity=true;  // parity of the resonance (1)
  double TimeWindow=10e20; //ns - gammas and electrons emitted after TimeWindow are not written in the output
  //-----------------------------
  // (1) true --> positive, false --> negative. If s-wave or d-wave, same parity as ZA in ground state. If p-wave, the other parity.
  // (2) The level scheme of the compound nucleus is created up to MaxExcEnergy
  //-------------------------------------------------------------------------------------


  //--------------------------------------------------------
  char* inputfname=0;
  char* outfname=argv[1];
  int exemode=1; //NuDEX_NCaptureCascadeGenerator01 outputfilebase LIBDIR ZA [keyname1] [val1] [keyname2] [val2] ...
  if((argc%2)==1){
    exemode=2; //NuDEX_NCaptureCascadeGenerator01 outputfilebase inputfile [keyname1] [val1] [keyname2] [val2] ...
  }
  if(exemode==2){
    inputfname=argv[2];
  }
  else{
    sprintf(LibDir,"%s",argv[2]);
    ZA=std::atoi(argv[3]);
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
    while(in>>word){
      if(word.c_str()[0]=='#'){in.ignore(10000,'\n');}
      if(word==string("END")){break;}
      else if(word==string("LIBDIR")){in>>LibDir;}
      else if(word==string("ZA")){in>>ZA;}
      else if(word==string("NCASCADES")){in>>NCascades;}
      else if(word==string("NEUTRONENERGY_MEV")){in>>NeutronEnergy;}
      else if(word==string("JSPIN")){in>>JSpin;}
      else if(word==string("PARITY")){double par; in>>par; if(par>0){parity=true;}else{parity=false;}}
      else if(word==string("DOTHERMALCASCADE")){double therm; in>>therm; if(therm<=0){DoThermal=false;}else{DoThermal=true;}}
      
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
      else if(word==string("SEED4")){in>>seed4;}

      else if(word==string("ELECTRONCONVERSIONFLAG")){in>>electronConversionFlag;}
      else if(word==string("TIMEWINDOW_NS")){in>>TimeWindow;}
      else if(word==string("PRIMARYTHCAPGAMNORM")){in>>primGamNormFactor;}
      else if(word==string("PRIMARYGAMMASECUT")){in>>primGamEcut;}
    }
    in.close();
  }
  //--------------------------------------------------------


  //--------------------------------------------------------
  //input parameters from the input line:
  int nparameters=0;
  int i_firstpar=4;
  if(exemode==1){
    nparameters=(argc-4)/2;
    i_firstpar=4;
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
    else if(string(parname)==string("NEUTRONENERGY_MEV")){NeutronEnergy=std::atof(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<NeutronEnergy<<std::endl;}
    else if(string(parname)==string("JSPIN")){JSpin=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<JSpin<<std::endl;}
    else if(string(parname)==string("PARITY")){double par=std::atof(argv[i_firstpar+2*i+1]); if(par<=0){parity=false;}else{parity=true;}  cout<<"      "<<parname<<"  "<<parity<<std::endl;}
    else if(string(parname)==string("DOTHERMALCASCADE")){ double dothermal=std::atof(argv[i_firstpar+2*i+1]); if(dothermal<=0){DoThermal=false;}else{DoThermal=true;}  cout<<"      "<<parname<<"  "<<DoThermal<<std::endl;}
    
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
    else if(string(parname)==string("SEED4")){seed4=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<seed4<<std::endl;}
    
    else if(string(parname)==string("ELECTRONCONVERSIONFLAG")){electronConversionFlag=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<electronConversionFlag<<std::endl;}
    else if(string(parname)==string("TIMEWINDOW_NS")){TimeWindow=std::atoi(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<TimeWindow<<std::endl;}
    else if(string(parname)==string("PRIMARYTHCAPGAMNORM")){primGamNormFactor=std::atof(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<primGamNormFactor<<std::endl;}
    else if(string(parname)==string("PRIMARYGAMMASECUT")){primGamEcut=std::atof(argv[i_firstpar+2*i+1]);  cout<<"      "<<parname<<"  "<<primGamEcut<<std::endl;}

    else{
      cout<<" ############ ERROR: Unknown input parameter ---> "<<parname<<"  ############"<<std::endl; return 1;
    }
  }
  //--------------------------------------------------------

  
  if(MaxExcEnergy==0){MaxExcEnergy=-1.0*std::max(NeutronEnergy*1.2,0.5);}

  
  //--------------------------------------------------------
  int Z=ZA/1000;
  int A=ZA-1000*Z+1;
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
  //Select the start level for the cascades:
  double Sn,I0;
  theStatisticalNucleus->GetSnAndI0(Sn,I0);
  double ExcitationEnergy=Sn+(A-1.)/(double)A*NeutronEnergy;
  int InitialLevel=-1;
  if(!DoThermal){
    if(JSpin<0){
      if(I0<=-100){
	std::cout<<" ############ Error in "<<__FILE__<<", line "<<__LINE__<<" ############"<<std::endl; exit(1);
      }
      JSpin=std::fabs(I0)+0.5;
      parity=true; if(I0<0){parity=false;}
    }
    InitialLevel=theStatisticalNucleus->GetClosestLevel(ExcitationEnergy,(int)(JSpin*2+0.1),parity);
    if(InitialLevel<0){ //there are no levels with the requested spin and parity
      InitialLevel=theStatisticalNucleus->GetClosestLevel(ExcitationEnergy,-1,parity);
    }
  }
  Level* iLevel=theStatisticalNucleus->GetLevel(InitialLevel);
  if(!DoThermal && (InitialLevel<0 || iLevel==0)){
    std::cout<<" ############ Error in "<<__FILE__<<", line "<<__LINE__<<" ############"<<std::endl; exit(1);
  }
  JSpin=iLevel->spinx2/2.;
  parity=iLevel->parity;
  theStatisticalNucleus->ChangeLevelSpinParityAndBR(InitialLevel,(int)(JSpin*2+0.1),parity,1,-1,seed4);
  
  int writepar=+1; if(!iLevel->parity){writepar=-1;}
  std::cout<<" Capture cascades starting in level i_level="<<InitialLevel<<", with E="<<iLevel->Energy<<" MeV, spin="<<iLevel->spinx2/2.<<", and parity="<<writepar<<std::endl;
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
  //Write input:
  std::ofstream outi(outfname_inp);
  outi<<std::endl;
  outi<<"LIBDIR "<<LibDir<<std::endl;
  outi<<"ZA "<<ZA<<std::endl;
  outi<<"NCASCADES "<<NCascades<<std::endl;
  outi<<"NEUTRONENERGY_MEV "<<NeutronEnergy<<std::endl;
  outi<<"JSPIN "<<JSpin<<std::endl;
  outi<<"PARITY "<<parity<<std::endl;
  outi<<"DOTHERMALCASCADE "<<DoThermal<<std::endl;
  outi<<"TIMEWINDOW_NS "<<TimeWindow<<std::endl;
  outi<<"SEED4 "<<seed4<<std::endl;
  outi<<std::endl;
  theStatisticalNucleus->PrintInput01(outi);
  if(!outi.good()){
    std::cout<<" ############# Error opening or writting in "<<outfname_inp<<" #############"<<std::endl; exit(1);
  }
  outi.close();
  //--------------------------------------------------------

     
  //--------------------------------------------------------
  //Generate the cascades:
  int Npar,Npar2;
  std::vector<char> pType;
  std::vector<double> pEnergy,pTime;

  std::ofstream out(outfname_cas);
  if(!out.good()){
    std::cout<<" ######## Error opening "<<outfname_cas<<" ########"<<std::endl; exit(1);
  }

  for(int i=0;i<NCascades;i++){
    Npar=theStatisticalNucleus->GenerateCascade(InitialLevel,ExcitationEnergy,pType,pEnergy,pTime);
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


















