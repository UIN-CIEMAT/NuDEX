
#include "NuDEXRandom.hh"

#if COMPILATIONTYPE == 1
//==============================================================================
NuDEXRandom::NuDEXRandom(unsigned int seed){
  theRandom=new TRandom2(seed);
}
NuDEXRandom::~NuDEXRandom(){
  delete theRandom;
}
void NuDEXRandom::SetSeed(unsigned int seed){
  theRandom->SetSeed(seed);
}
unsigned int NuDEXRandom::GetSeed(){
  return theRandom->GetSeed();
}
double NuDEXRandom::Uniform(double Xmin,double Xmax){
  return theRandom->Uniform(Xmin,Xmax);
}
unsigned int NuDEXRandom::Integer(unsigned int IntegerMax){
  return theRandom->Integer(IntegerMax);
}
double NuDEXRandom::Exp(double tau){
  return theRandom->Exp(tau);
}
double NuDEXRandom::Gaus(double mean,double sigma){
  return theRandom->Gaus(mean,sigma);
}
int NuDEXRandom::Poisson(double mean){
  return theRandom->Poisson(mean);
}
//==============================================================================
void NuDEXException(const char* originOfException, const char* exceptionCode,const char* ){
  std::cout<<" ############## Error in "<<originOfException<<", line "<<exceptionCode<<" ##############"<<std::endl; exit(1);
}
//==============================================================================

#elif COMPILATIONTYPE == 2
//==============================================================================
NuDEXRandom::NuDEXRandom(unsigned int seed){
  theEngine=new CLHEP::HepJamesRandom(seed);
  theRandFlat=new CLHEP::RandFlat(theEngine);
  theRandExponential=new CLHEP::RandExponential(theEngine);
  theRandGauss=new CLHEP::RandGauss(theEngine);
  theRandPoisson=new CLHEP::RandPoisson(theEngine);
}
NuDEXRandom::~NuDEXRandom(){

  //delete theRandFlat;
  //delete theRandExponential;
  //delete theRandGauss;
  //delete theRandPoisson;
  //delete theEngine;

}
void NuDEXRandom::SetSeed(unsigned int seed){
  theEngine->setSeed(seed);
  theRandGauss->setF(false);
}
unsigned int NuDEXRandom::GetSeed(){
  return (unsigned int)theEngine->getSeed();
}
double NuDEXRandom::Uniform(double Xmin,double Xmax){
  return theRandFlat->fire(Xmin,Xmax);
}
unsigned int NuDEXRandom::Integer(unsigned int IntegerMax){
  return theRandFlat->fireInt(IntegerMax); //bikerful!!!
}
double NuDEXRandom::Exp(double tau){
  return theRandExponential->fire(tau);
}
double NuDEXRandom::Gaus(double mean,double sigma){
  return theRandGauss->fire(mean,sigma);
}
int NuDEXRandom::Poisson(double mean){
  return theRandPoisson->fire(mean);
}
//==============================================================================
void NuDEXException(const char* originOfException, const char* exceptionCode,const char* description){
  G4Exception(originOfException,exceptionCode,FatalException,description);
}
//==============================================================================
#endif
