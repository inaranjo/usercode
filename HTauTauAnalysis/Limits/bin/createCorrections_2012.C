#include "TFile.h"
#include "TF1.h"
#include "TSystem.h"

#include "ratioEfficiencyElec.C"
#include "ratioEfficiencyMu.C"
#include "ratioEfficiencyTest.C"


///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
Double_t myFuncRatioElecIDBL(Double_t* x, Double_t *par) {
  double ratio = 1.;
  Float_t xx = x[0];
  //2012
  if( xx < 20 ) 
    ratio = 1.00;
  else if( xx >= 20 && xx< 30) 
    ratio = 0.922 ;
  else if( xx >= 30)
    ratio = 0.964;
 return ratio;
}
Double_t myFuncRatioElecIDEC(Double_t* x, Double_t *par) {
  double ratio = 1.;
  Float_t xx = x[0];
  //2012
  if( xx < 20 ) 
    ratio = 1.00;
  else if( xx >= 20 && xx< 30) 
    ratio = 0.944 ;
  else if( xx >= 30)
    ratio = 0.958;
  return ratio;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
Double_t myFuncRatioElecIsoBL(Double_t* x, Double_t *par) {
   double ratio = 1.;
  Float_t xx = x[0];
  //2012
  if( xx < 20 ) 
    ratio = 1.00;
  else if( xx >= 20 && xx< 30) 
    ratio = 0.974 ;
  else if( xx >= 30)
    ratio = 0.997;
  return ratio;

}

Double_t myFuncRatioElecIsoEC(Double_t* x, Double_t *par) {
  double ratio = 1.;
  Float_t xx = x[0];
  //2012
  if( xx < 20 ) 
    ratio = 1.00;
  else if( xx >= 20 && xx< 30) 
    ratio = 1.008 ;
  else if( xx >= 30)
    ratio = 0.983;
  return ratio;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
Double_t myFuncRatioElecIDIsoBL(Double_t* x, Double_t *par) {
  double ratio = 1.;
  Float_t xx = x[0];
  //2012
  if( xx < 20 ) 
    ratio = 1.00;
  else if( xx >= 20 && xx< 30) 
    ratio = 0.922*0.974 ;
  else if( xx >= 30)
    ratio = 0.964*0.997;
  return ratio;
}

Double_t myFuncRatioElecIDIsoEC(Double_t* x, Double_t *par) {
  double ratio = 1.;
  Float_t xx = x[0];
  //2012
  if( xx < 20 ) 
    ratio = 1.00;
  else if( xx >= 20 && xx< 30) 
    ratio = 0.944*1.008 ;
  else if( xx >= 30)
    ratio = 0.958*0.983;
  return ratio;
}


///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
Double_t myFuncRatioEle20BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle20EB = new ratioEfficiencyTest(20.97643939,1.15196354,2.27544602,1.01743868,2.04391816);

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec20MC = new ratioEfficiencyElec();

  Float_t xx = x[0];
  if (ratioEffElec20MC->mcEfficiency(xx, true)<=0) return -1;
  else return fitEffEle20EB->turnOn(xx)/ratioEffElec20MC->mcEfficiency(xx, true);
}

Double_t myFuncTurnOnEle20BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle20EB = new ratioEfficiencyTest(20.97643939,1.15196354,2.27544602,1.01743868,2.04391816 );

  Float_t xx = x[0];
  return fitEffEle20EB->turnOn(xx);
}

Double_t myFuncTurnOnEle20MCBL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec20MC = new ratioEfficiencyElec();

  Float_t xx = x[0];
  return ratioEffElec20MC->mcEfficiency(xx, true);
}

Double_t myFuncRatioEle20EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle20EC = new ratioEfficiencyTest(20.59874300,1.25425435,1.61098921,1.00146962,60.35067579);

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec20MC = new ratioEfficiencyElec();

  Float_t xx = x[0];
  if(ratioEffElec20MC->mcEfficiency(xx, false)<=0) return -1;
  else return fitEffEle20EC->turnOn(xx)/ratioEffElec20MC->mcEfficiency(xx, false);
}

Double_t myFuncTurnOnEle20EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle20EC = new ratioEfficiencyTest(20.59874300,1.25425435,1.61098921,1.00146962,60.35067579);

  Float_t xx = x[0];
  return fitEffEle20EC->turnOn(xx);
}

Double_t myFuncTurnOnEle20MCEC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec20MC = new ratioEfficiencyElec();

  Float_t xx = x[0];
  return ratioEffElec20MC->mcEfficiency(xx, false);
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

Double_t myFuncRatioEle22BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle22EB = new ratioEfficiencyTest(22.90752344,1.32376429,2.17813319,1.03674051,2.15454768);

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec20MC = new ratioEfficiencyElec();

  Float_t xx = x[0];
  if(ratioEffElec20MC->mcEfficiency(xx, true)<=0) return -1;
  else return fitEffEle22EB->turnOn(xx)/ratioEffElec20MC->mcEfficiency(xx, true);
}

Double_t myFuncTurnOnEle22BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle22EB = new ratioEfficiencyTest(22.90752344,1.32376429,2.17813319,1.03674051,2.15454768);

  Float_t xx = x[0];
  return fitEffEle22EB->turnOn(xx);
}

Double_t myFuncRatioEle22EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle22EC = new ratioEfficiencyTest(22.14553261,1.19913124,1.75642067,1.00826962,9.04331617);

  gSystem->Load("ratioEfficiencyElec_C.so");
  ratioEfficiencyElec* ratioEffElec20MC = new ratioEfficiencyElec();

  Float_t xx = x[0];
  if(ratioEffElec20MC->mcEfficiency(xx, false)) return -1;
  else return fitEffEle22EC->turnOn(xx)/ratioEffElec20MC->mcEfficiency(xx, false);
}

Double_t myFuncTurnOnEle22EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffEle22EC = new ratioEfficiencyTest(22.14553261,1.19913124,1.75642067,1.00826962,9.04331617);

  Float_t xx = x[0];
  return fitEffEle22EC->turnOn(xx);
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

Double_t myFuncRatioEleAllBL(Double_t* x, Double_t *par) {

  ratioEfficiencyTest* fitEffEle20BL = new ratioEfficiencyTest(20.97643939,1.15196354,2.27544602,1.01743868,2.04391816);
  ratioEfficiencyTest* fitEffEle22BL = new ratioEfficiencyTest(22.90752344,1.32376429,2.17813319,1.03674051,2.15454768);

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffElec20MCBL = new ratioEfficiencyTest(20.58604584,-1.89456806,3.69311772,1.05480046,1.28655181);

  float weightElec20 =  692.;
  float weightElec22 =  4390.;

  float total = weightElec20+weightElec22;
  weightElec20/=total;
  weightElec22/=total;

  Float_t xx = x[0];

  if(fitEffElec20MCBL->turnOn(xx)<=10e-3) return -1;
  else return (fitEffEle20BL->turnOn(xx) * weightElec20+
	       fitEffEle22BL->turnOn(xx) * weightElec22
	       )/fitEffElec20MCBL->turnOn(xx);
}

Double_t myFuncTurnOnEleAllBL(Double_t* x, Double_t *par) {

  ratioEfficiencyTest* fitEffEle20BL = new ratioEfficiencyTest(20.97643939,1.15196354,2.27544602,1.01743868,2.04391816);
  ratioEfficiencyTest* fitEffEle22BL = new ratioEfficiencyTest(22.90752344,1.32376429,2.17813319,1.03674051,2.15454768);

  float weightElec20 =  692.;
  float weightElec22 =  4390.;
  
  float total = weightElec20+weightElec22;
  weightElec20/=total;
  weightElec22/=total;

  Float_t xx = x[0];

  return (fitEffEle20BL->turnOn(xx) * weightElec20+
	  fitEffEle22BL->turnOn(xx) * weightElec22
	  );
}

Double_t myFuncRatioEleAllEC(Double_t* x, Double_t *par) {

  ratioEfficiencyTest* fitEffEle20EC = new ratioEfficiencyTest(20.59874300,1.25425435,1.61098921,1.00146962,60.35067579);
  ratioEfficiencyTest* fitEffEle22EC = new ratioEfficiencyTest(22.14553261,1.19913124,1.75642067,1.00826962,9.04331617);

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffElec20MCEC = new ratioEfficiencyTest(20.15425918,0.75449122,1.06027513,1.01106686,7.01956561);

  float weightElec20 =  692.;
  float weightElec22 =  4390.;

  float total = weightElec20+weightElec22;
  weightElec20/=total;
  weightElec22/=total;

  Float_t xx = x[0];
 
  if(fitEffElec20MCEC->turnOn(xx)<=10e-3) return -1;
  else return (fitEffEle20EC->turnOn(xx) * weightElec20 +
	       fitEffEle22EC->turnOn(xx) * weightElec22
	       )/fitEffElec20MCEC->turnOn(xx);
}

Double_t myFuncTurnOnEleAllEC(Double_t* x, Double_t *par) {

  ratioEfficiencyTest* fitEffEle20EC = new ratioEfficiencyTest(20.59874300,1.25425435,1.61098921,1.00146962,60.35067579);
  ratioEfficiencyTest* fitEffEle22EC = new ratioEfficiencyTest(22.14553261,1.19913124,1.75642067,1.00826962,9.04331617);

  float weightElec20 =  692.;
  float weightElec22 =  4390.;

  float total = weightElec20+weightElec22;
  weightElec20/=total;
  weightElec22/=total;

  Float_t xx = x[0];

  return (fitEffEle20EC->turnOn(xx) * weightElec20 +
	  fitEffEle22EC->turnOn(xx) * weightElec22
	  );
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
Double_t myFuncRatioMuIDBL(Double_t* x, Double_t *par) {
  double ratio = 1.;
  Float_t xx = x[0];
  //2012
  if( xx < 15 ) 
    ratio = 1.0;
  else if( xx >= 15 && xx< 20) 
    ratio = 0.989;
  else if( xx >= 20 && xx< 30)
    ratio = 0.991;
  else if( xx >= 30 )
    ratio = 0.989;
  return ratio;
}

Double_t myFuncRatioMuIDEC(Double_t* x, Double_t *par) {
  double ratio = 1.;
  Float_t xx = x[0];
  //2012
  if( xx < 15 ) 
    ratio = 1.0;
  else if( xx >= 15 && xx< 20) 
    ratio = 0.977;
  else if( xx >= 20 && xx< 30)
    ratio = 0.974;
  else if( xx >= 30 )
    ratio = 0.989;
  return ratio;
}


///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
Double_t myFuncRatioMuIsoBL(Double_t* x, Double_t *par) {
  double ratio = 1.;
  Float_t xx = x[0];
  //2012
  if( xx < 15 ) 
    ratio = 1.0;
  else if( xx >= 15 && xx< 20) 
    ratio = 0.945;
  else if( xx >= 20 && xx< 30)
    ratio = 1.005;
  else if( xx >= 30 )
    ratio = 0.993;
  return ratio;
}

Double_t myFuncRatioMuIsoEC(Double_t* x, Double_t *par) {
  double ratio = 1.;
  Float_t xx = x[0];
  //2012
  if( xx < 15 ) 
    ratio = 1.0;
  else if( xx >= 15 && xx< 20) 
    ratio = 1.047;
  else if( xx >= 20 && xx< 30)
    ratio = 0.992;
  else if( xx >= 30 )
    ratio = 1.005;
  return ratio;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
Double_t myFuncRatioMuIDIsoBL(Double_t* x, Double_t *par) {
  double ratio = 1.;
  Float_t xx = x[0];
  //2012
  if( xx < 15 ) 
    ratio = 1.0;
  else if( xx >= 15 && xx< 20) 
    ratio = 0.989*0.945;
  else if( xx >= 20 && xx< 30)
    ratio = 0.991*1.005;
  else if( xx >= 30 )
    ratio = 0.989*0.993;
  return ratio;
}

Double_t myFuncRatioMuIDIsoEC(Double_t* x, Double_t *par) {
  double ratio = 1.;
  Float_t xx = x[0];
  //2012
  if( xx < 15 ) 
    ratio = 1.0;
  else if( xx >= 15 && xx< 20) 
    ratio = 0.977*1.047;
  else if( xx >= 20 && xx< 30)
    ratio = 0.974*0.992;
  else if( xx >= 30 )
    ratio = 0.989*1.005;
  return ratio;
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
Double_t myFuncRatioMu18BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffMu18EB = new ratioEfficiencyTest(15.99983195,-0.39072829,0.28256338,1.72861719,0.95769408);

  gSystem->Load("ratioEfficiencyMu_C.so");
  ratioEfficiencyMu* ratioEffMu18MC = new ratioEfficiencyMu();

  Float_t xx = x[0];
  if(ratioEffMu18MC->mcEfficiency(xx, true)<=0)return -1;
  else return fitEffMu18EB->turnOn(xx)/ratioEffMu18MC->mcEfficiency(xx, true);
}

Double_t myFuncTurnOnMu18BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffMu18EB = new ratioEfficiencyTest(15.99983195,-0.39072829,0.28256338,1.72861719,0.95769408);

  Float_t xx = x[0];
  return fitEffMu18EB->turnOn(xx);
}

Double_t myFuncTurnOnMu18MCBL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyMu_C.so");
  ratioEfficiencyMu* ratioEffMu18MC = new ratioEfficiencyMu();

  Float_t xx = x[0];
  return ratioEffMu18MC->mcEfficiency(xx, true);
}

Double_t myFuncRatioMu18EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffMu18EC = new ratioEfficiencyTest(18.49754887,-0.16941614,0.26076717,1.05494469,1.53819978);

  gSystem->Load("ratioEfficiencyMu_C.so");
  ratioEfficiencyMu* ratioEffMu18MC = new ratioEfficiencyMu();

  Float_t xx = x[0];
  if(ratioEffMu18MC->mcEfficiency(xx, false)<=0)return -1;
  else return fitEffMu18EC->turnOn(xx)/ratioEffMu18MC->mcEfficiency(xx, false);
}

Double_t myFuncTurnOnMu18EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffMu18EC = new ratioEfficiencyTest(18.49754887,-0.16941614,0.26076717,1.05494469,1.53819978);

  Float_t xx = x[0];
  return fitEffMu18EC->turnOn(xx);
}

Double_t myFuncTurnOnMu18MCEC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyMu_C.so");
  ratioEfficiencyMu* ratioEffMu18MC = new ratioEfficiencyMu();

  Float_t xx = x[0];
  return ratioEffMu18MC->mcEfficiency(xx, false);
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

Double_t myFuncRatioMu17BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffMu17EB = new ratioEfficiencyTest(17.21270264,0.54997112,1.02874912,1.29646487,0.96724273);

  gSystem->Load("ratioEfficiencyMu_C.so");
  ratioEfficiencyMu* ratioEffMu17MC = new ratioEfficiencyMu();

  Float_t xx = x[0];
  if(ratioEffMu17MC->mcEfficiency(xx, true)<=0)return -1;
  else return fitEffMu17EB->turnOn(xx)/ratioEffMu17MC->mcEfficiency(xx, true);
}

Double_t myFuncTurnOnMu17BL(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffMu17EB = new ratioEfficiencyTest(17.21270264,0.54997112,1.02874912,1.29646487,0.96724273);

  Float_t xx = x[0];
  return fitEffMu17EB->turnOn(xx);
}

Double_t myFuncRatioMu17EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffMu17EC = new ratioEfficiencyTest(15.98037640,0.12062946,0.02183977,2.84751010,0.83985656);

  gSystem->Load("ratioEfficiencyMu_C.so");
  ratioEfficiencyMu* ratioEffMu17MC = new ratioEfficiencyMu();

  Float_t xx = x[0];
  if(ratioEffMu17MC->mcEfficiency(xx, false)<=0)return -1;
  else return fitEffMu17EC->turnOn(xx)/ratioEffMu17MC->mcEfficiency(xx, false);
}

Double_t myFuncTurnOnMu17EC(Double_t* x, Double_t *par) {

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffMu17EC = new ratioEfficiencyTest(15.98037640,0.12062946,0.02183977,2.84751010,0.83985656);

  Float_t xx = x[0];
  return fitEffMu17EC->turnOn(xx);
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

Double_t myFuncRatioMuAllBL(Double_t* x, Double_t *par) {

  ratioEfficiencyTest* fitEffMu18BL = new ratioEfficiencyTest(15.99983195,-0.39072829,0.28256338,1.72861719,0.95769408);
  ratioEfficiencyTest* fitEffMu17BL = new ratioEfficiencyTest(17.21270264,0.54997112,1.02874912,1.29646487,0.96724273);

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffMu18MCBL = new ratioEfficiencyTest(16.99389526,-0.04080190,0.00794730,1.60377906,0.99626161);

  float weightMu18 =  692.;
  float weightMu17 =  4390.;

  float total = weightMu18+weightMu17;
  weightMu18/=total;
  weightMu17/=total;

  Float_t xx = x[0];

  if(fitEffMu18MCBL->turnOn(xx)<=10e-3)return -1;
  else return (fitEffMu18BL->turnOn(xx) * weightMu18+
	       fitEffMu17BL->turnOn(xx) * weightMu17
	       )/fitEffMu18MCBL->turnOn(xx);
}

Double_t myFuncTurnOnMuAllBL(Double_t* x, Double_t *par) {
  
  ratioEfficiencyTest* fitEffMu18BL = new ratioEfficiencyTest(15.99983195,-0.39072829,0.28256338,1.72861719,0.95769408);
  ratioEfficiencyTest* fitEffMu17BL = new ratioEfficiencyTest(17.21270264,0.54997112,1.02874912,1.29646487,0.96724273);

  float weightMu18 =  692.;
  float weightMu17 =  4390.;

  float total = weightMu18+weightMu17;
  weightMu18/=total;
  weightMu17/=total;

  Float_t xx = x[0];

  return (fitEffMu18BL->turnOn(xx) * weightMu18+
	  fitEffMu17BL->turnOn(xx) * weightMu17
	  );
}


Double_t myFuncRatioMuAllEC(Double_t* x, Double_t *par) {

  ratioEfficiencyTest* fitEffMu18EC = new ratioEfficiencyTest(18.49754887,-0.16941614,0.26076717,1.05494469,1.53819978);
  ratioEfficiencyTest* fitEffMu17EC = new ratioEfficiencyTest(15.98037640,0.12062946,0.02183977,2.84751010,0.83985656);

  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* fitEffMu18MCEC = new ratioEfficiencyTest(16.99065795,-0.11993730,0.01384991,2.38867304,0.86552275);

  float weightMu18 =  692.;
  float weightMu17 =  4390.;

  float total = weightMu18+weightMu17;
  weightMu18/=total;
  weightMu17/=total;

  Float_t xx = x[0];

  if(fitEffMu18MCEC->turnOn(xx)<=10e-3)return -1;
  else return (fitEffMu18EC->turnOn(xx) * weightMu18+
	       fitEffMu17EC->turnOn(xx) * weightMu17
	       )/fitEffMu18MCEC->turnOn(xx);
}

Double_t myFuncTurnOnMuAllEC(Double_t* x, Double_t *par) {

  ratioEfficiencyTest* fitEffMu18EC = new ratioEfficiencyTest(18.49754887,-0.16941614,0.26076717,1.05494469,1.53819978);
  ratioEfficiencyTest* fitEffMu17EC = new ratioEfficiencyTest(15.98037640,0.12062946,0.02183977,2.84751010,0.83985656);

  float weightMu18 =  692.;
  float weightMu17 =  4390.;

  float total = weightMu18+weightMu17;
  weightMu18/=total;
  weightMu17/=total;

  Float_t xx = x[0];

  return (fitEffMu18EC->turnOn(xx) * weightMu18+
	  fitEffMu17EC->turnOn(xx) * weightMu17
	  );
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
// HLT PFTauLoose20 / MC Summer12 / elec+tau
Double_t myFuncTurnOnTauLoose20ElecTauMC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MC = new ratioEfficiencyTest(18.77448606,0.45765507,0.26077509,13.43372485,0.88037836);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20MC->turnOn(xx);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
// HLT PFTauLoose20 / Run2012A    / elec+tau

Double_t myFuncRatioTauLoose20ElecTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA = new ratioEfficiencyTest(18.84658959,0.25958704,0.17300958,2.43491208,0.85872017);
  
  TF1* turnOnTauLoose20ElecTauMC = new TF1("turnOnTauLoose20ElecTauMC", myFuncTurnOnTauLoose20ElecTauMC,0,400,0);

  Float_t xx = x[0];
  if(turnOnTauLoose20ElecTauMC->Eval(xx)<=0)return -1;
  else return ratioEffTauLoose20ElecTauRunA->turnOn(xx)/ turnOnTauLoose20ElecTauMC->Eval(xx);

}

Double_t myFuncTurnOnTauLoose20ElecTauRunA(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA = new ratioEfficiencyTest(18.84658959,0.25958704,0.17300958,2.43491208,0.85872017); 

  Float_t xx = x[0];

  return ratioEffTauLoose20ElecTauRunA->turnOn(xx);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
// HLT PFTauLoose20 / Run2012B    / elec+tau

Double_t myFuncRatioTauLoose20ElecTauRunB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunB = new ratioEfficiencyTest(18.48663118,1.63417147,20.25695815,138.55422224,0.89456038);
  
  TF1* turnOnTauLoose20ElecTauMC = new TF1("turnOnTauLoose20ElecTauMC", myFuncTurnOnTauLoose20ElecTauMC,0,400,0);

  Float_t xx = x[0];
  if(turnOnTauLoose20ElecTauMC->Eval(xx)<=0)return -1;
  else return ratioEffTauLoose20ElecTauRunB->turnOn(xx)/ turnOnTauLoose20ElecTauMC->Eval(xx);

}

Double_t myFuncTurnOnTauLoose20ElecTauRunB(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunB = new ratioEfficiencyTest(18.48663118,1.63417147,20.25695815,138.55422224,0.89456038); 

  Float_t xx = x[0];

  return ratioEffTauLoose20ElecTauRunB->turnOn(xx);

}

Double_t myFuncRatioTauElecTauAll(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
 
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA  = new ratioEfficiencyTest(18.84658959,0.25958704,0.17300958,2.43491208,0.85872017);
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunB  = new ratioEfficiencyTest(18.48663118,1.63417147,20.25695815,138.55422224,0.89456038);
  TF1* turnOnTauLoose20ElecTauMC = new TF1("turnOnTauLoose20ElecTauMC", myFuncTurnOnTauLoose20ElecTauMC,0,400,0);

  float weightA  =  692.;
  float weightB  =  4390.;

  float total = weightA+weightB;
  weightA /= total;
  weightB /= total;

  Float_t xx = x[0];

  float ratio = turnOnTauLoose20ElecTauMC->Eval(xx) > 0 ? 
    (ratioEffTauLoose20ElecTauRunA->turnOn(xx)  * weightA +
     ratioEffTauLoose20ElecTauRunB->turnOn(xx) * weightB
     ) / turnOnTauLoose20ElecTauMC->Eval(xx) : 0.;
  
  return ratio;
}

Double_t myFuncTurnOnTauElecTauAll(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
 
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunA  = new ratioEfficiencyTest(18.84658959,0.25958704,0.17300958,2.43491208,0.85872017);
  ratioEfficiencyTest* ratioEffTauLoose20ElecTauRunB  = new ratioEfficiencyTest(18.48663118,1.63417147,20.25695815,138.55422224,0.89456038);

  float weightA  =  692.;
  float weightB  =  4390.;

  float total = weightA+weightB;
  weightA /= total;
  weightB /= total;

  Float_t xx = x[0];

  return (ratioEffTauLoose20ElecTauRunA->turnOn(xx)  * weightA +
	  ratioEffTauLoose20ElecTauRunB->turnOn(xx) * weightB
	  ) ;
  }
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
////////////////////////////////////////////////////////
// HLT PFTauLoose20 / MC Summer12 / mu+tau
Double_t myFuncTurnOnTauLoose20MuTauMCBL(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MCBL = new ratioEfficiencyTest(18.86257072,0.25680380,0.16916101,2.42931257,0.89590264);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20MCBL->turnOn(xx);
}
Double_t myFuncTurnOnTauLoose20MuTauMCEC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MCEC = new ratioEfficiencyTest(18.74764561,1.82036845,701.46994969,101.57913480,0.82547043);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20MCEC->turnOn(xx);
}

////////////////////////////////////////////////////////
// HLT PFTauLoose20 / Run2012A    / mu+tau

Double_t myFuncRatioTauLoose20MuTauRunABL(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunABL = new ratioEfficiencyTest(18.52262128,1.85879597,3.48843815,1.15491294,1.02489024);
  
  TF1* turnOnTauLoose20MuTauMCBL = new TF1("turnOnTauLoose20MuTauMCBL", myFuncTurnOnTauLoose20MuTauMCBL,0,400,0);

  Float_t xx = x[0];
  if(turnOnTauLoose20MuTauMCBL->Eval(xx)<=0)return -1;
  else return ratioEffTauLoose20MuTauRunABL->turnOn(xx)/ turnOnTauLoose20MuTauMCBL->Eval(xx);
}

Double_t myFuncTurnOnTauLoose20MuTauRunABL(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunABL = new ratioEfficiencyTest(18.52262128,1.85879597,3.48843815,1.15491294,1.02489024); 

  Float_t xx = x[0];

  return ratioEffTauLoose20MuTauRunABL->turnOn(xx);
}

Double_t myFuncRatioTauLoose20MuTauRunAEC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunAEC = new ratioEfficiencyTest(18.90119559,0.14025596,0.14482632,1.56126508,0.81188198);
  
  TF1* turnOnTauLoose20MuTauMCEC = new TF1("turnOnTauLoose20MuTauMCEC", myFuncTurnOnTauLoose20MuTauMCEC,0,400,0);

  Float_t xx = x[0];
  if(turnOnTauLoose20MuTauMCEC->Eval(xx)<=0)return -1;
  else return ratioEffTauLoose20MuTauRunAEC->turnOn(xx)/ turnOnTauLoose20MuTauMCEC->Eval(xx);
}

Double_t myFuncTurnOnTauLoose20MuTauRunAEC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunAEC = new ratioEfficiencyTest(18.90119559,0.14025596,0.14482632,1.56126508,0.81188198);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20MuTauRunAEC->turnOn(xx);
}

////////////////////////////////////////////////////////
// HLT PFTauLoose20 / Run2012B    / mu+tau

Double_t myFuncRatioTauLoose20MuTauRunBBL(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBBL = new ratioEfficiencyTest(17.92648563,1.96846742,4.46406075,1.0223992,1.52260575);
  
  TF1* turnOnTauLoose20MuTauMCBL = new TF1("turnOnTauLoose20MuTauMCBL", myFuncTurnOnTauLoose20MuTauMCBL,0,400,0);

  Float_t xx = x[0];
  if(turnOnTauLoose20MuTauMCBL->Eval(xx)<=0)return -1;
  else return ratioEffTauLoose20MuTauRunBBL->turnOn(xx)/ turnOnTauLoose20MuTauMCBL->Eval(xx);
}

Double_t myFuncTurnOnTauLoose20MuTauRunBBL(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBBL = new ratioEfficiencyTest(17.92648563,1.96846742,4.46406075,1.0223992,1.52260575); 

  Float_t xx = x[0];

  return ratioEffTauLoose20MuTauRunBBL->turnOn(xx);

}

Double_t myFuncRatioTauLoose20MuTauRunBEC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBEC = new ratioEfficiencyTest(18.59856420,2.49132550,10.99643595,1.50651123,0.87952970);
  
  TF1* turnOnTauLoose20MuTauMCEC = new TF1("turnOnTauLoose20MuTauMCEC", myFuncTurnOnTauLoose20MuTauMCEC,0,400,0);

  Float_t xx = x[0];
  if(turnOnTauLoose20MuTauMCEC->Eval(xx)<=0)return -1;
  else return ratioEffTauLoose20MuTauRunBEC->turnOn(xx)/ turnOnTauLoose20MuTauMCEC->Eval(xx);
}

Double_t myFuncTurnOnTauLoose20MuTauRunBEC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBEC = new ratioEfficiencyTest(18.59856420,2.49132550,10.99643595,1.50651123,0.87952970);
  
  Float_t xx = x[0];

  return ratioEffTauLoose20MuTauRunBEC->turnOn(xx);
}

////////////////////////////////////////////////////////
// HLT PFTauLoose20 / Run2012    / mu+tau
Double_t myFuncRatioTauMuTauAllBL(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunABL  = new ratioEfficiencyTest(18.52262128,1.85879597,3.48843815,1.15491294,1.02489024);
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBBL  = new ratioEfficiencyTest(17.92648563,1.96846742,4.46406075,1.0223992,1.52260575);
  TF1* turnOnTauLoose20MuTauMCBL = new TF1("turnOnTauLoose20MuTauMC", myFuncTurnOnTauLoose20MuTauMCBL,0,400,0);

  float weightA  =  692.;
  float weightB  =  4390.;

  float total = weightA+weightB;
  weightA /= total;
  weightB /= total;

  Float_t xx = x[0];
  
  float ratio = turnOnTauLoose20MuTauMCBL->Eval(xx) > 0 ?
    (ratioEffTauLoose20MuTauRunABL->turnOn(xx)  * weightA +
     ratioEffTauLoose20MuTauRunBBL->turnOn(xx) * weightB
     )/turnOnTauLoose20MuTauMCBL->Eval(xx) : 0;
  return ratio;

}
Double_t myFuncTurnOnTauMuTauAllBL(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunABL  = new ratioEfficiencyTest(18.52262128,1.85879597,3.48843815,1.15491294,1.02489024);
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBBL  = new ratioEfficiencyTest(17.92648563,1.96846742,4.46406075,1.0223992,1.52260575);

  float weightA  =  692.;
  float weightB  =  4390.;

  float total = weightA+weightB;
  weightA /= total;
  weightB /= total;

  Float_t xx = x[0];

  return (ratioEffTauLoose20MuTauRunABL->turnOn(xx)  * weightA +
	  ratioEffTauLoose20MuTauRunBBL->turnOn(xx) * weightB
	  );
}

Double_t myFuncRatioTauMuTauAllEC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunAEC  = new ratioEfficiencyTest(18.90119559,0.14025596,0.14482632,1.56126508,0.81188198);
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBEC  = new ratioEfficiencyTest(18.59856420,2.49132550,10.99643595,1.50651123,0.87952970);
  TF1* turnOnTauLoose20MuTauMCEC = new TF1("turnOnTauLoose20MuTauMC", myFuncTurnOnTauLoose20MuTauMCEC,0,400,0);

  float weightA  =  692.;
  float weightB  =  4390.;

  float total = weightA+weightB;
  weightA /= total;
  weightB /= total;

  Float_t xx = x[0];

  float ratio = turnOnTauLoose20MuTauMCEC->Eval(xx) > 0 ?
    (ratioEffTauLoose20MuTauRunAEC->turnOn(xx)  * weightA +
     ratioEffTauLoose20MuTauRunBEC->turnOn(xx) * weightB
     )/turnOnTauLoose20MuTauMCEC->Eval(xx) : 0;
  return ratio;

}

Double_t myFuncTurnOnTauMuTauAllEC(Double_t* x, Double_t *par) {
  
  gSystem->Load("ratioEfficiencyTest_C.so");
  
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunAEC  = new ratioEfficiencyTest(18.90119559,0.14025596,0.14482632,1.56126508,0.81188198);
  ratioEfficiencyTest* ratioEffTauLoose20MuTauRunBEC  = new ratioEfficiencyTest(18.59856420,2.49132550,10.99643595,1.50651123,0.87952970);

  float weightA  =  692.;
  float weightB  =  4390.;

  float total = weightA+weightB;
  weightA /= total;
  weightB /= total;

  Float_t xx = x[0];

  return (ratioEffTauLoose20MuTauRunAEC->turnOn(xx)  * weightA +
	  ratioEffTauLoose20MuTauRunBEC->turnOn(xx) * weightB
	  );
}

/////////////////////////////////////////////////

void makeFile(){

  TFile* fout = new TFile("Corrections2012.root","RECREATE");

  TF1 *ratioElec20BL        = new TF1("ratioElec20BL",           myFuncRatioEle20BL,        15,800,0);
  TF1 *turnOnElec20BL       = new TF1("turnOnElec20BL",          myFuncTurnOnEle20BL,       15,800,0);
  TF1 *ratioElec20EC        = new TF1("ratioElec20EC",           myFuncRatioEle20EC,        15,800,0);
  TF1 *turnOnElec20EC       = new TF1("turnOnElec20EC",          myFuncTurnOnEle20EC,       15,800,0);
  TF1 *turnOnElec20MCBL     = new TF1("turnOnElec20MCBL",        myFuncTurnOnEle20MCBL,     15,800,0);
  TF1 *turnOnElec20MCEC     = new TF1("turnOnElec20MCEC",        myFuncTurnOnEle20MCEC,     15,800,0);

  TF1 *ratioElec22BL        = new TF1("ratioElec22BL",           myFuncRatioEle22BL,        15,800,0);
  TF1 *turnOnElec22BL       = new TF1("turnOnElec22BL",          myFuncTurnOnEle22BL,       15,800,0);
  TF1 *ratioElec22EC        = new TF1("ratioElec22EC",           myFuncRatioEle22EC,        15,800,0);
  TF1 *turnOnElec22EC       = new TF1("turnOnElec22EC",          myFuncTurnOnEle22EC,       15,800,0);

  TF1 *ratioElecAllBL       = new TF1("ratioElecAllBL",          myFuncRatioEleAllBL,       15,800,0);
  TF1 *turnOnElecAllBL      = new TF1("turnOnElecAllBL",         myFuncTurnOnEleAllBL,      15,800,0);
  TF1 *ratioElecAllEC       = new TF1("ratioElecAllEC",          myFuncRatioEleAllEC,       15,800,0);
  TF1 *turnOnElecAllEC      = new TF1("turnOnElecAllEC",         myFuncTurnOnEleAllEC,      15,800,0);

  TF1 *ratioElecIDBL        = new TF1("ratioElecIDBL",           myFuncRatioElecIDBL ,      15,800,0);
  TF1 *ratioElecIDEC        = new TF1("ratioElecIDEC",           myFuncRatioElecIDEC ,      15,800,0);
  TF1 *ratioElecIsoBL       = new TF1("ratioElecIsoBL",          myFuncRatioElecIsoBL ,     15,800,0);
  TF1 *ratioElecIsoEC       = new TF1("ratioElecIsoEC",          myFuncRatioElecIsoEC ,     15,800,0);
  TF1 *ratioElecIDIsoBL     = new TF1("ratioElecIDIsoBL",        myFuncRatioElecIDIsoBL ,   15,800,0);
  TF1 *ratioElecIDIsoEC     = new TF1("ratioElecIDIsoEC",        myFuncRatioElecIDIsoEC ,   15,800,0);

  TF1 *ratioMu18BL          = new TF1("ratioMu18BL",             myFuncRatioMu18BL,         15,800,0);
  TF1 *turnOnMu18BL         = new TF1("turnOnMu18BL",            myFuncTurnOnMu18BL,        15,800,0);
  TF1 *ratioMu18EC          = new TF1("ratioMu18EC",             myFuncRatioMu18EC,         15,800,0);
  TF1 *turnOnMu18EC         = new TF1("turnOnMu18EC",            myFuncTurnOnMu18EC,        15,800,0);
  TF1 *turnOnMu18MCBL       = new TF1("turnOnMu18MCBL",          myFuncTurnOnMu18MCBL,      15,800,0);
  TF1 *turnOnMu18MCEC       = new TF1("turnOnMu18MCEC",          myFuncTurnOnMu18MCEC,      15,800,0);

  TF1 *ratioMu17BL          = new TF1("ratioMu17BL",             myFuncRatioMu17BL,         15,800,0);
  TF1 *turnOnMu17BL         = new TF1("turnOnMu17BL",            myFuncTurnOnMu17BL,        15,800,0);
  TF1 *ratioMu17EC          = new TF1("ratioMu17EC",             myFuncRatioMu17EC,         15,800,0);
  TF1 *turnOnMu17EC         = new TF1("turnOnMu17EC",            myFuncTurnOnMu17EC,        15,800,0);

  TF1 *ratioMuAllBL         = new TF1("ratioMuAllBL",            myFuncRatioMuAllBL,        15,800,0);
  TF1 *turnOnMuAllBL        = new TF1("turnOnMuAllBL",           myFuncTurnOnMuAllBL,       15,800,0);
  TF1 *ratioMuAllEC         = new TF1("ratioMuAllEC",            myFuncRatioMuAllEC,        15,800,0);
  TF1 *turnOnMuAllEC        = new TF1("turnOnMuAllEC",           myFuncTurnOnMuAllEC,       15,800,0);
 
  TF1 *ratioMuIdBL          = new TF1("ratioMuIdBL",             myFuncRatioMuIDBL ,        15,800,0);
  TF1 *ratioMuIdEC          = new TF1("ratioMuIdEC",             myFuncRatioMuIDEC ,        15,800,0);
  TF1 *ratioMuIsoBL         = new TF1("ratioMuIsoBL",            myFuncRatioMuIsoBL ,       15,800,0);
  TF1 *ratioMuIsoEC         = new TF1("ratioMuIsoEC",            myFuncRatioMuIsoEC ,       15,800,0);
  TF1 *ratioMuIDIsoBL       = new TF1("ratioMuIDIsoBL",          myFuncRatioMuIDIsoBL ,     15,800,0);
  TF1 *ratioMuIDIsoEC       = new TF1("ratioMuIDIsoEC",          myFuncRatioMuIDIsoEC ,     15,800,0);

  TF1 *turnOnTauLoose20ElecTauMC      = new TF1("turnOnTauLoose20ElecTauMC",       myFuncTurnOnTauLoose20ElecTauMC,      18,800,0);
  TF1 *turnOnTauLoose20ElecTauRunA      = new TF1("turnOnTauLoose20ElecTauRunA",       myFuncTurnOnTauLoose20ElecTauRunA,      18,800,0);
  TF1 *turnOnTauLoose20ElecTauRunB      = new TF1("turnOnTauLoose20ElecTauRunB",       myFuncTurnOnTauLoose20ElecTauRunB,      18,800,0);
  TF1 *ratioTauLoose20ElecTauRunA      = new TF1("ratioTauLoose20ElecTauRunA",       myFuncRatioTauLoose20ElecTauRunA,      18,800,0);
  TF1 *ratioTauLoose20ElecTauRunB      = new TF1("ratioTauLoose20ElecTauRunB",       myFuncRatioTauLoose20ElecTauRunB,      18,800,0);

  TF1 *ratioTauElecTauAll    = new TF1("ratioTauElecTauAll", myFuncRatioTauElecTauAll ,19.5,800,0);
  TF1 *turnOnTauElecTauAll   = new TF1("turnOnTauElecTauAll",myFuncTurnOnTauElecTauAll,19.5,800,0);
 
  TF1 *turnOnTauLoose20MuTauMCBL      = new TF1("turnOnTauLoose20MuTauMCBL",       myFuncTurnOnTauLoose20MuTauMCBL,      18,800,0);
  TF1 *turnOnTauLoose20MuTauRunABL      = new TF1("turnOnTauLoose20MuTauRunABL",       myFuncTurnOnTauLoose20MuTauRunABL,      18,800,0);
  TF1 *turnOnTauLoose20MuTauRunBBL      = new TF1("turnOnTauLoose20MuTauRunBBL",       myFuncTurnOnTauLoose20MuTauRunBBL,      18,800,0);
  TF1 *turnOnTauLoose20MuTauMCEC      = new TF1("turnOnTauLoose20MuTauMCEC",       myFuncTurnOnTauLoose20MuTauMCEC,      18,800,0);
  TF1 *turnOnTauLoose20MuTauRunAEC      = new TF1("turnOnTauLoose20MuTauRunAEC",       myFuncTurnOnTauLoose20MuTauRunAEC,      18,800,0);
  TF1 *turnOnTauLoose20MuTauRunBEC      = new TF1("turnOnTauLoose20MuTauRunBEC",       myFuncTurnOnTauLoose20MuTauRunBEC,      18,800,0);
  TF1 *ratioTauLoose20MuTauRunABL      = new TF1("ratioTauLoose20MuTauRunABL",       myFuncRatioTauLoose20MuTauRunABL,      18,800,0);
  TF1 *ratioTauLoose20MuTauRunBBL      = new TF1("ratioTauLoose20MuTauRunBBL",       myFuncRatioTauLoose20MuTauRunBBL,      18,800,0);
  TF1 *ratioTauLoose20MuTauRunAEC      = new TF1("ratioTauLoose20MuTauRunAEC",       myFuncRatioTauLoose20MuTauRunAEC,      18,800,0);
  TF1 *ratioTauLoose20MuTauRunBEC      = new TF1("ratioTauLoose20MuTauRunBEC",       myFuncRatioTauLoose20MuTauRunBEC,      18,800,0);

  TF1 *ratioTauMuTauAllBL      = new TF1("ratioTauMuTauAllBL", myFuncRatioTauMuTauAllBL ,18,800,0);
  TF1 *turnOnTauMuTauAllBL     = new TF1("turnOnTauMuTauAllBL",myFuncTurnOnTauMuTauAllBL,18,800,0);
  TF1 *ratioTauMuTauAllEC      = new TF1("ratioTauMuTauAllEC", myFuncRatioTauMuTauAllEC ,18,800,0);
  TF1 *turnOnTauMuTauAllEC     = new TF1("turnOnTauMuTauAllEC",myFuncTurnOnTauMuTauAllEC,18,800,0);
 

  fout->cd();
  
  ratioElec20BL->SetNpx(3200);
  turnOnElec20BL->SetNpx(3200);
  ratioElec20EC->SetNpx(3200);
  turnOnElec20EC->SetNpx(3200);
  ratioElec22BL->SetNpx(3200);
  turnOnElec22BL->SetNpx(3200);
  ratioElec22EC->SetNpx(3200);
  turnOnElec22EC->SetNpx(3200);
  turnOnElec20MCEC->SetNpx(3200);
  turnOnElec20MCBL->SetNpx(3200);
  ratioElecAllBL->SetNpx(3200);
  turnOnElecAllBL->SetNpx(3200);
  ratioElecAllEC->SetNpx(3200);
  turnOnElecAllEC->SetNpx(3200);

  ratioElecIsoBL->SetNpx(3200);
  ratioElecIsoEC->SetNpx(3200);
  ratioMuIsoBL->SetNpx(3200);
  ratioMuIsoEC->SetNpx(3200);
  ratioElecIDIsoBL->SetNpx(6400);
  ratioElecIDIsoEC->SetNpx(6400);

  ratioMu18BL->SetNpx(3200);
  turnOnMu18BL->SetNpx(3200);
  ratioMu18EC->SetNpx(3200);
  turnOnMu18EC->SetNpx(3200);
  ratioMu17BL->SetNpx(3200);
  turnOnMu17BL->SetNpx(3200);
  ratioMu17EC->SetNpx(3200);
  turnOnMu17EC->SetNpx(3200);
  turnOnMu18MCBL->SetNpx(3200);
  turnOnMu18MCEC->SetNpx(3200);
  ratioMuAllBL->SetNpx(3200);
  turnOnMuAllBL->SetNpx(3200);
  ratioMuAllEC->SetNpx(3200);
  turnOnMuAllEC->SetNpx(3200);

  ratioMuIdBL->SetNpx(6400);
  ratioMuIdEC->SetNpx(6400);
  ratioMuIsoBL->SetNpx(6400);
  ratioMuIsoEC->SetNpx(6400);
  ratioMuIDIsoBL->SetNpx(6400);
  ratioMuIDIsoEC->SetNpx(6400);

  turnOnTauLoose20ElecTauMC->SetNpx(3200);     
  turnOnTauLoose20ElecTauRunA->SetNpx(3200);   
  turnOnTauLoose20ElecTauRunB->SetNpx(3200);   
  ratioTauLoose20ElecTauRunA->SetNpx(3200);    
  ratioTauLoose20ElecTauRunB->SetNpx(3200);  

  ratioTauElecTauAll->SetNpx(3200);
  turnOnTauElecTauAll->SetNpx(3200);

  turnOnTauLoose20MuTauMCBL->SetNpx(3200);     
  turnOnTauLoose20MuTauRunABL->SetNpx(3200);   
  turnOnTauLoose20MuTauRunBBL->SetNpx(3200);   
  turnOnTauLoose20MuTauMCEC->SetNpx(3200);     
  turnOnTauLoose20MuTauRunAEC->SetNpx(3200);   
  turnOnTauLoose20MuTauRunBEC->SetNpx(3200);  
  ratioTauLoose20MuTauRunABL->SetNpx(3200);   
  ratioTauLoose20MuTauRunBBL->SetNpx(3200);   
  ratioTauLoose20MuTauRunAEC->SetNpx(3200);    
  ratioTauLoose20MuTauRunBEC->SetNpx(3200);    

  ratioTauMuTauAllBL->SetNpx(3200);
  turnOnTauMuTauAllBL->SetNpx(3200);
  ratioTauMuTauAllEC->SetNpx(3200);
  turnOnTauMuTauAllEC->SetNpx(3200);



  ratioElec20BL->Write();
  turnOnElec20BL->Write();
  ratioElec20EC->Write();
  turnOnElec20EC->Write();
  turnOnElec20MCBL->Write();
  turnOnElec20MCEC->Write();
  ratioElec22BL->Write();
  turnOnElec22BL->Write();
  ratioElec22EC->Write();
  turnOnElec22EC->Write();
  ratioElecAllBL->Write();
  turnOnElecAllBL->Write();
  ratioElecAllEC->Write();
  turnOnElecAllEC->Write();

  ratioElecIDBL->Write();
  ratioElecIDEC->Write();
  ratioElecIsoBL->Write();
  ratioElecIsoEC->Write();
  ratioElecIDIsoBL->Write();
  ratioElecIDIsoEC->Write();

  ratioMu18BL->Write();
  turnOnMu18BL->Write();
  ratioMu18EC->Write();
  turnOnMu18EC->Write();
  ratioMu17BL->Write();
  turnOnMu17BL->Write();
  ratioMu17EC->Write();
  turnOnMu17EC->Write();
  turnOnMu18MCBL->Write();
  turnOnMu18MCEC->Write();
  ratioMuAllBL->Write();
  turnOnMuAllBL->Write();
  ratioMuAllEC->Write();
  turnOnMuAllEC->Write();

  ratioMuIdBL->Write();
  ratioMuIdEC->Write();
  ratioMuIsoBL->Write();
  ratioMuIsoEC->Write();
  ratioMuIDIsoBL->Write();
  ratioMuIDIsoEC->Write();

  turnOnTauLoose20ElecTauMC->Write();    
  turnOnTauLoose20ElecTauRunA->Write();  
  turnOnTauLoose20ElecTauRunB->Write();  
  ratioTauLoose20ElecTauRunA->Write();   
  ratioTauLoose20ElecTauRunB->Write(); 

  ratioTauElecTauAll->Write();
  turnOnTauElecTauAll->Write();

  turnOnTauLoose20MuTauMCBL->Write();    
  turnOnTauLoose20MuTauRunABL->Write();  
  turnOnTauLoose20MuTauRunBBL->Write();  
  turnOnTauLoose20MuTauMCEC->Write();    
  turnOnTauLoose20MuTauRunAEC->Write();  
  turnOnTauLoose20MuTauRunBEC->Write(); 
  ratioTauLoose20MuTauRunABL->Write();  
  ratioTauLoose20MuTauRunBBL->Write();  
  ratioTauLoose20MuTauRunAEC->Write();   
  ratioTauLoose20MuTauRunBEC->Write();   

  ratioTauMuTauAllBL->Write();
  turnOnTauMuTauAllBL->Write();
  ratioTauMuTauAllEC->Write();
  turnOnTauMuTauAllEC->Write();


  fout->Write();
  fout->Close();


}
