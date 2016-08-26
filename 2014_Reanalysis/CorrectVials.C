//This function will take the reference-corrected data and plot them, correcting
//them for the vial damage at the dose of interest.
#include <fstream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include "tinydir.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include "TBranch.h"
#include "TString.h"
#include "TH1D.h"
#include <fstream>
#include "TLeaf.h"
#include <math.h>
#include "TMath.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TSystem.h"
#include "TGraph.h"

using namespace std;

void CorrectVials(string ScintFname, string VialFname, int IDnum, int VialIDnum,
		  string Formulation, Long64_t IrradTime, string GraphOptions){
  //First get the damaged vial data.
  TFile* fvial = TFile::Open(VialFname.c_str());
  TChain* VialChain = (TChain*)fvial->Get("ChainResults");

  //I need to access the factor, the dose, the IDnum, (the timestamp later).
  double theVialFactor = -1;
  float theVialDose = -1;
  int theVialIDnum = -1;

  VialChain->SetBranchAddress("IDnum", &theVialIDnum);
  VialChain->SetBranchAddress("Dose", &theVialDose);
  VialChain->SetBranchAddress("Factor", &theVialFactor);

  Long64_t VialNumEntries = VialChain->GetEntries();

  //Make a map of factor as a function of dose
  map<double, float> DoseMapFactor;
  map<double, int> DoseMapNumEntries;
  for(int i = 0; i<VialNumEntries; i++){
    VialChain->GetEntry(i);
    if(VialIDnum==theVialIDnum){
      //Check if there is an entry at theVialDose, if there isn't add this, if
      //there is, take an average.
      if(DoseMapFactor.find(theVialDose)==DoseMapFactor.end()){
	DoseMapFactor[theVialDose] = theVialFactor;
	DoseMapNumEntries[theVialDose] = 1;
      }
      else{
	DoseMapFactor[theVialDose] =
	  ((DoseMapFactor[theVialDose]*DoseMapNumEntries[theVialDose])
	   + theVialFactor)/
	  (DoseMapNumEntries[theVialDose] + 1);
	DoseMapNumEntries[theVialDose]++;
      }
    }
  }

  //Get the scint file tree.
  TFile *fscint = TFile::Open(ScintFname.c_str());
  TChain* ScintChain = (TChain*)fscint->Get("ChainFinalResults");

  string* ScintName = NULL;
  Long64_t ScintTime = 0;
  int ScintIDnum = -100;
  float ScintDose = -1;
  string* ScintFormulation = NULL;
  double ScintCorrLY = -1;

  ScintChain->SetBranchAddress("TimeStamp", &ScintTime);
  ScintChain->SetBranchAddress("SampleName", &ScintName);
  ScintChain->SetBranchAddress("IDnum", &ScintIDnum);
  ScintChain->SetBranchAddress("Dose", &ScintDose);
  ScintChain->SetBranchAddress("Formulation", &ScintFormulation);
  ScintChain->SetBranchAddress("CorrLY", &ScintCorrLY);

  Long64_t ScintNumEntries = ScintChain->GetEntries();

  vector<double> Dose;
  vector<double> FactorFinal;
  vector<double> Time;

  for(int i = 0; i<ScintNumEntries; i++){
    ScintChain->GetEntry(i);
    if((ScintIDnum==IDnum)&&(*ScintFormulation==Formulation)&&
       (ScintTime>IrradTime)){
      //It's one of ours! Add it to the vector.
      Dose.push_back(ScintDose);
      Time.push_back(ScintTime);
      //Calculate the corrected LY (multiply by the vial factor at this dose).
      if(DoseMapFactor.find(ScintDose)!=DoseMapFactor.end()){
	//Get the corrected factor.
	if((ScintCorrLY>0.001)&&(ScintCorrLY<100)&&
	   (DoseMapFactor[ScintDose]>0.001)&&(DoseMapFactor[ScintDose]<100)){
	  FactorFinal.push_back(ScintCorrLY*DoseMapFactor[ScintDose]);
	  //  cout << "DoseMapFactor[" << ScintDose << "] = "
	  //   << DoseMapFactor[ScintDose] 
	  //   << ", ScintCorrLY = " << ScintCorrLY 
	  //   << ", SampleName = " << *ScintName << endl;
				//(DoseMapFactor[ScintDose]/
				//DoseMapNumEntries[ScintDose]));
	}
	else{//error state, don't add this point
	  Dose.pop_back();
	  Time.pop_back();
	}
      }
      else{
	//I'll throw an error, as I shouldn't get to this stage with my
	//measurements.
	cout << "CorrectVials ERROR: Couldn't find the following vial dose: "
	     << ScintDose << " Gy. Did you irradiate vials to this dose level?"
	     << endl;
	return;
      }
      if((FactorFinal.back()<0.001)||(FactorFinal.back()==DBL_MAX)){
	//error states, eliminate these points.
	FactorFinal.pop_back();
	Dose.pop_back();	
	Time.pop_back();
      }
    }
  }

  //Maybe get the mean and stdev at each dose level and plot those instead?
  //for(size_t i = 0; i<Dose.size(); i++){
  //cout << "Dose.at(" << i << ") = " << Dose.at(i)
  //	 << ", FactorFinal.at(" << i << ") = " << FactorFinal.at(i) << endl;
  
  //}

  //  for(map<char,int>::iterator it=mymap.begin(); it!=mymap.end(); ++it)
  //std::cout << it->first << " => " << it->second << '\n';

  TGraph* theGraph = new TGraph(Dose.size(), &Dose[0], &FactorFinal[0]);
  TGraph* theTimeGraph = new TGraph(Time.size(), &Time[0], &FactorFinal[0]);
  theGraph->Draw(GraphOptions.c_str());
  //theTimeGraph->Draw(GraphOptions.c_str());

  fvial->Close();
  fscint->Close();

}
