//This progam will run through every fsample name as specified in
//SampleNames_trunc.txt and, for each of these run CorrectRef.C.
//It will Chain the results files together.
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

using namespace std;

void Run_CorrectRef(){
  ifstream ifst("SampleNames_trunc.txt");
  string line;
  stringstream ss;
  string prefix;
  int IDnum;
  string command;
  int testval = 0;
  int linenum = 0;
  TChain* OutChain = new TChain("ChainFinalResults",
				"Chain of Corrected Shifting Results");
  while(std::getline(ifst,line)) {
    cout << "linenum = " << linenum << endl;
    command = ".x CorrectRef.C+(\"ChainedShiftingResults.root\", \"" 
      + line + "\")";
    gROOT->ProcessLine(command.c_str());
    command = "mv CorrectedShiftedResults.root Corrected_" + line + ".root";
    gSystem->Exec(command.c_str());
    command = "ls | grep Corrected_" + line + ".root";
    testval = gSystem->Exec(command.c_str());
    if(testval==0){
      command = "Corrected_" + line + ".root/FinalResults";
      OutChain->Add(command.c_str());
    }
    linenum++;
  }
  TFile* OutFile = new TFile("ChainedFinalResults.root", "RECREATE");
  OutChain->Write();
  OutFile->Close();
}
