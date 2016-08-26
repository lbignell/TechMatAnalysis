//This progam will run through every sample name as specified in
//the sample file and, for each of these run CorrectVial.C.
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

void Run_CorrectVial(string Scintdatafile, string Vialdatafile,
		    string samplefile){
  ifstream ifst(samplefile.c_str());
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
    command = ".x CorrectVial.C+(\"" + Scintdatafile + "\", \"" + Vialdatafile+ 
      "\", \"" + line + "\", 1428500000)";
    gROOT->ProcessLine(command.c_str());
    command = "mv VialCorrectedResults.root VialCorrected_" + line + ".root";
    gSystem->Exec(command.c_str());
    command = "ls | grep VialCorrected_" + line + ".root";
    testval = gSystem->Exec(command.c_str());
    if(testval==0){
      command = "VialCorrected_" + line + ".root/FinalResults";
      OutChain->Add(command.c_str());
    }
    linenum++;
  }
  TFile* OutFile = new TFile("FinalResults_Chained.root", "RECREATE");
  OutChain->Write();
  OutFile->Close();
}
