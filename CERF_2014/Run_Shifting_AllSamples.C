//This program will loop over all file names in UniqueSampleNames.txt;
//running ShiftingAlgorithm for each.
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

void Run_Shifting_AllSamples(string fname){
  TFile* InFile = TFile::Open(fname.c_str());
  //TChain* theChain = (TChain*)InFile->Get("LSSpec");
  TTree* theTree = (TTree*)InFile->Get("Results");
  ifstream ifst("SampleNames.txt");
  string line;
  stringstream ss;
  string prefix;
  int IDnum;
  string command;
  int testval = 0;
  TChain* OutChain = new TChain("ChainResults",
				"Chain of Shifting Results; for all Vials");
  while(std::getline(ifst,line)) {
    //cout << line << endl;
    //I need to check what type of vial is being used here; as they require
    //different thresholds...
    ss<<line;
    getline(ss, prefix, '_');
    if(prefix=="005A"){
      command = ".x ShiftingAlgorithm.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 20, 2)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Shifting_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="005B"){
      command = ".x ShiftingAlgorithm.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 20, 2)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Shifting_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="010A"){
      command = ".x ShiftingAlgorithm.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 25, 2)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Shifting_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="010B"){
      command = ".x ShiftingAlgorithm.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 30, 2)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Shifting_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="014A"){
      command = ".x ShiftingAlgorithm.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 35, 2)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Shifting_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="014B"){
      command = ".x ShiftingAlgorithm.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 40, 2)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Shifting_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="100"){
      command = ".x ShiftingAlgorithm.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 300, 2)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Shifting_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if((prefix=="EMPTY")||(prefix=="EMPTY\r")){
      command = ".x ShiftingAlgorithm.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 13, 2)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Shifting_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="STD"){
      command = ".x ShiftingAlgorithm.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 400, 2)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Shifting_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }  
    else{
      cout << "Run_Shifting_AllSamples::ERROR: Can't analyse sample "
	   << line << ", prefix = " << prefix << endl;
    }
    ss.str(std::string());
    ss.clear();
  }

  TFile* OutFile = new TFile("ChainedShiftingResults.root", "RECREATE");
  OutChain->Write();
  OutFile->Close();
}
