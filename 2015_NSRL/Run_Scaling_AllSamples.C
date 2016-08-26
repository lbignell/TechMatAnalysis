//This program will loop over all file names in UniqueSampleNames.txt;
//running ScalingAlgorithm for each.
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

void Run_Scaling_AllSamples(string fname, string SampleNamesFile){

  ifstream ifst(SampleNamesFile.c_str());
  string line;
  stringstream ss;
  string prefix;
  int IDnum;
  string command;
  int testval = 0;
  TChain* OutChain = new TChain("ChainResults",
				"Chain of Scaling Results; for all Vials");
  while(std::getline(ifst,line)) {
    //cout << line << endl;
    //I need to check what type of vial is being used here; as they require
    //different thresholds...
    ss<<line;
    getline(ss, prefix, '_');
    if(prefix=="005A"){
      command = ".x ScalingAlgorithmv2.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 20, 5)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ScalingResults.root Scaling_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Scaling_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="005B"){
      command = ".x ScalingAlgorithmv2.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 20, 5)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ScalingResults.root Scaling_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Scaling_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="010A"){
      command = ".x ScalingAlgorithmv2.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 25, 5)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ScalingResults.root Scaling_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Scaling_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="010B"){
      command = ".x ScalingAlgorithmv2.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 30, 5)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ScalingResults.root Scaling_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Scaling_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="014A"){
      command = ".x ScalingAlgorithmv2.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 35, 5)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ScalingResults.root Scaling_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Scaling_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="014B"){
      command = ".x ScalingAlgorithmv2.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 40, 5)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ScalingResults.root Scaling_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Scaling_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="100"){
      command = ".x ScalingAlgorithmv2.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 300, 5)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ScalingResults.root Scaling_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Scaling_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if((prefix=="EMPTY")||(prefix=="EMPTY\r")){
      command = ".x ScalingAlgorithmv2.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 13, 5)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ScalingResults.root Scaling_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Scaling_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="STD"){
      command = ".x ScalingAlgorithmv2.C+(\"" + fname + 
	"\", \"SampleName\", \"" + line + "\", 400, 5)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ScalingResults.root Scaling_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "Scaling_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }  
    else{
      cout << "Run_Scaling_AllSamples::ERROR: Can't analyse sample "
	   << line << ", prefix = " << prefix << endl;
    }
    ss.str(std::string());
    ss.clear();
  }

  TFile* OutFile = new TFile("ChainedScalingResults.root", "RECREATE");
  OutChain->Write();
  OutFile->Close();
}
