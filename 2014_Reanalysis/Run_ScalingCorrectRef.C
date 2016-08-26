//This program will loop over all file names in SampleNamesFile;
//running ScalingCorrectRef for each.
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

void Run_ScalingCorrectRef(string fname, string SampleNamesFile){

  ifstream ifst(SampleNamesFile.c_str());
  string line;
  stringstream ss;
  string prefix;
  string dose;
  string SampleNum;
  int IDnum;
  string command;
  int testval = 0;
  TChain* OutChain = new TChain("ChainResults",
				"Chain of Scaling Results; for all Vials");
  while(std::getline(ifst,line)) {
    cout << line << endl;
    //I need to check what type of vial is being used here; as they require
    //different thresholds...
    ss<<line;
    getline(ss, prefix, '_');
    if(prefix=="0050"){
      command = ".x ScalingCorrectRef.C+(\"" + fname + 
	"\", \"" + line + "\", \"0050_REF1\", 15)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ScalingCorrectRefResults.root ScalingCorrRef1_" + 
	line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "ScalingCorrRef1_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else if(prefix=="0100"){
      if(line!="0100_REF2"){
	command = ".x ScalingCorrectRef.C+(\"" + fname + 
	  "\", \"" + line + "\", \"0100_REF1\", 20)";
	cout << command << endl;
	gROOT->ProcessLine(command.c_str());
	command = ".! mv ScalingCorrectRefResults.root ScalingCorrRef1_" + 
	  line + ".root";
	gROOT->ProcessLine(command.c_str());
	command = "ScalingCorrRef1_" + line + ".root/Results";
	OutChain->Add(command.c_str());
      }
    }
    else if(prefix=="0140"){
      if(line!="0140_REF1"){
	command = ".x ScalingCorrectRef.C+(\"" + fname + 
	  "\", \"" + line + "\", \"0140_REF1\", 30)";
	cout << command << endl;
	gROOT->ProcessLine(command.c_str());
	command = ".! mv ScalingCorrectRefResults.root ScalingCorrRef1_"
	  + line + ".root";
	gROOT->ProcessLine(command.c_str());
	command = "ScalingCorrRef1_" + line + ".root/Results";
	OutChain->Add(command.c_str());
      }
    }
    else if(prefix=="1000"){
      command = ".x ScalingCorrectRef.C+(\"" + fname + 
	"\", \"" + line + "\", \"1000_REF1\", 300)";
      cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ScalingCorrectRefResults.root ScalingCorrRef1_" + 
	line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "ScalingCorrRef1_" + line + ".root/Results";
      OutChain->Add(command.c_str());
    }
    else{
      cout << "Run_ScalingCorrectRef_AllSamples::ERROR: Can't analyse sample "
	   << line << ", prefix = " << prefix << endl;
    }
    ss.str(std::string());
    ss.clear();
  }

  TFile* OutFile = new TFile("ChainedScalingCorrectRef.root", "RECREATE");
  OutChain->Write();
  OutFile->Close();
}
