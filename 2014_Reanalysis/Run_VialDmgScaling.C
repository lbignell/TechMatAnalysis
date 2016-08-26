//This program will loop over all file names in 'SampleNamesFile';
//running VialDmgScaling for each sample using the .root file 'fname'.
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

void Run_VialDmgScaling(string fname, string SampleNamesFile){

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
    //ss<<line;
    //getline(ss, prefix, '_');
    //if(prefix=="005A"){
    //command = ".x ScalingAlgorithmv2.C+(\"" + fname + 
    //	"\", \"SampleName\", \"" + line + "\", 20, 2)";
    //cout << command << endl;
    //gROOT->ProcessLine(command.c_str());
    //command = ".! mv ScalingResults.root Scaling_" + line + ".root";
    //gROOT->ProcessLine(command.c_str());
    //command = "Scaling_" + line + ".root/Results";
    //OutChain->Add(command.c_str());
    //}
    if((line.find("_REF")!=string::npos)||
       (line.find("_Leach")!=string::npos)||
       (line.find("STD")!=string::npos)||
       (line.find("EMPTY")!=string::npos)){
      //Found a reference, empty, leach, or std vial.
      //Do nothing for now...
    }
    //else if(prefix=="100"){
    else{//Assume it is one that we want.
      if(line.find("005")!=string::npos){
	command = ".x VialDmgScaling.C+(\"" + fname + 
	  "\", \"" + line + "\", 20)";
	cout << command << endl;
	gROOT->ProcessLine(command.c_str());
	command = ".! mv VialDmgScalingResults.root VialDmgScaling_" + line + 
	  ".root";
	gROOT->ProcessLine(command.c_str());
	command = "VialDmgScaling_" + line + ".root/Results";
	OutChain->Add(command.c_str());
      }
      else{
	command = ".x VialDmgScaling.C+(\"" + fname + 
	  "\", \"" + line + "\", 300)";
	cout << command << endl;
	gROOT->ProcessLine(command.c_str());
	command = ".! mv VialDmgScalingResults.root VialDmgScaling_" + line + 
	  ".root";
	gROOT->ProcessLine(command.c_str());
	command = "VialDmgScaling_" + line + ".root/Results";
	OutChain->Add(command.c_str());
      }
    }
    //else{
    //cout << "Run_Scaling_AllSamples::ERROR: Can't analyse sample "
    //	   << line << ", prefix = " << prefix << endl;
    //}
    ss.str(std::string());
    ss.clear();
  }

  TFile* OutFile = new TFile("ChainedVialDmgScalingResults.root", "RECREATE");
  OutChain->Write();
  OutFile->Close();
}
