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

void Run_Shifting_AllSamples_Chain(string fname){
  TFile* InFile = TFile::Open(fname.c_str());
  TChain* theChain = (TChain*)InFile->Get("LSChain");
  //TTree* theTree = (TTree*)InFile->Get("Results");
  ifstream ifst("SampleNames_trunc.txt");
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
    if(istringstream(prefix)>>IDnum){
      //it is a numeric prefix; process accordingly.
      if((IDnum==50)||(IDnum==52)){
	command = ".x ShiftingAlgorithm.C+(\"" + fname + "\", \"SampleName\", \"" + line + "\", 15, 5)";
	cout << command << endl;
	gROOT->ProcessLine(command.c_str());
	command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
	gROOT->ProcessLine(command.c_str());
	command = "ls | grep Shifting_" + line + ".root";
	testval = gSystem->Exec(command.c_str());
	if(testval==0){
	  command = "Shifting_" + line + ".root/Results";
	  OutChain->Add(command.c_str());
	}
	else{
	  command = "Shifting_" + line + "_.root/Results";
	  OutChain->Add(command.c_str());
	}
      }
      else if((IDnum==100)||(IDnum==102)){
	command = ".x ShiftingAlgorithm.C+(\"" + fname + "\", \"SampleName\", \"" + line + "\", 20, 5)";
	cout << command << endl;
	gROOT->ProcessLine(command.c_str());
	command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
	gROOT->ProcessLine(command.c_str());
	command = "ls | grep Shifting_" + line + ".root";
	testval = gSystem->Exec(command.c_str());
	if(testval==0){
	  command = "Shifting_" + line + ".root/Results";
	  OutChain->Add(command.c_str());
	}
	else{
	  command = "Shifting_" + line + "_.root/Results";
	  OutChain->Add(command.c_str());
	}
      }
      else if((IDnum==140)||(IDnum==142)){
	command = ".x ShiftingAlgorithm.C+(\"" + fname + "\", \"SampleName\", \"" + line + "\", 25, 5)";
	cout << command << endl;
	gROOT->ProcessLine(command.c_str());
	command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
	gROOT->ProcessLine(command.c_str());
	command = "ls | grep Shifting_" + line + ".root";
	testval = gSystem->Exec(command.c_str());
	if(testval==0){
	  command = "Shifting_" + line + ".root/Results";
	  OutChain->Add(command.c_str());
	}
	else{
	  command = "Shifting_" + line + "_.root/Results";
	  OutChain->Add(command.c_str());
	}
      }
      else if((IDnum==1000)||(IDnum==1002)){
	command = ".x ShiftingAlgorithm.C+(\"" + fname + "\", \"SampleName\", \"" + line + "\", 200, 5)";
	cout << command << endl;
	gROOT->ProcessLine(command.c_str());
	command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
	gROOT->ProcessLine(command.c_str());
	command = "ls | grep Shifting_" + line + ".root";
	testval = gSystem->Exec(command.c_str());
	if(testval==0){
	  command = "Shifting_" + line + ".root/Results";
	  OutChain->Add(command.c_str());
	}
	else{
	  command = "Shifting_" + line + "_.root/Results";
	  OutChain->Add(command.c_str());
	}
      }
      else{
	cout<< "ERROR: undefined sample! The name is: " << line << endl;
      }
    }
    //non-numeric prefix; process accordingly.
    else if((prefix=="EMPTY")||(prefix=="EMPTY\r")){
      command = ".x ShiftingAlgorithm.C+(\"" + fname + "\", \"SampleName\", \"" + line + "\", 12, 5)";
	cout << command << endl;
      gROOT->ProcessLine(command.c_str());
	command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
	gROOT->ProcessLine(command.c_str());
	command = "ls | grep Shifting_" + line + ".root";
	testval = gSystem->Exec(command.c_str());
	if(testval==0){
	  command = "Shifting_" + line + ".root/Results";
	  OutChain->Add(command.c_str());
	}
	else{
	  command = "Shifting_" + line + "_.root/Results";
	  OutChain->Add(command.c_str());
	}
    }
    else if(prefix=="STD"){
      command = ".x ShiftingAlgorithm.C+(\"" + fname + "\", \"SampleName\", \"" + line + "\", 200, 5)";
	cout << command << endl;
      gROOT->ProcessLine(command.c_str());
      command = ".! mv ShiftingResults.root Shifting_" + line + ".root";
      gROOT->ProcessLine(command.c_str());
      command = "ls | grep Shifting_" + line + ".root";
      testval = gSystem->Exec(command.c_str());
      if(testval==0){
	command = "Shifting_" + line + ".root/Results";
	OutChain->Add(command.c_str());
      }
      else{
	command = "Shifting_" + line + "_.root/Results";
	OutChain->Add(command.c_str());
      }
    }  
  }

  TFile* OutFile = new TFile("ChainedShiftingResults.root", "RECREATE");
  OutChain->Write();
  OutFile->Close();
  
}

