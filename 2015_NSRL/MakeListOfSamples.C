//This function just generates a file with a list of unique sample names, in
//a file called samplesfname.
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

using namespace std;

void MakeListOfSamples(string fname, string samplesfname){

  //Open the file.
  TFile* fdata = TFile::Open(fname.c_str());
  TTree* theTree = (TTree*)fdata->Get("LSSpec");

  Long64_t nentries = theTree->GetEntries();
  string* theName = NULL;
  theTree->SetBranchAddress("SampleName", &theName);

  vector< string > vNames;
  Long64_t numEntries = theTree->GetEntries();
  string thisName = "";
  bool isRpt = false;

  ofstream ofs(samplesfname.c_str(), ios::out|ios::trunc);

  for(Long64_t i = 0; i<nentries; i++){

    theTree->GetEntry(i);
    isRpt = false;

    if(vNames.size()!=0){
      for(int i = 0; i<vNames.size(); i++){
	if((*theName)==vNames.at(i)){isRpt = true;}
      }
      if(!isRpt){
	//Account for the windows '\r\n' newline
	if(theName->find("\r")!=string::npos)
	  {theName->replace(theName->find("\r"), string::npos, "");}
	ofs << *theName;
	ofs << "\n";
	vNames.push_back(*theName);
      }
    }
    else{
      //Account for the windows '\r\n' newline
      if(theName->find("\r")!=string::npos)
	{theName->replace(theName->find("\r"), string::npos, "");}
      ofs << *theName;
      ofs << "\n";
      vNames.push_back(*theName);
    }
  }

  ofs.close();
  fdata->Close();

}
