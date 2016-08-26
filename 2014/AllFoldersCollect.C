
//This macro will loop over all folders that contain the string '-irrad' and
//run CollectAllLSData.C for each. The output file will be moved to a new,
//unique name. All outputs will then be added to a TChain.
#include <fstream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include <string>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include "tinydir.h"
#include <vector>
#include <algorithm>

using namespace std;

void AllFoldersCollect(){

  //Initialize
  tinydir_dir dir;
  tinydir_file file;
  vector<string> dirList;

  TChain* theChain = new TChain("LSChain", "Chain of Collected LS Spec");

  //Get the list of files in the current path.
  tinydir_open(&dir, ".");
  while(dir.has_next){
    tinydir_readfile(&dir, &file);
    if(file.is_dir){
      //check if the directory is one of the ones to be added to the TChain.
      if(strstr(file.name,"-irrad_")){dirList.push_back(file.name);
	//cout << "Adding directory: " << file.name << endl;}
      }
    }
    tinydir_next(&dir);
  }

  cout << "Got directory names, running over each..." << endl;
  //Sort the folder names.
  std::sort(dirList.begin(), dirList.end());

  //For each folder; collect up the data; save; and add to the TChain.
  for(int i = 0; i< dirList.size(); i++){
    cout << "Running " << dirList.at(i) << "...";
    string command = ".x CollectAllLSData.C+(\"" + dirList.at(i) + "/\")";
    gROOT->ProcessLine(command.c_str());
    command = ".! mv CollectedLSSpecData.root Collected_" + 
      dirList.at(i) + ".root";
    gROOT->ProcessLine(command.c_str());
    command = "Collected_" + dirList.at(i) + ".root/LSSpec";
    theChain->Add(command.c_str());
    cout << "Done!" << endl;
  }

  //TFile* f = new TFile("ChainedLSData.root", "RECREATE");
  //theChain->Write();
  //f->Close();

  theChain->Merge("MergedLSData.root", "RECREATE");

}
