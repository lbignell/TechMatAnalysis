//This program will:
//-Get the .root file.
//-Make a vector of all the unique Sample Names in some file
//-Run ShiftingAlgorithm for each of these samples.
//-Combine all results of ShiftingAlgorithm into a Chain.
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
#include "TSelector.h"

using namespace std;

class AnalyseEverySample_Shifting : public TSelector{

public:

  TFile* fdata;
  TTree* DataTree;
  TBranch* NameBranch;
  TBranch* TimeBranch;
  TBranch* ChannelBranch;
  TBranch* CountsBranch;
  TBranch* PositionBranch;
  TBranch* SampleNumBranch;
  TBranch* CntTimeBranch;
  TBranch* HnumBranch;
  TBranch* TotCntsBranch;
  TBranch* IDnumBranch;
  TBranch* isRefBranch;

  //Create a vector of the sample names within the file.
  vector< string > vNames;
  //Long64_t nentries;
  string theName;
  Long64_t theTime;
  vector<double> theChannel;
  vector<double> theCounts;
  int thePosition;
  int theSampleNum;
  float theCntTime;
  float theHnum;
  int theTotCnts;
  int theIDnum;
  bool theisRef;

  TFile* OutFile;
  TTree* OutTree;
  string OutName;
  Long64_t OutTime;
  int OutAbsMeasNum;
  vector<double> OutChannel;
  vector<double> OutCounts;
  int OutPosition;
  int OutSampleNum;
  float OutCntTime;
  float OutHnum;
  int OutTotCnts;
  int OutIDnum;
  bool OutisRef;

  ofstream ofs;

  AnalyseEverySample_Shifting(TTree* = 0):
    theTime(0), OutTime(0), OutAbsMeasNum(0){
    //Do some initialisation here
    theChannel.assign(4096, 0.);
    theCounts.assign(4096, 0.);
    theName = "";
    thePosition = 0;
    theSampleNum = 0;
    theCntTime = 0;
    theHnum = 0;
    theTotCnts = 0;
    theIDnum = 0;
    theisRef = false;
  }
  virtual ~AnalyseEverySample_Shifting(){}
  virtual void    Init(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual Bool_t  Process(Long64_t entry);
  virtual void    Terminate();
  virtual Int_t   Version() const { return 2; }
 
  ClassDef(AnalyseEverySample_Shifting, 0);
};

void AnalyseEverySample_Shifting::Init(TTree* tree){
  cout << "In Init...";
  cout << "tree = " << tree << endl;
  //This function is for initializing a tree or chain.
  //nentries = tree->GetEntries;
  tree->SetMakeClass(1);
  //Assign the branches
  cout << "branchstatus:\t" << tree->GetBranchStatus("SampleName") << endl;
  tree->SetBranchAddress("SampleName", &theName, &NameBranch);
  tree->SetBranchAddress("TimeStamp", &theTime, &TimeBranch);
  tree->SetBranchAddress("Channel", &theChannel, &ChannelBranch);
  tree->SetBranchAddress("Counts",  &theCounts, &CountsBranch);
  tree->SetBranchAddress("Position", &thePosition, &PositionBranch);
  tree->SetBranchAddress("SampleNum", &theSampleNum, &SampleNumBranch);
  tree->SetBranchAddress("CntTime", &theCntTime, &CntTimeBranch);
  tree->SetBranchAddress("Hnum", &theHnum, &HnumBranch);
  tree->SetBranchAddress("TotCnts", &theTotCnts, &TotCntsBranch);
  tree->SetBranchAddress("IDnum", &theIDnum, &IDnumBranch);
  tree->SetBranchAddress("isRef", &theisRef, &isRefBranch);
  DataTree = tree;
  cout << "Done!" << endl;
}

void AnalyseEverySample_Shifting::SlaveBegin(TTree* tree){
  cout << "In SlaveBegin...";
  //SlaveBegin() is where you initialize trees/histograms.
  OutFile = new TFile("EverySample_Shifted.root", "RECREATE");
  OutTree = new TTree("Results", "The shifted results; all samples");

  //Create all of the branches.
  OutTree->Branch("SampleName", "string", &OutName);
  OutTree->Branch("TimeStamp", &OutTime, "TimeStamp/L");
  OutTree->Branch("AbsMeasNum", &OutAbsMeasNum, "AbsMeasNum/I");
  OutTree->Branch("Channel", "vector<double>", &OutChannel);
  OutTree->Branch("Counts", "vector<double>", &OutCounts);
  OutTree->Branch("Position", &OutPosition, "Position/I");
  OutTree->Branch("SampleNum", &OutSampleNum, "SampleNum/I");
  OutTree->Branch("CntTime", &OutCntTime, "CntTime/F");
  OutTree->Branch("Hnum", &OutHnum, "Hnum/F");
  OutTree->Branch("TotCnts", &OutTotCnts, "TotCnts/I");
  OutTree->Branch("IDnum", &OutIDnum, "IDnum/F");
  OutTree->Branch("isRef", &OutisRef, "isRef/O");

  ofs.open("SampleNames.txt", ios::out | ios::trunc);
  cout << "Done!" << endl;
}

Bool_t AnalyseEverySample_Shifting::Process(Long64_t entry){
  cout << "Processing, entry = " << entry << ", ";

  //This method is called for each entry in the tree (chain).
  OutAbsMeasNum = entry + DataTree->GetChainOffset();
  cout << "vNames.size() = " << vNames.size() << endl;

  cout << "AbsMeasNum = " << OutAbsMeasNum << endl;


  //Get the variables.
  cout << "Setting variables: OutChannel...";
  cout << "theChannel.size() = " << theChannel.size() << "...";
  cout << "theChannel = " << &theChannel << endl;
  OutChannel = theChannel;
  cout << "theCounts...";
  OutCounts = theCounts;
  cout << "thePosition...";
  OutPosition = thePosition;
  OutSampleNum = theSampleNum;
  OutCntTime = theCntTime;
  OutHnum = theHnum;
  OutTotCnts = theTotCnts;
  OutIDnum = theIDnum;
  OutisRef = theisRef;
  OutTime = theTime;
  OutName = theName;
  cout << "OutName = " << OutName << endl;

  bool isRpt = false;


  if(vNames.size()!=0){
    for(int i = 0; i<vNames.size(); i++){
      if(OutName==vNames.at(i)){isRpt = true;}
    }
    if(!isRpt){
      ofs << OutName;
      ofs << "\n";
      vNames.push_back(OutName);
    }
  }
  else{
    ofs << OutName;
    ofs << "\n";
    vNames.push_back(OutName);
  }

  //Fill the tree, if the timestamp is OK.
  if(OutTime>0){
    OutTree->Fill();
  }

  return kTRUE;
}

void AnalyseEverySample_Shifting::Terminate(){
  //Here's where to plot or save.
  ofs.close();
  OutTree->Write();
}
