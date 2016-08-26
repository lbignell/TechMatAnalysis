//This function will implement an algorithm to determine how much light yield
//gain/loss has been imparted upon a liquid scintillator upon being irradiated,
//by using unirradiated scintillator data as a reference.
//
//////////////////////////////////////////////////////////////////////////////
//Data are stored in ROOT files in subdirectories, so I want to have the file
//name as one of the arguments.
//
//The files store the sample names in a couple of branches; SampleName (a
//string) and IDnum (an integer). SampleName specifies individual vials, whereas
//IDnum specifies a type of sample composition (eg 5% WbLS).
//I want to be able to choose whether to access data by reference to the
//SampleName or IDnum, so I'll take a string that corresponds to the appropriate
//branch name as another argument.
//IDnum reference:
//-1 = unknown/indeteriminate
//0 = empty vial
//1 = irradiated (or ref) glass vial filled with fresh LS.
//2 = irradiated (or ref) PP vial filled with fresh LS.
//3 = irradiated (or ref) HDPE vial filled with fresh LS.
//10 = H-3 or C-14 (or other) standard vial.
//50 = 5% WbLS
//52 = 5% WbLS, 1% AA
//100 = 10% WbLS
//102 = 10% WbLS, 1% AA
//140 = 14% WbLS
//142 = 14% WbLS, 1% AA
//1000 = Pure LS
//1002 = Pure LS, 1% AA
//
//A third argument will need to specify which sample or type that I want.
//I'll make this a string, and parse it to an integer if IDnum was chosen.
//
//The final argument is a threshold index below which bins in the spectrum
//aren't examined to see if they depart from the reference data set in a
//statistically significant manner. The best value for the index depends upon
//the scintillator composition and should be set emprically.
//As a rule of thumb; the best values for % LS studied so far are:
//5% = 15
//10% = 20
//14% = 25
//Pure LS = 200
//
//This program will compare a single dataset (not repeated measurements across
//ROOT files, though in principle this could be possible if I stitch the 
//appropriate ROOT files together). The first measurement will be taken as the
//reference dataset against which all later datasets are compared. I will save
//the analysed data to a new ROOT file called ShiftingResults.

#include <string>
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"
#include "TBranch.h"
#include "TString.h"
#include "TH1D.h"
#include <fstream>
#include <iostream>
#include "TLeaf.h"
#include <math.h>
#include "TMath.h"
#include <sstream>

void* GetPtrToVal(TBranch* theBranch, int entry){
  const char* name = theBranch->GetName();
  theBranch->GetEntry(entry);
  TLeaf* theLeaf = theBranch->GetLeaf(name);
  return theLeaf->GetValuePointer();
}

void ShiftingAlgorithm(string DataFile, string BranchName, string SampleID,
		       int ThresholdIdx, int NumAvg){



  if(!((BranchName=="SampleName")||(BranchName=="IDnum"))){
    cout << "ERROR: this code only works with the second argument equal to"
	 << " 'SampleName' or 'IDnum'. You wrote: " << BranchName << "." << endl
	 << "Nothing has been done, exiting..." << endl;
    return;
  }

  cout << "Opening file " << DataFile << endl;
  cout << "SampleID = " << SampleID << endl;
  cout << "ThresholdIdx = " << ThresholdIdx << endl;
  cout << "NumAvg = " << NumAvg << endl;

  //Open the test/reference data.
  TFile* fdata = TFile::Open(DataFile.c_str());
  //TTree* theTree = (TTree*)fdata->Get("Results");
  cout << "Getting Chain... ";
  TChain* theTree = (TChain*)fdata->Get("LSChain");
  //cout << "Done! Getting Branches...";
  //TBranch* TypeBranch = (TBranch*)theTree->GetBranch(BranchName.c_str());
  //TBranch* CntBranch = (TBranch*)theTree->GetBranch("Counts");
  //TBranch* ChnBranch = (TBranch*)theTree->GetBranch("Channel");
  //TBranch* TimeBranch = (TBranch*)theTree->GetBranch("TimeStamp");
  //TBranch* MeasBranch = (TBranch*)theTree->GetBranch("AbsMeasNum");
  //TBranch* CntTimeBranch = (TBranch*)theTree->GetBranch("CntTime");
  //TBranch* HnumBranch = (TBranch*)theTree->GetBranch("Hnum");
  //TBranch* TotCntsBranch = (TBranch*)theTree->GetBranch("TotCnts");
  //TBranch* NameBranch = (TBranch*)theTree->GetBranch("SampleName");
  //TBranch* IDnumBranch = (TBranch*)theTree->GetBranch("IDnum");
  //TBranch* DoseBranch = (TBranch*)theTree->GetBranch("Dose");
  //TBranch* isRefBranch = (TBranch*)theTree->GetBranch("isRef");
  //cout << "Done!" << endl;

  string* theName = NULL;
  Long64_t theTime = 0;
  vector<double>* theChannel = NULL;
  vector<double>* theCounts = NULL;
  int thePosition = 0;
  int theSampleNum = 0;
  float theCntTime = 0;
  float theHnum = 0;
  int theTotCnts = 0;
  int theIDnum = -100;
  bool theisRef = false;
  int theAbsMeasNum = 0;

  theTree->SetBranchAddress("Counts", &theCounts);
  theTree->SetBranchAddress("Channel", &theChannel);
  theTree->SetBranchAddress("TimeStamp", &theTime);
  //theTree->SetBranchAddress("AbsMeasNum", &theAbsMeasNum);
  theTree->SetBranchAddress("CntTime", &theCntTime);
  theTree->SetBranchAddress("Hnum", &theHnum);
  theTree->SetBranchAddress("TotCnts", &theTotCnts);
  theTree->SetBranchAddress("SampleName", &theName);
  theTree->SetBranchAddress("IDnum", &theIDnum);
  theTree->SetBranchAddress("isRef", &theisRef);
  theTree->SetBranchAddress("Position", &thePosition);
  theTree->SetBranchAddress("SampleNum", &theSampleNum);
  cout << "Done!" << endl;

  cout << "Getting nentries...";
  Long64_t numEntries = theTree->GetEntries();
  cout << "Done! numEntries " << numEntries << endl;

  int IDnumber = 0;
  int thisIDnumber = 0;
  string thisSampleID = "";
  int loopCtr = 0;

  //Some vectors to store the ref data.
  vector<double> CntRef(4096,0);
  vector<double> Channel(4096,0);
  bool isIDnum = false;

  vector<double> thisCntRef(4096,0);
  int NumSamples = 0;

  //Now find the FIRST measurement's Cnt and Chn values.
  while(loopCtr<numEntries){
    //cout << "Getting entry " << loopCtr << endl;
    theTree->GetEntry(loopCtr);
    //IDnum case
    //DON'T do an average of the spectrum in the IDnum case!
    if(BranchName=="IDnum"){
      //We need to parse to SampleID to an integer.
      if(istringstream(SampleID)>>IDnumber){}//All good
      else{
	cout << "ERROR: your sample ID was rubbish (couldn't parse to int)."
	     << " SampleID = " << SampleID << endl
	     << "Doing nothing and exiting..." << endl;
	fdata->Close();
	return;
      }

      //cout << "Getting ID number...";
      //Now look at the IDnum for this event.
      //thisIDnumber = *(int*)GetPtrToVal(IDnumBranch, loopCtr);
      thisIDnumber = theIDnum;
      //cout << "Done! thisIDnumber = " << thisIDnumber << endl;

      if(thisIDnumber==IDnumber){
	//We've found the first one, fill the refCnt and refChn vectors.
	//CntRef = *(vector<double>*)GetPtrToVal(CntBranch, loopCtr);
	//Channel = *(vector<double>*)GetPtrToVal(ChnBranch, loopCtr);
	CntRef.swap(*theCounts);
	Channel.swap(*theChannel);
	isIDnum = true;
	break;
      }

    }
    //Now examine the SampleName case.
    else{

      //thisSampleID = *(string*)GetPtrToVal(NameBranch, loopCtr);
      thisSampleID = *theName;

      if(thisSampleID==SampleID){

	cout << "Found SampleID = " << thisSampleID << endl;

	//We've found the first one, fill the refCnt and refChn vectors.
	//thisCntRef = *(vector<double>*)GetPtrToVal(CntBranch, loopCtr);
	thisCntRef.swap(*theCounts);
	if(CntRef.size()!=0){
	  //Do an average
	  for(int i = 0; i<CntRef.size(); i++){
	    CntRef.at(i) = CntRef.at(i) + (thisCntRef.at(i)/NumAvg);
	  }
	  NumSamples++;
	  //cout << "NumSamples = " << NumSamples << ", 10th datum = "
	  //   << CntRef.at(9) << endl;
	}
	//else{
	//for(int i = 0; i<thisCntRef.size(); i++){
	//  CntRef.push_back(thisCntRef.at(i)/NumAvg);
	    //if(i==20){
	    //cout << "NumSamples = " << NumSamples << ", 10th datum = "
	    //	   << CntRef.at(9) << endl;
	    //}
	// }
	//}
	//Channel = *(vector<double>*)GetPtrToVal(ChnBranch, loopCtr);
	Channel.swap(*theChannel);
	isIDnum = false;
	if(NumSamples==NumAvg){
	  break;
	}
      }
    }
    loopCtr++;
  }

  //Define a vector for storing the test data
  vector<double> CntTest;
  Long64_t TestTime;
  float TestCntTime;
  float TestHnum;
  int TestTotCnts;
  int TestMeasNum;
  string TestSampleName = "";//default to skip state
  int TestIDnum = -2;//default to skip state
  //float TestDose;

  //Other variables needed for algorithm calculation.
  vector<double> NormDiff;
  vector<double> mean;
  vector<double> stdev;
  vector<double> tVals;
  vector<double> pVals;
  int TestIdx = 0;
  int TestChannel = 0;
  double RefStart = -1;
  double RefEnd = -1;
  double TestStart = -1;
  double TestEnd = -1;
  double MinFactor = DBL_MAX;
  double MinSumSq = DBL_MAX;  
  vector<double> MinScaledCntTest;
  vector<double> vecFactor;
  vector<double> vecSumSq;
  bool isRef = false;

  //Define a ROOT file (with tree) for saving the results.
  //Save results in a tree and save to disk.
  TFile* fout = new TFile("ShiftingResults.root", "RECREATE");
  TTree* outTree = new TTree("Results","Light yield shifting algo results");
  outTree->Branch("Channel", "vector<double>", &Channel);
  outTree->Branch("TimeStamp", &TestTime, "TimeStamp/L");
  outTree->Branch("CntTime", &TestCntTime, "CntTime/F");
  outTree->Branch("Hnum", &TestHnum, "Hnum/F");
  outTree->Branch("TotCnts", &TestTotCnts, "TotCnts/I");
  outTree->Branch("AbsMeasNum", &TestMeasNum, "AbsMeasNum/I");
  outTree->Branch("SampleName", "string", &TestSampleName);
  outTree->Branch("IDnum", &TestIDnum, "IDnum/I");
  //outTree->Branch("Dose", &TestDose, "Dose/F");
  outTree->Branch("TestCounts", "vector<double>", &CntTest);
  outTree->Branch("RefCounts", "vector<double>", &CntRef);
  outTree->Branch("NormDiff", "vector<double>", &NormDiff);
  outTree->Branch("Mean", "vector<double>", &mean);
  outTree->Branch("Stdev", "vector<double>", &stdev);
  outTree->Branch("tValues", "vector<double>", &tVals);
  outTree->Branch("pVals", "vector<double>", &pVals);
  outTree->Branch("TestIdx", &TestIdx, "TestIdx/I");
  outTree->Branch("TestChannel", &TestChannel, "TestChannel/I");
  outTree->Branch("TestStart", &TestStart, "TestStart/D");
  outTree->Branch("TestEnd", &TestEnd, "TestEnd/D");
  outTree->Branch("RefStart", &RefStart, "RefStart/D");
  outTree->Branch("RefEnd", &RefEnd, "RefEnd/D");
  outTree->Branch("Factor", &MinFactor, "Factor/D");
  outTree->Branch("SumSq", &MinSumSq, "SumSq/D");
  outTree->Branch("ScaledTestCounts", "vector<double>", &MinScaledCntTest);
  outTree->Branch("vecFactor", "vector<double>", &vecFactor);
  outTree->Branch("vecSumSq", "vector<double>", &vecSumSq);
  outTree->Branch("isRef", &isRef, "isRef/O");

  //loop over the entries, comparing with the reference data when a match is
  //found
  for(int k = 0; k<numEntries; k++){
    theTree->GetEntry(k);
    if(isIDnum){
      //Get next ID num and compare with the user input.
      //TestIDnum = *(int*)GetPtrToVal(IDnumBranch, k);
      TestIDnum = theIDnum;
      //TestSampleName = *(string*)GetPtrToVal(NameBranch, k);
      TestSampleName = *theName;
      if(TestIDnum == IDnumber){
	//Get the entry
	//CntTest = *(vector<double>*)GetPtrToVal(CntBranch, k);
	CntTest.swap(*theCounts);
      }
      else{TestIDnum = -2;}//Flag for unwanted entry.
    }
    else{
      //Case for SampleName
      //TestSampleName = *(string*)GetPtrToVal(NameBranch, k);
      TestSampleName = *theName;
      //      TestIDnum = *(int*)GetPtrToVal(IDnumBranch, k);
      TestIDnum = theIDnum;
      if(TestSampleName==SampleID){
	//Get the entry
	//CntTest = *(vector<double>*)GetPtrToVal(CntBranch, k);
	CntTest.swap(*theCounts);
      }
    }

    //Test whether to proceed or loop to the next entry
    if(((isIDnum)&&(TestIDnum==IDnumber))||
       ((!isIDnum)&&(TestSampleName==SampleID))){
      cout << "Analysing TestSampleName = " << TestSampleName << endl;

      if((CntTest.size()!=CntRef.size())){
	cout << "ERROR: Ref and Test data are not equal lengths" << endl
	     << "Ref size = " << CntRef.size() << ", Test size = " 
	     << CntTest.size() << endl
	     << "TestIDnum = " << TestIDnum << endl
	     << "TestSampleName = " << TestSampleName << endl
	     << endl;
	return;
      }

      //cout << "Processing, TestSampleName = " << TestSampleName << ". "
      //   << "TestIDnum = " << TestIDnum << endl;

      //Calculate the difference/sigma
      vector<double>::iterator testit = CntTest.begin();
      vector<double>::iterator refit = CntRef.begin();
      vector<double>::iterator chanit = Channel.begin();
      NormDiff.clear();
      vector<double> CumNormDiff;
      double tally = 0;
      RefStart = -1;
      RefEnd = -1;
      TestStart = -1;
      TestEnd = -1;
      int RefEndIdx = 0;
      int theIndex = 0;

      while((refit!=CntRef.end())&&(testit!=CntTest.end())){
	//get the 'normalised difference' values.
	if(*refit!=0){
	  NormDiff.push_back((*testit - *refit)/sqrt(*refit));
	  //Alternatively I can divide by sigma squared...
	  //NormDiff.push_back((*testit - *refit)/(*refit));
	  tally += (*testit - *refit)/sqrt(*refit);
	  CumNormDiff.push_back(tally);
	}    
	else{
	  NormDiff.push_back(0);
	  CumNormDiff.push_back(tally);  
	}
	
	//Also get the start and end point of each spectrum for renormalisation.
	if(*testit>10){
	  TestEnd = *chanit;
	  if(TestStart==-1){TestStart = *chanit;}
	}
	
	if(*refit>10){
	  RefEnd = *chanit;
	  RefEndIdx = theIndex;
	  if(RefStart==-1){RefStart = *chanit;}
	}
	
	testit++;
	refit++;
	chanit++;
	theIndex++;
      }
      
      double RefWidth = RefEnd - RefStart;
      double TestWidth = TestEnd - TestStart;  
      
      //I now need to step along NormDiff; getting the t-statistic.
      mean.clear();
      stdev.clear();
      double theMean;
      double theStdev;
      //Number of samples over which moving avg is taken.
      int avglen = 3;
      //Threshold, below which the test statistic isn't examined.
      int threshold = avglen + ThresholdIdx;
      //p-test criterion, make it 1%.
      double Criterion = 0.05;  
      bool IsStatSig = false;
      
      tVals.clear();
      pVals.clear();
      double thepVal;
      TestIdx = 0;
      TestChannel = 0;
      int numBelowCriterion = 0;
      int lengthBelowCriterion = 30;//semi-arbitrary, to prevent triggering on
      //pedestal
      
      for(int i = 0; i<NormDiff.size(); i++){
	if((i<threshold)||(i>(NormDiff.size()-avglen))){/*do nothing*/
	  mean.push_back(0);
	  stdev.push_back(1);
	  tVals.push_back(0);
	  pVals.push_back(-1);
	}
	else{
	  theMean = (CumNormDiff.at(i+floor(avglen/2)) - 
		     CumNormDiff.at(i-avglen+floor(avglen/2)))/avglen;

	  //Stdev should equal ~ 1.
	  //This is because when I took the difference between the ref and test
	  //I divided by the standard devation.
	  theStdev = 1;
	  tVals.push_back(sqrt(avglen)*(theMean/theStdev));
	  mean.push_back(theMean);
	  stdev.push_back(theStdev);
	  thepVal = TMath::StudentI(tVals.back(), (avglen-1));

	  if((thepVal!=-1)&&(thepVal<Criterion)&&(TestIdx==0)){
	    numBelowCriterion++;
	    if(numBelowCriterion>lengthBelowCriterion){
	      TestIdx = i-lengthBelowCriterion;
	      TestChannel = Channel.at(i-lengthBelowCriterion);
	      //cout << "TestIdx = " << TestIdx << endl
	      //   << "TestChannel = " << TestChannel << endl;
	      IsStatSig = true;
	    }
	  }
	  else{
	    numBelowCriterion = 0;
	  }
	  pVals.push_back(thepVal);
	}
      }

      if(IsStatSig==false){
    //	cout << "Data show no statistically significant difference!" << endl;
      }

      if(TestIdx<threshold){
	//cout << "TestIdx = " << TestIdx << "! Shifting to " << threshold
	//   << endl;
	TestIdx = threshold;
	TestChannel = Channel.at(TestIdx);
      }

  /////////////////////////////////////////////////////////////////////////////
  //Insert code to scale test data and perform minimisation here.
  /////////////////////////////////////////////////////////////////////////////
  //First guess of scale factor.
      double factor = 1.0;
      double SumSq = 0;
      double SumSqPrev = 0;
      double factorPrev = 1.03;
      //Loop infinitely until minimisation converges.
      vector<double> ScaledCntTest;
      MinScaledCntTest.clear();
      int NumBins = 0;
      
      int iterNum = 0;
      double testFactor;
      MinFactor = DBL_MAX;
      MinSumSq = DBL_MAX;  
      double NumZero = 0;
      
      vecFactor.clear();
      vecSumSq.clear();

      //Loop the iteration until convergence. -> Not used, I brute force instead
      //while((abs((SumSq-SumSqPrev)/SumSq)>0.0001)&&(iterNum<1000)){
      //Loop over all reasonable factor values.
      for(factor = -10.; factor<90.; factor = factor + 0.01){
	//cout << "Factor = " << factor << ", factor - floor(factor) = "
	//	 << factor - floor(factor) << endl;

	vecFactor.push_back(factor);
	ScaledCntTest.clear();
	for(int i = 0; i<4096; i++){
	  ScaledCntTest.push_back(0);
	}
	SumSq = 0;
	NumZero = 0;

	//Shift the spectrum by 'factor', between TestIdx and RefEndIdx
	for(int i = TestIdx; i<RefEndIdx; i++){
	  
	  if(factor!=floor(factor)){
	    ScaledCntTest.at(i+floor(factor)) +=
	      (1-(factor-floor(factor)))*CntTest.at(i)*(TestWidth/RefWidth);
	    
	    ScaledCntTest.at(i+floor(factor+1)) +=
	      (factor-floor(factor))*CntTest.at(i)*(TestWidth/RefWidth);
	  }
	  else{
	    ScaledCntTest.at(i+factor) = CntTest.at(i)*(TestWidth/RefWidth);
	  }
	}
	
	//Calculate sum of squared difference.
    //for(int i = TestIdx + floor(factor*TestIdx) + 1 + 1; i<RefEndIdx; i++){
	for(int i = TestIdx + floor(factor) + 1; i<RefEndIdx; i++){
	  if(ScaledCntTest.at(i)!=0){
	    //Weighted SumSq.
	    SumSq +=
	  pow((ScaledCntTest.at(i) - CntRef.at(i))/sqrt(ScaledCntTest.at(i)),2);

	    //if((ScaledCntTest.at(i)-CntRef.at(i))==0){
	    //cout << "Difference between Test and Ref = 0!" << endl;
	    //}
	  }
	  if(CntRef.at(i) == 0){
	    NumZero++;
	  }
	}

	//Reduced ChiSq
	//SumSq = SumSq/(RefEndIdx - TestIdx - NumZero - 3);
	SumSq = SumSq/(RefEndIdx - TestIdx - floor(factor) - 2);

    //Disqualify candidates where the test data is multiplied hugely and just
    //overlaps with some large amplitude noise
	if(NumZero>10){SumSq = DBL_MAX;}
	
	if(SumSq == 0){
	  SumSq = DBL_MAX;
	}
	
	if(SumSq<MinSumSq){
	  MinFactor = factor;
	  MinSumSq = SumSq;
	  MinScaledCntTest = ScaledCntTest;
	}
	
	
	vecSumSq.push_back(SumSq);
      }

      cout << "SampleID = " << SampleID << endl;
      cout << "MinFactor = " << MinFactor << endl;
      cout << "MinSumSq = " << MinSumSq << endl;

      //Set the other variables to be written out to the tree.
      //TestTime = *(Long64_t*)GetPtrToVal(TimeBranch, k);
      //TestCntTime = *(float*)GetPtrToVal(CntTimeBranch, k);
      //TestHnum = *(float*)GetPtrToVal(HnumBranch, k);
      //TestTotCnts = *(int*)GetPtrToVal(TotCntsBranch, k);
      //TestMeasNum = *(int*)GetPtrToVal(MeasBranch, k);
      //TestDose = *(float*)GetPtrToVal(DoseBranch, k);
      //TestSampleName = *(string*)GetPtrToVal(NameBranch, k);
      //TestIDnum = *(int*)GetPtrToVal(IDnumBranch, k);
      //isRef = *(bool*)GetPtrToVal(isRefBranch, k);
      TestTime = theTime;
      TestCntTime = theCntTime;
      TestHnum = theHnum;
      TestTotCnts = theTotCnts;
      TestMeasNum = theAbsMeasNum;
      cout << "*theName = " << *theName << endl;
      TestSampleName = SampleID;//*theName;
      TestIDnum = theIDnum;
      isRef = theisRef;
      outTree->Fill();
      cout << endl;
      
    }//End of contition that tests whether entry is wanted.
  }//End of the loop on entries....

  cout << "Writing tree to file." << endl << endl;
  outTree->Write();
  fout->Close();

}
