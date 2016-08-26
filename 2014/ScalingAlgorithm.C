//This function will implement an algorithm to determine how much light yield
//gain/loss has been imparted upon a liquid scintillator upon being irradiated,
//by using unirradiated scintillator data as a reference.

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

void ScalingAlgorithm(string DataFile, string BranchName, string SampleID, 
		      int ThresholdIdx, int NumAvg){

  if(!((BranchName=="SampleName")||(BranchName=="IDnum"))){
    cout << "ERROR: this code only works with the second argument equal to"
	 << " 'SampleName' or 'IDnum'. You wrote: " << BranchName << "." << endl
	 << "Nothing has been done, exiting..." << endl;
    return;
  }

  //Open the test/reference data.
  TFile* fdata = TFile::Open(DataFile.c_str());
  //TTree* TestTree = (TTree*)ftest->Get("LSspec");
  TChain* theTree = (TChain*)fdata->Get("LSChain");
  
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
  theTree->SetBranchAddress("AbsMeasNum", &theAbsMeasNum);
  theTree->SetBranchAddress("CntTime", &theCntTime);
  theTree->SetBranchAddress("Hnum", &theHnum);
  theTree->SetBranchAddress("TotCnts", &theTotCnts);
  theTree->SetBranchAddress("SampleName", &theName);
  theTree->SetBranchAddress("IDnum", &theIDnum);
  theTree->SetBranchAddress("isRef", &theisRef);
  theTree->SetBranchAddress("Position", &thePosition);
  theTree->SetBranchAddress("SampleNum", &theSampleNum);
 
  TBranch* TestBranch = (TBranch*)TestTree->GetBranch("Counts");
  vector<double> CntTest = *(vector<double>*)GetPtrToVal(TestBranch, 0);

  TFile* fref = TFile::Open(RefData.c_str());
  TTree* RefTree = (TTree*)fref->Get("LSspec");
  TBranch* RefBranch = (TBranch*)RefTree->GetBranch("Counts");
  TBranch* ChnBranch = (TBranch*)RefTree->GetBranch("Channel");
  vector<double> CntRef = *(vector<double>*)GetPtrToVal(RefBranch, 0);
  vector<double> Channel = *(vector<double>*)GetPtrToVal(ChnBranch, 0);

  //cout << "Read data in successfully!" << endl;

  if((CntTest.size()!=CntRef.size())){
    cout << "ERROR: Ref and Test data are not equal lengths" << endl
	 << "Ref size = " << CntRef.size() << ", Test size = " << CntTest.size()
	 << endl;
    return;
  }

  //Calculate the difference/sigma
  vector<double>::iterator testit = CntTest.begin();
  vector<double>::iterator refit = CntRef.begin();
  vector<double>::iterator chanit = Channel.begin();
  vector<double> NormDiff;
  vector<double> CumNormDiff;
  double tally = 0;
  double RefStart = -1;
  double RefEnd = -1;
  int RefEndIdx = 0;
  double TestStart = -1;
  double TestEnd = -1;
  int theIndex = 0;

  while((refit!=CntRef.end())&&(testit!=CntTest.end())){
    //get the 'normalised difference' values.
    //cout << "Getting 'Normalised Difference' cumsum" << endl;
    //cout << "*refit = " << *refit << ", *testit = " << *testit << endl;
    if(*refit!=0){
      NormDiff.push_back((*testit - *refit)/sqrt(*refit));
      //Alternatively I can divide by sigma squared...
      //NormDiff.push_back((*testit - *refit)/(*refit));
      //cout << "NormDiff = " << NormDiff.back() << endl;
      //cout << "testit = " << *testit << endl;
      //cout << "refit = " << *refit << endl;
      tally += (*testit - *refit)/sqrt(*refit);
      //Alternatively I can divide by sigma squared...
      //tally += (*testit - *refit)/(*refit);
      CumNormDiff.push_back(tally);
    }    
    else{
      NormDiff.push_back(0);
      CumNormDiff.push_back(tally);  
    }

    //Also get the start and end points of each spectrum for renormalisation.
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
  vector<double> mean;
  vector<double> stdev;
  double theMean;
  double theStdev;
  //Number of samples over which moving avg is taken.
  int avglen = 3;
  //Threshold, below which the test statistic isn't examined.
  int threshold = avglen + 0;//i'll set to min for now, may need changing later.
  //p-test criterion, make it 1%.
  double Criterion = 0.01;  

  vector<double> tVals;
  vector<double> pVals;
  double thepVal;
  int TestIdx = 0;
  int TestChannel = 0;
  int numBelowCriterion = 0;
  int lengthBelowCriterion = 30;

  for(int i = 0; i<NormDiff.size(); i++){
    if((i<threshold)||(i>(NormDiff.size()-avglen))){/*do nothing*/
      mean.push_back(0);
      stdev.push_back(1);
      tVals.push_back(0);
      pVals.push_back(-1);
      //cout << "doing nothing!" << endl;
    }
    else{
      theMean = 0;
      //theStdev = 0;
      //cout << "avglen = " << avglen << endl
      //   << "floor(avglen/2) = " << floor(avglen/2) << endl
      //   << "i = " << i << endl
      //   << "CumNormDiff.size() = " << CumNormDiff.size()
      //   << endl;
      theMean = (CumNormDiff.at(i+floor(avglen/2)) - 
		 CumNormDiff.at(i-avglen+floor(avglen/2)))/avglen;
      //cout << "theMean = " << theMean << ", i = " << i << endl;
      //for(int j = 0; j<avglen; j++){
      //theStdev += (1/(avglen-1))*pow((NormDiff.at(i-j)-theMean),2);
      //}
      //cout << "theStdev = " << theStdev << endl;
      
      //Stdev should equal ~ 1.
      theStdev = 1;
      tVals.push_back(sqrt(avglen)*(theMean/theStdev));
      mean.push_back(theMean);
      stdev.push_back(theStdev);
      thepVal = TMath::StudentI(tVals.back(), (avglen-1));
      //      if(i<100){
      //cout << "i = " << i << ", theMean = " << theMean << ", theStdev = "
      //     << theStdev << ", thepVal = " << thepVal << endl;
      //}
      if((thepVal!=-1)&&(thepVal<Criterion)&&(TestIdx==0)){
	numBelowCriterion++;
	//cout << "NumBelowCriterion = " << numBelowCriterion << endl
	//   << "thepVal = " << thepVal << endl
	//   << "i = " << i << endl
	//   << "Channel.at(i) = " << Channel.at(i) << endl;
	if(numBelowCriterion>lengthBelowCriterion){
	  TestIdx = i-lengthBelowCriterion;
	  TestChannel = Channel.at(i-lengthBelowCriterion);
	  cout << "TestIdx = " << TestIdx << endl
	       << "TestChannel = " << TestChannel << endl;
	}
      }
      else{
	numBelowCriterion = 0;
      }
      pVals.push_back(thepVal);
    }
  }

  if(TestIdx<10){
    cout << "TestIdx = " << TestIdx << "! Shifting to 10" << endl;
    TestIdx = 10;
    TestChannel = Channel.at(TestIdx);
  }

  //cout << "TestIdx = " << TestIdx << endl;
  //cout << "TestChannel = " << TestChannel << endl;
  
  /////////////////////////////////////////////////////////////////////////////
  //Insert code to scale test data and perform minimisation here.
  /////////////////////////////////////////////////////////////////////////////
  //First guess of scale factor.
  double factor = 1.05;
  double SumSq = 0;
  double SumSqPrev = 0;
  double factorPrev = 1.03;
  //Loop infinitely until minimisation converges.
  vector<double> ScaledCntTest;
  vector<double> MinScaledCntTest;
  int NumBins = 0;
  double SmallStepFactor = 0.15;
  double GradDescentFactor = 0.0000001;

  double derivative = 1;
  double derivativePrev;
  double dblderivative = -1000;

  int iterNum = 0;

  double testFactor;

  double MinFactor = DBL_MAX;
  double MinSumSq = DBL_MAX;  

  double NumZero = 0;

  int NumPtsUsed = 0;

  vector<double> vecFactor;
  vector<double> vecSumSq;

  vector<int> NumRpts;

  //Run through once to get the initial value of the SumSqPrev so the first
  //derivative estimate can work
  //  for(int k = 0; k<2; k++){
  //if(k==0){testFactor = factorPrev;}
  //else{testFactor = factor;}

  //ScaledCntTest.clear();
  //for(int i = 0; i<4096; i++){
  //ScaledCntTest.push_back(0);
  //}
  //for(int i = TestIdx; i<RefEndIdx; i++){
    //I can assume unit-spaced bins.
    //  NumBins = 1 + floor((i*testFactor - floor(i*testFactor))+testFactor);
  //for(int j = 0; j < NumBins; j++){
  //	//Calculate the increases in the bins.
  //	if(j==0){
  //	  ScaledCntTest.at(floor(testFactor*i + j)) += 
  //	    (1/testFactor)*(1 + floor(testFactor*i) - testFactor*i)*CntTest.at(i);
  //	}else if(j == NumBins-1){
  //	  ScaledCntTest.at(floor(testFactor*i + j)) += 
  //	    (1/testFactor)*(testFactor*(i + 1) - floor(testFactor*(i+1)))*
  //	    CntTest.at(i);
  //	}else{
  //	  ScaledCntTest.at(floor(testFactor*i + j)) += 
  //	    (1/testFactor)*CntTest.at(i);
  //	}
  //}
  //}
  //for(int i = TestIdx + 1; i<RefEndIdx; i++){
  //if(ScaledCntTest.at(i)!=0){
  //	if(k==0){SumSqPrev += pow((ScaledCntTest.at(i) - CntRef.at(i)),2);}
  //	else{SumSq += pow((ScaledCntTest.at(i) - CntRef.at(i)),2);}
  //}
  //}
  //}

  //Now determine an estimate of derivative, etc.
  //derivativePrev = (SumSq - SumSqPrev)/(factor - factorPrev);

  //  factor = factor - GradDescentFactor*derivativePrev;

  //cout << "Initial values: " << endl
    //<< "factor = " << factor << endl
    // << "factorPrev = " << factorPrev << endl
    //<< "SumSq = " << SumSq << endl
    //<< "SumSqPrev = " << SumSqPrev << endl
    // << "derivativePrev = " << derivativePrev << endl;
  
  //Loop the iteration until convergence.
  //while((abs((SumSq-SumSqPrev)/SumSq)>0.0001)&&(iterNum<1000)){
  //Loop over all reasonable factor values.
  for(factor = 0.9; factor<1.5; factor = factor + 0.00001){
    //cout << "Optimisation: iteration # " << iterNum << endl;
    //cout << "Initial Factor = " << factor << endl;

    vecFactor.push_back(factor);

    NumRpts.clear();
    ScaledCntTest.clear();
    for(int i = 0; i<4096; i++){
      ScaledCntTest.push_back(0);
      NumRpts.push_back(0);
    }
    //if(iterNum!=0){
    //SumSqPrev = SumSq;
    //derivativePrev = derivative;
    //}
    SumSq = 0;
    NumZero = 0;
    NumPtsUsed = 0;

    //if(factor<0){ factor = 1;}
    //multiply the spectrum by 'factor', between TestIdx and RefEndIdx
    for(int i = TestIdx; i<RefEndIdx; i++){
      //I can assume unit-spaced bins.
      //Old way
      //NumBins = ceil(factor*(i+1)) - floor(factor*i);
      //cout << "NumBins old = " << NumBins << endl;
      NumBins = 1 + floor((i*factor - floor(i*factor))+factor);
      for(int j = 0; j < NumBins; j++){
	//Calculate the increases in the bins.
	if(j==0){
	  ScaledCntTest.at(floor(factor*i + j)) += 
	    (1/factor)*(1 + floor(factor*i) - factor*i)*CntTest.at(i);
	  NumRpts.at(floor(factor*i + j))++;
	}else if(j == NumBins-1){
	  ScaledCntTest.at(floor(factor*i + j)) += 
	    (1/factor)*(factor*(i + 1) - floor(factor*(i+1)))*
	    CntTest.at(i);
	  NumRpts.at(floor(factor*i + j))++;
	}else{
	  ScaledCntTest.at(floor(factor*i + j)) += 
	    (1/factor)*CntTest.at(i);
	  NumRpts.at(floor(factor*i + j))++;
	}
      }
      //simplistic version for testing.
      //ScaledCntTest.at(floor(factor*i)) += CntTest.at(i);
      //cout << "i = " << i << ", factor = " << factor
      //   << ", floor(factor*i) = " << floor(factor*i) << endl;
      //cout << "CntTest.at(i) = " << CntTest.at(i) << endl;
    }
    //Calculate sum of squared difference.
    //for(int i = TestIdx + floor(factor*TestIdx) + 1 + 1; i<RefEndIdx; i++){
    for(int i = TestIdx + 1 + 1; i<RefEndIdx; i++){
      if(ScaledCntTest.at(i)!=0){
	if(factor>1){
	  if(NumRpts.at(i)>1){
	    SumSq +=
	      pow((ScaledCntTest.at(i) - CntRef.at(i))/sqrt(ScaledCntTest.at(i))
		  ,2);
	    NumPtsUsed++;
	  }
	}
	else{
	  SumSq +=
	    pow((ScaledCntTest.at(i) - CntRef.at(i))/sqrt(ScaledCntTest.at(i))
		,2);
	  NumPtsUsed++;
	}
	if((ScaledCntTest.at(i)-CntRef.at(i))==0){
	  cout << "Difference between Test and Ref = 0!" << endl;
	}
      }
      if(CntRef.at(i) == 0){
	NumZero++;
      }
    }

    //Reduced ChiSq
    SumSq = SumSq/(NumPtsUsed - 1);

    //Disqualify candidates where the test data is multiplied hugely and just
    //overlaps with some large amplitude noise
    if(NumZero>10){SumSq = DBL_MAX;}

    if(SumSq == 0){
      //cout << "SumSq = 0; setting SumSq to 1000" << endl;
      //SumSq = SumSqPrev*SumSqPrev;
      SumSq = DBL_MAX;
    }
    
    if(SumSq<MinSumSq){
      MinFactor = factor;
      MinSumSq = SumSq;
      MinScaledCntTest = ScaledCntTest;
    }

    vecSumSq.push_back(SumSq);
   
    //Determine whether the change was an improvement, using Newton's method.
    //derivative = (SumSq - SumSqPrev)/(factor - factorPrev);
    //cout << "SumSq = " << SumSq << endl;
    //cout << "Derivative = " << derivative << endl;
    //dblderivative = (derivative - derivativePrev)/(factor - factorPrev);
    //cout << "dblderivative = " << dblderivative << endl;
    //factorPrev = factor;
    //Netwon's Method.
    //factor = factor - SmallStepFactor*(derivative/dblderivative);
    //Gradient Descent.
    //factor = factor - GradDescentFactor*derivative;
    //cout << "New Factor = " << factor << endl;
    //iterNum++;
    //cout << "SumSq = " << SumSq << endl;
    //cout << "iterNum = " << iterNum << endl;
  }
  cout << "Factor = " << MinFactor << endl;
  cout << "MinSumSq = " << MinSumSq << endl;
  //  cout << "Number of Iterations = " << iterNum << endl;

  //cout << "Optimisation complete. Factor = " << factor << endl;

  //Save results in a tree and save to disk.
  TFile* fout = new TFile("ScalingResults.root", "RECREATE");
  TTree* outTree = new TTree("Results","Light yield scaling algorithm results");
  outTree->Branch("Channel", "vector<double>", &Channel);
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
  outTree->Fill();
  outTree->Write();
  fout->Close();

}
