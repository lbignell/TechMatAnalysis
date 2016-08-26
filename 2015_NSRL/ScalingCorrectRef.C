//This function will implement an algorithm to determine how much light yield
//gain/loss occurs in a scintillator upon being irradiated,
//using an unirradiated identical scintillator as a reference.
//Note that I'm assuming that both (or at least the irradiated scintillator)
//have been transferred to unirradiated vials.
//
//I will use the *next* unirradiated scintillator as the reference data set.
//(This has the disadvantage of not having a very large reference data set, but
//the advantage that instrument drift is automatically accounted for).
//
//I can assume that the sample and reference data are stored in the same ROOT
//file, and that the reference sample name is just <SampleName>_REF.
//
//To simplify the analysis, I will also impose a requirement that the data come
//from a single measurement run (seperate runs can be chained together after
//this analysis).

//Arguments:
//The root data file containing the spectral data as part of a tree (string).
//
//The sample name (string).
//
//The threshold index, which controls the minium channel the pedestal threshold
//can be set as (int).
//
//
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

void ScalingCorrectRef(string DataFile, string SampleName,
		       int ThresholdIdx){

  cout << "Opening file " << DataFile << endl;
  cout << "SampleName = " << SampleName << endl;
  cout << "ThresholdIdx = " << ThresholdIdx << endl;

  //Open the test/reference data.
  TFile* fdata = TFile::Open(DataFile.c_str());
  TTree* theTree = (TTree*)fdata->Get("LSSpec");

  string* theName = NULL;
  Long64_t theTime = 0;
  vector<double>* theChannel = NULL;
  vector<double>* theCounts = NULL;
  int thePosition = 0;
  int theSampleNum = 0;
  float theCntTime = 0;
  int theTotCnts = 0;
  int theIDnum = -100;
  bool theisRef = false;
  int theMeasNum = 0;
  float theDose = -1;
  string* theFormulation = NULL;
  int theAbsPos = -1;

  theTree->SetBranchAddress("Counts", &theCounts);
  theTree->SetBranchAddress("Channel", &theChannel);
  theTree->SetBranchAddress("Position", &thePosition);
  theTree->SetBranchAddress("TimeStamp", &theTime);
  theTree->SetBranchAddress("SampleNum", &theSampleNum);
  theTree->SetBranchAddress("CntTime", &theCntTime);
  theTree->SetBranchAddress("TotCnts", &theTotCnts);
  theTree->SetBranchAddress("MeasNum", &theMeasNum);
  theTree->SetBranchAddress("SampleName", &theName);
  theTree->SetBranchAddress("IDnum", &theIDnum);
  theTree->SetBranchAddress("Dose", &theDose);
  theTree->SetBranchAddress("isRef", &theisRef);
  theTree->SetBranchAddress("Formulation", &theFormulation);
  theTree->SetBranchAddress("AbsPos", &theAbsPos);

  Long64_t numEntries = theTree->GetEntries();

  //Now I need to go through and get the event numbers of the count and ref
  //measurements.
  string thisSampleName = "";
  int loopCtr = 0;
  vector<int> SampleEvts;
  vector<int> RefEvts;
  int NumSamples = 0;
  int NumRefs = 0;

  while(loopCtr<numEntries){
    //
    theTree->GetEntry(loopCtr);

    thisSampleName = *theName;
    //cout << "thisSampleName = " << thisSampleName << endl;
    if(thisSampleName==SampleName){
      //cout << "Found SampleName = " << thisSampleName << endl;
      NumSamples++;
      SampleEvts.push_back(loopCtr);
    }
    else if(thisSampleName==(SampleName+"_REF")){
      NumRefs++;
      RefEvts.push_back(loopCtr);
    }
    loopCtr++;
  }

  //NumRefs should equal NumSamples, or NumSamples - 1 (if the run finished on
  //the sample).
  if(NumRefs==NumSamples){cout << "NumRefs = NumSamples = " << NumRefs << endl;}
  else if(NumRefs==(NumSamples-1))
    {cout << "NumRefs = " << NumRefs 
	  << ", NumSamples = " << NumSamples << endl;}
  else{cout << "WARNING!!! NumRefs = " << NumRefs 
	    << ", NumSamples = " << NumSamples << endl;}

  //Some vectors to store the ref data.
  vector<double> CntRef(4096,0);
  vector<double> Channel(4096,0);
  double RefTotCnts = 0;

  //Define a vector for storing the test data
  vector<double> CntTest(4096,0);
  Long64_t TestTime = -1;
  float TestCntTime = -1;
  int TestTotCnts = -1;
  int TestMeasNum = -1;
  string TestSampleName = "";//default to skip state
  int TestIDnum = -2;//default to skip state
  bool isRef = false;

  //Added these
  float TestDose = -1;
  int TestPosition = -1;
  int TestSampleNum = -1;
  string TestFormulation = "";
  int TestAbsPos = -1;

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

  //Define a ROOT file (with tree) for saving the results.
  //Save results in a tree and save to disk.
  TFile* fout = new TFile("ScalingCorrectRefResults.root", "RECREATE");
  TTree* outTree = new TTree("Results","Light yield scaling algo results");
  outTree->Branch("Channel", "vector<double>", &Channel);
  outTree->Branch("Position", &TestPosition, "Position/I");
  outTree->Branch("TimeStamp", &TestTime, "TimeStamp/L");
  outTree->Branch("SampleNum", &TestSampleNum, "SampleNum/I");
  outTree->Branch("CntTime", &TestCntTime, "CntTime/F");
  outTree->Branch("TotCnts", &TestTotCnts, "TotCnts/I");
  outTree->Branch("MeasNum", &TestMeasNum, "MeasNum/I");
  outTree->Branch("SampleName", "string", &TestSampleName);
  outTree->Branch("IDnum", &TestIDnum, "IDnum/I");
  outTree->Branch("Dose", &TestDose, "Dose/F");
  outTree->Branch("isRef", &isRef, "isRef/O");
  outTree->Branch("Formulation", "string", &TestFormulation);
  outTree->Branch("AbsPos", &TestAbsPos, "AbsPos/I");

  outTree->Branch("TestCounts", "vector<double>", &CntTest);
  outTree->Branch("RefCounts", "vector<double>", &CntRef);
  outTree->Branch("RefTotCnts", &RefTotCnts, "RefTotCnts/I");
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

  //Loop over the number of reference measurements.
  for(size_t k = 0; k<RefEvts.size(); k++){

    //Get the reference event spectrum.
    theTree->GetEntry(RefEvts.at(k));
    CntRef.swap(*theCounts);
    RefTotCnts = theTotCnts;

    //Now get the test event spectrum and fill some tree variables.
    theTree->GetEntry(SampleEvts.at(k));
    CntTest.swap(*theCounts);
    Channel.swap(*theChannel);

    TestTotCnts = theTotCnts;
    TestSampleName = *theName;
    TestIDnum = theIDnum;
    TestPosition = thePosition;
    TestSampleNum = theSampleNum;
    TestTime = theTime;
    TestCntTime = theCntTime;
    TestMeasNum = theMeasNum;
    TestDose = theDose;
    isRef = theisRef;
    TestFormulation = *theFormulation;
    TestAbsPos = theAbsPos;
    
    //Test whether to proceed or loop to the next entry
    if(TestSampleName==SampleName){
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


      //Calculate the difference/sigma
      vector<double>::iterator testit = CntTest.begin();
      vector<double>::iterator refit = CntRef.begin();
      vector<double>::iterator chanit = Channel.begin();
      int counter = 0;
      tVals.clear();
      pVals.clear();

      //significance criterion
      double Criterion = 0.05;  
      bool IsStatSig = false;
      int numBelowCriterion = 0;
      int lengthBelowCriterion = 5;//semi-arbitrary parameter, to help prevent
      //triggering on a random fluctuation.
      int AvgLength = 3;

      int RefEndIdx = -1;
      RefStart = -1;
      RefEnd = -1;
      TestStart = -1;
      TestEnd = -1;
      TestIdx = 0;

      double integralRatio = RefTotCnts/TestTotCnts;
      //scale the test vector by the ratio of integrated spectra:
      std::transform(CntTest.begin(), CntTest.end(), CntTest.begin(),
		     std::bind1st(std::multiplies<double>(),integralRatio));

      while((refit!=CntRef.end())&&(testit!=CntTest.end())){
	stdev.push_back(sqrt(*testit));
	if((counter<ThresholdIdx)||(*testit<10)||(*refit<10)){/*do nothing*/
	  tVals.push_back(0);
	  pVals.push_back(-1);//an unrealistic value to reset the NumBelow
	}
	else{
	  //Calculate the t-value at this point.
	  //sigma = sqrt(N)
	  //tVals.push_back((*refit - *testit)/sqrt(*testit));
	  tVals.push_back(((*refit-*testit)*(*refit - *testit))/(*testit));

	  //double thepVal = TMath::StudentI(tVals.back(),*testit);//technically, we
	  //double thepVal = TMath::Prob(tVals.back(), 1);
	  double SumChi2 = 0;
	  //vector<double>::reverse_iterator tvalsit = tVals.rbegin();
	  for(int i = -AvgLength; i<0; i++){
	    SumChi2 += tVals.at(counter+i+1);
	  }
	  double thepVal = TMath::Prob(SumChi2, AvgLength-1);
	  //have 0 DoF...

	  //cout << "Channel = " << counter << ", Counts = " << *testit
	  //   << ", difference = " << (*testit - *refit)
	  //   << ", tval = " << tVals.back()
	  //   << ", pval = " << thepVal
	  //   << ", SumChi2 = " << SumChi2 << endl;

	  if((thepVal!=-1)&&(thepVal<Criterion)&&(TestIdx==0)){
	    numBelowCriterion++;
	    //cout << "pval below criterion!, numBelowCriterion = "
	    //	 << numBelowCriterion << endl;
	    if(numBelowCriterion>lengthBelowCriterion){
	      TestIdx = counter-lengthBelowCriterion;
	      TestChannel = Channel.at(counter-lengthBelowCriterion);
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

	if(*testit>100){
	  TestEnd = Channel.at(counter);
	  if(TestStart==-1){TestStart = Channel.at(counter);}
	}

	if(*refit>100){
	  RefEnd = Channel.at(counter);
	  RefEndIdx = counter;
	  if(RefStart==-1){RefStart = Channel.at(counter);}
	}

	testit++;
	refit++;
	counter++;
      }

      if(IsStatSig==false){
	cout << "Data show no statistically significant difference!" << endl;
      }

      //Here's where I forcibly re-set the obtained threshold if it is too close
      //to the edge estimate.
      if(TestIdx<ThresholdIdx){
	cout << "TestIdx = " << TestIdx << "! Shifting to " << ThresholdIdx
	     << endl;
	TestIdx = ThresholdIdx;
	TestChannel = Channel.at(TestIdx);
      }
      //Check that there are 'enough' bins to do the Chi2 over.
      else if(((TestEnd - TestIdx)<40)&&(theIDnum<6)&&(theIDnum!=4)){
	cout << "Fit range is not wide enough! TestIdx = " << TestIdx 
	     << ", TestEnd = " << TestEnd << ". Shifting TestIdx to ";
	TestIdx = TestIdx - (40 - (TestEnd-TestIdx));
	if(TestIdx<ThresholdIdx){TestIdx=ThresholdIdx;}
	cout << TestIdx << "." << endl;
      }
      else if(((TestEnd - TestIdx)<60)&&(theIDnum==10)){
	cout << "Fit range is not wide enough! TestIdx = " << TestIdx 
	     << ", TestEnd = " << TestEnd << ". Shifting TestIdx to ";
	TestIdx = TestIdx - (60 - (TestEnd-TestIdx));
	cout << TestIdx << "." << endl;
      }
      else if(((TestEnd - TestIdx)<80)&&(theIDnum==10)||(theIDnum==14)){
	cout << "Fit range is not wide enough! TestIdx = " << TestIdx 
	     << ", TestEnd = " << TestEnd << ". Shifting TestIdx to ";
	TestIdx = TestIdx - (80 - (TestEnd-TestIdx));
	cout << TestIdx << "." << endl;
      }
      else if(((TestEnd - TestIdx)<425)&&(theIDnum==100)||(theIDnum==4)){
	//else if(((theIDnum==100)||(theIDnum==4))){
	cout << "Fit range is not wide enough! TestIdx = " << TestIdx 
	     << ", TestEnd = " << TestEnd << ". Shifting TestIdx to ";
	TestIdx = 300;
	//TestIdx = (425 - (TestEnd-TestIdx));
	if(TestIdx<ThresholdIdx){TestIdx=ThresholdIdx;}
	cout << TestIdx << "." << endl;
      }
      else{ cout << "TestIdx found!" << endl
		 << "TestIdx = " << TestIdx << endl; }

      //Protect against TestIdx going negative!
      if(TestIdx<ThresholdIdx){
	cout << "Whoops! TestIdx = " << TestIdx 
	     << "! Shifting to " << ThresholdIdx << endl;
	TestIdx = ThresholdIdx;
	TestChannel = Channel.at(TestIdx);
      }


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
      MinScaledCntTest.clear();
      int NumBins = 0;
      //double SmallStepFactor = 0.15;
      //double GradDescentFactor = 0.0000001;
      
      //double derivative = 1;
      //double derivativePrev;
      //double dblderivative = -1000;
      
      //int iterNum = 0;
      
      //double testFactor;
      
      MinFactor = DBL_MAX;
      MinSumSq = DBL_MAX;  
      
      double NumZero = 0;
      
      int NumPtsUsed = 0;
      
      vecFactor.clear();
      vecSumSq.clear();
      
      //vector<int> NumRpts;
      
      //Loop the iteration until convergence.
      //while((abs((SumSq-SumSqPrev)/SumSq)>0.0001)&&(iterNum<1000)){
      //Loop over all reasonable factor values.
      for(factor = 0.9; factor<1.5; factor = factor + 0.0001){

	vecFactor.push_back(factor);
	//NumRpts.clear();
	ScaledCntTest.clear();
	for(int i = 0; i<4096; i++){
	  ScaledCntTest.push_back(0);
	  //NumRpts.push_back(0);
	}

	SumSq = 0;
	NumZero = 0;
	NumPtsUsed = 0;
	
	//multiply the spectrum by 'factor', between TestIdx and RefEndIdx
	for(int i = TestIdx; i<RefEndIdx; i++){
	  //I can assume unit-spaced bins.
	  NumBins = 1 + (floor((i+1)*factor) - floor(i*factor));
	  for(int j = 0; j < NumBins; j++){
	    //Calculate the increases in the bins.
	    if(j==0){
	      //Fill 1st bin
	      ScaledCntTest.at(floor(factor*i)) += 
		(1/factor)*(1 + floor(factor*i) - factor*i)*CntTest.at(i);
	      //NumRpts.at(floor(factor*i))++;
	    }else if(j == NumBins-1){
	      ScaledCntTest.at(floor(factor*(i+1))) += 
		(1/factor)*(factor*(i + 1) - floor(factor*(i+1)))*
		CntTest.at(i);
	      //NumRpts.at(floor(factor*i + j))++;
	    }else{
	      ScaledCntTest.at(floor(factor*i + j)) += 
		(1/factor)*CntTest.at(i);
	      //NumRpts.at(floor(factor*i + j))++;
	    }
	  }
	}

	//Calculate sum of squared difference.
	for(int i=TestIdx + floor(factor*TestIdx) + 1 + 1; i<RefEndIdx; i++){
	//for(int i = TestIdx + 5; i<RefEndIdx; i++){
	  if(ScaledCntTest.at(i)!=0){
	    //if(factor>1){
	      //if(NumRpts.at(i)>1){//This is to exclude points that have only
	      //been added to once (the 'turn on' of the algorithm).
		SumSq += pow((ScaledCntTest.at(i) - CntRef.at(i))
			     /sqrt(ScaledCntTest.at(i)),2);
		NumPtsUsed++;
		//cout << "SumSq = " << SumSq << ", ScaledCntTest.at(" << i
		//   << ") = " << ScaledCntTest.at(i) << ", CntRef.at(" << i 
		//   << ") = " << CntRef.at(i) << endl;
		//}
		//}
		//else{
		//SumSq += pow((ScaledCntTest.at(i) - 
		//	    CntRef.at(i))/sqrt(ScaledCntTest.at(i)),2);
		//NumPtsUsed++;
		//}
	    if((ScaledCntTest.at(i)-CntRef.at(i))==0){
	      cout << "Difference between Test and Ref = 0!" << endl;
	    }
	  }
	  if(CntRef.at(i) == 0){NumZero++;}
	}
	
	//Reduced ChiSq
	SumSq = SumSq/(NumPtsUsed - 1);
	
	//Disqualify candidates where the test data is multiplied hugely and
	//overlaps with some large amplitude noise
	if(NumZero>10){
	  //cout << "NumZero = " << NumZero << ">10!" << endl;
	  SumSq = DBL_MAX;
	}
	
	if(SumSq == 0){
	  //cout << "SumSq = 0; setting SumSq to DBL_MAX" << endl;
	  //SumSq = SumSqPrev*SumSqPrev;
	  SumSq = DBL_MAX;
	}
	
	if(SumSq<MinSumSq){
	  MinFactor = factor;
	  MinSumSq = SumSq;
	  MinScaledCntTest = ScaledCntTest;
	}
	
	vecSumSq.push_back(SumSq);
	
	//Determine whether the change was an improvement using Newton's method.
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
 
      //Implement the factor correction here:
      double CorrectedFactor = 0;
      if(theDose<15){
	//5% WbLS
	CorrectedFactor = MinFactor/
	  (1.32705 - 0.518823*MinFactor + 0.191817*MinFactor*MinFactor);
      }
      else if(theDose>=15){
	//Pure LS
	CorrectedFactor = MinFactor/
	  (1.11252 - 0.18939*MinFactor + 0.0768342*MinFactor*MinFactor);
      }
      else{CorrectedFactor=MinFactor;}
     
      cout << "SampleName = " << SampleName << endl;
      cout << "MinFactor = " << MinFactor << endl;
      cout << "CorrectedFactor = " << CorrectedFactor << endl;
      cout << "Reduced ChiSq = " << MinSumSq << endl;

      //MinFactor is what's actually written to the tree, so overwrite this.
      MinFactor = CorrectedFactor;

      outTree->Fill();
      cout << endl;     
 
    }//End of contition that tests whether entry is wanted.
  }//End of the loop on entries....

  outTree->Write();
  fout->Close();

}
