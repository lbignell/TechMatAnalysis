//This function will implement an algorithm to determine how much light yield
//gain/loss has been imparted upon a liquid scintillator upon being irradiated,
//by using unirradiated scintillator data as a reference.
//
//////////////////////////////////////////////////////////////////////////////
//Data are stored in ROOT files in subdirectories, so I want to have the file
//name as the first argument.
//
//The files store the sample names in a couple of branches; SampleName (a
//string) and IDnum (an integer). SampleName specifies individual vials, whereas
//IDnum specifies a type of sample composition (eg 5% WbLS).
//I want to be able to choose whether to access data by reference to the
//SampleName or IDnum, so I'll take a string that corresponds to the appropriate
//branch name as the second argument.
//IDnum reference:
//-1 = unknown/indeteriminate
//0 = empty vial
//1 = irradiated (or ref) glass vial filled with fresh LS.
//2 = irradiated (or ref) PP vial filled with fresh LS.
//3 = irradiated (or ref) HDPE vial filled with fresh LS.
//4 = H-3 or C-14 (or other) standard vial.
//5 = 5% WbLS
//10 = 10% WbLS
//14 = 14% WbLS
//100 = Pure LS
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

void ScalingAlgorithmBenchmark(string DataFile, int ThresholdIdx){

  cout << "Opening file " << DataFile << endl;
  cout << "ThresholdIdx = " << ThresholdIdx << endl;

  //Open the test/reference data.
  TFile* fdata = TFile::Open(DataFile.c_str());
  TTree* theTree = (TTree*)fdata->Get("Results");

  vector<double>* theChannel = NULL;
  vector<double>* theCounts = NULL;
  double theTrueFactor = -1;

  theTree->SetBranchAddress("Spectrum", &theCounts);
  theTree->SetBranchAddress("Energy", &theChannel);
  theTree->SetBranchAddress("Factor", &theTrueFactor);

  Long64_t numEntries = theTree->GetEntries();

  //Some vectors to store the ref data.
  vector<double> CntRef;
  vector<double> Channel;
  bool isIDnum = false;

  int NumSamples = 0;
  double RefTotCnts = 0;

  Long64_t loopCtr = 0;

  //Now find the Cnt and Chn values for theTrueFactor==1.
  while(loopCtr<numEntries){
    theTree->GetEntry(loopCtr);
    cout << "theTrueFactor = " << theTrueFactor << endl;
    if((theTrueFactor<1.000000001)&&(0.99999999999<theTrueFactor)){
	CntRef.swap(*theCounts);
	Channel.swap(*theChannel);
	break;
    }
    loopCtr++;
  }

  //Define a vector for storing the test data
  vector<double> CntTest;
  double TrueFactorTest;

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

  double ThresholdLY = -1;

  //Define a ROOT file (with tree) for saving the results.
  //Save results in a tree and save to disk.
  TFile* fout = new TFile("ScalingResults.root", "RECREATE");
  TTree* outTree = new TTree("Results","Light yield scaling algo results");
  outTree->Branch("Channel", "vector<double>", &Channel);
  outTree->Branch("TrueFactor", &TrueFactorTest, "TrueFactor/D");
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
  outTree->Branch("ThresholdLY", &ThresholdLY, "ThresholdLY/D");

  //loop over the entries, comparing with the reference data when a match is
  //found
  for(int k = 0; k<numEntries; k++){
    theTree->GetEntry(k);

    CntTest.swap(*theCounts);
    
    if((CntTest.size()!=CntRef.size())){
      cout << "ERROR: Ref and Test data are not equal lengths" << endl
	   << "Ref size = " << CntRef.size() << ", Test size = " 
	   << CntTest.size() << endl
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
    
    //double integralRatio = 1;//RefTotCnts/TestTotCnts;
    //scale the test vector by the ratio of integrated spectra:
    //std::transform(CntTest.begin(), CntTest.end(), CntTest.begin(),
    //		   std::bind1st(std::multiplies<double>(),integralRatio));

    while((refit!=CntRef.end())&&(testit!=CntTest.end())){
      stdev.push_back(sqrt(*testit));
      if((counter<ThresholdIdx)||(*testit<10)||(*refit<10)){/*do nothing*/
	tVals.push_back(0);
	pVals.push_back(-1);//an unrealistic value to reset the NumBelow
      }
      else{
	//Calculate the t-value at this point.
	tVals.push_back(((*refit-*testit)*(*refit - *testit))/(*testit));
	double SumChi2 = 0;
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
      
      if(*testit>1000){
	TestEnd = Channel.at(counter);
	if(TestStart==-1){TestStart = Channel.at(counter);}
      }
      
      if(*refit>1000){
	RefEnd = Channel.at(counter);
	RefEndIdx = counter;
	if(RefStart==-1){RefStart = Channel.at(counter);}
      }
      
      testit++;
      refit++;
      counter++;
    }
    
    double TestWidth = TestEnd - TestStart;
    double RefWidth = RefEnd - RefStart;
    
    if(IsStatSig==false){
      cout << "Data show no statistically significant difference!" << endl;
    }
    
    if(TestIdx<ThresholdIdx){
      cout << "TestIdx = " << TestIdx << "! Shifting to " << ThresholdIdx
	   << endl;
      TestIdx = ThresholdIdx;
      TestChannel = Channel.at(TestIdx);
    }
    //Check that there are 'enough' bins to do the Chi2 over.
    else if((TestEnd - TestIdx)<400){//&&(theIDnum<6)&&(theIDnum!=4)){
      TestIdx = TestIdx - (400 - (TestEnd-TestIdx));
      if(TestIdx<ThresholdIdx){TestIdx=ThresholdIdx;}
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
    double SmallStepFactor = 0.15;
    double GradDescentFactor = 0.0000001;
    
    double derivative = 1;
    double derivativePrev;
    double dblderivative = -1000;
    
    int iterNum = 0;
    
    double testFactor;
    
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
    for(factor = 0.9; factor<2; factor = factor + 0.0001){
      
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
	//Old way
	//NumBins = ceil(factor*(i+1)) - floor(factor*i);
	//cout << "NumBins old = " << NumBins << endl;
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
	//simplistic version for testing.
	//ScaledCntTest.at(floor(factor*i)) += CntTest.at(i);
	//cout << "i = " << i << ", factor = " << factor
	//   << ", floor(factor*i) = " << floor(factor*i) << endl;
	//cout << "CntTest.at(i) = " << CntTest.at(i) << endl;
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
    
    cout << "MinFactor = " << MinFactor << endl;
    cout << "Reduced ChiSq = " << MinSumSq << endl;
    
    //Set the other variables to be written out to the tree.
    TrueFactorTest = theTrueFactor;
    ThresholdLY = TestEnd/RefEnd;
    cout << "Scaling LY = " << 1/MinFactor 
      	 << ", ThresholdLY = " << ThresholdLY << endl;

    outTree->Fill();
    cout << endl;     

}//End of the loop on entries....

  //cout << "Writing tree to file." << endl << endl;
  outTree->Write();
  fout->Close();

}
