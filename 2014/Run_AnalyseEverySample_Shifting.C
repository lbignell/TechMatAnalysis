//THIS IS A MACRO, DON'T COMPILE
{
  TFile* f = TFile::Open("Collected_Pre-irrad_01_141014.root");
  //TChain* theChain = (TChain*)f->Get("LSSpec");
  TTree* theChain = (TTree*)f->Get("LSSpec");
  cout << "theChain = " << theChain << " has " << theChain->GetEntries()
       << " entries." << endl;
  theChain->Process("AnalyseEverySample_Shifting.C+");
}
