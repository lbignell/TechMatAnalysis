//This macro will call all of the necessary functions for analysis.
#include "TROOT.h"
#include <string>

using namespace std;

void RunAll()
{
  string command = ".x Run_Scaling_AllSamples.C+(\"CollectedNSRLAll2.root\"," +
    std::string("\"SampleNamesNSRL.txt\")");
  gROOT->ProcessLine(command.c_str());
  command = ".x Run_CorrectRef.C+(\"ChainedScalingResults.root\"," +
    std::string("\"SampleNamesNSRL.txt\")");
  gROOT->ProcessLine(command.c_str());
  command = ".x Run_VialDmgScaling.C+(\"CollectedNSRLAll2.root\"," +
    std::string("\"SampleNamesVials.txt\")");
  gROOT->ProcessLine(command.c_str());
  command = ".x Run_CorrectVial.C+(\"ChainedFinalResults.root\"," +
    std::string("\"ChainedVialDmgScalingResults.root\",") +
    std::string("\"SampleNamesNSRL.txt\")");
  gROOT->ProcessLine(command.c_str());

}
