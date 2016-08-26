#include <vector>
#include "TF1.h"

using namespace std;

void ConvGaus(vector<double>& SimEnergy,
	      vector<double>& SimData, double GausWidth,
	      vector<double>& Result){

  vector<double> ConvResult(SimEnergy.size(), 0);

  if(SimEnergy.size()!=SimData.size()){
    cout << "Energy and Spectrum vectors aren't the same size!" << endl;
    return;
  }

  //TF1* fgauss = new TF1("fgauss", "TMath::Gaus(x,[0],[1],1)");
  //fgauss->SetParameter(1,GausWidth);

  //double ConvWidth = 10;

  for(size_t n = 0; n<SimData.size(); n++){
    //Do the convolution
    //fgauss->SetParameter(0, SimEnergy.at(n));
    //cout << "n  = " << n << endl;
    for(int m = -SimData.size(); m<int(SimData.size());
	m++){
      //cout << "m = " << m << endl;
      //Conv(n) = SimData(m)*gauss(n-m)
      double GaussVal = (1/(sqrt(2*3.1415)*GausWidth))*
	exp(-0.5*pow(((m)/GausWidth),2));
      if(((n-m)>0)&&((n-m)<SimData.size())){
	ConvResult.at(n) += SimData.at(n-m)*GaussVal;
	//if(isnan(ConvResult.at(n))){
	//cout<< "n = " << n << ", m = " << m
	//    << ", SimData.at(n-m) = " << SimData.at(n-m) 
	//    << ", GaussVal = " << GaussVal << endl;
	//}
      }
    }
  }

  //Normalise?
  double DataSum = 0;
  double ConvSum = 0;
  for(size_t i = 0; i<SimData.size(); i++){
    if(ConvResult.at(i)>0){
      DataSum+=SimData.at(i);
      ConvSum+=ConvResult.at(i);
    }
  }

  //cout << "DataSum = " << DataSum << ", ConvSum = " << ConvSum << endl;

  for(size_t j = 0; j<SimData.size(); j++)
    {ConvResult.at(j) = (DataSum/ConvSum)*ConvResult.at(j);}

  Result.swap(ConvResult);

  return;
}
