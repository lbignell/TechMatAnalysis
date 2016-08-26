#include <iostream>    
#include <fstream>     
#include <sstream>
#include <string>
#include <vector>

using namespace std;

void GetTimeInfo(string fname, vector<string>& TimeInfo){
  std::ifstream ifs;
  ifs.open (fname.c_str(), std::ifstream::in);
  string theLine;
  getline(ifs, theLine);
  stringstream ss(theLine);
  string val;
  while(ss >> val){
    if(val.find(":")==string::npos){
      TimeInfo.push_back(val);
    }
    else{
      stringstream ss2(val);
      string tok;
      while(getline(ss2, tok, ':')){
	TimeInfo.push_back(tok);
      }
    }
  }
}

int ConvertToInt(string str){

  stringstream ss(str);
  int val;
  ss >> val;
  if(ss.fail()){
    cout << "Couldn't convert " << str << " to an integer!" << endl;
    return 0;
  }
  return val;
}

void TestRead() {
  cout << "Getting Time info...";
  string theFile = "Pre-irrad_01_141014/MeasurementNo00175.csv";
  vector <string> Tinfo;
  GetTimeInfo(theFile, Tinfo);
  cout << "done!" << endl;
  for(size_t i = 0; i<Tinfo.size(); i++){
    cout << "[" << i << "] = " << Tinfo.at(i) << endl;
  }

  cout << "testing ConvertToInt..." << endl;

  int year = ConvertToInt(Tinfo.at(5));
  int day = ConvertToInt(Tinfo.at(3));
  year = year + 5;
  cout << "in 5 years it will be: " << year << endl;

}
