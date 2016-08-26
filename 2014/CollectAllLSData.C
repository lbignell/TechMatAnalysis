//This macro will collect up all of the fluorescence spectrum data in the
//directory and write them to a single root file.
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
#include <vector>

using namespace std;

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

int GetTimeInfo(string fname, vector<string>& TimeInfo){
  std::ifstream ifs;
  ifs.open (fname.c_str(), std::ifstream::in);
  string theLine;
  getline(ifs, theLine);
  ifs.close();
  if((theLine.at(0) == 'I')||(theLine.length()==0)){
    //This is the first character of "Instrument Type" (in the summary files,
    // but not in the data files). Skip this one and make a note.
    return -1;
  }
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
  return 0;
}

//Define months
string months[] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
		   "Sep", "Oct", "Nov", "Dec"};

int GetMonthNum(string MonthName){
  for(int MonthNum = 1; MonthNum<13; MonthNum++){
    if(months[MonthNum-1] ==  MonthName){
      return MonthNum;
    }
  }
  //Error past this point, return 0
  cout << "Error: Month name '" << MonthName << "' is incorrect format" << endl;
  return 0;

}

int GrabFile(string path, string fname, vector<double>& ch,
	      vector<double>& cnt, int& pos, fstream& dump, 
	      time_t& theTimestamp, int& theSampleNum, float& theCntTime,
	      float& theHnum, int& theTotCnt){

  char theLine[256];
  FILE *thefile = fopen((path + fname).c_str(), "r");

  if(thefile == NULL){
    cout << "The file: " << fname << " does not exist or could not be opened"
	 << endl;
    fclose(thefile);
    return -2;
  }
  //else{ cout << "Reading file: " << fname << endl;} 

  vector<string> Tinfo;
  int rval =  GetTimeInfo(path + fname, Tinfo);
  if(rval==-1){
    cout << "Skipping file: " << fname << ". This is a summary file." << endl;
    fclose(thefile);
    return -1;
  }

  //Tinfo now contains the time info from the first line; get that info...
  if(Tinfo.size()!=9){
    cout << "Error! something went wrong while reading the file!" << endl;
    for(size_t i =0; i<Tinfo.size(); i++){
      cout << "Tinfo.at(" << i << ") = " << Tinfo.at(i) << endl;
    }
    fclose(thefile);
    return -3;
  }
  int year = ConvertToInt(Tinfo.at(5));
  string month = Tinfo.at(4);
  int MonthNumber = GetMonthNum(month);
  int day = ConvertToInt(Tinfo.at(3));
  int hr = ConvertToInt(Tinfo.at(6));
  int min = ConvertToInt(Tinfo.at(7));
  int sec = ConvertToInt(Tinfo.at(8));

  struct tm StructTime;
  StructTime.tm_year = year - 1900;    //Defined as # years since 1900
  StructTime.tm_mon = MonthNumber - 1; //Defined as # months since Jan
  StructTime.tm_mday = day;
  StructTime.tm_hour = hr;
  StructTime.tm_min = min;
  StructTime.tm_sec = sec;

  theTimestamp = mktime(&StructTime);

  //////////////////////////////////////////////////////////////////////////////
  //Read using C. This code has issues sometimes getting the correct time.
  //Get the next line
  //fgets(theLine, sizeof(theLine), thefile);
  //if((theLine[0] == 'I')||(strlen(theLine)==0)){
    //This is the first character of "Instrument Type" (in the summary files,
    // but not in the data files). Skip this one and make a note.
    //cout << "Skipping file: " << fname << ". This is a summary file." << endl;
    //fclose(thefile);
    //return -1;
  //}
  //struct tm StructTime;
  //int year = 0;
  //char month[256];
  //int day = 0;
  //int hr = 0;
  //int min = 0;
  //int sec = 0;
  //sscanf(theLine, "%*s %*s %*s %i %s %i %i:%i:%i", &day, month, &year,
  //	 &hr, &min, &sec);
  //cout << "year = " << year << endl
  //   << "month = " << month << endl
  //   << "day = " << day << endl
  //   << "hr = " << hr << endl
  //   << "min = " << min << endl
  //   << "sec = " << sec << endl;
  //string strMonth = month;
  //int MonthNumber = GetMonthNum(strMonth);
  //if(MonthNumber == 0){
  //cout << "FileName = " << fname << endl;
  //}
  //StructTime.tm_year = year - 1900;    //Defined as # years since 1900
  //StructTime.tm_mon = MonthNumber - 1; //Defined as # months since Jan
  //StructTime.tm_mday = day;
  //StructTime.tm_hour = hr;
  //StructTime.tm_min = min;
  //StructTime.tm_sec = sec;

  //theTimestamp = mktime(&StructTime);
  //cout << "StructTime.tm_year = " << StructTime.tm_year << endl
  //   << "StructTime.tm_mon = " << StructTime.tm_mon << endl
  //   << "StructTime.tm_mday = " << StructTime.tm_mday << endl
  //   << "StructTime.tm_hour = " << StructTime.tm_hour << endl
  //   << "StructTime.tm_min = " << StructTime.tm_min << endl
  //   << "StructTime.tm_sec = " << StructTime.tm_sec << endl
  //   << "Timestamp = " << theTimestamp << endl;

  //skip next 8 lines
  //for(int i = 0; i<8; i++){fgets(theLine, sizeof(theLine), thefile);}
  //If using c++ method to get time, skip next 9 lines...
  for(int i = 0; i<9; i++){fgets(theLine, sizeof(theLine), thefile);}


  //The next line has the sample, rack-pos, time
  fgets(theLine, sizeof(theLine), thefile);
  //Format the line into values.
  sscanf(theLine, "\"Sample, Rack-Pos, Time:\" %i %*c%*c%*c %i %f",
	 &theSampleNum, &pos, &theCntTime);

  //cout << "SampleNum = " << theSampleNum << endl
  //   << "Position = " << pos << endl
  //   << "Count Time = " << theCntTime << endl;

  //The Next line has teh total counts:
  fgets(theLine, sizeof(theLine), thefile);
  //Format the line into values
  sscanf(theLine, "\"H#, Total Counts:\" %f %i", &theHnum, &theTotCnt);
  //  cout << "Horrock's Number = " << theHnum << endl
  //   << "Total Cnts = " << theTotCnt << endl;

  //skip next 4 lines
  for(int i = 0; i<4; i++){fgets(theLine, sizeof(theLine), thefile);}

  int theChan = 0;
  int theCnt = 0;

  //loop over remaining lines.
  while(fgets(theLine, sizeof(theLine), thefile) != NULL){
    sscanf(theLine, "%i %i", &theChan, &theCnt);
    ch.push_back(theChan);
    cnt.push_back(theCnt);
    dump << ch.back() << "," << pos << "," << cnt.back() << endl;
  }

  fclose(thefile);
  return 0;

}


void CollectAllLSData(string path){

  string theFile; 
  tinydir_dir dir;
  tinydir_file file;
  fstream filedump;
  string MapFileName = path;
  string suffix = "_Info.csv";
  if(path.find("_csv")!=string::npos){
    //The new style folder name
    MapFileName.replace(MapFileName.find("_csv/"), 5, suffix);
  }    
  else{
    //the old-style folder name
    MapFileName.replace(MapFileName.find("/"), 1, suffix);
  }
  cout << "MapFileName = " << MapFileName << endl;
  ifstream theMapFile((path + MapFileName).c_str());
  //FILE *theMapFile = fopen((path + MapFileName).c_str(), "r");
  //Read the MapFile data into a string vector.
  //char theLine[256];
  //Skip the first line.
  //fgets(theLine, sizeof(theLine), theMapFile);
  //Get the remaining lines
  string theLine;
  //skip the first line
  getline(theMapFile, theLine);
  vector<string> vSampleNames;
  string substr;
  //Get the remaining lines
  while(getline(theMapFile, theLine)){
    std::stringstream ss;
    ss << theLine;
    //skip the header
    getline(ss, substr, ',');
    //while(ss.good()){
    //getline(ss, substr, ',');
    while(getline(ss, substr, ',')){  
    if((substr!="")&&(substr!="\r")&&(substr.size()>2)){
	vSampleNames.push_back(substr);
      }
    }
  }

  for(int j = 0; j<vSampleNames.size(); j++){
    cout << "vSampleNames.at(" << j << ") = " << vSampleNames.at(j) << endl;
  }

  filedump.open("DataDump.csv",
		ios_base::in | ios_base::out | ios_base::trunc);
  cout << "Processing files in " << path << "... " << endl;

  TFile* Collected = TFile::Open("CollectedLSSpecData.root", "RECREATE");
  TTree* theTree = new TTree("LSSpec", "Collected LS Spectra");
  vector<double> channel;
  vector<double> counts;
  time_t timestamp;
  int position;
  int sampleNum;
  float cntTime;
  float HorrockNum;
  int TotCnts;
  int MeasNum;
  string SampleName;
  int IDnum = 0;
  float Dose = 0;
  //int RptNum = 0;
  bool isRef = false;

  theTree->Branch("Channel", "vector<double>", &channel);
  theTree->Branch("Counts", "vector<double>", &counts);
  theTree->Branch("Position", &position, "Position/I");
  theTree->Branch("TimeStamp" , &timestamp, "TimeStamp/L");
  theTree->Branch("SampleNum", &sampleNum, "SampleNum/I");
  theTree->Branch("CntTime", &cntTime, "CntTime/F");
  theTree->Branch("Hnum", &HorrockNum, "Hnum/F");
  theTree->Branch("TotCnts", &TotCnts, "TotCnts/I");
  theTree->Branch("MeasNum", &MeasNum, "MeasNum/I");
  theTree->Branch("SampleName", "string", &SampleName);
  theTree->Branch("IDnum", &IDnum, "IDnum/I");
  theTree->Branch("Dose", &Dose, "Dose/F");
  theTree->Branch("isRef", &isRef, "isRef/O");
  //theTree->Branch("RptNum", &RptNum, "RptNum/I");

  //make a loop over the file names:
  tinydir_open(&dir, path.c_str());
  //Make a vector to contain the sorted list
  vector< string > dirList;

  while (dir.has_next)
    {
      tinydir_readfile(&dir, &file);
      if(!file.is_dir){
	//Add to the vector
	dirList.push_back(string(file.name));
      }
      tinydir_next(&dir);
    }

  //This makes sure that the list is in order
  std::sort(dirList.begin(), dirList.end());

  int NumSamples = vSampleNames.size();

  int SampleNum = 0;
  bool justFailed;

  cout << "Number of files in directory = " << dirList.size() << endl;

  //for(std::vector<string>::iterator it = dirList.begin();
  //  it != dirList.end(); it++){
  for(int k = 0; k<dirList.size(); k++){
    //if(strstr(*it->c_str(), ".csv")
    // &&(!(strstr(*it->c_str(), "DataDump")))
    // &&(!(strstr(*it->c_str(), "Info")))){
    if(dirList.at(k).find("MeasurementNo")!=string::npos){  
      //printf("Processing in file %s \n", dirList.at(k).c_str());
      int returnval = GrabFile(path, dirList.at(k), channel, counts, position,
			       filedump, timestamp, sampleNum, cntTime,
			       HorrockNum, TotCnts);
	if(returnval == 0){
	  if(justFailed==true){
	    //investigate timestamp
	    char buff[20];
	    strftime(buff, 20, "%Y-%m-%d %H:%M:%S", localtime(&timestamp));
	    cout << "Time after summary set: " << buff << endl;
	  }
	  justFailed = false;
	  sscanf(dirList.at(k).c_str(), "MeasurementNo%d.csv", &MeasNum);
	  //Get sample name, composition code, and dose
	  SampleName = vSampleNames.at(SampleNum);
	  stringstream ss;
	  ss << SampleName;
	  string dummy;
	  getline(ss, dummy, '_');
	  if(istringstream(dummy)>>IDnum){
	    // cout << "IDnum = " << IDnum << endl;
	  }//All good
	  else if((dummy == "EMPTY")||(dummy == " EMPTY")||(dummy == "EMPTY ")
		  ||(dummy == "EMPTY\r")){
	    IDnum = 0;
	  }
	  else if((dummy == "GLASS")||(dummy == " GLASS")||(dummy == "GLASS ")
		  ||(dummy == "GLASS\r")){
	    IDnum = 1;
	  }
	  else if((dummy == "PP")||(dummy == " PP")||(dummy == "PP ")
		  ||(dummy == "PP\r")){
	    IDnum = 2;
	  }
	  else if((dummy == "HDPE")||(dummy == " HDPE")||(dummy == "HDPE ")
		  ||(dummy == "HDPE\r")){
	    IDnum = 3;
	  }
	  else if((dummy == "STD")||(dummy == " STD")||(dummy == "STD ")
		  ||(dummy == "STD\r")){
	    IDnum = 10;
	  }
	  else{
	    cout << "IDnum couldn't be set, string was " << dummy << endl;
	    IDnum = -1;
	  }
	  //Later on I can come up with a way of getting the dose
	  //get()
	  //Dose
	  //increment/reset the SampleNum
	  if(SampleNum == NumSamples - 1){
	    cout << "Resetting SampleNum, SampleNum = " << SampleNum << endl;
	    cout << "Last vial = " << SampleName << endl;
	    SampleNum = 0;
	  }
	  else{SampleNum++;}
	  
	  //Determine if the sample was a reference:
	  if(SampleName.find("REF")!=string::npos){
	    isRef = true;\
	    //cout << "Found REF file; name = " << SampleName << endl;
	  }
	  else{isRef = false;}

	  theTree->Fill();
	  channel.clear();
	  counts.clear();
	}
	else{
	  sscanf(dirList.at(k).c_str(), "MeasurementNo%d.csv", &MeasNum);
	  cout << "Error Returned; MeasNum = " << MeasNum << endl;
	  justFailed = true;
	}
    }
  }

  theTree->Write();
  Collected->Close();

}

