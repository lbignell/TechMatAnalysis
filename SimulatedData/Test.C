//This is a macro.
//void Test()
{
  vector<double> Data(100,1);
  Data.at(50) = 5;
  
  vector<double> xvals;
  for(int i = 0; i<100; i++){xvals.push_back(i);}
  
  TGraph* Orig = new TGraph(100, &xvals[0], &Data[0]);
  Orig->Draw("");
  Orig->SetLineWidth(2);
  
  vector<double> Convolution;
  
  gROOT->ProcessLine(".x ConvGaus.C+(xvals, Data, 1, Convolution);");
 
  TGraph* Conv1 = new TGraph(100, &xvals[0], &Convolution[0]);
  Conv1->Draw("SAME");
  Conv1->SetLineColor(2);
}
