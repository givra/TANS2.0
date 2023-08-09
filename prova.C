#include "TFile.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "Riostream.h"
TH1F* maniphist();

void prova(){
  TH1F* histog = maniphist();
  new TCanvas();
  histog->SetLineColor(kBlue);
  histog->SetMinimum(0);
  //histog->Draw();
  
  TAxis *xa=histog->GetXaxis();
  float step = xa->GetBinWidth(1);
  int b1=xa->FindBin(-2.);
  int b2=xa->FindBin(2.);
  int nobins=b2-b1+1;
 
  float Xlow=xa->GetBinLowEdge(b1); //-2.04
  float Xhig=xa->GetBinUpEdge(b2); //2.04
 
  
  float Ymax, Xmax;
  
  for(int r = 0; r < nobins ;r++){
        	       
        	       if(Ymax < histog->GetBinContent(r)){
        	                  Ymax = histog->GetBinContent(r);
        	                  
        	                  Xmax = r;
        	                  }
        	        }
  
  float u1, u2;

  float fmax = Ymax;

  float Xt;
  float Yt;
  float f;
  
  int conto = 0;
  
  cout<<"Xmax "<<Xmax<<"nbins "<<nobins<<"Xlow "<<Xlow<<"Xhig "<<Xhig<<endl;
  TH1D* hist = new TH1D ("h","distribuzione eta" , nobins, Xlow, Xhig);
  hist->SetMinimum(0);
  
  int nentries = 50000;
  
  for(int entry = 0; entry<nentries; entry++){
  do{
  
  conto++;
  u1 = gRandom->Rndm();
  u2 = gRandom->Rndm();
  
  Yt = fmax*u2;
  Xt = (nobins*u1)/1;
  f = histog->GetBinContent(Xt);}
  while(f<=Yt);
  
  //cout<<conto<<endl;
  //*N+=conto;
  
  hist->Fill(((Xhig-Xlow)*(Xt)/(nobins))+Xlow-((Xhig-Xlow)/(2*nobins)));//-((Xhig-Xlow)/2*nobins));
  //cout<<((Xhig-Xlow)*Xt/nobins)+Xlow<<endl;
  }
  hist->Draw();
}
  TH1F* maniphist(){
  TFile F("kinem.root");
  TH1F *disteta = (TH1F*)F.Get("heta2");
  disteta->SetDirectory(0);
  disteta->SetMinimum(0);
  F.Close();
  TAxis *xa=disteta->GetXaxis();
  float step = xa->GetBinWidth(1);
  int b1=xa->FindBin(-2.);
  int b2=xa->FindBin(2.);
 
  float xlow=xa->GetBinLowEdge(b1);
  float xhig=xa->GetBinUpEdge(b2);
  int nobins=b2-b1+1;
  float step2 = (xhig-xlow)/nobins;
  cout << "Check: "<<step<<"; "<<step2<<endl;
  TH1F* heta2 = new TH1F("heta2","#eta distribution 2",nobins,xlow,xhig);
  int j=1;
  for(int i=b1;i<=b2;i++)heta2->SetBinContent(j++,disteta->GetBinContent(i));
  //  heta2->Draw();
  //new TCanvas();
  //disteta->Draw();
  heta2->SetLineColor(2);
 // heta2->Draw("same");
  return heta2;
}

/*float reiezione (int* N){

  float u1;
  float u2;
  float fmax = pow(a, -1);

  float Xt;
  float Yt;
  float f;
  
  int conto = 0;

  do{
  
  conto++;
  u1 = gRandom->Rndm();
  u2 = gRandom->Rndm();
  
  Yt = fmax*u2;
  Xt = Xmin + (Xmax - Xmin)*u1;
  f = pow((sin(Xt)*sin(Xt)) + a*cos(Xt)*cos(Xt),-1);}
  while(f<=Yt);
  
  //cout<<conto<<endl;
  *N+=conto;
  
  return Xt;
  
}
*/

