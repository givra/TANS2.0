#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TAxis.h"

#include "Punto2.h"
#include "Vertex.h"
#include "TracciaMC.h"
#include "Multis.h"
#include "Smearing.h"

using namespace std;

void Eff(){
	
	TStopwatch time;	
	 
    int numeroeventi;
	  
	const int size = 10; 			// dimensione array per TGraph  
	float moltep[size] = {0,0,0,0,0,0,0,0,0,0}; 		// array contenente le moltep di ogni evento da passare al TGraph
	float errmoltep[size] = {0,0,0,0,0,0,0,0,0,0};	// binomiale??
	float eff[size] = {0,0,0,0,0,0,0,0,0,0};			// array contenente le efficienze
	float erreff[size] = {0,0,0,0,0,0,0,0,0,0};
	float differenza = 0;
	
	// histo per efficienza
	static TH1D* histM3 = new TH1D("hM3","2.5 < molteplicità < 3.5", 150, -0.15, 0.15);// fatto per molteplicità 3
	static TH1D* histM5 = new TH1D("hM5","4.5 < molteplicità < 5.5", 150, -0.15, 0.15);
	static TH1D* histM6 = new TH1D("hM6","5.5 < molteplicità < 6.5", 150, -0.15, 0.15);
	static TH1D* histM7 = new TH1D("hM7","6.5 < molteplicità < 7.5", 150, -0.15, 0.15);
	static TH1D* histM8 = new TH1D("hM8","7.5 < molteplicità < 8.5", 150, -0.15, 0.15);
	static TH1D* histM12 = new TH1D("hM12","11.5 < molteplicità < 12.5", 250, -0.15, 0.15);
	static TH1D* histM22 = new TH1D("hM22","21.5 < molteplicità < 22.5", 250, -1.5, 1.5);
	static TH1D* histM32 = new TH1D("hM32","31.5 < molteplicità < 32.5", 250, -1.5, 1.5);
	static TH1D* histM42 = new TH1D("hM42","41.5 < molteplicità < 42.5", 250, -1.5, 1.5);
	static TH1D* histM52 = new TH1D("hM52","51.5 < molteplicità < 52.5", 250, -1.5, 1.5);
	
	
	typedef struct {
           float X,Y,Z;
           int mult;} VTX;
    static VTX point;
	static VTX pointRec;
	
	typedef struct{
			float x0,x1,x2,x3,x4,x5,x6,x7,x8,x9;		// Z simulate con un certo valore di moltep, sarà denominatore dell'efficienza
			} vett;
		static vett denEff;
		
	TFile hfile2("htree2.root");
	TTree *tree2 = (TTree*)hfile2.Get("T2");
    TBranch *b1=tree2->GetBranch("VertMult");
	b1->SetAddress(&point.X);
		
	TFile hfile3("htree3.root");
	TTree *tree3 = (TTree*)hfile3.Get("T3");
	TBranch *branch3 = tree3->GetBranch("denEff");
    branch3->SetAddress(&denEff.x0);
	
	TFile hfile4("htree4.root");
	TTree *tree4 = (TTree*)hfile4.Get("T4");
	TBranch *branch4 = tree4->GetBranch("Zrec");
    branch4->SetAddress(&pointRec.X);
	
	numeroeventi = tree2->GetEntries(); //acquisico informazione sul numero di eventi nel mio detector
	
	float Zrec[numeroeventi];
    float Zsim[numeroeventi];
	
       
	// loop sugli ingressi nel TTree
    for(int ev=0;ev<numeroeventi;ev++){
                
		tree3->GetEvent(ev);
        tree2->GetEvent(ev);
		tree4->GetEvent(ev);
		
		Zsim[ev] = point.Z; 
		Zrec[ev] = pointRec.Z;
		
		if(pointRec.Z != 100.){	
			differenza = pointRec.Z	- Zsim[ev];
			
			switch(point.mult)
				{case 3:
					histM3->Fill(differenza);
					moltep[0] = point.mult;
					//cout << " moltep " << point.mult << " diff " << differenza << endl;
					break;
				 case 5:
					histM5->Fill(differenza);
					moltep[1] = point.mult;
					//cout << " moltep " << point.mult << " diff " << differenza << endl;
					break;
				 case 6:
					histM6->Fill(differenza);
					moltep[2] = point.mult;
					//cout << " moltep " << point.mult << " diff " << differenza << endl;
					break;
				 case 7:
					histM7->Fill(differenza);
					moltep[3] = point.mult;
					//cout << " moltep " << point.mult << " diff " << differenza << endl;
					break;
				 case 8:
					histM8->Fill(differenza);
					moltep[4] = point.mult;
					//cout << " moltep " << point.mult << " diff " << differenza << endl;
					break;
				 case 12:
					histM12->Fill(differenza);
					moltep[5] = point.mult;
					//cout << " moltep " << point.mult << " diff " << differenza << endl;
					break;
				 case 22:
					histM22->Fill(differenza);
					
					moltep[6] = point.mult;
					break;
				 case 32:
					histM32->Fill(differenza);
					
					moltep[7] = point.mult;
					break;
				 case 42:
					histM42->Fill(differenza);
					
					moltep[8] = point.mult;
					break;
				 case 52:
					histM52->Fill(differenza);
					moltep[9] = point.mult;
					//cout << " moltep " << point.mult << " diff " << differenza << endl;
					break;
				}
		}
	} 		// fine ciclo eventi
	
	// # di Z ricostruite con una certa moltep su quelle simulate
		eff[0] = histM3->GetEntries()/denEff.x0;		
		eff[1] = histM5->GetEntries()/denEff.x1;
		eff[2] = histM6->GetEntries()/denEff.x2;
		eff[3] = histM7->GetEntries()/denEff.x3;
		eff[4] = histM8->GetEntries()/denEff.x4;
		eff[5] = histM12->GetEntries()/denEff.x5;
		eff[6] = histM22->GetEntries()/denEff.x6;
		eff[7] = histM32->GetEntries()/denEff.x7;
		eff[8] = histM42->GetEntries()/denEff.x8;
		eff[9] = histM52->GetEntries()/denEff.x9;
		
	// fit gaussiani
		//TCanvas *cM3=new TCanvas("cM3","cM3",800,600);
		TF1 *fM3 = new TF1("fM3","gaus",-0.06,0.06);
		fM3->SetLineColor(kRed);
		//fM3->SetParameters(1,60);	// media
		//fM3->SetParameters(2,7);	// sigma
		histM3->Fit(fM3,"NR+");
		gStyle->SetOptFit(0);
		//fM3->Draw("same");
		//histM3->Draw("same");
		
		//TCanvas *cM5=new TCanvas("cM5","cM5",800,600);
		TF1 *fM5 = new TF1("fM5","gaus",-0.06,0.06);
		fM5->SetLineColor(kRed);
		//fM5->SetParameters(1,60);	// media
		//fM5->SetParameters(2,7);	// sigma
		histM5->Fit(fM5,"NR+");
		gStyle->SetOptFit(0);
		//fM5->Draw("same");
		//histM5->Draw("same");
		
		//TCanvas *M6=new TCanvas("M6","M6",800,600);
		TF1 *fM6 = new TF1("fM6","gaus",-0.07,0.07);
		fM6->SetLineColor(kRed);
		//fM6->SetParameters(1,60);	// media
		//fM6->SetParameters(2,7);	// sigma
		histM6->Fit(fM6,"NR+");
		gStyle->SetOptFit(0);
		//fM6->Draw("same");
		//histM6->Draw("same");
		
		//TCanvas *M7=new TCanvas("M7","M7",800,600);
		TF1 *fM7 = new TF1("fM7","gaus",-0.07,0.07);
		fM7->SetLineColor(kRed);
		//fM7->SetParameters(1,60);	// media
		//fM7->SetParameters(2,7);	// sigma
		histM7->Fit(fM7,"NR+");
		gStyle->SetOptFit(0);
		//fM7->Draw("same");
		//histM7->Draw("same");
		
		//TCanvas *M8=new TCanvas("M8","M8",800,600);
		TF1 *fM8 = new TF1("fM8","gaus",-0.05,0.05);
		fM8->SetLineColor(kRed);
		//fM8->SetParameters(1,60);	// media
		//fM8->SetParameters(2,7);	// sigma
		histM8->Fit(fM8,"NR+");
		gStyle->SetOptFit(0);
		//fM8->Draw("same");
		//histM8->Draw("same");
		
		//TCanvas *M12=new TCanvas("M12","M12",800,600);
		TF1 *fM12 = new TF1("fM12","gaus",-0.05,0.05);
		fM12->SetLineColor(kRed);
		//fM12->SetParameters(1,60);	// media
		//fM12->SetParameters(2,7);	// sigma
		histM12->Fit(fM12,"NR+");
		gStyle->SetOptFit(0);
		//fM12->Draw("same");
		//histM12->Draw("same");
		
		//TCanvas *M22=new TCanvas("M22","M22",800,600);
		TF1 *fM22 = new TF1("fM22","gaus",-0.04,0.04);
		fM22->SetLineColor(kRed);
		//fM22->SetParameters(1,60);	// media
		//fM22->SetParameters(2,7);	// sigma
		histM22->Fit(fM22,"NR+");
		gStyle->SetOptFit(0);
		//fM22->Draw("same");
		//histM22->Draw("same");
		
		//TCanvas *M32=new TCanvas("M32","M32",800,600);
		TF1 *fM32 = new TF1("fM32","gaus",-0.04,0.04);
		fM32->SetLineColor(kRed);
		//fM32->SetParameters(1,60);	// media
		//fM32->SetParameters(2,7);	// sigma
		histM32->Fit(fM32,"NR+");
		gStyle->SetOptFit(0);
		//fM32->Draw("same");
		//histM32->Draw("same");
		
		//TCanvas *M42=new TCanvas("M42","M42",800,600);
		TF1 *fM42 = new TF1("fM42","gaus",-0.03,0.03);
		fM42->SetLineColor(kRed);
		//fM42->SetParameters(1,60);	// media
		//fM42->SetParameters(2,7);	// sigma
		histM42->Fit(fM42,"NR+");
		gStyle->SetOptFit(0);
		//fM42->Draw("same");
		//histM42->Draw("same");
		
		//TCanvas *M52=new TCanvas("M52","M52",800,600);
		TF1 *fM52 = new TF1("fM52","gaus",-0.03,0.03);
		fM52->SetLineColor(kRed);
		//fM52->SetParameters(1,60);	// media
		//fM52->SetParameters(2,7);	// sigma
		histM52->Fit(fM52,"NR+");
		gStyle->SetOptFit(0);
		//fM52->Draw("same");
		//histM52->Draw("same");
		
		//TCanvas *c1=new TCanvas("c1","c1",800,600);
		TGraphErrors *graphE= new TGraphErrors(size,moltep,eff,errmoltep,erreff);
		graphE->SetMarkerSize(0.9);//https://root.cern.ch/doc/master/classTAttMarker.html
		graphE->SetMarkerStyle(43);
		graphE->SetTitle("Efficiency vs molteplicity");
		graphE->GetXaxis()->SetTitle("molteplicity[]");
		graphE->GetYaxis()->SetTitle("efficiency[]");
		//c1->cd();
		graphE->Draw("ap");
		
		TFile file1("efficienza.root", "recreate");
		histM3->Write(); 
	    histM5->Write(); 
        histM6->Write(); 
        histM7->Write(); 
        histM8->Write(); 
        histM12->Write();
        histM22->Write();
        histM32->Write();
        histM42->Write();
        histM52->Write();
        graphE->Write();
        file1.Close();
        delete histM3;
        delete histM5;
        delete histM6;
        delete histM7;
        delete histM8;
        delete histM12;
        delete histM22;
        delete histM32;
        delete histM42;
        delete histM52;
        delete graphE;
        graphE = nullptr;
		
		double TT = time.CpuTime();	
         cout<<"Il tempo impiegato dalla CPU è "<<TT<<" s"<<endl;  
		
}