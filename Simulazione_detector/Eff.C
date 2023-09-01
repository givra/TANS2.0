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
	float risol[size] = {0,0,0,0,0,0,0,0,0,0};			// array contenente le risoluzioni
	float errrisol[size] = {0,0,0,0,0,0,0,0,0,0};
	float differenza = 0;
	
	// histo per efficienza
	static TH1D* histM3 = new TH1D("hM3","2.5 < molteplicità < 3.5", 15, -0.05, 0.05);// fatto per molteplicità 3
	static TH1D* histM5 = new TH1D("hM5","4.5 < molteplicità < 5.5", 25, -0.05, 0.05);
	static TH1D* histM6 = new TH1D("hM6","5.5 < molteplicità < 6.5", 25, -0.05, 0.05);
	static TH1D* histM7 = new TH1D("hM7","6.5 < molteplicità < 7.5", 25, -0.05, 0.05);
	static TH1D* histM8 = new TH1D("hM8","7.5 < molteplicità < 8.5", 25, -0.05, 0.05);
	static TH1D* histM12 = new TH1D("hM12","11.5 < molteplicità < 12.5", 50, -0.05, 0.05);
	static TH1D* histM22 = new TH1D("hM22","21.5 < molteplicità < 22.5", 50, -0.04, 0.04);
	static TH1D* histM32 = new TH1D("hM32","31.5 < molteplicità < 32.5", 50, -0.04, 0.04);
	static TH1D* histM42 = new TH1D("hM42","41.5 < molteplicità < 42.5", 50, -0.04, 0.04);
	static TH1D* histM52 = new TH1D("hM52","51.5 < molteplicità < 52.5", 50, -0.04, 0.04);
	
	
	  //definizione struct
        typedef struct {
           float X,Y,Z;
           int mult;} VTX1;
        static VTX1 point;
        
        typedef struct {
           float Z;} VTX;
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
    branch4->SetAddress(&pointRec.Z);
	
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
		
		// errore binomiale
		erreff[0] = pow((eff[0]*(1-eff[0])/denEff.x0),0.5);
		erreff[1] = pow((eff[1]*(1-eff[1])/denEff.x1),0.5);
		erreff[2] = pow((eff[2]*(1-eff[2])/denEff.x2),0.5);
		erreff[3] = pow((eff[3]*(1-eff[3])/denEff.x3),0.5);
		erreff[4] = pow((eff[4]*(1-eff[4])/denEff.x4),0.5);
		erreff[5] = pow((eff[5]*(1-eff[5])/denEff.x5),0.5);
		erreff[6] = pow((eff[6]*(1-eff[6])/denEff.x6),0.5);
		erreff[7] = pow((eff[7]*(1-eff[7])/denEff.x7),0.5);
		erreff[8] = pow((eff[8]*(1-eff[8])/denEff.x8),0.5);
		erreff[9] = pow((eff[9]*(1-eff[9])/denEff.x9),0.5);
		
		// binomiale anche per molteplicità????
		
		
	// fit gaussiani
		//TCanvas *cM3=new TCanvas("cM3","cM3",800,600);
		//cM3->SetFillColor(0);
		//cM3->cd();
		TF1 *fM3 = new TF1("fM3","gaus",-0.05,0.05);
		fM3->SetLineColor(kRed);
		//fM3->SetParameters(1,60);	// media
		//fM3->SetParameters(2,7);	// sigma
		histM3->Fit(fM3,"NR+");
		gStyle->SetOptFit(0);
		double mediaM3 = abs(fM3->GetParameter(1));
		double errmediaM3 = abs(fM3->GetParError(1));
		double sigmaM3 = fM3->GetParameter(2);
		double errsigmaM3 = fM3->GetParError(2);
		risol[0] = sigmaM3;
		errrisol[0] = errsigmaM3;
		fM3->Draw("same");
		histM3->DrawCopy("same");
		
		//TCanvas *cM5=new TCanvas("cM5","cM5",800,600);
		TF1 *fM5 = new TF1("fM5","gaus",-0.04,0.04);
		fM5->SetLineColor(kRed);
		//fM5->SetParameters(1,60);	// media
		//fM5->SetParameters(2,7);	// sigma
		histM5->Fit(fM5,"NR+");
		gStyle->SetOptFit(0);
		double mediaM5 = abs(fM5->GetParameter(1));
		double errmediaM5 = abs(fM5->GetParError(1));
		double sigmaM5 = fM5->GetParameter(2);
		double errsigmaM5 = fM5->GetParError(2);
		risol[1] = sigmaM5;
		errrisol[1] = errsigmaM5;
		//fM5->Draw("same");
		//histM5->Draw("same");
		
		//TCanvas *M6=new TCanvas("M6","M6",800,600);
		TF1 *fM6 = new TF1("fM6","gaus",-0.04,0.04);
		fM6->SetLineColor(kRed);
		//fM6->SetParameters(1,60);	// media
		//fM6->SetParameters(2,7);	// sigma
		histM6->Fit(fM6,"NR+");
		gStyle->SetOptFit(0);
		double mediaM6 = abs(fM6->GetParameter(1));
		double errmediaM6 = abs(fM6->GetParError(1));
		double sigmaM6 = fM6->GetParameter(2);
		double errsigmaM6 = fM6->GetParError(2);
		risol[2] = sigmaM6;
		errrisol[2] = errsigmaM6;
		//fM6->Draw("same");
		//histM6->Draw("same");
		
		//TCanvas *M7=new TCanvas("M7","M7",800,600);
		TF1 *fM7 = new TF1("fM7","gaus",-0.04,0.04);
		fM7->SetLineColor(kRed);
		//fM7->SetParameters(1,60);	// media
		//fM7->SetParameters(2,7);	// sigma
		histM7->Fit(fM7,"NR+");
		gStyle->SetOptFit(0);
		double mediaM7 = abs(fM7->GetParameter(1));
		double errmediaM7 = abs(fM7->GetParError(1));
		double sigmaM7 = fM7->GetParameter(2);
		double errsigmaM7 = fM7->GetParError(2);
		risol[3] = sigmaM7;
		errrisol[3] = errsigmaM7;
		//fM7->Draw("same");
		//histM7->Draw("same");
		
		
		//TCanvas *M8=new TCanvas("M8","M8",800,600);
		TF1 *fM8 = new TF1("fM8","gaus",-0.04,0.04);
		fM8->SetLineColor(kRed);
		//fM8->SetParameters(1,60);	// media
		//fM8->SetParameters(2,7);	// sigma
		histM8->Fit(fM8,"NR+");
		gStyle->SetOptFit(0);
		double mediaM8 = abs(fM8->GetParameter(1));
		double errmediaM8 = abs(fM8->GetParError(1));
		double sigmaM8 = fM8->GetParameter(2);
		double errsigmaM8 = fM8->GetParError(2);
		risol[4] = sigmaM8;
		errrisol[4] = errsigmaM8;
		//fM8->Draw("same");
		//histM8->Draw("same");
		
		//TCanvas *M12=new TCanvas("M12","M12",800,600);
		TF1 *fM12 = new TF1("fM12","gaus",-0.03,0.03);
		fM12->SetLineColor(kRed);
		//fM12->SetParameters(1,60);	// media
		//fM12->SetParameters(2,7);	// sigma
		histM12->Fit(fM12,"NR+");
		gStyle->SetOptFit(0);
		double mediaM12 = abs(fM12->GetParameter(1));
		double errmediaM12 = abs(fM12->GetParError(1));
		double sigmaM12 = fM12->GetParameter(2);
		double errsigmaM12 = fM12->GetParError(2);
		risol[5] = sigmaM12;
		errrisol[5] = errsigmaM12;
		//fM12->Draw("same");
		//histM12->Draw("same");
		
		//TCanvas *M22=new TCanvas("M22","M22",800,600);
		TF1 *fM22 = new TF1("fM22","gaus",-0.025,0.025);
		fM22->SetLineColor(kRed);
		//fM22->SetParameters(1,60);	// media
		//fM22->SetParameters(2,7);	// sigma
		histM22->Fit(fM22,"NR+");
		gStyle->SetOptFit(0);
		double mediaM22 = abs(fM22->GetParameter(1));
		double errmediaM22 = abs(fM22->GetParError(1));
		double sigmaM22 = fM22->GetParameter(2);
		double errsigmaM22 = fM22->GetParError(2);
		risol[6] = sigmaM22;
		errrisol[6] = errsigmaM22;
		//fM22->Draw("same");
		//histM22->Draw("same");
		
		//TCanvas *M32=new TCanvas("M32","M32",800,600);
		TF1 *fM32 = new TF1("fM32","gaus",-0.02,0.02);
		fM32->SetLineColor(kRed);
		//fM32->SetParameters(1,60);	// media
		//fM32->SetParameters(2,7);	// sigma
		histM32->Fit(fM32,"NR+");
		gStyle->SetOptFit(0);
		double mediaM32 = abs(fM32->GetParameter(1));
		double errmediaM32 = abs(fM32->GetParError(1));
		double sigmaM32 = fM32->GetParameter(2);
		double errsigmaM32 = fM32->GetParError(2);
		risol[7] = sigmaM32;
		errrisol[7] = errsigmaM32;
		//fM32->Draw("same");
		//histM32->Draw("same");
		
		//TCanvas *M42=new TCanvas("M42","M42",800,600);
		TF1 *fM42 = new TF1("fM42","gaus",-0.015,0.015);
		fM42->SetLineColor(kRed);
		//fM42->SetParameters(1,60);	// media
		//fM42->SetParameters(2,7);	// sigma
		histM42->Fit(fM42,"NR+");
		gStyle->SetOptFit(0);
		double mediaM42 = abs(fM42->GetParameter(1));
		double errmediaM42 = abs(fM42->GetParError(1));
		double sigmaM42 = fM42->GetParameter(2);
		double errsigmaM42 = fM42->GetParError(2);
		risol[8] = sigmaM42;
		errrisol[8] = errsigmaM42;
		//fM42->Draw("same");
		//histM42->Draw("same");
		
		//TCanvas *M52=new TCanvas("M52","M52",800,600);
		TF1 *fM52 = new TF1("fM52","gaus",-0.015,0.015);
		fM52->SetLineColor(kRed);
		//fM52->SetParameters(1,60);	// media
		//fM52->SetParameters(2,7);	// sigma
		histM52->Fit(fM52,"NR+");
		gStyle->SetOptFit(0);
		double mediaM52 = abs(fM52->GetParameter(1));
		double errmediaM52 = abs(fM52->GetParError(1));
		double sigmaM52 = fM52->GetParameter(2);
		double errsigmaM52 = fM52->GetParError(2);
		risol[9] = sigmaM52;
		errrisol[9] = errsigmaM52;
		//fM52->Draw("same");
		//histM52->Draw("same");
		
		for(int l=0; l<=9; l++) cout<<risol[l]<<endl;
		
// __________________________________________ efficienza _________________________________________
		//TCanvas *c1=new TCanvas("c1","c1",800,600);
		TGraphErrors *graphE= new TGraphErrors(size,moltep,eff,errmoltep,erreff);
		graphE->SetMarkerSize(1);//https://root.cern.ch/doc/master/classTAttMarker.html
		graphE->SetMarkerStyle(33);
		graphE->SetTitle("Efficiency vs Molteplicity");
		graphE->GetXaxis()->SetTitle("Molteplicity[]");
		graphE->GetYaxis()->SetTitle("Efficiency[]");
		//c1->cd();
		graphE->Draw("ap");

	
// ____________________________________________ risoluzione ________________________________________		
		TGraphErrors *graphR= new TGraphErrors(size,moltep,risol,errmoltep,errrisol);
		graphR->SetMarkerSize(1);//https://root.cern.ch/doc/master/classTAttMarker.html
		graphR->SetMarkerStyle(33);
		graphR->SetTitle("Std Deviation vs Molteplicity");
		graphR->GetXaxis()->SetTitle("Molteplicity[]");
		graphR->GetYaxis()->SetTitle("Std Deviation[]");
		//c1->cd();
		graphR->Draw("ap");
		
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
		graphR->Write();
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
       // delete graphE;
        //graphE = nullptr;
	//	delete graphR;
       // graphR = nullptr;
		
		double TT = time.CpuTime();	
         cout<<"Il tempo impiegato dalla CPU è "<<TT<<" s"<<endl;  
		
}
