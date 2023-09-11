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

#include "Punto.h"
#include "Vertex.h"
#include "Traccia.h"
#include "Multis.h"
#include "Smearing.h"

using namespace std;

void Eff(){
	
	TStopwatch time;	
	 
    int numeroeventi;
	  
	const int size = 10; 			// dimensione array per TGraph  
	float moltep[size] = {0,0,0,0,0,0,0,0,0,0}; 		// array contenente le moltep di ogni evento da passare al TGraph
	float errmoltep[size] = {0,0,0,0,0,0,0,0,0,0};	// binomiale
	float eff[size] = {0,0,0,0,0,0,0,0,0,0};			// array contenente le efficienze in funz di moltep
	float erreff[size] = {0,0,0,0,0,0,0,0,0,0};
	float eff2[size] = {0,0,0,0,0,0,0,0,0,0};			// array contenente le efficienze in funz di ztrue
	float erreff2[size] = {0,0,0,0,0,0,0,0,0,0};
	float risol[size] = {0,0,0,0,0,0,0,0,0,0};			// array contenente le risoluzioni in funz di moltep
	float errrisol[size] = {0,0,0,0,0,0,0,0,0,0};
	float risol2[size] = {0,0,0,0,0,0,0,0,0,0};			// array contenente le risoluzioni in funz di ztrue
	float errrisol2[size] = {0,0,0,0,0,0,0,0,0,0};
	float binZsim[size] = {-13,-9.5,-6.5,-3.75,-1.25,1.25,3.75,6.5,9.5,13};
        float errbinZsim[size] = {2,1.5,1.5,1.25,1.25,1.25,1.25,1.5,1.5,2};
        float effZsim[size] = {0,0,0,0,0,0,0,0,0,0};
        float effZrec[size] = {0,0,0,0,0,0,0,0,0,0};
        float sigma2[size] = {0,0,0,0,0,0,0,0,0,0};
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
	
	static TH1D* histR0 = new TH1D("hR0","-15 < Ztrue < -11", 50, -0.05, 0.05);
	static TH1D* histR1 = new TH1D("hR1","-11 < Ztrue < -8", 50, -0.05, 0.05);
	static TH1D* histR2 = new TH1D("hR2","-8 < Ztrue < -5", 50, -0.05, 0.05);
	static TH1D* histR3 = new TH1D("hR3","-5 < Ztrue < -2.5", 50, -0.05, 0.05);
	static TH1D* histR4 = new TH1D("hR4","-2.5 < Ztrue < -0", 50, -0.05, 0.05);
	static TH1D* histR5 = new TH1D("hR5","0 < Ztrue < 2.5", 50, -0.05, 0.05);
	static TH1D* histR6 = new TH1D("hR6","2.5 < Ztrue < 5", 50, -0.04, 0.04);
	static TH1D* histR7 = new TH1D("hR7","5 < Ztrue < 8", 50, -0.04, 0.04);
	static TH1D* histR8 = new TH1D("hR8","8 < Ztrue < 11", 50, -0.04, 0.04);
	static TH1D* histR9 = new TH1D("hR9","11 < Ztrue < 15", 50, -0.04, 0.04);
	
	
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
	
	float Zrec[numeroeventi]; //!!!!FORSE NON SERVONO
    float Zsim[numeroeventi];
    int conto = 0;
	// loop sugli ingressi nel TTree
    for(int ev=0;ev<numeroeventi;ev++){
                
		tree3->GetEvent(ev);
        tree2->GetEvent(ev);
		tree4->GetEvent(ev);
		
		Zsim[ev] = point.Z; 
		
		
		if(pointRec.Z != 100.){	
			differenza = pointRec.Z	- Zsim[ev];
			if(differenza>1) conto++;
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
		
		if((point.Z>=-15.)&&(point.Z<=-11.)){
			effZsim[0]++; 
			histR0->Fill(differenza);
		}
		                                     
		if((point.Z>-11.)&&(point.Z<=-8.)){
			effZsim[1]++; 
			histR1->Fill(differenza);
		}
		                                    
		if((point.Z>-8.)&&(point.Z<=-5.)){
			effZsim[2]++; 
			histR2->Fill(differenza);
		}
		                                    
		if((point.Z>-5.)&&(point.Z<=-2.5)){
			effZsim[3]++; 
			histR3->Fill(differenza);
		}
		                                    
		if((point.Z>-2.5)&&(point.Z<=0.)){
			effZsim[4]++; 
			histR4->Fill(differenza);
		}
		                                     
		if((point.Z>0.)&&(point.Z<=2.5)){
			effZsim[5]++; 
			histR5->Fill(differenza);
		}		    
			
		if((point.Z>2.5)&&(point.Z<=5.)){
			effZsim[6]++; 
			histR6->Fill(differenza);
		}
		                                     
		if((point.Z>5.)&&(point.Z<=8.)){
			effZsim[7]++; 
			histR7->Fill(differenza);
		}
		                                     
		if((point.Z>8.)&&(point.Z<=11.)){
			effZsim[8]++; 
			histR8->Fill(differenza);
		}
		                                   
		if((point.Z>11.)&&(point.Z<=15.)){
			effZsim[9]++; 
			histR9->Fill(differenza);
		}
		
		
		
		if((pointRec.Z>=-15.)&&(pointRec.Z<=-11.))effZrec[0]++; 
		                                     
		if((pointRec.Z>-11.)&&(pointRec.Z<=-8.))effZrec[1]++; 
		                                    
		if((pointRec.Z>-8.)&&(pointRec.Z<=-5.))effZrec[2]++; 
		                                    
		if((pointRec.Z>-5.)&&(pointRec.Z<=-2.5))effZrec[3]++; 
		                                    
		if((pointRec.Z>-2.5)&&(pointRec.Z<=0.))effZrec[4]++; 
		                                     
		if((pointRec.Z>0.)&&(pointRec.Z<=2.5))effZrec[5]++; 
		                                     
		if((pointRec.Z>2.5)&&(pointRec.Z<=5.))effZrec[6]++; 
		                                     
		if((pointRec.Z>5.)&&(pointRec.Z<=8.))effZrec[7]++; 
		                                     
		if((pointRec.Z>8.)&&(pointRec.Z<=11.))effZrec[8]++; 
		                                   
		if((pointRec.Z>11.)&&(pointRec.Z<=15.))effZrec[9]++; 
	} 		// fine ciclo eventi
	
	for(int ii=0; ii<size; ii++){ eff2[ii] = effZrec[ii] / effZsim[ii];}
	cout<<"Quante differenze > 1: "<<conto<<endl;
	
	 // errore binomiale
		erreff2[0] = pow((eff2[0]*(1-eff2[0])/effZsim[0]),0.5);
		erreff2[1] = pow((eff2[1]*(1-eff2[1])/effZsim[1]),0.5);
		erreff2[2] = pow((eff2[2]*(1-eff2[2])/effZsim[2]),0.5);
		erreff2[3] = pow((eff2[3]*(1-eff2[3])/effZsim[3]),0.5);
		erreff2[4] = pow((eff2[4]*(1-eff2[4])/effZsim[4]),0.5);
		erreff2[5] = pow((eff2[5]*(1-eff2[5])/effZsim[5]),0.5);
		erreff2[6] = pow((eff2[6]*(1-eff2[6])/effZsim[6]),0.5);
		erreff2[7] = pow((eff2[7]*(1-eff2[7])/effZsim[7]),0.5);
		erreff2[8] = pow((eff2[8]*(1-eff2[8])/effZsim[8]),0.5);
		erreff2[9] = pow((eff2[9]*(1-eff2[9])/effZsim[9]),0.5);
	
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
		
		
		// rad(N) per moltep
		for(int i = 0; i < 9; i++){
			errmoltep[i] = pow(moltep[i],0.5);
		}
		
		
	// fit gaussiani risoluzione vs moltep
	
	cout << "************************** fit gaus risoluzione vs moltep ************************" << endl;
		TF1 *fM3 = new TF1("fM3","gaus",-0.05,0.05);
		fM3->SetLineColor(kRed);
		histM3->Fit(fM3,"NR+");
		gStyle->SetOptFit(0);
		double sigmaM3 = fM3->GetParameter(2);
		double errsigmaM3 = fM3->GetParError(2);
		risol[0] = sigmaM3;
		errrisol[0] = errsigmaM3;
		
		TF1 *fM5 = new TF1("fM5","gaus",-0.04,0.04);
		fM5->SetLineColor(kRed);
		histM5->Fit(fM5,"NR+");
		gStyle->SetOptFit(0);
		double sigmaM5 = fM5->GetParameter(2);
		double errsigmaM5 = fM5->GetParError(2);
		risol[1] = sigmaM5;
		errrisol[1] = errsigmaM5;
		
		TF1 *fM6 = new TF1("fM6","gaus",-0.04,0.04);
		fM6->SetLineColor(kRed);
		histM6->Fit(fM6,"NR+");
		gStyle->SetOptFit(0);
		double sigmaM6 = fM6->GetParameter(2);
		double errsigmaM6 = fM6->GetParError(2);
		risol[2] = sigmaM6;
		errrisol[2] = errsigmaM6;
		
		TF1 *fM7 = new TF1("fM7","gaus",-0.04,0.04);
		fM7->SetLineColor(kRed);
		histM7->Fit(fM7,"NR+");
		gStyle->SetOptFit(0);
		double sigmaM7 = fM7->GetParameter(2);
		double errsigmaM7 = fM7->GetParError(2);
		risol[3] = sigmaM7;
		errrisol[3] = errsigmaM7;
		
		TF1 *fM8 = new TF1("fM8","gaus",-0.04,0.04);
		fM8->SetLineColor(kRed);
		histM8->Fit(fM8,"NR+");
		gStyle->SetOptFit(0);
		double sigmaM8 = fM8->GetParameter(2);
		double errsigmaM8 = fM8->GetParError(2);
		risol[4] = sigmaM8;
		errrisol[4] = errsigmaM8;
		
		TF1 *fM12 = new TF1("fM12","gaus",-0.03,0.03);
		fM12->SetLineColor(kRed);
		histM12->Fit(fM12,"NR+");
		gStyle->SetOptFit(0);
		double sigmaM12 = fM12->GetParameter(2);
		double errsigmaM12 = fM12->GetParError(2);
		risol[5] = sigmaM12;
		errrisol[5] = errsigmaM12;
		
		TF1 *fM22 = new TF1("fM22","gaus",-0.025,0.025);
		fM22->SetLineColor(kRed);
		histM22->Fit(fM22,"NR+");
		gStyle->SetOptFit(0);
		double sigmaM22 = fM22->GetParameter(2);
		double errsigmaM22 = fM22->GetParError(2);
		risol[6] = sigmaM22;
		errrisol[6] = errsigmaM22;
		
		TF1 *fM32 = new TF1("fM32","gaus",-0.02,0.02);
		fM32->SetLineColor(kRed);
		histM32->Fit(fM32,"NR+");
		gStyle->SetOptFit(0);
		double sigmaM32 = fM32->GetParameter(2);
		double errsigmaM32 = fM32->GetParError(2);
		risol[7] = sigmaM32;
		errrisol[7] = errsigmaM32;
		
		TF1 *fM42 = new TF1("fM42","gaus",-0.015,0.015);
		fM42->SetLineColor(kRed);
		histM42->Fit(fM42,"NR+");
		gStyle->SetOptFit(0);
		double sigmaM42 = fM42->GetParameter(2);
		double errsigmaM42 = fM42->GetParError(2);
		risol[8] = sigmaM42;
		errrisol[8] = errsigmaM42;
		
		TF1 *fM52 = new TF1("fM52","gaus",-0.015,0.015);
		fM52->SetLineColor(kRed);
		histM52->Fit(fM52,"NR+");
		gStyle->SetOptFit(0);
		double sigmaM52 = fM52->GetParameter(2);
		double errsigmaM52 = fM52->GetParError(2);
		risol[9] = sigmaM52;
		errrisol[9] = errsigmaM52;
		
//______________________________ fit gaussiani risoluzione vs zsim____________________________
	
	cout << "******************* fit gaus risoluzione vs zsim ********************************" << endl;
		TF1 *fR0 = new TF1("fR0","gaus",-0.03,0.03);
		fR0->SetLineColor(kRed);
		histR0->Fit(fR0,"NR+");
		gStyle->SetOptFit(0);
		double sigmaR0 = fR0->GetParameter(2);
		double mediaR0 = TMath::Abs(fR0->GetParameter(1));
		double errsigmaR0 = fR0->GetParError(2);
		risol2[0] = sigmaR0/mediaR0;
		sigma2[0] = sigmaR0;
		errrisol2[0] = errsigmaR0;
		
		TF1 *fR1 = new TF1("fR1","gaus",-0.03,0.03);
		fR1->SetLineColor(kRed);
		histR1->Fit(fR1,"NR+");
		gStyle->SetOptFit(0);
		double sigmaR1 = fR1->GetParameter(2);
		double mediaR1 = TMath::Abs(fR1->GetParameter(1));
		double errsigmaR1 = fR1->GetParError(2);
		risol2[1] = sigmaR1/mediaR1;
		sigma2[1] = sigmaR1;
		errrisol2[1] = errsigmaR1;
		
		TF1 *fR2 = new TF1("fR2","gaus",-0.03,0.03);
		fR2->SetLineColor(kRed);
		histR2->Fit(fR2,"NR+");
		gStyle->SetOptFit(0);
		double sigmaR2 = fR2->GetParameter(2);
		double mediaR2 = TMath::Abs(fR2->GetParameter(1));
		double errsigmaR2 = fR2->GetParError(2);
		risol2[2] = sigmaR2/mediaR2;
		sigma2[2] = sigmaR2;
		errrisol2[2] = errsigmaR2;
		
		TF1 *fR3 = new TF1("fR3","gaus",-0.03,0.03);
		fR3->SetLineColor(kRed);
		histR3->Fit(fR3,"NR+");
		gStyle->SetOptFit(0);
		double sigmaR3 = fR3->GetParameter(2);
		double mediaR3 = TMath::Abs(fR3->GetParameter(1));
		double errsigmaR3 = fR3->GetParError(2);
		risol2[3] = sigmaR3/mediaR3;
		sigma2[3] = sigmaR3;
		errrisol2[3] = errsigmaR3;
		
		TF1 *fR4 = new TF1("fR4","gaus",-0.03,0.03);
		fR4->SetLineColor(kRed);
		histR4->Fit(fR4,"NR+");
		gStyle->SetOptFit(0);
		double sigmaR4 = fR4->GetParameter(2);
		double mediaR4 = TMath::Abs(fR4->GetParameter(1));
		double errsigmaR4 = fR4->GetParError(2);
		risol2[4] = sigmaR4/mediaR4;
		sigma2[4] = sigmaR4;
		errrisol2[4] = errsigmaR4;
		
		TF1 *fR5 = new TF1("fR5","gaus",-0.03,0.03);
		fR5->SetLineColor(kRed);
		histR5->Fit(fR5,"NR+");
		gStyle->SetOptFit(0);
		double sigmaR5 = fR5->GetParameter(2);
		double mediaR5 = TMath::Abs(fR5->GetParameter(1));
		double errsigmaR5 = fR5->GetParError(2);
		risol2[5] = sigmaR5/mediaR5;
		sigma2[5] = sigmaR5;
		errrisol2[5] = errsigmaR5;
		
		TF1 *fR6 = new TF1("fR6","gaus",-0.03,0.03);
		fR6->SetLineColor(kRed);
		histR6->Fit(fR6,"NR+");
		gStyle->SetOptFit(0);
		double sigmaR6 = fR6->GetParameter(2);
		double mediaR6 = TMath::Abs(fR6->GetParameter(1));
		double errsigmaR6 = fR6->GetParError(2);
		risol2[6] = sigmaR6/mediaR6;
		sigma2[6] = sigmaR6;
		errrisol2[6] = errsigmaR6;
		
		TF1 *fR7 = new TF1("fR7","gaus",-0.03,0.03);
		fR7->SetLineColor(kRed);
		histR7->Fit(fR7,"NR+");
		gStyle->SetOptFit(0);
		double sigmaR7 = fR7->GetParameter(2);
		double mediaR7 = TMath::Abs(fR7->GetParameter(1));
		double errsigmaR7 = fR7->GetParError(2);
		risol2[7] = sigmaR7/mediaR7;
		sigma2[7] = sigmaR7;
		errrisol2[7] = errsigmaR7;
		
		TF1 *fR8 = new TF1("fR8","gaus",-0.03,0.03);
		fR8->SetLineColor(kRed);
		histR8->Fit(fR8,"NR+");
		gStyle->SetOptFit(0);
		double sigmaR8 = fR8->GetParameter(2);
		double mediaR8 = TMath::Abs(fR8->GetParameter(1));
		double errsigmaR8 = fR8->GetParError(2);
		risol2[8] = sigmaR8/mediaR8;
		sigma2[8] = sigmaR8;
		errrisol2[8] = errsigmaR8;
		
		TF1 *fR9 = new TF1("fR9","gaus",-0.03,0.03);
		fR9->SetLineColor(kRed);
		histR9->Fit(fR9,"NR+");
		gStyle->SetOptFit(0);
		double sigmaR9 = fR9->GetParameter(2);
		double mediaR9 = TMath::Abs(fR9->GetParameter(1));
		double errsigmaR9 = fR9->GetParError(2);
		risol2[9] = sigmaR9/mediaR9;
		sigma2[9] = sigmaR9;
		errrisol2[9] = errsigmaR9;
		
		
// __________________________________________ efficienza vs moltep _________________________________________
		//TCanvas *c1=new TCanvas("c1","c1",800,600);
		TGraphErrors *graphE= new TGraphErrors(size,moltep,eff,errmoltep,erreff);
		graphE->SetMarkerSize(1);//https://root.cern.ch/doc/master/classTAttMarker.html
		graphE->SetMarkerStyle(33);
		graphE->SetTitle("Efficiency Vs Molteplicity");
		graphE->GetXaxis()->SetTitle("Molteplicity[]");
		graphE->GetYaxis()->SetTitle("Efficiency[]");
		
		//graphE->Draw("ap");
// __________________________________________ efficienza vs Ztrue _________________________________________		
		//TCanvas *c1=new TCanvas("c1","c1",800,600);
		TGraphErrors *graphE2= new TGraphErrors(size,binZsim,eff2,errbinZsim,erreff2);
		graphE2->SetMarkerSize(1);//https://root.cern.ch/doc/master/classTAttMarker.html
		graphE2->SetMarkerStyle(33);
		graphE2->SetTitle("Efficiency Vs Ztrue");
		graphE2->GetXaxis()->SetTitle("Ztrue[cm]");
		graphE2->GetYaxis()->SetTitle("Efficiency[]");
		//c1->cd();
		//graphE2->Draw("ap");
// ____________________________________________ risoluzione vs moltep ________________________________________		
		TGraphErrors *graphR= new TGraphErrors(size,moltep,risol,errmoltep,errrisol);
		graphR->SetMarkerSize(1);//https://root.cern.ch/doc/master/classTAttMarker.html
		graphR->SetMarkerStyle(33);
		graphR->SetTitle("Resolution Vs Molteplicity");
		graphR->GetXaxis()->SetTitle("Molteplicity[]");
		graphR->GetYaxis()->SetTitle("Resolution[cm]");
		//c1->cd();
		//graphR->Draw("ap");
// ____________________________________________ risoluzione vs Ztrue ________________________________________		
		TGraphErrors *graphR2= new TGraphErrors(size,binZsim,sigma2,errbinZsim,errrisol2);
		graphR2->SetMarkerSize(1);//https://root.cern.ch/doc/master/classTAttMarker.html
		graphR2->SetMarkerStyle(33);
		graphR2->SetTitle("Resolution Vs Ztrue");
		graphR2->GetXaxis()->SetTitle("Ztrue[cm]");
		graphR2->GetYaxis()->SetTitle("Resolution[cm]");
		//c1->cd();
		//graphR2->Draw("ap");
		
		
		
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
		histR0->Write();
		histR1->Write();
		histR2->Write();
		histR3->Write();
		histR4->Write();
		histR5->Write();
		histR6->Write();
		histR7->Write();
		histR8->Write();
		histR9->Write();
        graphE->Write();
		graphE2->Write();
		graphR->Write();
		graphR2->Write();
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
		delete histR0;
		delete histR1;
		delete histR2;
		delete histR3;
		delete histR4;
		delete histR5;
		delete histR6;
		delete histR7;
		delete histR8;
		delete histR9;
		
        delete graphE;
        graphE = nullptr;
		delete graphE2;
        graphE2 = nullptr;
		delete graphR;
        graphR = nullptr;
		delete graphR2;
        graphR2 = nullptr;
		
		
		double TT = time.CpuTime();	
         cout<<"Il tempo impiegato dalla CPU è "<<TT<<" s"<<endl;  
         
         MemInfo_t memInfo;
         gSystem->GetMemInfo(&memInfo);
cout << "Mem Used = " << memInfo.fMemUsed << " MB"<<endl; //returning value in MB
		
}
