#include <iostream>
#include <fstream>
#include <math.h>

//#include "Punto.h"
#include "Punto2.h"
#include "TClonesArray.h"
#include "Vertex.h"
#include "TracciaMC.h"
#include "Multis.h"
#include "Smearing.h"
#include "TRandom3.h"

#include <vector>

using namespace std;

void Ricostruzione2(){
	
	int numeroeventi;
	
	TStopwatch time;
	
	//TCanvas *c1=new TCanvas("c1","c1",800,600);
    static TH1D* hist1 = new TH1D ("histofin","Zrec-Ztrue" , 250, -0.05, 0.05);
	static TH1D* hist2 = new TH1D ("h2","Zreconstructed versione2" , 54, -13.5, 13.5);
	// histo per efficienza
/*	static TH1D* histM3 = new TH1D("hM3","2.5 < molteplicità < 3.5", 150, -0.15, 0.15);// fatto per molteplicità 3
	static TH1D* histM5 = new TH1D("hM5","4.5 < molteplicità < 5.5", 150, -0.15, 0.15);
	static TH1D* histM6 = new TH1D("hM6","5.5 < molteplicità < 6.5", 150, -0.15, 0.15);
	static TH1D* histM7 = new TH1D("hM7","6.5 < molteplicità < 7.5", 150, -0.15, 0.15);
	static TH1D* histM8 = new TH1D("hM8","7.5 < molteplicità < 8.5", 150, -0.15, 0.15);
	static TH1D* histM12 = new TH1D("hM12","11.5 < molteplicità < 12.5", 250, -0.15, 0.15);
	static TH1D* histM22 = new TH1D("hM22","21.5 < molteplicità < 22.5", 250, -1.5, 1.5);
	static TH1D* histM32 = new TH1D("hM32","31.5 < molteplicità < 32.5", 250, -1.5, 1.5);
	static TH1D* histM42 = new TH1D("hM42","41.5 < molteplicità < 42.5", 250, -1.5, 1.5);
	static TH1D* histM52 = new TH1D("hM52","51.5 < molteplicità < 52.5", 250, -1.5, 1.5);
*/	
   
	static int iter; //iteratore
	float deltaPhi;
	
    float Z3, Z2, TanTheta;
			
	float q;      //termine noto retta
				
	static float R2 = 4.;
	static float R3 = 7.;
			
	float differenza = 0;
	//float Zsim1;
	//float Zrec1;
	float Zrec2;
	//float nZrec = 0;		// float altrimenti approssima l'eff all'int 0
	const int size = 10; 			// dimensione array per TGraph
	float moltep[size] = {0,0,0,0,0,0,0,0,0,0}; 		// array contenente le moltep di ogni evento da passare al TGraph
	float errmoltep[size] = {0,0,0,0,0,0,0,0,0,0};	// binomiale??
	float eff[size] = {0,0,0,0,0,0,0,0,0,0};			// array contenente le efficienze
	float erreff[size] = {0,0,0,0,0,0,0,0,0,0};
	
	
        //definizione struct
        typedef struct {
           float X,Y,Z;
           int mult;} VTX;
        static VTX point;
		static VTX pointRec;
		
		typedef struct{
			float x0,x1,x2,x3,x4,x5,x6,x7,x8,x9;		// Z simulate con un certo valore di moltep, sarà denominatore dell'efficienza
			} vett;
		static vett denEff;
  
        //apertura file di input
        TFile hfile2("htree2.root");
  
        //lettura TTree  e branch
        TTree *tree2 = (TTree*)hfile2.Get("T2");
        TBranch *b1=tree2->GetBranch("VertMult");
        TBranch *b2=tree2->GetBranch("HitsPrimo");
        TBranch *b3=tree2->GetBranch("HitsSecondo");
        TBranch *b4=tree2->GetBranch("HitsTerzo");
		
		TFile hfile3("htree3.root");
		TTree *tree3 = (TTree*)hfile3.Get("T3");
		TBranch *branch3 = tree3->GetBranch("denEff");
		
		TFile hfile4("htree4.root", "RECREATE");
		TTree *tree4 = new TTree("T4","TTree dati efficienza");
		tree4->Branch("Zrec",&pointRec.X,"X/F:Y:Z:mult/I");
		
		
		
        numeroeventi = tree2->GetEntries(); //acquisico informazione sul numero di eventi nel mio detector
        	
	//arrays contenenti Z ricostruite e simulate
  //float Zrec[numeroeventi];
        float Zsim[numeroeventi];	
		
				
  
        //dichiarazione TClonesArray
        //TClonesArray *hits1 = new TClonesArray("Punto2",numeroeventi); //per la ricostruzione non usiamo il primo layer
        TClonesArray *hits2 = new TClonesArray("Punto2",100);
        TClonesArray &ptrhits2 = *hits2; 
      
        TClonesArray *hits3 = new TClonesArray("Punto2",100);
        TClonesArray &ptrhits3 = *hits3;
		
		//float denEff[10];
  
        //definizione degli indirizzi per la lettura dei dati su ttree
        b1->SetAddress(&point.X);
		
        //b2->SetAddress(&hits1); //relativo al primo layer
        b3->SetAddress(&hits2);
        b4->SetAddress(&hits3);
		branch3->SetAddress(&denEff.x0);
  
        //arrays in cui salvo le buone combinazioni di hits
    /*    TClonesArray *hitsgood2 = new TClonesArray("Punto2",100);
        TClonesArray &ptrhitsgood2 = *hitsgood2;
      
        TClonesArray *hitsgood3 = new TClonesArray("Punto2",100);
        TClonesArray &ptrhitsgood3 = *hitsgood3;
    */    
        //provo invece con vectors di float
        vector<float> Zgood2;
        vector<float> Zgood3;
		
		
		
	// loop sugli ingressi nel TTree
        for(int ev=0;ev<numeroeventi;ev++){
                
                iter = 0;
				pointRec.X = 0;
				pointRec.Y = 0;
				pointRec.mult = 0;
				
				tree3->GetEvent(ev);
                tree2->GetEvent(ev);
                //cout<<"Evento "<<ev<<"; Molteplicita= "<<point.mult<<endl; 
                    
                Zsim[ev] = point.Z; //riempio l'array coi dati simulati

                for(int a = 0; a < point.mult; a++){		// loop su hits layer 2
			     
					Punto2 *tst2;
					tst2=(Punto2*)hits2->At(a);
						
					if(tst2->Getnum()!=-1){			// escludiamo hits fuori dal rivelatore
						for(int b = 0; b < point.mult; b++){		// loop su hits layer 3
				  
						Punto2 *tst3;
						tst3=(Punto2*)hits3->At(b);
						if(tst3->Getnum()!=-1){
							deltaPhi = tst2->GetPhi() - tst3->GetPhi();
								if(abs(deltaPhi) < 0.01){						
									//new(ptrhitsgood2[iter])Punto2(0., tst2->GetZ(), tst2->Getnum());
									//new(ptrhitsgood3[iter])Punto2(0., tst3->GetZ(), tst3->Getnum());
									Zgood2.push_back(tst2->GetZ());
									Zgood3.push_back(tst3->GetZ());
									iter++;
								}
					//else continue;
							}
				   //else continue;
						}
			       }
			   //else continue;
			   }
			
				
	         // __________ricostruzione vertice__________
			if(Zgood2.size() != Zgood3.size()) cout << " Zgood2 e Zgood3 hanno dimensioni diverse: " << Zgood2.size() << " " << Zgood3.size() << endl;	
			
		// for(int e=0; e<ptrhitsgood2.GetEntries(); e++){	
			for(int e=0; e<Zgood2.size(); e++){		// ciclo sui punti buoni
			
		/*	Punto2 *tstgood2;
			Punto2 *tstgood3;
			tstgood2=(Punto2*)hitsgood2->At(e);
			tstgood3=(Punto2*)hitsgood3->At(e);
				
			Z2 = tstgood2->GetZ();
				
			Z3 = tstgood3->GetZ();
		*/	
				Z2 = Zgood2[e];
					
				Z3 = Zgood3[e];
				
				TanTheta = (R3-R2)/(Z3-Z2);
				
				q = R3 - Z3*TanTheta; //termine noto retta
				
				Zrec2 = -q/TanTheta;
					
				hist2 -> Fill(Zrec2);		// riempio l'histo con gli Zrec dati dalle coppie "buone"
			
			}
					
	        //Zrec[ev] = hist2->GetMean(); //salvo la media di tutto il mio istogramma, di tutte le Zrec2 (anche quelle derivate da coppie sbagliate) nel vector Zrec (solo una per evento)	
	        
	        int Ymax = 0;
			int Xmax = 0;	
		
        	for(int r = 0; r < 54 ;r++){		// 54 bin decisi arbitrariamente
        	       //cout<<"Bin numero  "<<r<<" valore associato di Z "<<-13.75 + r*0.5<<" numero di Zrec "<<hist2->GetBinContent(r)<<endl;
        	       if(Ymax < hist2->GetBinContent(r)){
        	                  Ymax = hist2->GetBinContent(r);
        	                  
        	                  Xmax = r;
        	                  }
        	        }
        	// sotto le 2 tracce non ricostruiamo -> associamo Z=100 per poi scartare il dato
        	if(Ymax > 2){
				//if(hist2->GetMean() < 5.3*3 && hist2->GetMean() > 5.3*3){		// controllo su 3sigma
	//Zrec[ev] = hist2->GetMean(); //Zrec[ev] = -13.75 + Xmax*0.5;
				pointRec.Z = hist2->GetMean();
				tree4->Fill();
				//nZrec++;		// conto Z buone
				//}
			}
        	else{pointRec.Z = 100.;
				tree4->Fill();
			}		// Zrec[ev] = 100.;

               
                hist2->Reset("ICESM");
                //ptrhitsgood2.Clear();
				//ptrhitsgood3.Clear();
				Zgood2.clear();
				Zgood3.clear();
				
			
				
			if(pointRec.Z != 100.){		// Zrec[ev] != 100.
			differenza = pointRec.Z	- Zsim[ev];		// Zrec[ev] - Zsim[ev];
			hist1->Fill(differenza);
			
			/*
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
				*/
			}
			
		
		}		// fine ciclo eventi
		
		// # di Z ricostruite con una certa moltep su quelle simulate
	/*	eff[0] = histM3->GetEntries()/denEff.x0;		
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
		TF1 *fM3 = new TF1("fM3","gaus",-0.06,0.06);
		fM3->SetLineColor(kRed);
		//fM3->SetParameters(1,60);	// media
		//fM3->SetParameters(2,7);	// sigma
		histM3->Fit(fM3,"NR+");
		//histM3->Draw("same");
		
		TF1 *fM5 = new TF1("fM5","gaus",-0.06,0.06);
		fM5->SetLineColor(kRed);
		//fM5->SetParameters(1,60);	// media
		//fM5->SetParameters(2,7);	// sigma
		histM5->Fit(fM5,"NR+");
		//histM5->Draw("same");
		
		TF1 *fM6 = new TF1("fM6","gaus",-0.07,0.07);
		fM6->SetLineColor(kRed);
		//fM6->SetParameters(1,60);	// media
		//fM6->SetParameters(2,7);	// sigma
		histM6->Fit(fM6,"NR+");
		//histM6->Draw("same");
		
		TF1 *fM7 = new TF1("fM7","gaus",-0.07,0.07);
		fM7->SetLineColor(kRed);
		//fM7->SetParameters(1,60);	// media
		//fM7->SetParameters(2,7);	// sigma
		histM7->Fit(fM7,"NR+");
		//histM7->Draw("same");
		
		TF1 *fM8 = new TF1("fM8","gaus",-0.05,0.05);
		fM8->SetLineColor(kRed);
		//fM8->SetParameters(1,60);	// media
		//fM8->SetParameters(2,7);	// sigma
		histM8->Fit(fM8,"NR+");
		//histM8->Draw("same");
		
		TF1 *fM12 = new TF1("fM12","gaus",-0.05,0.05);
		fM12->SetLineColor(kRed);
		//fM12->SetParameters(1,60);	// media
		//fM12->SetParameters(2,7);	// sigma
		histM12->Fit(fM12,"NR+");
		//histM12->Draw("same");
		
		TF1 *fM22 = new TF1("fM22","gaus",-0.04,0.04);
		fM22->SetLineColor(kRed);
		//fM22->SetParameters(1,60);	// media
		//fM22->SetParameters(2,7);	// sigma
		histM22->Fit(fM22,"NR+");
		//histM22->Draw("same");
		
		TF1 *fM32 = new TF1("fM32","gaus",-0.04,0.04);
		fM32->SetLineColor(kRed);
		//fM32->SetParameters(1,60);	// media
		//fM32->SetParameters(2,7);	// sigma
		histM32->Fit(fM32,"NR+");
		//histM32->Draw("same");
		
		TF1 *fM42 = new TF1("fM42","gaus",-0.03,0.03);
		fM42->SetLineColor(kRed);
		//fM42->SetParameters(1,60);	// media
		//fM42->SetParameters(2,7);	// sigma
		histM42->Fit(fM42,"NR+");
		//histM42->Draw("same");
		
		TF1 *fM52 = new TF1("fM52","gaus",-0.03,0.03);
		fM52->SetLineColor(kRed);
		//fM52->SetParameters(1,60);	// media
		//fM52->SetParameters(2,7);	// sigma
		histM52->Fit(fM52,"NR+");
		//histM52->Draw("same");
				
		for(int i = 0; i < 10; i++){
			 cout << i << " " << eff[i] << " " << moltep[i] << endl;
			
		 }
		
		
		//TCanvas *c1=new TCanvas("c1","c1",800,600);
		TGraphErrors *graphE= new TGraphErrors(size,moltep,eff,errmoltep,erreff);
		graphE->SetMarkerSize(0.9);//https://root.cern.ch/doc/master/classTAttMarker.html
		graphE->SetMarkerStyle(43);
		graphE->SetTitle("Efficiency vs molteplicity");
		graphE->GetXaxis()->SetTitle("molteplicity[]");
		graphE->GetYaxis()->SetTitle("efficiency[]");
		//c1->cd();
		graphE->Draw("ap");
	*/	
		delete hist2;
	
        hfile2.Close();
			
	//creo grafico Ztrue-Zrec, avrò numeroeventi dati
/*	TCanvas *c1=new TCanvas("c1","c1",800,600);
        TH1D* hist1 = new TH1D ("histofin","Ztrue-Zrec" , 250, -0.05, 0.05);	
	         		 
	for(int u = 0; u < numeroeventi; u++){
		 
	      Zsim1 = Zsim[u];
	      Zrec1 = Zrec[u];
	    
	     if(Zrec1 != 100.){
	      differenza = Zsim1 - Zrec1;
		 
	      hist1->Fill(differenza); //grafico la differenza
	      }}
*/		 		 
	 hist1->Draw("ep");
	 
	 TFile file1("zrec-ztrue-tree.root", "recreate");
	 hist1->Write();
	/* histM3->Write(); 
	 histM5->Write(); 
	 histM6->Write(); 
	 histM7->Write(); 
	 histM8->Write(); 
	 histM12->Write();
	 histM22->Write();
	 histM32->Write();
	 histM42->Write();
	 histM52->Write();
	 graphE->Write();*/
	 file1.Close();
	 delete hist1;
	/* delete histM3;
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
	*/ 
	 hfile4.Write();              
     hfile4.Close();
				
         double TT = time.CpuTime();	
         cout<<"Il tempo impiegato dalla CPU è "<<TT<<" s"<<endl;  
  
         }
 
