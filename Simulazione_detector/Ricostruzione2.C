#include <iostream>
#include <fstream>
#include <math.h>

#include "Punto.h"
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
	TH1D* hist2 = new TH1D ("h2","Zreconstructed versione2" , 54, -13.5, 13.5);
	
	int iter; //iteratore
	float deltaPhi;
	
        float Z3, Z2, TanTheta;
			
	float q;      //termine noto retta
				
	float R2 = 4.;
	float R3 = 7.;
			
	float differenza = 0;
	float Zsim1;
	float Zrec1, Zrec2;
			
        //definizione struct
        typedef struct {
           float X,Y,Z;
           int mult;} VTX;
        static VTX point;
  
        //apertura file di input
        TFile hfile2("htree2.root");
  
        //lettura TTree  e branch
        TTree *tree2 = (TTree*)hfile2.Get("T2");
        TBranch *b1=tree2->GetBranch("VertMult");
        TBranch *b2=tree2->GetBranch("HitsPrimo");
        TBranch *b3=tree2->GetBranch("HitsSecondo");
        TBranch *b4=tree2->GetBranch("HitsTerzo");
  
        numeroeventi = tree2->GetEntries(); //acquisico informazione sul numero di eventi nel mio detector
        	
	//arrays contenenti Z ricostruite e simulate
        float Zrec[numeroeventi];
        float Zsim[numeroeventi];	
  
        //dichiarazione TClonesArray
        //TClonesArray *hits1 = new TClonesArray("Punto2",numeroeventi); //per la ricostruzione non usiamo il primo layer
        TClonesArray *hits2 = new TClonesArray("Punto2",100);
        TClonesArray &ptrhits2 = *hits2; 
      
        TClonesArray *hits3 = new TClonesArray("Punto2",100);
        TClonesArray &ptrhits3 = *hits3;
  
        //definizione degli indirizzi per la lettura dei dati su ttree
        b1->SetAddress(&point.X);
        //b2->SetAddress(&hits1); //relativo al primo layer
        b3->SetAddress(&hits2);
        b4->SetAddress(&hits3);
  
        //arrays in cui salvo le buone combinazioni di hits
        TClonesArray *hitsgood2 = new TClonesArray("Punto2",100);
        TClonesArray &ptrhitsgood2 = *hitsgood2;
      
        TClonesArray *hitsgood3 = new TClonesArray("Punto2",100);
        TClonesArray &ptrhitsgood3 = *hitsgood3;
        
        //provo invece con vectors di float
        //vector<float> Zgood2;
        //vector<float> Zgood3;
	
  	// loop sugli ingressi nel TTree
        for(int ev=0;ev<numeroeventi;ev++){
                
                iter = 0;
                
                tree2->GetEvent(ev);
                //cout<<"Evento "<<ev<<"; Molteplicita= "<<point.mult<<endl; 
                // cout<<"X,Y,Z = "<<point.X<<"; "<<point.Y<<"; "<<point.Z<<endl;
                    
                Zsim[ev] = point.Z; //riempio il vector coi dati simulati

                for(int a = 0; a < point.mult; a++){
			     
			  Punto2 *tst2;
			  tst2=(Punto2*)hits2->At(a);
				
			  if(tst2->Getnum()!=-1){
			      for(int b = 0; b < point.mult; b++){
				  
				  Punto2 *tst3;
				  tst3=(Punto2*)hits3->At(b);
				  if(tst3->Getnum()!=-1){
					deltaPhi = tst2->GetPhi() - tst3->GetPhi();
					if(abs(deltaPhi) < 0.01){						
						new(ptrhitsgood2[iter])Punto2(0., tst2->GetZ(), tst2->Getnum());
						new(ptrhitsgood3[iter])Punto2(0., tst3->GetZ(), tst3->Getnum());
						//Zgood2.push_back(tst2->GetZ());
						//Zgood3.push_back(tst3->GetZ());
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
			
		 for(int e=0; e<ptrhitsgood2.GetEntries(); e++){
		 //for(int e=0; e<Zgood2.size(); e++){
			
			Punto2 *tstgood2;
			Punto2 *tstgood3;
			tstgood2=(Punto2*)hitsgood2->At(e);
			tstgood3=(Punto2*)hitsgood3->At(e);
				
			Z2 = tstgood2->GetZ();
				
			Z3 = tstgood3->GetZ();
			
			//Z2 = Zgood2[e];
				
			//Z3 = Zgood3[e];
																
		 	TanTheta = (R3-R2)/(Z3-Z2);
		 	
		 	q = R3 - Z3*TanTheta; //termine noto retta
			
			Zrec2 = -q/TanTheta;
				 
			hist2 -> Fill(Zrec2);
							
			}
					
	        //Zrec[ev] = hist2->GetMean(); //salvo la media di tutto il mio istogramma, di tutte le Zrec2 (anche quelle derivate da coppie sbagliate) nel vector Zrec (solo una per evento)	
	        
	        int Ymax = 0;
		int Xmax = 0;	
		
        	for(int r = 0; r < 54 ;r++){
        	       //cout<<"Bin numero  "<<r<<" valore associato di Z "<<-13.75 + r*0.5<<" numero di Zrec "<<hist2->GetBinContent(r)<<endl;
        	       if(Ymax < hist2->GetBinContent(r)){
        	                  Ymax = hist2->GetBinContent(r);
        	                  
        	                  Xmax = r;
        	                  }
        	        }
        	
        	if(Ymax > 4) Zrec[ev] = hist2->GetMean();//Zrec[ev] = -13.75 + Xmax*0.5;
        	else Zrec[ev] = 100.;
        	
               
                hist2->Reset("ICESM");
                ptrhitsgood2.Clear();
		ptrhitsgood3.Clear();
		//Zgood2.clear();
		//Zgood3.clear();
		
	}
		
	delete hist2;
	
        hfile2.Close();
			
	//creo grafico Ztrue-Zrec, avrò numeroeventi dati
	TCanvas *c1=new TCanvas("c1","c1",800,600);
        TH1D* hist1 = new TH1D ("histofin","Ztrue-Zrec" , 250, -0.05, 0.05);	
	         		 
	for(int u = 0; u < numeroeventi; u++){
		 
	      Zsim1 = Zsim[u];
	      Zrec1 = Zrec[u];
	    
	     if(Zrec1 != 100.){
	      differenza = Zsim1 - Zrec1;
		 
	      hist1->Fill(differenza); //grafico la differenza
	      }}
		 		 
	 hist1->Draw("ep");
	 
	 TFile file1("ztrue-zrec-tree.root", "recreate");
	 hist1->Write();
	 file1.Close();
				
         double TT = time.CpuTime();	
         cout<<"Il tempo impiegato dalla CPU è "<<TT<<" s"<<endl;  
  
         }
 
