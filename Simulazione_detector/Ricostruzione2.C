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
	
	TH1D* hist2 = new TH1D ("h2","Zreconstructed versione2" , 80, -14., 14.);
	
	int iter; //iteratore
	float deltaPhi;
	
        float X2, Y2, Z2, Phi2; //coordinate intersezione secondo layer
        float X3, Y3, Z3, Phi3;   //coordinate intersezione terzo layer
			
	float C[3];   //coefficienti retta
	float t;      //parametro retta
				
	float R2 = 4.;
	float R3 = 7.;
			
	float Phi, Theta;
			
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
						new(ptrhitsgood2[iter])Punto2(tst2->GetPhi(), tst2->GetZ(), tst2->Getnum());
						new(ptrhitsgood3[iter])Punto2(tst3->GetPhi(), tst3->GetZ(), tst3->Getnum());
						iter++;
						}
					//else continue;
					}
				   //else continue;
				   }
			       }
			   //else continue;
			   }
			
				
	         // __________tentativo ricostruzione vertice__________
		 // eq retta passante per due punti A(x2,y2,z2) e B(x3,y3,z3)
		 // x0 = x3 + c1t		c1 = sin(theta)*cos(phi)
		 // y0 = y3 + c2t		c2 = sin(theta)*sin(phi)
		 // z0 = z3 + c3t		c3 = cos(phi)
		 // noi vogliamo solo la z
		 // impongo y=0 e x=0, trovo quindi due valori di t e li medio tra loro (x3/c1 e y3/c2)
			
				
		 for(int e=0; e<ptrhitsgood2.GetEntries(); e++){
			
			Punto2 *tstgood2;
			Punto2 *tstgood3;
			tstgood2=(Punto2*)hitsgood2->At(e);
			tstgood3=(Punto2*)hitsgood3->At(e);
				
			Phi2 = tstgood2->GetPhi();
			X2 = R2*TMath::Cos(Phi2);
			Y2 = R2*TMath::Sin(Phi2);
			Z2 = tstgood2->GetZ();
				
			Phi3 = tstgood3->GetPhi();
			X3 = R3*TMath::Cos(Phi3);
			Y3 = R3*TMath::Sin(Phi3);
			Z3 = tstgood3->GetZ();
																
			//calcolo Phi e Theta della retta passante tra i due punti
		        if(X3 < X2){
			    Phi = TMath::ATan((Y3-Y2)/(X3-X2)) + TMath::Pi();
			    }
			else Phi = TMath::ATan((Y3-Y2)/(X3-X2));
								
			if(Z3 < Z2){
			    Theta = TMath::ATan((R3-R2)/(Z3-Z2)) + TMath::Pi();
			    }
			else Theta = TMath::ATan((R3-R2)/(Z3-Z2));
				
		        C[0] =  TMath::Sin(Theta)*TMath::Cos(Phi);		
			C[1] =  TMath::Sin(Theta)*TMath::Sin(Phi);		        
			C[2] =  TMath::Cos(Theta);	
			
		        t = (X3/C[0] + Y3/C[1])/2;
				 
		        //l'intersezione con la beam line sarà il mio Zrec2
			Zrec2 = Z3 - C[2]*t;
				 
			hist2 -> Fill(Zrec2);				
			}
						
		Zrec[ev] = hist2->GetMean(); //salvo la media di tutto il mio istogramma, di tutte le Zrec2 (anche quelle derivate da coppie sbagliate) nel vector Zrec (solo una per evento)	
        		
                hist2->Reset("ICESM");
                ptrhitsgood2.Clear();
		ptrhitsgood3.Clear();
		}
	
		
	delete hist2;
        hfile2.Close();
			
	//creo grafico Ztrue-Zrec, avrò numeroeventi dati
	TCanvas *c1=new TCanvas("c1","c1",800,600);
        TH1D* hist1 = new TH1D ("h2","Ztrue-Zrec" , 100, -1, 1);	
	         		 
	for(int u = 0; u < numeroeventi; u++){
		 
	      Zsim1 = Zsim[u];
	      Zrec1 = Zrec[u];
		 
	      differenza = Zsim1 - Zrec1;
		 
	      hist1->Fill(differenza); //grafico la differenza
	      }
		 		 
	 hist1->Draw();
	 
	 TFile file1("ztrue-zrec-tree.root", "recreate");
	 hist1->Write();
	 file1.Close();
				
         double TT = time.CpuTime();	
         cout<<"Il tempo impiegato dalla CPU è "<<TT<<" s"<<endl;
  
         }
 
