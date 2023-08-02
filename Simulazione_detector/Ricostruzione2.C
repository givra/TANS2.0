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

//float T (float x3, float y3, float C[3]);



void Ricostruzione2(){
	
	int numeroeventi;
	TStopwatch time;
	TH1D* hist2 = new TH1D ("h2","zreconstructed versione2" , 80, -14., 14.);
	
	
	        float X2, Y2, Z2, Phi2; //coordinate intersezione secondo layer
			
		float X3, Y3, Z3, Phi3;   //coordinate intersezione terzo layer
			
		float C[3];   //coefficienti retta
		float t;      //parametro retta
				
		float R2 = 4.;
		float R3 = 7.;
			
		float Phi, Theta;
			
		float differenza = 0;
		float Zsim1;
	        float Zrec1, Zrec;
			
			
		 
 //----------------------------------------lettura del tree		

  // definizione struct
  typedef struct {
    float X,Y,Z;
    int mult;} VTX;
  static VTX point;
  //Apertura file di input
  TFile hfile2("htree2.root");
  
  //Lettura TTree  e branch
  TTree *tree2 = (TTree*)hfile2.Get("T2");
  TBranch *b1=tree2->GetBranch("VertMult");
  TBranch *b2=tree2->GetBranch("HitsPrimo");
  TBranch *b3=tree2->GetBranch("HitsSecondo");
  TBranch *b4=tree2->GetBranch("HitsTerzo");
  
  // Dichiarazione TClonesArray
  //TClonesArray *hits1 = new TClonesArray("Punto2",numeroeventi);
  TClonesArray *hits2 = new TClonesArray("Punto2",100);
  TClonesArray &ptrhits2 = *hits2;									// per qualche motivo serve l'alias per poi riempirlo
  TClonesArray *hits3 = new TClonesArray("Punto2",100);
  TClonesArray &ptrhits3 = *hits3;
  
  // Definizione degli indirizzi per la lettura dei dati su ttree
  b1->SetAddress(&point.X);
 //b2->SetAddress(&hits1);
  b3->SetAddress(&hits2);
  b4->SetAddress(&hits3);

  
  TClonesArray *hitsgood2 = new TClonesArray("Punto2",100);
  TClonesArray &ptrhitsgood2 = *hitsgood2;
  TClonesArray *hitsgood3 = new TClonesArray("Punto2",100);
  TClonesArray &ptrhitsgood3 = *hitsgood3;

	
  numeroeventi = tree2->GetEntries();		
  // arrays contenenti Z ricostruite e simulate
  float Zvec[numeroeventi];
  float Zsim[numeroeventi];


  
  // loop sugli ingressi nel TTree
  for(int ev=0;ev<tree2->GetEntries();ev++){
  
    tree2->GetEvent(ev);
    cout<<"Evento "<<ev<<"; Molteplicita= "<<point.mult<<endl;
   // cout<<"X,Y,Z = "<<point.X<<"; "<<point.Y<<"; "<<point.Z<<endl;
    
   Zsim[ev] = point.Z;		// riempio il vector coi dati simulati

         //int numerocoppie = 0;
        
			for(int a = 0; a < point.mult; a++){
			     Punto2 *tst2;
			     tst2=(Punto2*)hits2->At(a);
				 
				 //hits2->Clear();
				 //hits3->Clear();
				for(int b = 0; b < point.mult; b++){
				Punto2 *tst3;
				tst3=(Punto2*)hits3->At(b);
				
				//cout<<tst2->GetPhi()<<"        "<<tst3->GetPhi()<<endl;

					float deltaPhi = tst2->GetPhi() - tst3->GetPhi();
					if(abs(deltaPhi) < 0.01){
						
						
						new(ptrhitsgood2[a])Punto2(tst2->GetPhi(), tst2->GetZ(), tst2->Getnum());
						new(ptrhitsgood3[b])Punto2(tst3->GetPhi(), tst3->GetZ(), tst3->Getnum());
						
					}
					else continue;
					
				}
				
			}
			
				
			// __________tentativo ricostruzione vertice__________
			// eq retta passante per due punti A(x2,y2,z2) e B(x3,y3,z3)
			// x0 = x3 + c1t		c1 = sin(theta)*cos(phi)
			// y0 = y3 + c2t		c2 = sin(theta)*sin(phi)
			// z0 = z3 + c3t		c3 = cos(phi)
			// noi vogliamo solo la z
			// impongo y=0 e x=0, trovo quindi due valori di t e li medio tra loro (x3/c1 e y3/c2)
			
				cout<<ptrhitsgood2.GetEntries()<<endl;
			for(int e=0; e<ptrhitsgood2.GetEntries(); e++)
			{
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
				
								
				//cout<<" label 2 "<<tstgood2->Getnum()<<"  "<<"label 3 "<<tstgood3->Getnum()<<endl;
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
				 
				 //l'intersezione con la beam line sarà il mio Zrec
				 Zrec = Z3 - C[2]*t;
				 
				 hist2 -> Fill(Zrec);
				
		}

			
			Zvec[ev] = hist2->GetMean(); //salvo la media di tutto il mio istogramma, di tutte le Zrec (anche quelle derivate da coppie sbagliate) nel vector Zvec (solo una per evento)	
                
		
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
		 Zrec1 = Zvec[u];
		 
		 differenza = Zsim1 - Zrec1;
		 
		 hist1 -> Fill(differenza); //grafico la differenza
		 }
		 
		
		 
		hist1->Draw();
		TFile file1("ztrue-zrec-tree.root", "recreate");
		hist1->Write();
		file1.Close();
		
		
  double TT = time.CpuTime();	
  cout<<"Il tempo impiegato dalla CPU è "<<TT<<" s"<<endl;
  
 }
 
 
