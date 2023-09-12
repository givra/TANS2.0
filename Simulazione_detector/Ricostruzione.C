#include <iostream>
#include <fstream>
#include <math.h>

#include "Punto.h"
#include "TClonesArray.h"
#include "TRandom3.h"

#include <vector>

using namespace std;

void Ricostruzione(){
	
	int numeroeventi;
	
	TStopwatch time;
		
        TH1D* hist1 = new TH1D ("histofin","Zrec-Ztrue" , 250, -0.05, 0.05);
	TH1D* hist2 = new TH1D ("h2","Zreconstructed Tracklets" , 54, -13.5, 13.5);
   
	float deltaPhi;
	
        float Z3, Z2, TanTheta;
			
	float q;      //termine noto retta
	int Ymax;			
				
	static float R2 = 4.;
	static float R3 = 7.;
			
	float differenza = 0;
			
        //definizione struct
        typedef struct {
           float X,Y,Z;
           int mult;} VTX1;
        static VTX1 point;
		
		typedef struct {
        float Z;} VTX;
		static VTX pointRec;
  
        //apertura file di input
        TFile hfile2("htree2.root");
  
        //lettura TTree  e branch
        TTree *tree2 = (TTree*)hfile2.Get("T2");
        TBranch *b1=tree2->GetBranch("VertMult");     
        TBranch *b3=tree2->GetBranch("HitsSecondo");
        TBranch *b4=tree2->GetBranch("HitsTerzo");
		
	TFile hfile4("htree4.root", "RECREATE");
	TTree *tree4 = new TTree("T4","TTree dati efficienza");
	tree4->Branch("Zrec",&pointRec.Z,"Z/I");
		
        numeroeventi = tree2->GetEntries(); //acquisico informazione sul numero di eventi nel mio detector
        	      		
        //dichiarazione TClonesArray
        TClonesArray *hits2 = new TClonesArray("Punto",100);
        TClonesArray *hits3 = new TClonesArray("Punto",100);
  
        //definizione degli indirizzi per la lettura dei dati su ttree
        b1->SetAddress(&point.X);
        b3->SetAddress(&hits2);
        b4->SetAddress(&hits3);
  
        // salvo le buone combinazioni di hit nei vectors di float
        vector<float> Zgood2;
        vector<float> Zgood3;
	float percento = 0;
		
	// loop sugli ingressi nel TTree
        for(int ev=0;ev<numeroeventi;ev++){
            if((ev+1)%(numeroeventi/20)==0){
				percento = percento + 5;
				cout<< percento << "%" << endl; //controllo su come procede la simulazione
			        }
				
             tree2->GetEvent(ev); 
                              
             for(int a = 0; a < hits2->GetEntries(); a++){		// loop su hits layer 2
			     
					Punto *tst2;
					tst2=(Punto*)hits2->At(a);
						
					if(tst2->Getnum()!=-1){			// escludiamo hits fuori dal rivelatore
						for(int b = 0; b < hits2->GetEntries(); b++){		// loop su hits layer 3
				  
						Punto *tst3;
						tst3=(Punto*)hits3->At(b);
						if(tst3->Getnum()!=-1){
							deltaPhi = tst2->GetPhi() - tst3->GetPhi();
							
								if(abs(deltaPhi) < 0.01){
									Zgood2.push_back(tst2->GetZ());
									Zgood3.push_back(tst3->GetZ());
								}
							}
							tst3->Clear();
						}
			       }
				   tst2->Clear();
				   
			   }
			

	    // __________ricostruzione vertice__________
			
		for(int e=0; e<Zgood2.size(); e++){		// ciclo sui punti buoni
		
		
				TanTheta = (R3-R2)/(Zgood3[e]-Zgood2[e]);
		
				q = R3 - Zgood3[e]*TanTheta; //termine noto retta
											
				hist2 -> Fill(-q/TanTheta);		// riempio l'histo con gli Zrec dati dalle coppie "buone"
		
			}
			 
	        Ymax = hist2->GetBinContent(hist2->GetMaximumBin()); //valore del bin più alto
					       				
        	// sotto le 2 tracce non ricostruiamo -> associamo Z=100 per poi scartare il dato
        	if(Ymax > 2){
				//if(hist2->GetMean() < 5.3*3 && hist2->GetMean() > 5.3*3){		// controllo su 3sigma
				pointRec.Z = hist2->GetMean();
				tree4->Fill();
				
				differenza = pointRec.Z	- point.Z;		
			        hist1->Fill(differenza);
				
			}
        	else{pointRec.Z = 100.;
				tree4->Fill();
			}		

               
                hist2->Reset("ICESM");
                
                Zgood2.clear();
		Zgood3.clear();
									
                hits2->Clear();
                hits3->Clear();
				 
				 
		}		// fine ciclo eventi
				
		
		
	 delete hist2;
	 hfile2.Close();
						 	
	 TFile file1("zrec-ztrue.root", "recreate");
	 hist1->Write();	
	 delete hist1;
	 file1.Close();
	
	 hfile4.Write();              
         hfile4.Close();
				
         double TT = time.CpuTime();			
         cout<< endl <<"Il tempo impiegato dalla CPU è "<<TT<<" s"<<endl;
         
         MemInfo_t memInfo;
         gSystem->GetMemInfo(&memInfo);
         cout << "Mem Used = " << memInfo.fMemUsed << " MB"<<endl; //returning value in MB
         }
