//Progetto esame - Simulazione Montecarlo di una collisione protone protone

#include <iostream>
#include <fstream>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "Punto.h"
#include "Punto2.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "Vertex.h"
#include "TracciaMC.h"
#include "Multis.h"
#include "TRandom3.h"

#include <vector>

using namespace std;

void Simulazione2(bool multipleScatt){
	
  TStopwatch time;
	
 //bool multipleScatt = 1;		// 0 = false, 1 = true
 
  // ______________grafico 3d intersezioni_______________
		
	TFile hfile2("htree2.root","RECREATE");	  
        TTree *tree2 = new TTree("T2","TTree con 2 branches");
        
        int numeroeventi = 1000;
        
        TClonesArray *ptrhits1 = new TClonesArray("Punto2",100);
        TClonesArray &hits1 = *ptrhits1;
        
        TClonesArray *ptrhits2 = new TClonesArray("Punto2",100);
        TClonesArray &hits2 = *ptrhits2;
        
        TClonesArray *ptrhits3 = new TClonesArray("Punto2",100);
        TClonesArray &hits3 = *ptrhits3;
        
      
        
        //ofstream fileout1("Zsim.txt"); //forse questo diventa inutile perchè ho già una foglia del tree con la Z del vertice
		 
		   typedef struct{
                   float X,Y,Z;
                   int mult;} VTX;
    
    
                   static VTX point;
		 
		   tree2->Branch("VertMult",&point.X,"X/F:Y:Z:mult/I");
                   tree2->Branch("HitsPrimo",&ptrhits1);
                   tree2->Branch("HitsSecondo",&ptrhits2);
                   tree2->Branch("HitsTerzo",&ptrhits3);
                   
                  
                   std::array<float,2> intersez;
                  
                   Vertex vertice0;
                   
              
	for(int evento = 0; evento < numeroeventi; evento++){
	

		    vertice0.NewVertex();
		 		 
		     point.mult = vertice0.GetMolteplicity();
                     point.X = vertice0.GetX();
                     point.Y = vertice0.GetY();
                     point.Z = vertice0.GetZ();

		       		 
		  // __________creo oggetto Multis per il Multiple Scattering____________
		  
		  Multis catter;
		 
		 // __________calcolo intersezioni su ogni layer____________
		
		 int i = 0;
		 int numParticella = 0;		// label per identificare le tracce simulate
		 
		 int numSpuri = 5;
		 
		 do{		 		
			 TracciaMC tr(0.,0.,point.X,point.Y,point.Z);
			 
			 tr.SetDistribEta(); 
			 tr.SetDistribPhi();
			 tr.Theta();
			 
			 tr.CalcCoeff();
			 
			 intersez = tr.intersezione(3);
			 if(intersez[1] < -27/2 || intersez[1] > 27/2) continue;
			 
			 
			 for(int t=0; t<3; t++){
			 
				intersez = tr.intersezione(t+1);										
				catter.NuoviAngoli(multipleScatt);
				catter.VarioAngolo(tr);
				tr.SetOrigine(intersez, t+1);
				
				Smearing ringo;
				ringo.smearZ(intersez[1]);
                                ringo.smearPhi(intersez[0], t+1);
				
				if(t==0) new(hits1[i])Punto2(ringo.GetPhirec(),ringo.GetZrec(), numParticella);
				if(t==1) new(hits2[i])Punto2(ringo.GetPhirec(),ringo.GetZrec(), numParticella);
				if(t==2) new(hits3[i])Punto2(ringo.GetPhirec(),ringo.GetZrec(), numParticella);
				
				}			
			
			 i++;
			 numParticella++;
			 
		 }while( i < point.mult);
               
                // aggiungo ora numSpuri punti spuri su ciascun layer
		 for(int u=0; u<numSpuri; u++){
		    new(hits1[point.mult + u])Punto2(gRandom->Rndm()*2*TMath::Pi(),gRandom->Rndm()*27 - (27./2),-20);
		    new(hits2[point.mult + u])Punto2(gRandom->Rndm()*2*TMath::Pi(),gRandom->Rndm()*27 - (27./2),-20);
		    new(hits3[point.mult + u])Punto2(gRandom->Rndm()*2*TMath::Pi(),gRandom->Rndm()*27 - (27./2),-20);
		 }
		 tree2->Fill();
                 ptrhits1->Clear();
                 ptrhits2->Clear();
                 ptrhits3->Clear();
		 }
		
		
		 
		 // Save all objects in this file
                 hfile2.Write();

                 // Close the file. 
                 hfile2.Close();
                 
                 
                 
                  double TT = time.CpuTime();	
  cout<<"Il tempo impiegato dalla CPU è "<<TT<<" s"<<endl;
	
 }
