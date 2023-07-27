//Progetto esame - Simulazione Montecarlo di una collisione protone protone

#include <iostream>
#include <fstream>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "Punto.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "Vertex.h"
#include "TracciaMC.h"
#include "Multis.h"
#include "TRandom3.h"

#include <vector>

using namespace std;

void Simulazione(){
	
  TStopwatch time;
	
 bool multipleScatt = 1;		// 0 = false, 1 = true
 
  // ______________grafico 3d intersezioni_______________
		
	TFile hfile("htree.root","RECREATE");	  
        TTree *tree = new TTree("T","TTree con 2 branches");
        
        int numeroeventi = 1000;
        
        TClonesArray *ptrhits1 = new TClonesArray("Punto",numeroeventi);
        TClonesArray &hits1 = *ptrhits1;
        
        TClonesArray *ptrhits2 = new TClonesArray("Punto",numeroeventi);
        TClonesArray &hits2 = *ptrhits2;
        
        TClonesArray *ptrhits3 = new TClonesArray("Punto",numeroeventi);
        TClonesArray &hits3 = *ptrhits3;
        
      
        
        ofstream fileout1("Zsim.txt"); //forse questo diventa inutile perchè ho già una foglia del tree con la Z del vertice
		 
		   typedef struct{
                   float X,Y,Z;
                   int mult;} VTX;
    
    
                   static VTX point;
		 
		   tree->Branch("VertMult",&point.X,"X/F:Y:Z:mult/I");
                   tree->Branch("HitsPrimo",&ptrhits1,"HitsPrimo/F");
                   tree->Branch("HitsSecondo",&ptrhits2,"HitsSecondo/F");
                   tree->Branch("HitsTerzo",&ptrhits3,"HitsTerzo/F");
                   
                   vector<float> coordinate, intersez;
                   Vertex vertice0;
                   
              
	for(int evento = 0; evento < numeroeventi; evento++){
	

		    vertice0.NewVertex();
		 		 
		     point.mult = vertice0.GetMolteplicity();
                     point.X = vertice0.GetX();
                     point.Y = vertice0.GetY();
                     point.Z = vertice0.GetZ();

		 
		
		 fileout1 <<point.Z<<endl ;
		 
		 
		  // __________creo oggetto Multis per il Multiple Scattering____________
		  
		  Multis catter;
		 
		 // __________calcolo intersezioni su ogni layer____________
		
		 int i = 0;
		 do{
		 
		
			 TracciaMC tr(0.,0.,point.X,point.Y,point.Z);
			 
			 tr.SetDistribEta(); 
			 tr.SetDistribPhi();
			 tr.Theta();
			 
			 tr.CalcCoeff();
			 
			 intersez = tr.intersezione(3);
			 if(intersez[2] < -27/2 || intersez[2] > 27/2) continue;
			 
			 
			 
			 for(int t=0; t<3; t++){
			 
				intersez = tr.intersezione(t+1);										
				catter.NuoviAngoli(multipleScatt);
				catter.VarioAngolo(tr);
				tr.SetOrigine(intersez);
				
				if(t==0) new(hits1[i])Punto(intersez[0],intersez[1],intersez[2]);
				if(t==1) new(hits2[i])Punto(intersez[0],intersez[1],intersez[2]);
				if(t==2) new(hits3[i])Punto(intersez[0],intersez[1],intersez[2]);
				}			
			
			 i++;
			 
		 }while( i < point.mult);
               
               /*  // Debug
    printf("Evento %d - moltepl: %d\n",evento,point.mult);
    printf("x= %f ; y= %f; z= %f \n",point.X,point.Y,point.Z);
    printf("Entries nel TClonesArray: %d\n",ptrhits1->GetEntries());
    for (int j=0; j<hits1.GetEntries(); j++){
      Punto *tst=(Punto*)ptrhits1->At(j);
      cout<<"int primo layer "<<j<<") x, y, z = "<<tst->GetX()<<"; "<<tst->GetY()<<"; "<<tst->GetZ()<<endl;}*/
    
    // fine del debug
		 
		 tree->Fill();
                 ptrhits1->Clear();
                 ptrhits2->Clear();
                 ptrhits3->Clear();
		 }
		
		
		 fileout1.close();
		   // Save all objects in this file
                 hfile.Write();

                // Close the file. 
                hfile.Close();
                 
                  double TT = time.CpuTime();	
  cout<<"Il tempo impiegato dalla CPU è "<<TT<<" s"<<endl;
	
 }
