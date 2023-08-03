//Progetto esame - Simulazione Montecarlo di una collisione protone protone

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

#include "Punto2.h"
#include "Vertex.h"
#include "TracciaMC.h"
#include "Multis.h"
#include "Smearing.h"

using namespace std;

void Simulazione2(bool multipleScatt){ // 0 = false, 1 = true
	
        TStopwatch time;	
	 
        //creo Tree dove intendo salvare i miei dati
	TFile hfile2("htree2.root","RECREATE");	  
        TTree *tree2 = new TTree("T2","TTree con 3 branches");
        
        int numeroeventi = 1000;
        
        TClonesArray *ptrhits1 = new TClonesArray("Punto2",100);
        TClonesArray &hits1 = *ptrhits1;
        
        TClonesArray *ptrhits2 = new TClonesArray("Punto2",100);
        TClonesArray &hits2 = *ptrhits2;
        
        TClonesArray *ptrhits3 = new TClonesArray("Punto2",100);
        TClonesArray &hits3 = *ptrhits3;
        
        //definisco una struct che contiene l'informazione sulla posizione e la molteplicità del vertice
		 
        typedef struct{
               float X,Y,Z;
               int mult;} VTX;
        static VTX point;
		 
        //creo i branches del tree
        
        tree2->Branch("VertMult",&point.X,"X/F:Y:Z:mult/I");
        tree2->Branch("HitsPrimo",&ptrhits1);
        tree2->Branch("HitsSecondo",&ptrhits2);
        tree2->Branch("HitsTerzo",&ptrhits3);
                   
        //creo gli oggetti da usare per ogni evento della simulazione         
        
        std::array<float,2> intersez; //array contenente le coordinate delle intersezioni
                  
        Vertex vertice0; //vertice
        Multis catter; //multiple scattering
        
        int numSpuri = 5; //numero di eventi spuri simulati per ogni evento
        int numParticella; // label per identificare le tracce simulate
        int count=0;
              
	for(int evento = 0; evento < numeroeventi; evento++){
	         //if(evento%50000==0) cout<<"siamo arrivati al numero "<<evento<<endl; //controllo su come procede la simulazione
		 vertice0.NewVertex(); //estraggo nuove coordinate casuali del vertice
		 	
		 //salvo informazioni nell'oggetto point	 
		 point.mult = vertice0.GetMolteplicity();
                 point.X = vertice0.GetX();
                 point.Y = vertice0.GetY();
                 point.Z = vertice0.GetZ();
		       		 		 
		 //calcolo intersezioni su ogni layer		
		 numParticella = 0;		
		 		 
		 do{		 		
			 TracciaMC tr(0.,0.,point.X,point.Y,point.Z); //creo oggetto traccia inserendo le coordinate del Vertice
			 
			 //imposto la direzione iniziale della particella, estratta casualmente
			 tr.SetDistribEta(); 
			 tr.SetDistribPhi();
			 tr.Theta();			 
			 tr.CalcCoeff();
			 
			 //controllo che la particella non esca dal terzo layer del rivelatore
			/* intersez = tr.intersezione(3);
			 if(intersez[1] < -27/2 || intersez[1] > 27/2) {
			 count++;
			 cout<<"andata fuori numero "<<count<<endl;
			 continue;
			 }
			 */
			 
			 for(int t=0; t<3; t++){
			 
		          intersez = tr.intersezione(t+1); //calcolo intersezione con il layer
		          
		          if((intersez[1] >= -27/2) && (intersez[1] <= 27/2)){									
			    
			    catter.NuoviAngoli(multipleScatt); //estraggo nuovi angoli del multiple scattering casuali (se multipleScatt = 1)
		            catter.VarioAngolo(tr); //conseguentemente cambiano gli angoli della direzione della particella 
			    tr.SetOrigine(intersez, t+1); //l'intersezione è un punto della nuova retta descritta dalla particella
				
			    Smearing ringo; //estratti gaussianamente un dZ e dPhi che si andranno a sommare alle variabili vere
			    ringo.smearZ(intersez[1]);
                            ringo.smearPhi(intersez[0], t+1);
                                
                            //inserisco le hits dentro ai branches del tree   				
			    if(t==0) new(hits1[numParticella])Punto2(ringo.GetPhirec(),ringo.GetZrec(), numParticella);
			    if(t==1) new(hits2[numParticella])Punto2(ringo.GetPhirec(),ringo.GetZrec(), numParticella);
			    if(t==2) new(hits3[numParticella])Punto2(ringo.GetPhirec(),ringo.GetZrec(), numParticella);
			    }
			  
			    else{
			    if(t==0) new(hits1[numParticella])Punto2(intersez[0],intersez[1], -1);
			    if(t==1) new(hits2[numParticella])Punto2(intersez[0],intersez[1], -1);
			    if(t==2) new(hits3[numParticella])Punto2(intersez[0],intersez[1], -1);
			    }
				
			  }			
		
			  numParticella++;
			 
		   }while(numParticella < point.mult);
               
                 // aggiungo ora numSpuri punti spuri su ciascun layer
		 for(int u=0; u<numSpuri; u++){
		    new(hits1[point.mult + u])Punto2(gRandom->Rndm()*2*TMath::Pi(),gRandom->Rndm()*27 - (27./2),-20);
		    new(hits2[point.mult + u])Punto2(gRandom->Rndm()*2*TMath::Pi(),gRandom->Rndm()*27 - (27./2),-20);
		    new(hits3[point.mult + u])Punto2(gRandom->Rndm()*2*TMath::Pi(),gRandom->Rndm()*27 - (27./2),-20);
		    }
		    
		 tree2->Fill(); //riempie il tree
		 
		 //ripristino i TClonesArrays
                 ptrhits1->Clear();
                 ptrhits2->Clear();
                 ptrhits3->Clear();
		 }
				 
		 //salvo il tree sul file e lo chiudo
                 hfile2.Write();              
                 hfile2.Close();
                 
                 double TT = time.CpuTime();	
                 cout<<"Il tempo impiegato dalla CPU è "<<TT<<" s"<<endl;
	
 }
