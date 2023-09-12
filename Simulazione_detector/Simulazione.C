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
#include "TH1F.h"
#include "TAxis.h"

#include "Punto.h"
#include "Vertex.h"
#include "Traccia.h"
#include "Multis.h"
#include "Smearing.h"

using namespace std;

void Simulazione(bool multipleScatt, bool molteplicity, bool eta){ // 0 = false, 1 = true
	
        TStopwatch time;

        //creo i Tree 
		TFile hfile("htree.root","RECREATE");	
        TTree *tree = new TTree("T","TTree dati simulati");
        
        int numeroeventi = 10000;
        
        TClonesArray *ptrhits1 = new TClonesArray("Punto",100);
        TClonesArray &hits1 = *ptrhits1;
        
        TClonesArray *ptrhits2 = new TClonesArray("Punto",100);
        TClonesArray &hits2 = *ptrhits2;
        
        TClonesArray *ptrhits3 = new TClonesArray("Punto",100);
        TClonesArray &hits3 = *ptrhits3;
		
        //definisco una struct che contiene l'informazione sulla posizione e la molteplicità del vertice
		 
        typedef struct{
               float X,Y,Z;
               int mult;} VTX;
        static VTX point;
		 
        //creo i branches del tree
        
        tree->Branch("VertMult",&point.X,"X/F:Y:Z:mult/I");
        tree->Branch("HitsPrimo",&ptrhits1);
        tree->Branch("HitsSecondo",&ptrhits2);
        tree->Branch("HitsTerzo",&ptrhits3);
	              
        //creo gli oggetti da usare per ogni evento della simulazione         
        
        std::array<float,2> intersez; //array contenente le coordinate delle intersezioni
                  
        Vertex vertice0; //vertice
        Multis catter; //multiple scattering
        Traccia tr; //traccia
        
		//Scelta della distribuzione della molteplicità
        TH1F* hist1;
		 //distribuzione di molteplicita' 0->uniforme 1->kinem.root
        
		if(molteplicity==1){
               TFile F("kinem.root");
               TH1F *disteta = (TH1F*)F.Get("hm");
               disteta->SetDirectory(0);
               disteta->SetMinimum(0);
               F.Close();
               TAxis *xa=disteta->GetXaxis();
               float step = xa->GetBinWidth(1);
               int b1=xa->FindBin(2.6);
               int b2=xa->FindBin(70.4821);
 
               float xlow=xa->GetBinLowEdge(b1);
               float xhig=xa->GetBinUpEdge(b2);
               int nobins=b2-b1+1;
               float step2 = (xhig-xlow)/nobins;
  
               hist1 = new TH1F("hist1","molteplicity distribution",nobins,xlow,xhig);
               int j=1;
			   for(int i=b1;i<=b2+1;i++)hist1->SetBinContent(j++,disteta->GetBinContent(i));
               
          }
		  
		 //Scelta della distribuzione di Eta
        TH1F* hist2;
        //distribuzione di eta 0->uniforme 1->kinem.root
        
		if(eta==1){
               TFile F("kinem.root");
               TH1F *disteta = (TH1F*)F.Get("heta2");
               disteta->SetDirectory(0);
               disteta->SetMinimum(0);
               F.Close();
               TAxis *xa=disteta->GetXaxis();
               float step = xa->GetBinWidth(1);
               int b1=xa->FindBin(-2.);
               int b2=xa->FindBin(2.);
 
               float xlow=xa->GetBinLowEdge(b1);
               float xhig=xa->GetBinUpEdge(b2);
               int nobins=b2-b1+1;
               float step2 = (xhig-xlow)/nobins;
  
               hist2 = new TH1F("hist2","#eta distribution 2",nobins,xlow,xhig);
               int j=1;
               for(int i=b1;i<=b2+1;i++)hist2->SetBinContent(j++,disteta->GetBinContent(i));
               
          }
        
        int numSpuri = 0; //numero di eventi spuri simulati per ogni evento
        int numParticella; // label per identificare le tracce simulate
        
        float u1, u2;
        int fmax1 = 27201;
        int fmax2 = 15344;
        float Xt, Yt, f;
		float percento = 0;
             	
	for(int evento = 0; evento < numeroeventi; evento++){
	     if((evento+1)%(numeroeventi/20)==0){
			 percento = percento + 5;
			 cout<< percento << "%" << endl; //controllo su come procede la simulazione
		 }
		 
		 vertice0.NewVertex(); //estraggo nuove coordinate casuali del vertice
		 numSpuri = (gRandom->Rndm()*3)/1;
		 
		 if(molteplicity==0) vertice0.SetMoltUniform();
		 
		 else {//distribuzione data del file kinem.root
        
                //il numero di bins è 141
                //il range va da 2.6 a 70.4821
                //Xmax è 1 e Ymax è 27201
         
                do{  
                    u1 = gRandom->Rndm();
                    u2 = gRandom->Rndm();
  
                    Yt = fmax1*u2;
                    Xt = (141*u1)/1;
                    f = hist1->GetBinContent(Xt);}
                while(f<=Yt);
       
                vertice0.SetMolt(((67.9821)*(Xt-1)/(141))+2.5);  
         }
		 
		 
		 //salvo informazioni nell'oggetto point	 
		 point.mult = vertice0.GetMolteplicity();
		 
                 point.X = vertice0.GetX();
                 point.Y = vertice0.GetY();
                 point.Z = vertice0.GetZ();
		       		 		 
		 //calcolo intersezioni su ogni layer		
		 numParticella = 0;		
		 		 	
		 do{		 		
			 tr.SetOrigine(point.X,point.Y,point.Z); //setto come origine della traccia le coordinate del Vertice
			 
			 if(eta==0) tr.SetEtaUni();
			 
			 else {//distribuzione data del file kinem.root
        
                         //il numero di bins è 34
                         //il range va da -2.04 a 2.04
                         //Xmax è 33 e Ymax è 15344
         
                         do{  
                             u1 = gRandom->Rndm();
                             u2 = gRandom->Rndm();
  
                             Yt = fmax2*u2;
                             Xt = (34*u1)/1;
                             f = hist2->GetBinContent(Xt);}
                         while(f<=Yt);
       
                         tr.SetEta(((4.08)*(Xt-1)/(34))-2.04);
                         }
                        
			 //imposto la direzione iniziale della particella, estratta casualmente			 
			 tr.SetPhi();
			 tr.Theta();			 
			 tr.CalcCoeff();
			 
			 
			 for(int t=0; t<3; t++){
			 
		          intersez = tr.intersezione(t+1); //calcolo intersezione con il layer
		          
		        if((intersez[1] >= -27/2) && (intersez[1] <= 27/2)){									
			    
					catter.NuoviAngoli(multipleScatt); // estraggo nuovi angoli del multiple scattering casuali (se multipleScatt = 1)
		            catter.VarioAngolo(tr); // conseguentemente cambiano gli angoli della direzione della particella 
					tr.SetHit(intersez, t+1); // l'intersezione è un punto della nuova retta descritta dalla particella
				
					Smearing ringo; //estratti gaussianamente un dZ e dPhi che si andranno a sommare alle variabili vere
					ringo.smearZ(intersez[1]);
					ringo.smearPhi(intersez[0], t+1);
									
								//inserisco le hits dentro ai branches del tree   				
					if(t==0) new(hits1[numParticella])Punto(ringo.GetPhirec(),ringo.GetZrec(), numParticella);
					if(t==1) new(hits2[numParticella])Punto(ringo.GetPhirec(),ringo.GetZrec(), numParticella);
					if(t==2) new(hits3[numParticella])Punto(ringo.GetPhirec(),ringo.GetZrec(), numParticella);
			    }
				// se sta fuori dal rivelatore assegno label -1
			    else{
					if(t==0) new(hits1[numParticella])Punto(intersez[0],intersez[1], -1);
					if(t==1) new(hits2[numParticella])Punto(intersez[0],intersez[1], -1);
					if(t==2) new(hits3[numParticella])Punto(intersez[0],intersez[1], -1);
			    }
				
			  }			
		
			  numParticella++;
			 
		   }while(numParticella < point.mult);
               
                 // aggiungo ora numSpuri punti spuri su ciascun layer
		 for(int u=0; u<numSpuri; u++){
		    new(hits2[point.mult + u])Punto(gRandom->Rndm()*2*TMath::Pi(),gRandom->Rndm()*27 - (27./2),-20);
		    new(hits3[point.mult + u])Punto(gRandom->Rndm()*2*TMath::Pi(),gRandom->Rndm()*27 - (27./2),-20);
		    }
		    
		 tree->Fill(); //riempie il tree
		 	
                 
		 //ripristino i TClonesArrays
                 ptrhits1->Clear();
                 ptrhits2->Clear();
                 ptrhits3->Clear();
		 }		// fine ciclo eventi
		 
		 hfile.Write();              
         hfile.Close();		
		 
		 if(molteplicity==1) delete hist1;
		 if(eta==1) delete hist2;
		 
         double TT = time.CpuTime();	
         cout<< endl <<"Il tempo impiegato dalla CPU è "<<TT<<" s"<<endl;
         
         MemInfo_t memInfo;
         gSystem->GetMemInfo(&memInfo);
cout << "Mem Used = " << memInfo.fMemUsed << " MB"<<endl; //returning value in MB
	
 }
 
 

  
  
  
 

  
  
  
