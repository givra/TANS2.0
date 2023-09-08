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

void Simulazione2(bool multipleScatt, bool molteplicity, bool eta){ // 0 = false, 1 = true
	
        TStopwatch time;

        //creo Tree dove intendo salvare i miei dati
		TFile hfile2("htree2.root","RECREATE");	
        TTree *tree2 = new TTree("T2","TTree con 3 branches");
		// tree2->SetDirectory(&hfile2);
		
		// se provo a salvare tutto sul tree2 mi dà segmentation violation :)
		TFile hfile3("htree3.root","RECREATE");			
		TTree *tree3 = new TTree("T3","TTree dati efficienza");
		// tree3->SetDirectory(&hfile3);
        
        int numeroeventi = 1000000;
        
        TClonesArray *ptrhits1 = new TClonesArray("Punto",100);
        TClonesArray &hits1 = *ptrhits1;
        
        TClonesArray *ptrhits2 = new TClonesArray("Punto",100);
        TClonesArray &hits2 = *ptrhits2;
        
        TClonesArray *ptrhits3 = new TClonesArray("Punto",100);
        TClonesArray &hits3 = *ptrhits3;
		
		typedef struct{
			float x0,x1,x2,x3,x4,x5,x6,x7,x8,x9;		// Z simulate con un certo valore di moltep, sarà denominatore dell'efficienza
			} vett;
		static vett denEff;
		
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
		
		tree3->Branch("denEff",&denEff.x0, "x0:x1:x2:x3:x4:x5:x6:x7:x8:x9/F");
                   
        //creo gli oggetti da usare per ogni evento della simulazione         
        
        std::array<float,2> intersez; //array contenente le coordinate delle intersezioni
                  
        Vertex vertice0; //vertice
        Multis catter; //multiple scattering
        Traccia tr; //traccia
        
        //Scelta della distribuzione della molteplicità
        TH1F* hist1;
        //bool molteplicity = 1; //distribuzione di eta 0->uniforme 1->kinem.root
        
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
        //bool eta = 1; //distribuzione di eta 0->uniforme 1->kinem.root
        
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
        
        int numSpuri = 5; //numero di eventi spuri simulati per ogni evento
        int numParticella; // label per identificare le tracce simulate
        int count=0;
        
        float u1, u2;
        int fmax1 = 27201;
        int fmax2 = 15344;
        float Xt, Yt, f;
		
		float count3 = 0;
		float count5 = 0;
		float count6 = 0;
		float count7 = 0;
		float count8 = 0;
		float count12 = 0;
		float count22 = 0;
		float count32 = 0;
		float count42 = 0;
		float count52 = 0;
		
		int percento = 0;
              
	for(int evento = 0; evento < numeroeventi; evento++){
	         if((evento+1)%(numeroeventi/20)==0){
	         percento = percento + 5;
	         cout<<percento<<" %"<<endl;} //controllo su come procede la simulazione
	        
		 vertice0.NewVertex(); //estraggo nuove coordinate casuali del vertice
		 //calcolo la molteplicità del vertice in base alla mia scelta	
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
       
                         vertice0.SetMolt(((67.9821)*(Xt-1)/(141))+2.5);  //-((4.08)/(2*34)));
                         }
		 //salvo informazioni nell'oggetto point	 
		 point.mult = vertice0.GetMolteplicity();
		 
                 point.X = vertice0.GetX();
                 point.Y = vertice0.GetY();
                 point.Z = vertice0.GetZ();
		       		 		 
		 // conto quante Z ho per ogni valore di molteplicità
		 switch(point.mult)
		 {case 3:
				count3++;
				//cout << " n evento " << evento << " count3 " << count3 << endl;
				break;
		  case 5:
				count5++;
				//cout << " n evento " << evento << " count5 " << count5 << endl;
				break;
		  case 6:
				count6++;
				//cout << " n evento " << evento << " count6 " << count6 << endl;
				break;
		  case 7:
				count7++;
				//cout << " n evento " << evento << " count7 " << count7 << endl;
		        break;
		  case 8:
				count8++;
				//cout << " n evento " << evento << " count8 " << count8 << endl;
		        break;
		  case 12:
				count12++;
				//cout << " n evento " << evento << " count12 " << count12 << endl;
		        break;
		  case 22:
				count22++;
				//cout << " n evento " << evento << " count22 " << count22 << endl;
		        break;
		  case 32:
				count32++;
				//cout << " n evento " << evento << " count32 " << count32 << endl;
		        break;
		  case 42:
				count42++;
				//cout << " n evento " << evento << " count42 " << count42 << endl;
		        break;
		  case 52:
				count52++;
				//cout << " n evento " << evento << " count52 " << count52 << endl;
		        break;
			 
		 }
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
       
                         tr.SetEta(((4.08)*(Xt-1)/(34))-2.04); //-((4.08)/(2*34)));
                         }
                        
			 //imposto la direzione iniziale della particella, estratta casualmente			 
			 tr.SetPhi();
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
		    //new(hits1[point.mult + u])Punto(gRandom->Rndm()*2*TMath::Pi(),gRandom->Rndm()*27 - (27./2),-20);
		    new(hits2[point.mult + u])Punto(gRandom->Rndm()*2*TMath::Pi(),gRandom->Rndm()*27 - (27./2),-20);
		    new(hits3[point.mult + u])Punto(gRandom->Rndm()*2*TMath::Pi(),gRandom->Rndm()*27 - (27./2),-20);
		    }
		    
		 tree2->Fill(); //riempie il tree
		 
		 //ripristino i TClonesArrays
                 ptrhits1->Clear();
                 ptrhits2->Clear();
                 ptrhits3->Clear();
		 }		// fine ciclo eventi
		 
		 hfile2.Write();              
         hfile2.Close();
		 
		 denEff.x0 = count3;
		 denEff.x1 = count5;
		 denEff.x2 = count6;
		 denEff.x3 = count7;
		 denEff.x4 = count8;
		 denEff.x5 = count12;
		 denEff.x6 = count22;
		 denEff.x7 = count32;
		 denEff.x8 = count42;
		 denEff.x9 = count52;
		 tree3->Fill();
		 

			// cout << denEff.x0 << " " ;
		 
		
		 		 
		 //salvo il tree sul file e lo chiudo
                 
                 hfile3.Write();
		 hfile3.Close();
		 
		 delete hist1;
		 delete hist2;
				 
                 double TT = time.CpuTime();	
                 cout<<endl<<"Il tempo impiegato dalla CPU è "<<TT<<" s"<<endl;
	
 }
 
 

  
  
  
