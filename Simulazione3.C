//Progetto esame - Simulazione Montecarlo di una collisione protone protone

#include <iostream>
#include <fstream>
#include <math.h>

#include "Vertex.h"
#include "TracciaMC.h"
#include "Multis.h"
#include "TRandom3.h"

#include <vector>
#include <chrono>

using namespace std;

void Simulazione3(){
	
//Inserimento del seed manuale

 /*unsigned int seed;
	
   cout << " inserire un seed: " << endl;
   cin >> seed;	
 gRandom->SetSeed(seed);*/
 
 TStopwatch time; //Fa partire il tempo per cronometrare il CPU time
 std::chrono::steady_clock::time_point tcount_start, tcount_end;
	
//Chiede se attivare o meno il multiple scattering ---- non lo chiede di volta in volta, come mai??
 bool multipleScatt = 0;	// 0 = false, 1 = true
 cout << " attivare multiple scattering? " << endl;
 cin >> multipleScatt;
 
        
		
	ofstream fileout0("Zsim.txt"); //andremo a segnarci la Z simulata (la nostra Ztrue) del vertice dell'evento (solo una per evento)
	ofstream fileout1("Zrec.txt"); //ci segnamo la Z ricostruita per poi poter fare un confronto con la simulazione		 

        int numeroeventi = 100; //decidiamo il numero di eventi da simulare
        
        //conviene inizializzare le variabili che devo usare, una volta sola, fuori dal for!
        
        vector<float> coordinate;
        int molteplicità;
        
        int i;
        int ii;
        int iii;
        int j;
        
        vector<float> intersez;
        
        float Xtrue;
	float Ytrue;
	float Ztrue;
	
	int l;
	int m;
	
	int numerocoppie;
	
	float X2;
	float Y2;   //coordinate intersezione secondo layer
	float Z2;
			
	float X3;
	float Y3;   //coordinate intersezione terzo layer
	float Z3;
				
	float C[3];   //coefficienti retta
	float t;      //parametro retta
				
	float R2 = 4.;
	float R3 = 7.;
			
	float Phi;
	float Theta;
			
	float Zrec;
	float Xrec;
	float Yrec;
			
	float Zsim;
	
	float differenza;
	float Zsim1;
	float Zrec1;
	
	float mat[4][81][2];

        double tim;
        double tot = 0;
        
        TH1D* hist = new TH1D ("h1","zreconstructed" , 54, -13.5, 13.5);
	
	for(int evento = 0; evento < numeroeventi; evento++){ //ciclo sul numero di eventi
	
	tcount_start = std::chrono::steady_clock::now();

//File relativi al singolo evento, andranno cancellati e sovrascritti ad ogni ciclo di for
	
ofstream fileout("intersez_simulate.txt");
ifstream fin("intersez_simulate.txt");  //file dove scrivo le intersezioni coi 3 layer
ifstream fin1("intersez_simulate.txt");	//ho un secondo ifstream per trovare la molteplicità	
        
//NON SERVE ? ofstream fout("intersez_rec.txt"); //intersezioni misurate coi 3 layer (dopo lo smearing)
          
ofstream fout1("dati_tracklets.txt"); //ci scrivo le coppie di intersezioni 2° e 3° layer che stanno dentro a un angolo deltaPhi 0.01 rad   
ifstream fin2("dati_tracklets.txt"); 

  	
//TIME	@100 eventi 38.85 ms
	
	//SIMULAZIONE
		 
		 Vertex vertice0; //creo vertice con le sue 3 coordinate
		 vertice0.NewVertex();
		
		 coordinate = vertice0.GetCoordinate();
		 molteplicità = vertice0.GetMolteplicity();
		 		 		
		 fileout0 << coordinate[2] << endl ; //mi salvo su Zsim la Z del vertice
		 
		  // __________creo oggetto Multis per il Multiple Scattering____________
		  
		  Multis catter;
		 
		 // __________calcolo intersezioni su ogni layer____________
		 
		 ii = 0;
		 do{
		 
		 
			 TracciaMC tr(0.,0.,coordinate[0],coordinate[1],coordinate[2]);
			 
			 tr.SetDistribEta(); 
			 tr.SetDistribPhi();
			 tr.Theta();
			 
			 tr.CalcCoeff();
			 
			 intersez = tr.intersezione(3);
			 if(intersez[2] < -27/2 || intersez[2] > 27/2) continue; //Ragionare se questa condizione va bene
			 
			
			 
			 for(int t=0; t<3; t++){
			 
				intersez = tr.intersezione(t+1);
				
				catter.NuoviAngoli(multipleScatt);
				catter.VarioAngolo(tr);
				tr.SetOrigine(intersez);
				fileout << intersez[0] << " " << intersez[1] << " " << intersez[2] << "		";}
				
				fileout<<endl;
			 
			 
			
			 ii++;
			 
		 }while( ii < molteplicità);
                
		 fileout.close();
		 
		 
		 
		 
  
//TIME @100 eventi 200 ms
		 
	//tcount_end_sim = std::chrono::steady_clock::now();
	//cout<< "Tempo simulazione: "<<std::chrono::duration_cast<std::chrono::microseconds>(tcount_end-tcount_start).count()<<endl;
	//tcount_start_ric = std::chrono::steady_clock::now();
		 
   //RICOSTRUZIONE
                        
			
			//per contare la molteplicità
			l = 0;
			do{
				fin1 >> Xtrue;
				l++;
			}while(!fin1.eof());
			fin1.close();
			
			m = l/9; //molteplicità
			
		
			
//TIME @100 eventi 230 ms		

			//ci serve salvare i dati del file rintersez_rec su una matrice per poter fare i confronti a coppie per ricostruire le tracklets
			//float mat[4][m][2];
			
			/*
			4 -> numero coordinate, X - Y - Z - Phi
			m -> molteplicità evento
			3 -> numero intersezioni coi layers*/	
			//usare vector o tree !!!
			
			i = 0;		// indice colonna(layer)
			j = 0;		// indice riga(molteplicità)
			
			do{				
				fin >> Xtrue;
				fin >> Ytrue; //salvo le coordinate delle intersezioni simulate
				fin >> Ztrue;
				
				Smearing ringo; //faccio lo smearing
				
				ringo.smearZ(Ztrue);
				ringo.smearPhi(Xtrue, Ytrue);
				
				//x,y,z,phi per ogni layer
				if (i%3 == 0) {};
				if ((i%3 == 1)||(i%3 == 2)) {
				mat[0][j][i%3] = ringo.GetXrec();
				mat[1][j][i%3] = ringo.GetYrec();  //salvo le coordinate smearate e quindi misurate nella matrice
				mat[2][j][i%3] = ringo.GetZrec();
				mat[3][j][i%3] = ringo.GetPhirec();}
								
				
				i++;
				
				if(i%3 == 0){
					//fout << endl;  relativa al file intersez_rec MA CHE NON SERVE!
					j++;
				}
				
			}while(i < l/3);
			
			fin.close();
			//fout.close();
	
				
			
//TIME @100 eventi 280 ms
			
//trovo il numero di coppie che stanno dentro l'angolo di Phi 0.01rad e le salvo nel file dati_tracklets
			numerocoppie = 0;
			
			for(int a = 0; a < m; a++){
				for(int b = 0; b < m; b++){
					float deltaPhi = mat[3][a][0] - mat[3][b][1];
					if(abs(deltaPhi) < 0.01){
						
						
						
						fout1 << mat[0][a][0] << " " << mat[1][a][0] << " " << mat[2][a][0] << " " << mat[0][b][1] << " " << mat[1][b][1] << " " << mat[2][b][1] << endl;
						numerocoppie++;
					}
				}
			}
	
//TIME	@ 100 eventi 360 ms		
			
			
			// __________tentativo ricostruzione vertice__________
			// eq retta passante per due punti A(x2,y2,z2) e B(x3,y3,z3)
			// x0 = x3 + c1t		c1 = sin(theta)*cos(phi)
			// y0 = y3 + c2t		c2 = sin(theta)*sin(phi)
			// z0 = z3 + c3t		c3 = cos(phi)
			// noi vogliamo solo la z
			// impongo y=0 e x=0, trovo quindi due valori di t e li medio tra loro (x3/c1 e y3/c2)
			
			
			
			//TH1D* hist = new TH1D ("h1","zreconstructed" , 54, -13.5, 13.5);
				
			for(int e=0; e<numerocoppie; e++)
			{
			        //prendo le coppie buone leggendo il file dati_tracklets
				 fin2 >> X2;
				 fin2 >> Y2;
				 fin2 >> Z2;
				 
				 fin2 >> X3;
				 fin2 >> Y3;
				 fin2 >> Z3;
				 
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
				 
				 hist -> Fill(Zrec);
				 
				
		}

		fileout1 << hist->GetMean() << endl; //salvo la media di tutto il mio istogramma, di tutte le Zrec (anche quelle derivate da coppie sbagliate) nel file Zrec (solo una per evento)	
                
                hist->Reset("ICESM");
                //delete hist;
		fin2.close();
		fout1.close();
		}
		
		delete hist;
		
		tcount_end = std::chrono::steady_clock::now();
			tim = std::chrono::duration_cast<std::chrono::microseconds>(tcount_end-tcount_start).count();
  tot = tim + tot;
  cout<<"TIME:"<< tot <<endl;		
			
//TIME	@100 eventi circa 400-500ms (grande variabilità)	
		
				

		/*hist->Draw();
		TFile file("zreconstructed.root", "recreate");
		hist->Write();
		file.Close();*/
             //    tcount_end_ric = std::chrono::steady_clock::now();
	//cout<<"Tempo ricostruzione:"<< std::chrono::duration_cast<std::chrono::microseconds>(tcount_end-tcount_start).count()<<endl;
	//cout<<"Delta(sim-ric):"<< std::chrono::duration_cast<std::chrono::microseconds>((tcount_end_sim-tcount_start_sim)-(tcount_end_ric-tcount_start_ric)).count()<<"evento con molt: "<<molteplicità<<endl;
	//tcount_start_sim = std::chrono::steady_clock::now();
		 
		 
		 fileout1.close();
		 fileout0.close();
		 
		 //creo grafico Ztrue-Zrec, avrò numeroeventi dati
		 TCanvas *c1=new TCanvas("c1","c1",800,600);
		 TH1D* hist1 = new TH1D ("h2","Ztrue-Zrec" , 100, -0.5, 0.5);
		 
		 ifstream filein0("Zsim.txt");
	         ifstream filein1("Zrec.txt");	
	         
	         
		 
		 for(int u = 0; u < numeroeventi; u++){
		 
		 filein0 >> Zsim1;
		 filein1 >> Zrec1;
		 
		 differenza = Zsim1 - Zrec1;
		 hist1 -> Fill(differenza); //grafico la differenza
		 }
		 
		hist1->Draw();
		TFile file1("ztrue-zrec.root", "recreate");
		hist1->Write();
		file1.Close();
		//delete hist1;
		
  double TT = time.CpuTime();	
  cout<<"Il tempo impiegato dalla CPU è "<<TT<<" s"<<endl;
  
 
		
	
 }
 
