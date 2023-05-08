//Progetto esame - Simulazione Montecarlo di una collisione protone protone

#include <iostream>
#include <fstream>
#include <math.h>

#include "Vertex.h"
#include "TracciaMC.h"
#include "Multis.h"
#include "TRandom3.h"

#include <vector>

using namespace std;

void Simulazione3(){
	
	/*unsigned int seed;
	
	cout << " inserire un seed: " << endl;
	cin >> seed;
	
 gRandom->SetSeed(seed);*/
	
//NON FUNZIONA...PERCHÈ???
 bool multipleScatt = 0;	// 0 = false, 1 = true
 cout << " attivare multiple scattering? " << endl;
 cin >> multipleScatt;
 
 
 
  // ______________grafico 3d intersezioni_______________
		  
        
		
	ofstream fileout0("Zsim.txt");
	ofstream fileout1("Zrec.txt");		 

        int numeroeventi = 500;
        
	for(int evento = 0; evento < numeroeventi; evento++){
	
	fstream file("intersez_simulate.txt", ios::out | ios::in);	// in teoria così dovrebbe essere in modalità scrittura e lettura
        //ifstream fin("intersez_simulate.txt");
        //ifstream fin1("intersez_simulate.txt");		
        
        ofstream fout("intersez_rec.txt");
                	
	fstream fin2("dati_tracklets.txt", ios::in | ios::out);
	//ofstream fout1("dati_tracklets.txt");
	
	//SIMULAZIONE
		 Vertex vertice0;
		 
		 vector<float> coordinate = vertice0.GetCoordinate();
		 
		 int molteplicità = vertice0.GetMolteplicity();
		 		 
		
		 fileout0 << coordinate[2] << endl ;
		 
		  // __________creo oggetto Multis per il Multiple Scattering____________
		  
		  Multis catter;
		 
		 // __________calcolo intersezioni su ogni layer____________
		 
		 int ii = 0;
		 do{
		 
		 
			 TracciaMC tr(0.,0.,coordinate);
			 
			 tr.SetDistribEta(); 
			 tr.SetDistribPhi();
			 tr.Theta();
			 
			 tr.CalcCoeff();
			 
			 vector<float> intersez = tr.intersezione(3);
			 if(intersez[2] < -27/2 || intersez[2] > 27/2) continue;
			 
			
			 
			 for(int t=0; t<3; t++){
			 
				intersez = tr.intersezione(t+1);
				
				catter.NuoviAngoli(multipleScatt);
				catter.VarioAngolo(tr);
				tr.SetOrigine(intersez);
				file << intersez[0] << " " << intersez[1] << " " << intersez[2] << "		";}
				
				file<<endl;
			 
			 
			
			 ii++;
			 
		 }while( ii < molteplicità);
                
		 file.close();
		 //file.clear();
		 
		 
   //RICOSTRUZIONE
                        float Xtrue;
			float Ytrue;
			float Ztrue;
			

			// per contare la molteplicità
			int l = 0;
			do{
				file >> Xtrue;
				l++;
			}while(!file.eof());
			file.close();
			//file.clear();
			file.seekg(0,ios::beg); 		// sintassi copiata da "esercizio1.C" dovrebbe ritornare all'inizio del file senza chiudere e riaprire
			
			int m = l/9; 		// molteplicità
			
			float mat[4][m][3];		// contiene le coordinate delle intersezioni
			
			int i = 0;		// indice colonna(layer)
			int j = 0;		// indice riga(molteplicità)
			
			do{
				
				file >> Xtrue;
				file >> Ytrue;
				file >> Ztrue;
				

				Smearing ringo;
				
				ringo.smearZ(Ztrue);
				ringo.smearPhi(Xtrue, Ytrue);
				
				//x,y,z,phi per ogni layer
				
				mat[0][j][i%3] = ringo.GetXrec();
				mat[1][j][i%3] = ringo.GetYrec();
				mat[2][j][i%3] = ringo.GetZrec();
				mat[3][j][i%3] = ringo.GetPhirec();
				
				fout << mat[0][j][i%3] << " " << mat[1][j][i%3] << " " << mat[2][j][i%3] << " " << mat[3][j][i%3] << "		";
				
				i++;
				
				if(i%3 == 0){
					fout << endl;
					j++;
				}
				
			}while(i < l/3);
			
			file.close();
			file.clear();
			fout.close();
			fout.clear();
			
			
			
			int numerocoppie = 0;
			
			for(int a = 0; a < m; a++){
				for(int b = 0; b < m; b++){
					float deltaPhi = mat[3][a][1] - mat[3][b][2];
					if(abs(deltaPhi) < 0.01){
						
						
						
						fin2 << mat[0][a][1] << " " << mat[1][a][1] << " " << mat[2][a][1] << " " << mat[0][b][2] << " " << mat[1][b][2] << " " << mat[2][b][2] << endl;
						numerocoppie++;
					}
				}
			}
			
			
			
			// __________tentativo ricostruzione vertice__________
			// eq retta passante per due punti A(x2,y2,z2) e B(x3,y3,z3)
			// x0 = x3 + c1t		c1 = sin(theta)*cos(phi)
			// y0 = y3 + c2t		c2 = sin(theta)*sin(phi)
			// z0 = z3 + c3t		c3 = cos(phi)
			// noi vogliamo solo la z
			// impongo y=0 e x=0, trovo quindi due valori di t e li medio tra loro (x3/c1 e y3/c2)
			
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
			
			float Zmedio = 0.;
			float Zsim;
			
			TH1D* hist = new TH1D ("h1","zreconstructed" , 54, -13.5, 13.5);
				
			for(int e=0; e<numerocoppie; e++)
			{
			
				 fin2 >> X2;
				 fin2 >> Y2;
				 fin2 >> Z2;
				 
				 fin2 >> X3;
				 fin2 >> Y3;
				 fin2 >> Z3;
			
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
				 
				 Zrec = Z3 - C[2]*t;
				 
				 hist -> Fill(Zrec);
				
		}

		fileout1 << hist->GetMean() << endl;	
                
                delete hist;
		fin2.close();
		fin2.clear();
		//fout1.close();
		//fout1.clear();
		

		/*hist->Draw();
		TFile file("zreconstructed.root", "recreate");
		hist->Write();
		file.Close();*/
   
		 }
		 
		 fileout1.close();
		 //fileout1.clear();
		 fileout0.close();
		 //fileout0.clear();
		 
		 
		 TCanvas *c1=new TCanvas("c1","c1",800,600);
		 TH1D* hist1 = new TH1D ("h2","ztrue-zrec" , 100, -0.5, 0.5);
		 
		 ifstream filein0("Zsim.txt");
	         ifstream filein1("Zrec.txt");	
	         
	         float differenza;
	         float Zsim1;
	         float Zrec1;
		 
		 for(int u = 0; u < numeroeventi; u++){
		 
		 filein0 >> Zsim1;
		 filein1 >> Zrec1;
		 
		 differenza = Zsim1 - Zrec1;
		 hist1 -> Fill(differenza);}
		 
		hist1->Draw();
		TFile file1("ztrue-zrec.root", "recreate");
		hist1->Write();
		
		//file1.close();
		/*fileout1.close();
		fileout1.clear();
		fileout0.close();
		fileout0.clear();*/
		filein0.close();
		//filein0.clear();
		filein1.close();
		//filein1.clear();
	
 }
