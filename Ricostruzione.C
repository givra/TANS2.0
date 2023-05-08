
#include <iostream>
#include <fstream>
#include <math.h>

#include "Vertex.h"
#include "TracciaMC.h"
#include "Multis.h"
#include "Smearing.h"
#include "TRandom3.h"

#include <vector>

using namespace std;

//float T (float x3, float y3, float C[3]);

void Ricostruzione(){
	
	TH1D* hist = new TH1D ("h1","zreconstructed" , 80, -14., 14.);
	
	ifstream fin("intersez_simulate.txt");
	ifstream fin1("intersez_simulate.txt");
	ifstream fin2("dati_tracklets.txt");
	//ifstream filein1("Zvertex.txt");
	
	ofstream fout("intersez_rec.txt");
	ofstream fout1("dati_tracklets.txt");
	ofstream fileout1("Zrec.txt");
			
	
	for(int evento = 0; evento < 2; evento++){
	
			float Xtrue;
			float Ytrue;
			float Ztrue;
			

			// per contare la molteplicità
			int l = 0;
			do{
				fin1 >> Xtrue;
				l++;
			}while(!fin1.eof());
			fin1.close();
			
			int m = l/9; 		// molteplicità
			
			float mat[4][m][3];		// contiene le coordinate delle intersezioni
			
			int i = 0;		// indice colonna(layer)
			int j = 0;		// indice riga(molteplicità)
			
			do{
				
				fin >> Xtrue;
				fin >> Ytrue;
				fin >> Ztrue;
				

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
			
			fin.close();
			fout.close();
			
			/*
			for(int k = 0; k < 4; k++){
				for(int j = 0; j < l/9; j++){
					for(int i = 0; i < 3; i++){
						cout << mat[k][j][i] << " ";
					}
					cout << endl;
				}
				cout << endl;
			}*/
			
			
			int numerocoppie = 0;
			
			for(int a = 0; a < m; a++){
				for(int b = 0; b < m; b++){
					float deltaPhi = mat[3][a][1] - mat[3][b][2];
					if(abs(deltaPhi) < 0.01){
						
						
						//		X2	|	Y2	|	Z2	|	X3	|	Y3	|	Z3
						fout1 << mat[0][a][1] << " " << mat[1][a][1] << " " << mat[2][a][1] << " " << mat[0][b][2] << " " << mat[1][b][2] << " " << mat[2][b][2] << endl;
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
				 Zmedio += Zrec;
				
				 hist -> Fill(Zrec);
				
		}

			//filein1 >> Zsim;
			fileout1 << Zmedio/numerocoppie << endl;

		fin2.close();
		fout1.close();
		fileout1.close();
		//filein1.close();
                
                cout<<hist->GetMean() <<endl;
		hist->Draw();
		TFile file("zreconstructed.root", "recreate");
		hist->Write();
		file.Close();
		
		
	}
}

/*
float T (float x, float y, float C[3])
{
       float r = 4.;
       float delta = (x*C[0] + y*C[1])*(x*C[0] + y*C[1]) - (C[0]*C[0] + C[1]*C[1])*(x*x + y*y - r*r);
       float num = -(x*C[0] + y*C[1]) + TMath::Sqrt(delta);
       float den = C[0]*C[0] + C[1]*C[1];
       return num/den;
       
       }*/
