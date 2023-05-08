//Progetto esame - Simulazione Montecarlo di una collisione protone protone

#include <iostream>
#include <fstream>
#include <math.h>
#include "TTree.h"
#include "TBranch.h"
#include "Vertex.h"
#include "TracciaMC.h"
#include "Multis.h"
#include "TRandom3.h"

#include <vector>

using namespace std;

void Simulazione2(){
	
	/*unsigned int seed;
	
	cout << " inserire un seed: " << endl;
	cin >> seed;
	
 gRandom->SetSeed(seed);*/
	
 bool multipleScatt = 0;		// 0 = false, 1 = true
 cout << " attivare multiple scattering? " << endl;
 cin >> multipleScatt;
 
  /* ______________grafico 3d intersezioni_______________
		  TCanvas *c1=new TCanvas("c1","c1",800,600);
		  
		  TH3F* h3 = new TH3F("h3", "Simulazione Detector", 16, -8, 8, 16, -8, 8, 28, -14, 14);
		  h3->Draw(); 
		  h3->GetXaxis()->SetTitle("X");
		  h3->GetYaxis()->SetTitle("Y");
		  h3->GetZaxis()->SetTitle("Z");
		*/  
    TFile hfile("htree.root","RECREATE");

	for(int evento = 0; evento < 1; evento++){
		
		 Vertex vertice0;
		 
		 vector<float> coordinate = vertice0.GetCoordinate();
		 
		 int molteplicità = vertice0.GetMolteplicity();
		 vector<float> intersez;
		 
		 float x2,y2,z2,x3,y3,z3;
		 
        //ofstream fileout("intersez_simulate.txt");		// creo file output
        //ofstream fileout1("Zsim.txt");		 
		  
		  
		 //TGraph *g = new TGraph(1 + molteplicità*3);
		 //TPolyMarker3D *g = new TPolyMarker3D(1 + molteplicità*3);
		 
		 //TPolyLine3D *traccie[molteplicità];
		 //for(int s=0; s<=molteplicità; s++) traccie[s] = new TPolyLine3D(4);
		 
		 
		 //g->SetPoint(0, coordinate[0], coordinate[1], coordinate[2]);
		 
		TTree* tree(NULL);
	    tree = new TTree("T","TREEPROVA");
	    
		//tree->Branch("Zsim", &coordinate[2], "coordinate[2]/D");
	    tree->Branch("x2", &x2, "x2/D");
	    tree->Branch("y2", &y2, "y2/D");
	    tree->Branch("z2", &z2, "z2/D");
	    tree->Branch("x3", &x3, "x3/D");
	    tree->Branch("y3", &y3, "y3/D");
	    tree->Branch("z3", &z3, "z3/D");
	    
		
		 cout<<"Molteplicità: "<<molteplicità<<endl;
		 cout << " Z vertice "<< coordinate[2] << endl;
		 //fileout1 << coordinate[2] << endl ;
		 
		  // __________creo oggetto Multis per il Multiple Scattering____________
		  
		  Multis catter;
		 
		 // __________calcolo intersezioni su ogni layer____________
		 
		 int i = 0;
		 do{
		 //for(int i=0; i<=molteplicità; i++){ 
		 
			 //traccie[i] -> SetPoint(0, coordinate[0], coordinate[1], coordinate[2]);
			 TracciaMC tr(0.,0.,coordinate);
			 
			 tr.SetDistribEta(); 
			 tr.SetDistribPhi();
			 tr.Theta();
			 
			 tr.CalcCoeff();
			 
			 intersez = tr.intersezione(3);
			 if(intersez[2] < -27/2 || intersez[2] > 27/2) continue;
			 
			 /*
			 // ________intersez layer 1 + multiple scatt ___________
			 intersez = tr.intersezione(1);
			 traccie[i]->SetPoint(1, intersez[0], intersez[1], intersez[2]);
			 
			 Multis catter1;
			 catter1.VarioAngolo(tr);
			 tr.SetOrigine(intersez);
			 
			 fout << intersez[0] << " " << intersez[1] << " " << intersez[2] << "	|	";
			 
			 // _________intersez layer 2 + multiple scatt___________
			 intersez = tr.intersezione(2);
			 traccie[i]->SetPoint(2, intersez[0], intersez[1], intersez[2]);

			 Multis catter2;
			 catter2.VarioAngolo(tr);
			 tr.SetOrigine(intersez);
			 
			 fout << intersez[0] << " " << intersez[1] << " " << intersez[2] << "	|	";
			 
			 // _________intersez layer 3 + multiple scatt___________
			 intersez = tr.intersezione(3);
			 traccie[i]->SetPoint(3, intersez[0], intersez[1], intersez[2]);
			 fout << intersez[0] << " " << intersez[1] << " " << intersez[2] << "\n";
			 */
			 
			 for(int t=0; t<3; t++){
			 
				intersez = tr.intersezione(t+1);
				//traccie[i]->SetPoint(t+1, intersez[0], intersez[1], intersez[2]);
				
				catter.NuoviAngoli(multipleScatt);
				catter.VarioAngolo(tr);
				tr.SetOrigine(intersez);
				//fileout << intersez[0] << " " << intersez[1] << " " << intersez[2] << "		";
					if(t == 1){
						x2 = intersez[0];
						y2 = intersez[1];
						z2 = intersez[2];
					
					}
					else if(t == 2){
						x3 = intersez[0];
						y3 = intersez[1];
						z3 = intersez[2];
					
				    }
				}
				
				//fileout<<endl;
			 cout << "x2: " << x2 << " y2: " << y2 << " z2: " << z2 << endl;
			 cout << "x3: " << x3 << " y3: " << y3 << " z3: " << z3 << endl;
			// traccie[i]->Draw();
			tree->Fill();
			
			 i++;
			 
		 }while( i < molteplicità);
		 
         //fileout.close();
		 //fileout1.close();      

		 }
		 //}
		hfile.Write();		// svuota il buffer di qualsiasi cosa contenga
		hfile.Close();
		 //c1->Update(); 
	//}
 }
