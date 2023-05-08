//Progetto esame - Simulazione Montecarlo di una collisione protone protone

#include <iostream>
#include <math.h>

#include "Vertex.h"
#include "TracciaMC.h"
#include "Multis.h"

#include <vector>

void Simulazione(){

 Vertex vertice0;
 vertice0.SetMoltUniform();
 
 vector<float> coordinate = vertice0.GetCoordinate();
 
 cout<<"VERTICE"<<endl;
 cout<<"Coordinata x: "<<coordinate[0]<<"; Coordinata y: "<<coordinate[1]<<"; Coordinata z: "<<coordinate[2]<<endl;
 
 TGraph *g = new TGraph(3);
 g->SetPoint(0, coordinate[0], coordinate[1]);
 

 int molteplicità = vertice0.GetMolteplicity();
 
 cout<<"Molteplicità: "<<molteplicità<<endl;
 
 
 TracciaMC tr(0.,0.,coordinate);
 
 tr.SetDistribEta(); 
 tr.SetDistribPhi();
 tr.Theta();
 
 tr.CalcCoeff();
 
 vector<float> intersez = tr.intersezione(1);
 
 cout<<"INTERSEZIONE 1"<<endl;
 cout<<"Coordinata x: "<<intersez[0]<<"; Coordinata y: "<<intersez[1]<<"; Coordinata z: "<<intersez[2]<<endl;
 g->SetPoint(1, intersez[0], intersez[1]);
 
 vector<float> cfactor = tr.GetC();
 cout<<"c1 prima del MS: "<<cfactor[0]<<"c2 prima del MS: "<<cfactor[1]<<endl;
 
 Multis catter1;
 catter1.NuoviAngoli(1);
 catter1.VarioAngolo(tr);
 
 cfactor = tr.GetC();
 cout<<"c1 dopo il MS: "<<cfactor[0]<<"c2 dopo il MS: "<<cfactor[1]<<endl;
 
 tr.SetOrigine(intersez);
 
 intersez = tr.intersezione(2);

 cout<<"INTERSEZIONE 2"<<endl;
 cout<<"Coordinata x: "<<intersez[0]<<"; Coordinata y: "<<intersez[1]<<"; Coordinata z: "<<intersez[2]<<endl;
 g->SetPoint(2, intersez[0], intersez[1]);

 Multis catter2;
 catter2.NuoviAngoli(1);
 catter2.VarioAngolo(tr);
 tr.SetOrigine(intersez);
 
 intersez = tr.intersezione(3);
 
 cout<<"INTERSEZIONE 3"<<endl;
 cout<<"Coordinata x: "<<intersez[0]<<"; Coordinata y: "<<intersez[1]<<"; Coordinata z: "<<intersez[2]<<endl;
 g->SetPoint(3, intersez[0], intersez[1]);
 g->Draw("A*");
 
 TLine *l = new TLine(coordinate[0], coordinate[1], intersez[0], intersez[1]);
 l -> Draw();
 }
