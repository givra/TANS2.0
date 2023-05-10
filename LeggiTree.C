#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "Punto.h"
#include "TMath.h"
#include "TClonesArray.h"

void LeggiTree(){

  int numeroeventi = 100;
  // definizione struct
  typedef struct {
    float X,Y,Z;
    int mult;} VTX;
  static VTX point;
  // Dichiarazione TClonesArray
  TClonesArray *hits1 = new TClonesArray("Punto",numeroeventi);
  TClonesArray *hits2 = new TClonesArray("Punto",numeroeventi);
  TClonesArray *hits3 = new TClonesArray("Punto",numeroeventi);
  //Apertura file di input
  TFile hfile("htree.root");
  //Lettura TTree  e branch
  TTree *tree = (TTree*)hfile.Get("T");
  TBranch *b1=tree->GetBranch("VertMult");
  TBranch *b2=tree->GetBranch("HitsPrimo");
  TBranch *b3=tree->GetBranch("HitsSecondo");
  TBranch *b4=tree->GetBranch("HitsTerzo");

  // Definizione degli indirizzi per la lettura dei dati su ttree
  b1->SetAddress(&point.X);
  b2->SetAddress(&hits1);
  b3->SetAddress(&hits2);
  b4->SetAddress(&hits3);

  // loop sugli ingressi nel TTree
  for(int ev=0;ev<tree->GetEntries();ev++){
    tree->GetEvent(ev);
    cout<<"Evento "<<ev<<"; Molteplicita= "<<point.mult<<endl;
    cout<<"X,Y,Z = "<<point.X<<"; "<<point.Y<<"; "<<point.Z<<endl;
    int num=hits1->GetEntries();
    cout<<"Numero di elementi nel TClonesArray "<<num<<endl; //viene a tutti zero (nel senso di solo un dato)
    //è giusto?? Forse vuol dire che nel buffer ce ne sta solo uno che prende e da 
    for (int j=0; j<num; j++){
      Punto *tst=(Punto*)hits1->At(j);
      cout<<"Intersezione primo layer"<<j<<") x, y, z = "<<tst->GetX()<<"; "<<tst->GetY()<<"; "<<tst->GetZ()<<endl;
    }
  }
 
}
/*
#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "Punto.h"
#include "TMath.h"
#include "TClonesArray.h"

void LeggiTree(){
  // definizione struct
  typedef struct {
    double X,Y,Z;
    int mult;} VTX;
  static VTX point;
  // Dichiarazione TClonesArray
  TClonesArray *hits1 = new TClonesArray("Punto",100);
  //Apertura file di input
  TFile hfile("htree.root");
  //Lettura TTree  e branch
  TTree *tree = (TTree*)hfile.Get("T");
  TBranch *b1=tree->GetBranch("VertMult");
  TBranch *b2=tree->GetBranch("HitsPrimo");

  // Definizione degli indirizzi per la lettura dei dati su ttree
  b1->SetAddress(&point.X);
  b2->SetAddress(&hits1);

  // loop sugli ingressi nel TTree
  for(int ev=0;ev<tree->GetEntries();ev++){
    tree->GetEvent(ev);
    cout<<"Evento "<<ev<<"; Molteplicita= "<<point.mult<<endl;
    cout<<"X,Y,Z = "<<point.X<<"; "<<point.Y<<"; "<<point.Z<<endl;
    int num=hits1->GetEntries();
    cout<<"Numero di elementi nel TClonesArray "<<num<<endl;
    for (int j=0; j<num; j++){
      Punto *tst=(Punto*)hits1->At(j);
      cout<<"Punto "<<j<<") x, y, z = "<<tst->GetX()<<"; "<<tst->GetY()<<"; "<<tst->GetZ()<<endl;
    }
  }
 
}*/
