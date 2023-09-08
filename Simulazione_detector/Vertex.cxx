#include <Riostream.h>
#include "Vertex.h" 
#include "TObject.h" 
#include "TRandom3.h"
#include "TMath.h"

ClassImp(Vertex) 

Vertex::Vertex():TRandom3(){ 
      //unità di misura: cm

      X = 0.;
      Y = 0.;
      Z = 0.;

      moltep = 1;
      }
 
 
void Vertex::NewVertex(){

      X = gRandom->Gaus(0,0.01); //coordinata x
      Y = gRandom->Gaus(0,0.01); //coordinata y
      Z = gRandom->Gaus(0,5.3); //coordinata z
      }
 
void Vertex::SetMolt(float moltepl){

      moltep = moltepl;     
      }
      
void Vertex::SetMoltUniform(){
  
      moltep = 1 + (gRandom->Rndm()*80)/1;
      }

int Vertex::GetMolteplicity(){return moltep;}

float Vertex::GetX(){return X;}

float Vertex::GetY(){return Y;}

float Vertex::GetZ(){return Z;}

Vertex::~Vertex() {
	//distruttore di default
}




