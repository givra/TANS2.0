#include "TObject.h"
#include "Punto.h"

ClassImp(Punto)


Punto::Punto():TObject(),
 fnum(0),
 fZ(0.),
 fPhi(0.){
   // default constructor
 }

Punto::Punto(float Phi, float Z, int num):TObject(),
 fnum(num),
 fZ(Z),
 fPhi(Phi){
	//standard constructor 
}	     

Punto::~Punto()	 {
  // destructor
}
