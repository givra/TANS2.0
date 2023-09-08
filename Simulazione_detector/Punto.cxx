#include "TObject.h"
#include "Punto.h"



ClassImp(Punto)

//________________________________________________________________________
Punto::Punto():TObject(),
 fnum(0),
 fZ(0.),
 fPhi(0.){
   // default constructor
 }


//___________________________________________________________________________
Punto::Punto(float Phi, float Z, int num):TObject(),
 fnum(num),
 fZ(Z),
 fPhi(Phi){
	//standard constructor 
}	     

//___________________________________________________________________________
Punto::~Punto()	 {
  // destructor
}
