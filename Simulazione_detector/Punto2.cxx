#include "TObject.h"
#include "Punto2.h"



ClassImp(Punto2)

//________________________________________________________________________
Punto2::Punto2():TObject(),
 fX(0.),
 fY(0.),
 fZ(0.),
 fPhi(0.){
   // default constructor
 }


//___________________________________________________________________________
Punto2::Punto2(float X, float Y, float Z, float Phi):TObject(),
 fX(X),
 fY(Y),
 fZ(Z),
 fPhi(Phi){
	//standard constructor 
}	     

//___________________________________________________________________________
Punto2::~Punto2()	 {
  // destructor
}
