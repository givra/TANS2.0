#include "TObject.h"
#include "Punto2.h"



ClassImp(Punto2)

//________________________________________________________________________
Punto2::Punto2():TObject(),
 fnum(0),
 fZ(0.),
 fPhi(0.){
   // default constructor
 }


//___________________________________________________________________________
Punto2::Punto2(float Phi, float Z, int num):TObject(),
 fnum(num),
 fZ(Z),
 fPhi(Phi){
	//standard constructor 
}	     

//___________________________________________________________________________
Punto2::~Punto2()	 {
  // destructor
}
