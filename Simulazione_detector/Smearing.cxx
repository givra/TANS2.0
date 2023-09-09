#include <Riostream.h>
#include "Smearing.h" 
#include "TObject.h" 
#include "TRandom3.h"
#include "TMath.h"

ClassImp(Smearing)

Smearing::Smearing():TRandom3(),
	Zrec(0.),
	Phirec(0.)
	{
		//////////////////////////////
	}
	
void Smearing::smearZ(float Ztrue){
		
        Zrec = Ztrue + gRandom->Gaus(0,0.012);
	}
	
void Smearing::smearPhi(float Phitrue, int layer){
		
       float R = 1;
       
       switch(layer)
	    {case 1:
	     R=3;
	     break;
	     
             case 2:
	     R=4;
	     break;
	     
	     case 3:
	     R=7;
	     break;
	     
             default:
	     cout<<"Ci sono solo 3 layer, inserisci un numero compreso tra 1 e 3"<<endl;
	     break;
	     }
				
       Phirec = Phitrue + (gRandom->Gaus(0,0.003))/R;		
       }
	
float Smearing::GetZrec(){

       return Zrec;
       }
	
float Smearing::GetPhirec(){

       return Phirec;
       }
	
Smearing::~Smearing() {
	//distruttore di default
}
