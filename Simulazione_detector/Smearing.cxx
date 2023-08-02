#include <Riostream.h>
#include "Smearing.h" 
#include "TObject.h" 
#include "TRandom3.h"
#include "TMath.h"

ClassImp(Smearing)

Smearing::Smearing():TRandom3(),
	Zrec(0.),
	Phirec(0.)
	//Xrec(0.),
	//Yrec(0.)
	{
		//////////////////////////////
	}
	
	void Smearing::smearZ(float Ztrue){
		
		Zrec = Ztrue + gRandom->Gaus(0,0.012);
	}
	
	void Smearing::smearPhi(float Phitrue, int layer){
		
		float R = 1;
       
       switch(layer)
	   { case 1:
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
	/*vecchia versione con le tre coordinate X, Y e Z
	void Smearing::smearPhi(float Xtrue, float Ytrue){
		
		float r = TMath::Sqrt(Xtrue*Xtrue + Ytrue*Ytrue);
		float Phitrue;
		if(Xtrue < 0){
			Phitrue = TMath::ATan(Ytrue/Xtrue) + TMath::Pi();
		}
		else Phitrue = TMath::ATan(Ytrue/Xtrue);
		
		Phirec = Phitrue + (gRandom->Gaus(0,0.003))/r;
		Xrec = r*TMath::Cos(Phirec);
		Yrec = r*TMath::Sin(Phirec);
	}*/
	
	/*float Smearing::GetXrec(){
		return Xrec;
	}
	
	float Smearing::GetYrec(){
		return Yrec;
	}*/
	
	float Smearing::GetZrec(){
		return Zrec;
	}
	
	float Smearing::GetPhirec(){
		return Phirec;
	}
	
Smearing::~Smearing() {
	//distruttore di default
}
