#include <Riostream.h>
#include "Smearing.h" 
#include "TObject.h" 
#include "TRandom3.h"
#include "TMath.h"

ClassImp(Smearing)

Smearing::Smearing():TRandom3(),
	Zrec(0.),
	Phirec(0.),
	Xrec(0.),
	Yrec(0.)
	{
		//////////////////////////////
	}
	
	void Smearing::smearZ(float Ztrue){
		
		Zrec = Ztrue + gRandom->Gaus(0,0.012);
	}
	
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
	}
	
	float Smearing::GetXrec(){
		return Xrec;
	}
	
	float Smearing::GetYrec(){
		return Yrec;
	}
	
	float Smearing::GetZrec(){
		return Zrec;
	}
	
	float Smearing::GetPhirec(){
		return Phirec;
	}