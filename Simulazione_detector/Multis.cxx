#include <Riostream.h>
#include "Multis.h"
#include "TracciaMC.h" 
#include "TObject.h" 
#include "TRandom3.h"
#include "TMath.h"

ClassImp(Multis) //serve a root 


Multis::Multis():TRandom3(),
 Thetap(0.),
 Phip(0.)
{ 
 //////////////////// 
}
 
void Multis::NuoviAngoli(bool m){
	if(m){
		Phip = gRandom->Rndm()*2*TMath::Pi();
		Thetap = gRandom->Gaus(0,0.001);
	}
	else{
		Phip = 0.;
		Thetap = 0.;
	}
}

void Multis::VarioAngolo(TracciaMC track){

 float Phi = track.GetPhi();
 float Theta = track.GetTheta();

 std::array<float,3> cd = track.GetC();

 float mr[3][3];

 mr[0][0] = -TMath::Sin(Phi);
 mr[1][0] = TMath::Cos(Phi);
 mr[2][0] = 0.;
 mr[0][1] = -TMath::Cos(Phi)*TMath::Cos(Theta);
 mr[1][1] = -TMath::Cos(Theta)*TMath::Sin(Phi);
 mr[2][1] = TMath::Sin(Theta);
 mr[0][2] = TMath::Sin(Theta)*TMath::Cos(Phi);
 mr[1][2] = TMath::Sin(Theta)*TMath::Sin(Phi);
 mr[2][2] = TMath::Cos(Theta);
 
 float cdp[3];
 
 cdp[0] = TMath::Sin(Thetap)*TMath::Cos(Phip);
 cdp[1] = TMath::Sin(Thetap)*TMath::Sin(Phip);
 cdp[2] = TMath::Cos(Thetap);
 
 for(int i=0; i<3; i++){
   cd[i]=0.;
   for(int j=0; j<3; j++){
     cd[i]+=mr[i][j]*cdp[j];
     }
   }
 
 track.SetCoeff(cd);
}

Multis::~Multis() {
	//distruttore di default
}
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      

