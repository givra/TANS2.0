#include <Riostream.h>
#include "Traccia.h" 
#include "TObject.h" 
#include "TRandom3.h"
#include "TMath.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TFile.h"
#include <fstream>

ClassImp(Traccia)

Traccia::Traccia():TRandom3(),
fEta(0.),
fPhi(0.),
fTheta(0.),
fT(0.),
fH(27.)
{
 
  // origine = new float[3];
  // fC = new float[3];
   origine[0] = 0.;	
   origine[1] = 0.;
   origine[2] = 0.;
   fC[0] = 0.;
   fC[1] = 0.;
   fC[2] = 0.;
}

// MEMBER FUNCTIONS DEFINITIONS
void Traccia::SetEtaUni(){
              fEta = (gRandom->Rndm()*2) - 1;
}

void Traccia::SetEta(float Eta){
              fEta = Eta;
}

void Traccia::SetPhi(){
	fPhi = gRandom->Rndm()*2*TMath::Pi();		// distrib unif [0,2Pi]
}

void Traccia::Theta(){
	float argo = TMath::Exp(-fEta);
        fTheta = 2*TMath::ATan(argo);
}

void Traccia::SetOrigine(float Xo, float Yo, float Zo){
 
   //origine = new float[3];
   origine[0] = Xo;	
   origine[1] = Yo;
   origine[2] = Zo;
   
}
void Traccia::CalcCoeff(){
	fC[0] = TMath::Sin(fTheta)*TMath::Cos(fPhi);		// c1
	fC[1] = TMath::Sin(fTheta)*TMath::Sin(fPhi);		        // c2
	fC[2] = TMath::Cos(fTheta);				// c3
}

void Traccia::SetCoeff(std::array<float,3> C){	
       for(int i=0; i<3; i++){
         fC[i]=C[i];}
}

void Traccia::SetHit(std::array<float,2> inter, int lay){ 

         float R = 1;
       
       switch(lay)
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

         origine[0] = R*TMath::Cos(inter[0]);
         origine[1] = R*TMath::Sin(inter[0]);
         origine[2] = inter[1];
}

//funzione per il calcolo del parametro t
float t (float r, float x0, float y0, float C[3])
{
 
       float delta = (x0*C[0] + y0*C[1])*(x0*C[0] + y0*C[1]) - (C[0]*C[0] + C[1]*C[1])*(x0*x0 + y0*y0 - r*r);
       float num = -(x0*C[0] + y0*C[1]) + TMath::Sqrt(delta);
       float den = C[0]*C[0] + C[1]*C[1];
       return num/den;
       
       }
       
std::array<float,2> Traccia::intersezione(int layer){

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
   
       fT = t(R,origine[0],origine[1],fC);
       
       std::array<float,2> coord_int;
      // coord_int.reserve(2);
      if((origine[1] + fC[1]*fT)>=0) coord_int[0] = TMath::ACos((origine[0] + fC[0]*fT)/R);
      if((origine[1] + fC[1]*fT)<0) coord_int[0] = 2*TMath::Pi() - TMath::ACos((origine[0] + fC[0]*fT)/R);
      
      if(coord_int[0]<0) coord_int[0] = coord_int[0] + 2*TMath::Pi();
      if(coord_int[0]>2*TMath::Pi()) coord_int[0] = coord_int[0] - 2*TMath::Pi();
      
      coord_int[1] = origine[2] + fC[2]*fT;
      // coord_int.push_back(TMath::ACos((origine[0] + fC[0]*fT)/R));
      // coord_int.push_back(origine[2] + fC[2]*fT);
       /*for(int i = 0; i < 3; i++){
		coord_int.push_back(origine[i] + fC[i]*fT);
	}*/
	
       return coord_int;}
       
float Traccia::GetPhi(){
       return fPhi;				
}

float Traccia::GetTheta(){
	return fTheta;				
}

std::array<float,3> Traccia::GetC(){
	
	/*vector<float> C;

        C.push_back(fC[0]);
        C.push_back(fC[1]);
        C.push_back(fC[2]);*/
        
        std::array<float,3> C;
        for(int i = 0; i < 3; i++){
		C[i] = fC[i];
	}

        return C;}
        
/*vector<float> Traccia::GetO(){
	
	vector<float> O;

        O.push_back(origine[0]);
        O.push_back(origine[1]);
        O.push_back(origine[2]);

        return O;}*/
       
float Traccia::GetT(){
	return fT;				
}

/*int Traccia::GetLabel(){
	return fLabel;
}
*/
Traccia::~Traccia() {
	// distruttore di default
}
