#include <Riostream.h>
#include "TracciaMC.h" 
#include "TObject.h" 
#include "TRandom3.h"
#include "TMath.h"

ClassImp(TracciaMC)
/*
TracciaMC::TracciaMC():TRandom3(),
fEta(0.),
fPhi(0.),
fTheta(0.),
fT(0.),
fH(27.)
{

   origine = new float[3];
   fC = new float[3];
   origine[0] = 0.;	
   origine[1] = 0.;
   origine[2] = 0.;
   fC[0] = 0.;
   fC[1] = 0.;
   fC[2] = 0.;
}
*/

TracciaMC::TracciaMC(float theta, float phi, vector<float> orig):TRandom3(),
fEta(0.),
fPhi(phi),
fTheta(theta),
fT(0.),
fH(27.)
{ 
   origine = new float[3];
   origine[0] = orig[0];	
   origine[1] = orig[1];
   origine[2] = orig[2];
   
   fC = new float[3];
}

// MEMBER FUNCTIONS DEFINITIONS

void TracciaMC::SetDistribEta(){
	//bisogna capire
	//fEta = Gaus(0,1); 
	fEta = (gRandom->Rndm()*2) - 1;
}

void TracciaMC::SetDistribPhi(){
	fPhi = gRandom->Rndm()*2*TMath::Pi();		// distrib unif [0,2Pi]
}

void TracciaMC::Theta(){
	float argo = TMath::Exp(-fEta);
        fTheta = 2*TMath::ATan(argo);
}

void TracciaMC::CalcCoeff(){
	fC[0] = TMath::Sin(fTheta)*TMath::Cos(fPhi);		// c1
	fC[1] = TMath::Sin(fTheta)*TMath::Sin(fPhi);        // c2
	fC[2] = TMath::Cos(fTheta);							// c3
}

void TracciaMC::SetCoeff(vector<float> C){
       for(int i=0; i<3; i++){
         fC[i]=C[i];}
}

void TracciaMC::SetOrigine(vector<float> O){
       for(int i=0; i<3; i++){
         origine[i]=O[i];}
}

//funzione per il calcolo del parametro t
float t (float r, float x0, float y0, float C[3])
{
 
       float delta = (x0*C[0] + y0*C[1])*(x0*C[0] + y0*C[1]) - (C[0]*C[0] + C[1]*C[1])*(x0*x0 + y0*y0 - r*r);
       float num = -(x0*C[0] + y0*C[1]) + TMath::Sqrt(delta);
       float den = C[0]*C[0] + C[1]*C[1];
       return num/den;
       
       }
       
vector<float> TracciaMC::intersezione(int layer){

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
	     cout<< "Ci sono solo 3 layer, inserisci un numero compreso tra 1 e 3" <<endl;
	     break;
		}
   
       fT = t(R,origine[0],origine[1],fC);
       
       vector<float> coord_int;
       
       for(int i = 0; i < 3; i++){
		coord_int.push_back(origine[i] + fC[i]*fT);
	}
	
       return coord_int;}
       
float TracciaMC::GetPhi(){
       return fPhi;				
}

float TracciaMC::GetTheta(){
	return fTheta;				
}

vector<float> TracciaMC::GetC(){
	
	vector<float> C;

        C.push_back(fC[0]);
        C.push_back(fC[1]);
        C.push_back(fC[2]);

        return C;}
        
vector<float> TracciaMC::GetO(){
	
	vector<float> O;

        O.push_back(origine[0]);
        O.push_back(origine[1]);
        O.push_back(origine[2]);

        return O;}
       
float TracciaMC::GetT(){
	return fT;				
}
