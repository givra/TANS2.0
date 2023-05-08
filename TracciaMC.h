// Classe per direzione particella post scattering

#ifndef TRACCIAMC_H
#define TRACCIAMC_H

#include <vector>
#include "TObject.h" 
#include "TRandom3.h"

class TracciaMC : public TRandom3 {
	
	public:
	
	//TracciaMC();	//costruttore di default
	TracciaMC(float theta, float phi, vector<float> orig);
		
	void SetDistribEta();			// distribuzione da vedere, per ora è gaussiana
	void SetDistribPhi();			// uniforme
	void Theta();		       // funzie che da Eta (Pseudorapidità) calcola Theta
	
	void CalcCoeff();	// calcola c1,c2,c3 a partire dagli angoli
	void SetCoeff(vector<float> C);
	
	void SetOrigine(vector<float> O);
	vector<float> GetO();
	
	float GetPhi();
	float GetTheta();
	
	float GetT();
	
	vector<float> GetC();
	
	vector<float> intersezione(int layer);
	
	private:
	float fEta;
	float fPhi;
	float fTheta;
	
	float* origine;
	
	float* fC; 		// vettore con c1, c2, c3
	float fT;
	float fH;
	
};
#endif
