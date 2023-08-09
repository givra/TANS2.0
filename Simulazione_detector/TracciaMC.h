// Classe per direzione particella post scattering

#ifndef TRACCIAMC_H
#define TRACCIAMC_H

#include <vector>
#include "TObject.h" 
#include "TRandom3.h"
#include "TH1F.h"


using namespace std;

class TracciaMC : public TRandom3 {
	
	public:
	
	TracciaMC();	//costruttore di default
	//TracciaMC(float theta, float phi, float Xo, float Yo, float Zo);
	
	virtual ~TracciaMC();
	
	void SetOrigine(float Xo, float Yo, float Zo);	
	
	void SetEtaUni();
	void SetEta(float Eta);   //scegli se estrarre con uniforme (1) o con kinem.root (2)
	
	void SetPhi();			// uniforme
	void Theta();		       // funzie che da Eta (Pseudorapidità) calcola Theta
	
	void CalcCoeff();	// calcola c1,c2,c3 a partire dagli angoli
	void SetCoeff(std::array<float,3> C);	 
	
	void SetHit(std::array<float,2> inter, int lay);	
	//vector<float> GetO();
	
	float GetPhi();
	float GetTheta();
	
	float GetT();
	//int GetLabel();
	
	std::array<float,3> GetC();
	
	std::array<float,2> intersezione(int layer);
	
	private:
	//int fLabel;
	float fEta;
	float fPhi;
	float fTheta;
	
	float origine[3];
	
	float fC[3]; 		// vettore con c1, c2, c3
	float fT;
	float fH;
	
	
};
#endif
