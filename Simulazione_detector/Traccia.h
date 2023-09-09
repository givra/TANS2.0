#ifndef Traccia_H
#define Traccia_H

#include <vector>
#include "TObject.h" 
#include "TRandom3.h"
#include "TH1F.h"

using namespace std;

class Traccia : public TRandom3 {
	
	public:
	
	Traccia();	//costruttore di default
	
	virtual ~Traccia();
	
	void SetOrigine(float Xo, float Yo, float Zo);	
	
	void SetEtaUni();
	void SetEta(bool Eta);   //scegli se estrarre con uniforme (0) o con kinem.root (1)
	
	void SetPhi();	// uniforme
	void Theta();  // funzione che da Eta (Pseudorapidità) calcola Theta
	
	void CalcCoeff();  // calcola c1,c2,c3 a partire dagli angoli
	void SetCoeff(std::array<float,3> C);	 
	
	void SetHit(std::array<float,2> inter, int lay);	
	
	float GetPhi();
	float GetTheta();
	
	float GetT();
	
	std::array<float,3> GetC();
	
	std::array<float,2> intersezione(int layer);
	
	private:
	
	bool fEta;
	float fPhi;
	float fTheta;
	
	float origine[3];
	
	float fC[3];  //vettore con c1, c2, c3
	float fT;
	float fH;
	
	
};
#endif
