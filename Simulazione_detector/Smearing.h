 // Classe per smearare le intersezioni
 
 #ifndef SMEARING_H
 #define SMEARING_H
 
 #include <vector>
 #include "TObject.h" 
 #include "TRandom3.h"
 
 class Smearing : public TRandom3 {

public:

	Smearing();
	virtual ~Smearing();
	void smearZ(float Ztrue);
	void smearPhi(float Phitrue, int layer);
	//float GetXrec();
	//float GetYrec();
	float GetZrec();
	float GetPhirec();
	
private:

	float Zrec;
	float Phirec;
	//float Xrec;
	//float Yrec;
	
 };
 
 #endif
