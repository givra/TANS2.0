#ifndef Punto_H
#define Punto_H

#include "TObject.h"

class Punto : public TObject
{

public:

Punto();
Punto(float Phi, float Z, int num);

virtual ~Punto();

 float GetZ() const {return fZ;}
 float GetPhi() const {return fPhi;}
 int Getnum() const {return fnum;}

private:


int fnum;
float fZ;
float fPhi;

//ClassDef(Punto,1)
};


#endif
