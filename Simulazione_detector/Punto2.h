#ifndef PUNTO2_H
#define PUNTO2_H

#include "TObject.h"

class Punto2 : public TObject
{

public:

Punto2();
Punto2(float Phi, float Z, int num);

virtual ~Punto2();

 float GetZ() const {return fZ;}
 float GetPhi() const {return fPhi;}
 int Getnum() const {return fnum;}

private:


int fnum;
float fZ;
float fPhi;

//ClassDef(Punto2,1)
};


#endif
