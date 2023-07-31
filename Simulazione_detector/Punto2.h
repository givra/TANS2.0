#ifndef PUNTO2_H
#define PUNTO2_H

#include "TObject.h"

class Punto2 : public TObject
{

public:

Punto2();
Punto2(float X, float Y, float Z, float Phi);

virtual ~Punto2();

 float GetX() const {return fX;} 
 float GetY() const {return fY;}
 float GetZ() const {return fZ;}
 float GetPhi() const {return fPhi;}

private:


float fX;
float fY;
float fZ;
float fPhi;

//ClassDef(Punto2,1)
};


#endif
