#ifndef PUNTO_H
#define PUNTO_H

#include "TObject.h"

class Punto : public TObject
{

public:

Punto();
Punto(float X, float Y, float Z);

virtual ~Punto();

 float GetX() const {return fX;} 
 float GetY() const {return fY;}
 float GetZ() const {return fZ;}


private:


float fX;
float fY;
float fZ;

ClassDef(Punto,1)
};


#endif
