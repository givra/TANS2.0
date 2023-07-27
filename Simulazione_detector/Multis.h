//Classe per simulare gli effetti del Multiple Scattering  

#ifndef MULTIS_H
#define MULTIS_H

#include <vector>
#include "TObject.h" 
#include "TRandom3.h"

#include "TracciaMC.h"

class Multis : public TRandom3 { 

public: 
  
  Multis();
  
 
  void VarioAngolo(TracciaMC track); //vado a modificare fC di track
  void NuoviAngoli(bool m);
    
private: 
  
  float Phip;
  float Thetap;
 
  };
  
  #endif
