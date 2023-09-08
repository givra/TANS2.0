//Classe per simulare gli effetti del Multiple Scattering  

#ifndef MULTIS_H
#define MULTIS_H

#include <vector>
#include "TObject.h" 
#include "TRandom3.h"

#include "Traccia.h"

class Multis : public TRandom3 { 

public: 
  
  Multis();
  virtual ~Multis();
   
  void VarioAngolo(Traccia track); //vado a modificare fC di track
  void NuoviAngoli(bool m);
    
private: 
  
  float Phip;
  float Thetap;
 
  };
  
  #endif
