//Classe per Simulare il vertice primario della collisione

#ifndef VERTEX_H
#define VERTEX_H

#include <vector>
#include "TObject.h" 
#include "TRandom3.h"

class Vertex : public TRandom3 { //Lo si fa ereditare da TRandom??


public: 
  
  Vertex(); //costruttore di default, coordinate gaussiane , molteplicità fissa a 40
 
  void SetZUniform(); //impostare Z con distribuzione uniforme tra -20 e 20
  void SetMoltUniform(); //settare distribuzione uniforme per molteplicità tra 1 e 80
  
  vector<float> GetCoordinate();
  int GetMolteplicity();
  

  
private: 
  
  
  float X;
  float Y;
  float Z;
  
  int moltep;
  

  };
  
  #endif
