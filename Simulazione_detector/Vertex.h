//Classe per Simulare il vertice primario della collisione

#ifndef VERTEX_H
#define VERTEX_H

#include <vector>
#include "TObject.h" 
#include "TRandom3.h"

class Vertex : public TRandom3 { //Lo si fa ereditare da TRandom??


public: 
  
  Vertex(); //costruttore di default, x,y,z=0 e molt=1
  virtual ~Vertex();
  
  void NewVertex(); //coordinate gaussiane , distribuzione uniforme per molteplicità tra 1 e 80
 
  void SetZUniform(); //impostare Z con distribuzione uniforme tra -20 e 20
  
  
  vector<float> GetCoordinate();
  int GetMolteplicity();
  
  float GetX();
  float GetY();
  float GetZ();
  
  
private: 
  
  
  float X;
  float Y;
  float Z;
  
  int moltep;
  

  };
  
  #endif
