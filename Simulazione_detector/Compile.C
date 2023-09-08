#include "TString.h"
#include "Traccia.h"
#include "Vertex.h"
#include "Multis.h"
#include "Smearing.h"
#include "Punto.h"

#include "TMath.h"

void Compile(TString myopt="fast"){
  TString opt;
  if(myopt.Contains("force")){
    opt = "kfg";
  }
  else {
    opt = "kg";
  }
  gSystem->CompileMacro("Traccia.cxx",opt.Data());
  gSystem->CompileMacro("Vertex.cxx",opt.Data());
  gSystem->CompileMacro("Multis.cxx",opt.Data());
  gSystem->CompileMacro("Smearing.cxx",opt.Data());
  gSystem->CompileMacro("Punto.cxx",opt.Data());
  //gSystem->CompileMacro("Punto2.cxx",opt.Data());
 // gSystem->CompileMacro("Simulazione2.C",opt.Data());
 // gSystem->CompileMacro("Ricostruzione2.C",opt.Data());
}
