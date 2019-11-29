// Authors: Spencer Behling, Benjamin Fenker, and Melissa Anholm - 2013

#ifndef K37ComplexGammaFunction_h
#define K37ComplexGammaFunction_h 1

#define kComplexGammeFunctionCoefficientsSize (10)

#include "globals.hh"

class K37ComplexGammaFunction 
{
 public:
  K37ComplexGammaFunction();
  ~K37ComplexGammaFunction();
  void computeK37ComplexGammaFunction(char, G4double, G4double);
  G4double realPart(char, G4double, G4double);
  G4double absSquaredComplexGamma(char, G4double, G4double);
  G4double realGamma(G4double);
  G4double squaredRealGamma(G4double);

 private:
  G4double coefficients[kComplexGammeFunctionCoefficientsSize];
  G4double Real_1;
  G4double Real_0;
  G4double Imaginary_1;
  G4double NA;
  G4double Z1;
  G4double TH;
  G4double GR;
  G4double GI;
  G4double T;
  G4double GR1;
  G4double GI1;
  G4double TH1;
  G4double SR;
  G4double SI;
  G4double Z2;
  G4double TH2;
  G4double G0;
};

#endif
