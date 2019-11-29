// Authors: Spencer Behling and Benjamin Fenker 2013

#include <iomanip>
#include <cmath>
#include "K37ComplexGammaFunction.hh"
#include "globals.hh"



K37ComplexGammaFunction::K37ComplexGammaFunction()
  : Real_1(0), Real_0(0), Imaginary_1(0), NA(0), Z1(0), TH(0), GR(0), GI(0),
    T(0), GR1(0), GI1(0), TH1(0), SR(0), SI(0), Z2(0), TH2(0), G0(0) {
  coefficients[0] = 8.333333333333333e-2;
  coefficients[1] = -2.777777777777778e-3;
  coefficients[2] = 7.936507936507937e-4;
  coefficients[3] = -5.952380952380952e-4;
  coefficients[4] = 8.417508417508418e-4;
  coefficients[5] = -1.917526917526918e-3;
  coefficients[6] = 6.410256410256410e-3;
  coefficients[7] = -2.955065359477124e-2;
  coefficients[8] = 1.796443723688307e-1;
  coefficients[9] = -1.39243221690590;
}

K37ComplexGammaFunction::~K37ComplexGammaFunction() {
}

void K37ComplexGammaFunction::computeK37ComplexGammaFunction(char choice,
                                                             double Real,
                                                             double Imaginary) {
  if (choice != 'l' && choice != 'g') {
    G4cout << "The choices are 'l' for log gamma and 'g' for gamma" << G4endl;
    G4cout << "!!!!!!! You put: '" << choice << "' try again!!!!!!" << G4endl;
    return;
  }

  // G4cout<< std::boolalpha<< (Imaginary == 0.0) << G4endl;
  // G4cout<< std::boolalpha<< (Real == int(Real)) << G4endl;
  // G4cout<< std::boolalpha<< (Real < 0.0) << G4endl;
  if (Imaginary == 0.0 && (Real == static_cast<int>(Real)) && Real < 0.0) {
    G4cout << "there is a first order pole at negative integers" << G4endl;
    return;
  } else if (Real < 0.0) {
    Real_1 = Real;
    Imaginary_1 = Imaginary;
    Real = -Real;
    Imaginary = -Imaginary;
  }
  Real_0 = Real;
  if (Real < 7.0) {
    NA = static_cast<int>(7.0 - Real);
    Real_0 = Real + NA;
  }
  Z1 = sqrt(std::pow(Real_0, 2.0)+std::pow(Imaginary, 2.0));
  TH = std::atan(Imaginary/Real_0);

  GR = (Real_0-0.5)*log(Z1)-TH*Imaginary-Real_0+0.5*log(2.0*M_PI);
  GI = TH*(Real_0-0.5)+Imaginary*log(Z1)-Imaginary;

  for (int K =1; K <=10; ++K) {
    T = pow(Z1, (1.0-2.0*K));
    GR = GR+coefficients[K-1]*T*cos((2.0*K-1.0)*TH);
    GI = GI-coefficients[K-1]*T*sin((2.0*K-1.0)*TH);
  }

  if (Real < 7.0) {
    GR1 = 0.0;
    GI1 = 0.0;
    for (int J = 0; J <= (NA-1); ++J) {
      // G4cout<< J << G4endl;
      GR1 += 0.5*log(pow((Real+J), 2.0)+pow(Imaginary, 2.0));
      GI1+=atan(Imaginary/(Real+J));
    }
    GR-=GR1;
    GI-=GI1;
    //     G4cout<< GR1<<"   "<< GI1<<"   "<< GR<<"   " << GI<< G4endl;
  }


  if (Real_1 < 0.0) {
    Z1 = sqrt(Real*Real+Imaginary*Imaginary);
    TH1 = atan(Imaginary/Real);
    SR = -sin(M_PI*Real)*cosh(M_PI*Imaginary);
    SI = -cos(M_PI*Real)*sinh(M_PI*Imaginary);
    Z2 = sqrt(SR*SR+SI*SI);
    TH2 = atan(SI/SR);
    // G4cout<< "It is here" << G4endl;
    // G4cout<< Z1<<"  "<< TH1<< "  "<< SR << G4endl;
    if (SR < 0.0) {
      TH2 += M_PI;
    }

    GR = log(M_PI/(Z1*Z2))-GR;
    GI = -TH1-TH2-GI;
    Real = Real_1;
    Imaginary = Imaginary_1;
  }

  if (choice == 'l') {
    G4cout << choice << " " << Real << " " << Imaginary << " "<< GR <<" "
           << GI << G4endl;
  }

  if (choice == 'g') {
    G0 = exp(GR);
    GR = G0*cos(GI);
    GI = G0*sin(GI);
    G4cout << choice << " " << Real << " " << Imaginary << " " << GR << " "
           << GI << G4endl;
  }
}

double K37ComplexGammaFunction::realPart(char choice, double Real,
                                         double Imaginary) {
  if (choice != 'l' && choice != 'g') {
    return 4.23e17;
  }

  if (Imaginary == 0.0 && (Real == static_cast<int>(Real)) && Real < 0.0) {
    return 4.23e17;
  } else if (Real < 0.0) {
    Real_1 = Real;
    Imaginary_1 = Imaginary;
    Real = -Real;
    Imaginary = -Imaginary;
  }
  Real_0 = Real;
  if (Real < 7.0) {
    NA = static_cast<int>(7.0 - Real);
    Real_0 = Real + NA;
  }
  Z1 = sqrt(std::pow(Real_0, 2.0)+std::pow(Imaginary, 2.0));
  TH  = std::atan(Imaginary/Real_0);

  GR = (Real_0-0.5)*log(Z1)-TH*Imaginary-Real_0+0.5*log(2.0*M_PI);
  GI = TH*(Real_0-0.5)+Imaginary*log(Z1)-Imaginary;

  for (int K =1; K <=10; ++K) {
    T = pow(Z1, (1.0-2.0*K));
    GR = GR+coefficients[K-1]*T*cos((2.0*K-1.0)*TH);
    GI = GI-coefficients[K-1]*T*sin((2.0*K-1.0)*TH);
  }

  if (Real < 7.0) {
    GR1 = 0.0;
    GI1 = 0.0;
    for (int J = 0; J <= (NA-1); ++J) {
      GR1 += 0.5*log(pow((Real+J), 2.0)+pow(Imaginary, 2.0));
      GI1 += atan(Imaginary/(Real+J));
    }
    GR-=GR1;
    GI-=GI1;
  }

  if (Real_1 < 0.0) {
    Z1 = sqrt(Real*Real+Imaginary*Imaginary);
    TH1 = atan(Imaginary/Real);
    SR = -sin(M_PI*Real)*cosh(M_PI*Imaginary);
    SI = -cos(M_PI*Real)*sinh(M_PI*Imaginary);
    Z2 = sqrt(SR*SR+SI*SI);
    TH2 = atan(SI/SR);

    if (SR < 0.0) {
      TH2 += M_PI;
    }

    GR = log(M_PI/(Z1*Z2))-GR;
    GI = -TH1-TH2-GI;
    Real = Real_1;
    Imaginary = Imaginary_1;
  }

  if (choice == 'g') {
    G0 = exp(GR);
    GR = G0*cos(GI);
    GI = G0*sin(GI);
  }

  return GR;
}

double K37ComplexGammaFunction::absSquaredComplexGamma(char choice, double Real,
                                                       double Imaginary) {
  if ( (choice != 'l') && (choice != 'g') ) {
    return 4.23e17;
  }

  if (Imaginary == 0.0 && (Real == static_cast<int>(Real)) && Real < 0.0) {
    return 4.23e17;
  } else if (Real < 0.0) {
    Real_1 = Real;
    Imaginary_1 = Imaginary;
    Real = -Real;
    Imaginary = -Imaginary;
  }

  Real_0 = Real;
  if (Real < 7.0) {
    NA = static_cast<int>(7.0 - Real);
    Real_0 = Real + NA;
  }
  Z1 = sqrt(std::pow(Real_0, 2.0)+std::pow(Imaginary, 2.0));
  TH  = std::atan(Imaginary/Real_0);

  GR = (Real_0-0.5)*log(Z1)-TH*Imaginary-Real_0+0.5*log(2.0*M_PI);
  GI = TH*(Real_0-0.5)+Imaginary*log(Z1)-Imaginary;

  for (int K =1; K <=10; ++K) {
    T = pow(Z1, (1.0-2.0*K));
    GR = GR+coefficients[K-1]*T*cos((2.0*K-1.0)*TH);
    GI = GI-coefficients[K-1]*T*sin((2.0*K-1.0)*TH);
  }

  if (Real < 7.0) {
    GR1 = 0.0;
    GI1 = 0.0;
    for (int J = 0; J <= (NA-1); ++J) {
      GR1 += 0.5*log(pow((Real+J), 2.0)+pow(Imaginary, 2.0));
      GI1 += atan(Imaginary/(Real+J));
    }
    GR-=GR1;
    GI-=GI1;
  }


  if (Real_1 < 0.0) {
    Z1 = sqrt(Real*Real+Imaginary*Imaginary);
    TH1 = atan(Imaginary/Real);
    SR = -sin(M_PI*Real)*cosh(M_PI*Imaginary);
    SI = -cos(M_PI*Real)*sinh(M_PI*Imaginary);
    Z2 = sqrt(SR*SR+SI*SI);
    TH2 = atan(SI/SR);

    if (SR < 0.0) {
      TH2+=M_PI;
    }

    GR = log(M_PI/(Z1*Z2))-GR;
    GI = -TH1-TH2-GI;
    Real = Real_1;
    Imaginary = Imaginary_1;
  }

  if (choice == 'g') {
    G0 = exp(GR);
    GR = G0*cos(GI);
    GI = G0*sin(GI);
  }

  return (GR+GI)*(GR-GI);
}

G4double K37ComplexGammaFunction::realGamma(G4double x) {
/* Gamma function in double precision */
  G4int k, n;
  G4double w, y;

  n = x < 1.5 ? -((G4int) (2.5 - x)) : (G4int) (x - 1.5);
  w = x - (n + 2);
  y = ((((((((((((-1.99542863674e-7 * w + 1.337767384067e-6) * w -
                 2.591225267689e-6) * w - 1.7545539395205e-5) * w +
               1.45596568617526e-4) * w - 3.60837876648255e-4) * w -
             8.04329819255744e-4) * w + 0.008023273027855346) * w -
           0.017645244547851414) * w - 0.024552490005641278) * w +
         0.19109110138763841) * w - 0.233093736421782878) * w -
       0.422784335098466784) * w + 0.99999999999999999;
  if (n > 0) {
    w = x - 1;

    for (k= 2; k <= n; k++) {
      w *= x - k;
    }
  } else {
    w = 1;
    for (k = 0; k > n; k--) {
      y *= x - k;
    }
  }
  return w / y;
}

G4double K37ComplexGammaFunction::squaredRealGamma(G4double x) {
/* Gamma function in double precision */
  G4int k, n;
  G4double w, y;

  n = x < 1.5 ? -((G4int) (2.5 - x)) : (G4int) (x - 1.5);
  w = x - (n + 2);
  y = ((((((((((((-1.99542863674e-7 * w + 1.337767384067e-6) * w -
                 2.591225267689e-6) * w - 1.7545539395205e-5) * w +
               1.45596568617526e-4) * w - 3.60837876648255e-4) * w -
             8.04329819255744e-4) * w + 0.008023273027855346) * w -
           0.017645244547851414) * w - 0.024552490005641278) * w +
         0.19109110138763841) * w - 0.233093736421782878) * w -
       0.422784335098466784) * w + 0.99999999999999999;
  if (n > 0) {
    w = x - 1;

    for (k= 2; k <= n; k++) {
      w *= x - k;
    }
  } else {
    w = 1;
    for (k = 0; k > n; k--) {
      y *= x - k;
    }
  }
  return (w/y)*(w/y);
}

