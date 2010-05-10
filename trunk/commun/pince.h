#pragma once

#include "matrice.h"
#include "fichier.h"
#include "patternsproject.h"

//#include <afx.h>		//class CString
//#include <afxdlgs.h>	//class CFileDialog

class WindPatternsProject;
class Pince
{
    public:
        Pince(int n1, int n2);
        Pince()
        {
                diffAmps=false; debug=false;
                glue = false;
        }
        virtual ~Pince();

        void SetFunction1 ( Matrice *func );
        void SetFunction2 ( Matrice *func );
        void SetP1 ( Matrice *p );
        void SetP2 ( Matrice *p );
        void SetPos (double pa, double pf);

        void SetAmps (double aa, int _ia,  double af, int _if);
        
        void SetDiffAmps (double a1, double f1, double a2, double f2 );
        void SetDiffAmps0 (double a1, double a2);
        int nerv1, nerv2, ipa1, ipf1, i01, ipa2, ipf2, i02;

        double AmpA0, AmpF0, AmpA, AmpF, PosA, PosF;

        double AmpA01, AmpA02, AmpA1, AmpF1, AmpA2, AmpF2;
        bool diffAmps;
        bool debug;
        bool glue;
        Matrice *function1, *function2;
        //Matrice *amp1, *amp2;
        Matrice *P1, *P2;
        void makeSrez(double deb1, double fin1, double deb2, double fin2) ;
        void print();
};

Matrice* GetFunctionSrez(Matrice*X, Matrice *Y, double pa, double pf, int *ia);

Matrice* GetFunctionSrezDeb(Matrice*X, Matrice *Y, double pa);

double pinceFunctionA(WindPatternsProject* gfd, double x);

double pinceFunctionF(WindPatternsProject* gfd, double x);

void CalculWidth(Matrice *Xd, Matrice *Yd, Matrice *Xd0, Matrice *Yd0, double *width);

double CalculW0byHW1(double h, double w1);

double CalculHbyw0w1(double w0, double w1);

double functionPincePower(double x, double power);

double functionPinceArctanSym(double x, double k1, double k2);

double functionPinceArctan(double x, double k1, double k2, double k3);

void makePointPince(Matrice *Xd, Matrice *Yd, Matrice *P,
        double PosPinceBA, double PosPinceBF, Matrice **Xr, Matrice **Yr, Matrice **Pr, int *ipair, int *ipfir);

void CalculPinceAloneAbs(WindPatternsProject* gfd, Matrice *Xd, Matrice *Yd,
        Matrice *Xd0, Matrice *Yd0,
        Matrice *X, Matrice *Y, Matrice *Z, Matrice *P,
        Matrice *X0, Matrice *Y0, Matrice *Z0, Matrice *P0,
        double PosPinceBA, double AmpPinceBA, double PosPinceBF, double AmpPinceBF, double modePincesT,
        Matrice **Xdp, Matrice **Ydp, Matrice **Pp);

void CalculPinceAlone(WindPatternsProject* gfd,Matrice *Xd, Matrice *Yd, Matrice *Xd0, Matrice *Yd0, Matrice *X, Matrice *Y, Matrice *Z, Matrice *P, Matrice *X0, Matrice *Y0, Matrice *Z0, Matrice *P0,
        double PosPinceBA, double AmpPinceBA, double PosPinceBF, double AmpPinceBF, double modePincesT,
        Matrice **Xdp, Matrice **Ydp, Matrice **Pp);


void CalculPinceAloneAbsFunc(Matrice *Xd, Matrice *Yd,
                             Matrice *Xd0, Matrice *Yd0,
                            Matrice *X, Matrice *Y, Matrice *Z, Matrice *P,
        Matrice *X0, Matrice *Y0, Matrice *Z0, Matrice *P0, Pince* p, double ampa, double ampf, int no,
        Matrice **Xdp, Matrice **Ydp, Matrice **Pp);

void CalculPinceNew(Matrice *Xd1, Matrice *Yd1, Matrice *Xd2, Matrice *Yd2,
        Matrice *X1, Matrice *Y1, Matrice *Z1, Matrice *P1,
        Matrice *X2, Matrice *Y2, Matrice *Z2, Matrice *P2, Pince* p,
        Matrice **Xdp1, Matrice **Ydp1, Matrice **Pp1,
        Matrice **Xdp2, Matrice **Ydp2, Matrice **Pp2);

void CalculPince(WindPatternsProject* gfd, Matrice *Xd1, Matrice *Yd1, Matrice *Xd2, Matrice *Yd2,
        Matrice *X1, Matrice *Y1, Matrice *Z1, Matrice *P1,
        Matrice *X2, Matrice *Y2, Matrice *Z2, Matrice *P2,
        double PosPinceBA, double AmpPinceBA, double PosPinceBF, double AmpPinceBF, double modePincesT,
        Matrice **Xdp1, Matrice **Ydp1, Matrice **Pp1, Matrice **Xdp2, Matrice **Ydp2, Matrice **Pp2);

void CalculPinceDiff(WindPatternsProject* gfd, Matrice *Xd1, Matrice *Yd1, Matrice *Xd2, Matrice *Yd2,
        Matrice *X1, Matrice *Y1, Matrice *Z1, Matrice *P1,
        Matrice *X2, Matrice *Y2, Matrice *Z2, Matrice *P2,
        double PosPinceBA1, double AmpPinceBA1, double PosPinceBF1, double AmpPinceBF1,
        double PosPinceBA2, double AmpPinceBA2, double PosPinceBF2, double AmpPinceBF2,
        double mode,
        Matrice **Xdp1, Matrice **Ydp1, Matrice **Pp1, Matrice **Xdp2, Matrice **Ydp2, Matrice **Pp2);

Pince* getPince(WindPatternsProject* gfd, int nerv1, int nerv2, int face);

void getPinceFunctions(WindPatternsProject* gfd, int nerv1, int nerv2, int face, int* ipa, int* ipf, Matrice** Pout,
                        Matrice** FuncFunction1,Matrice** FuncFunction2, double* ampFfA, double* ampFfF,
                        Matrice** RadiusFunction, double* ampRfA, double* ampRfF);

void getPinceRadiusFunction(WindPatternsProject* gfd, int nerv1, int nerv2, int face);

void goCalcPinceLen(WindPatternsProject* gfd, int noNerv, int face, double *len);

Pince* getPincePlus(WindPatternsProject* gfd,
                        int nerv1, int nerv2,
                        float deb1, float fin1, float deb2, float fin2, int faceDeb, int faceFin );

void CalculPincePlusNew(Matrice *Xd1, Matrice *Yd1, Matrice *Xd2, Matrice *Yd2,
                            Pince* p,
                            Matrice **Xdp1, Matrice **Ydp1,
                            Matrice **Xdp2, Matrice **Ydp2);

void CalculPincePlusAloneAbsFunc(Matrice *Xd, Matrice *Yd,
                                Matrice *Xd0, Matrice *Yd0, Pince* p, int no,
                                Matrice **Xdp, Matrice **Ydp);

