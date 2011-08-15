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

        void SetFunction1 ( Matrix *func );
        void SetFunction2 ( Matrix *func );
        void SetP1 ( Matrix *p );
        void SetP2 ( Matrix *p );
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
        Matrix *function1, *function2;
        //Matrix *amp1, *amp2;
        Matrix *P1, *P2;
        void makeSrez(double deb1, double fin1, double deb2, double fin2) ;
        void print();
};

Matrix* GetFunctionSrez(Matrix*X, Matrix *Y, double pa, double pf, int *ia);

Matrix* GetFunctionSrezDeb(Matrix*X, Matrix *Y, double pa);

Matrix* GetFunctionSrezDebFin(Matrix*X, Matrix *Y, double pa, double pf);

double pinceFunctionA(WindPatternsProject* gfd, double x);

double pinceFunctionF(WindPatternsProject* gfd, double x);

void calcWidth(Matrix *Xd, Matrix *Yd, Matrix *Xd0, Matrix *Yd0, double *width);

double calcW0byHW1(double h, double w1);

double calcHbyw0w1(double w0, double w1);

double functionPincePower(double x, double power);

double functionPinceArctanSym(double x, double k1, double k2);

double functionPinceArctan(double x, double k1, double k2, double k3);

void makePointPince(Matrix *Xd, Matrix *Yd, Matrix *P,
        double PosPinceBA, double PosPinceBF, Matrix **Xr, Matrix **Yr, Matrix **Pr, int *ipair, int *ipfir);

void calcPinceAloneAbs(WindPatternsProject* gfd, Matrix *Xd, Matrix *Yd,
        Matrix *Xd0, Matrix *Yd0,
        Matrix *X, Matrix *Y, Matrix *Z, Matrix *P,
        Matrix *X0, Matrix *Y0, Matrix *Z0, Matrix *P0,
        double PosPinceBA, double AmpPinceBA, double PosPinceBF, double AmpPinceBF, double modePincesT,
        Matrix **Xdp, Matrix **Ydp, Matrix **Pp);

void calcPinceAlone(WindPatternsProject* gfd,Matrix *Xd, Matrix *Yd, Matrix *Xd0, Matrix *Yd0, Matrix *X, Matrix *Y, Matrix *Z, Matrix *P, Matrix *X0, Matrix *Y0, Matrix *Z0, Matrix *P0,
        double PosPinceBA, double AmpPinceBA, double PosPinceBF, double AmpPinceBF, double modePincesT,
        Matrix **Xdp, Matrix **Ydp, Matrix **Pp);


void calcPinceAloneAbsFunc(Matrix *Xd, Matrix *Yd,
                             Matrix *Xd0, Matrix *Yd0,
                            Matrix *X, Matrix *Y, Matrix *Z, Matrix *P,
        Matrix *X0, Matrix *Y0, Matrix *Z0, Matrix *P0, Pince* p, double ampa, double ampf, int no,
        Matrix **Xdp, Matrix **Ydp, Matrix **Pp);

void calcPinceNew(Matrix *Xd1, Matrix *Yd1, Matrix *Xd2, Matrix *Yd2,
        Matrix *X1, Matrix *Y1, Matrix *Z1, Matrix *P1,
        Matrix *X2, Matrix *Y2, Matrix *Z2, Matrix *P2, Pince* p,
        Matrix **Xdp1, Matrix **Ydp1, Matrix **Pp1,
        Matrix **Xdp2, Matrix **Ydp2, Matrix **Pp2);

void calcPince(WindPatternsProject* gfd, Matrix *Xd1, Matrix *Yd1, Matrix *Xd2, Matrix *Yd2,
        Matrix *X1, Matrix *Y1, Matrix *Z1, Matrix *P1,
        Matrix *X2, Matrix *Y2, Matrix *Z2, Matrix *P2,
        double PosPinceBA, double AmpPinceBA, double PosPinceBF, double AmpPinceBF, double modePincesT,
        Matrix **Xdp1, Matrix **Ydp1, Matrix **Pp1, Matrix **Xdp2, Matrix **Ydp2, Matrix **Pp2);

void calcPinceDiff(WindPatternsProject* gfd, Matrix *Xd1, Matrix *Yd1, Matrix *Xd2, Matrix *Yd2,
        Matrix *X1, Matrix *Y1, Matrix *Z1, Matrix *P1,
        Matrix *X2, Matrix *Y2, Matrix *Z2, Matrix *P2,
        double PosPinceBA1, double AmpPinceBA1, double PosPinceBF1, double AmpPinceBF1,
        double PosPinceBA2, double AmpPinceBA2, double PosPinceBF2, double AmpPinceBF2,
        double mode,
        Matrix **Xdp1, Matrix **Ydp1, Matrix **Pp1, Matrix **Xdp2, Matrix **Ydp2, Matrix **Pp2);

Pince* getPince(WindPatternsProject* gfd, int nerv1, int nerv2, int face);

void getPinceFunctions(WindPatternsProject* gfd, int nerv1, int nerv2, int face, int* ipa, int* ipf, Matrix** Pout,
                        Matrix** FuncFunction1,Matrix** FuncFunction2, double* ampFfA, double* ampFfF,
                        Matrix** RadiusFunction, double* ampRfA, double* ampRfF);

void getPinceRadiusFunction(WindPatternsProject* gfd, int nerv1, int nerv2, int face);

void goCalcPinceLen(WindPatternsProject* gfd, int noNerv, int face, double *len);

Pince* getPincePlus(WindPatternsProject* gfd,
                        int nerv1, int nerv2,
                        float deb1, float fin1, float deb2, float fin2, int faceDeb, int faceFin );

void calcPincePlusNew(Matrix *Xd1, Matrix *Yd1, Matrix *Xd2, Matrix *Yd2,
                            Pince* p,
                            Matrix **Xdp1, Matrix **Ydp1,
                            Matrix **Xdp2, Matrix **Ydp2);

void calcPincePlusAloneAbsFunc(Matrix *Xd, Matrix *Yd,
                                Matrix *Xd0, Matrix *Yd0, Pince* p, int no,
                                Matrix **Xdp, Matrix **Ydp);

