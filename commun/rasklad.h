#pragma once

#include "afx.h"		//class CString
#include "afxdlgs.h"	//class CFileDialog

#include "profil.h"
#include "pince.h"
#include "matrice.h"
#include "fichier.h"

#define HACCUR 0.005f

class WindPatternsProject;

class Rasklad
{
    public:
        Rasklad ()
        {

        }
        Rasklad (int n);
        virtual ~Rasklad();

        double *pinceLA0Amp1;
        double *pinceLAAmp1;
        double *pinceLFAmp1;

        double *pinceRA0Amp1;
        double *pinceRAAmp1;
        double *pinceRFAmp1;

        Matrice** funcL1;
        Matrice** funcL2;
        Matrice** funcR1;
        Matrice** funcR2;

        double *lenp1;

        double *pinceLA0Amp2;
        double *pinceLAAmp2;
        double *pinceLFAmp2;

        double *pinceRA0Amp2;
        double *pinceRAAmp2;
        double *pinceRFAmp2;
        double *lenp2;
        double *coeffn;
        double *coeffd;
        int quantDiag;
        int *noNervD;
        int *noNervD1;
        int *noNervD2;
        int nrp;
        bool isCenterPanel;
};

const int REP_TRIANGLE = 1;
const int REP_LINE = 2;
const int REP_CROSS = 3;
const int REP_0 = 4;
const int REP_VENTDEB = 5;
const int REP_VENTFIN = 6;
const int REP_KLAPAN_FIN = 7;
const int REP_V = 8;
const int REP_MIDDLE_LINE = 9;

Rasklad* CalculIndepPinceRasklad(WindPatternsProject* gfd, Forme* F);

void SaveRasklad2(WindPatternsProject* gfd, Rasklad* rasklad);

void CalculPatronPosNerv(WindPatternsProject* gfd, int noNerv1, bool sym1, int FaceDeb1, int FaceFin1, double Deb1, double Fin1,
        int noNerv2, bool sym2, int FaceDeb2, int FaceFin2, double Deb2, double Fin2,
        Matrice **Xdev1, Matrice **Ydev1, Matrice **Xdev2, Matrice **Ydev2,
        Matrice **X0, Matrice **Y0, Matrice **Z0, Matrice **P0,
        Matrice **X1, Matrice **Y1, Matrice **Z1, Matrice **P1);
void CalculPatronKlapan(WindPatternsProject* gfd, int noNerv1, bool sym1, int noNerv2, bool sym2, double posKlapanInt, double posKlapanFin,
        Matrice **Xdev1, Matrice **Ydev1, Matrice **Xdev2, Matrice **Ydev2,
        Matrice **X0, Matrice **Y0, Matrice **Z0, Matrice **P0,
        Matrice **X1, Matrice **Y1, Matrice **Z1, Matrice **P1);

void CalculPatron(WindPatternsProject* gfd, int noNerv1, bool sym1, int FaceDeb1, int FaceFin1, double Deb1, double Fin1,
        int noNerv2, bool sym2, int FaceDeb2, int FaceFin2, double Deb2, double Fin2,
        Matrice **Xdev1, Matrice **Ydev1, Matrice **Xdev2, Matrice **Ydev2,
        Matrice **X0, Matrice **Y0, Matrice **Z0, Matrice **P0,
        Matrice **X1, Matrice **Y1, Matrice **Z1, Matrice **P1);

double CalculWidthNervs(WindPatternsProject* gfd, int noNerv1, int noNerv2, int face);

void GetMiddleProfile (WindPatternsProject* gfd, Forme* F, int nerv1, int nerv2, int face, Matrice** XProf, Matrice** YProf);

void GetMiddleProfileBal(WindPatternsProject* gfd, Forme* F, int nerv1, int nerv2, int face, Matrice** XProf, Matrice** YProf);

ProfilGeom* getProfile(WindPatternsProject* gfd, Forme* F, int nerv);

void goCalcIndepPinceNew(WindPatternsProject* gfd, int noNerv, int face, double *pLA, double *pLF, Matrice** fl, double *pRA, double *pRF, Matrice** fr, double *len);

Matrice* getReperPoints(WindPatternsProject* gfd, Matrice* Xd, Matrice* Yd, Matrice* P, 
						int nerv, float deb, int faceDeb, float fin, int faceFin, bool klapanShov);


Matrice* getReperPointsPince(WindPatternsProject* gfd,  Matrice* Xd0, Matrice* Yd0, double coeff, Matrice* Xd, Matrice* Yd, Matrice* P, 
								int nerv, float deb, int faceDeb, float fin, int faceFin, bool klapanShov);


void makePosProfile(Matrice* Ext, Matrice* Int, double percent, Matrice** ExtRes);

void goCalcNervureWithPince(WindPatternsProject* gfd, int noNerv1, int face1, int noNerv2, int face2, double lp1, double lp2, double *rcoeff);

void goCalcIndepPince(int noNerv, int face, double *pLA, double *pLF, double *pRA, double *pRF, double *len);

double calculPolyLength(Matrice* X, Matrice* Y, double xrp, double yrp);

int calcRepPointByLength(Matrice* Xd, Matrice* Yd, double l, double* x, double* y);

void printM(Matrice* m);