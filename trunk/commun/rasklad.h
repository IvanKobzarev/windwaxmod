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

        Matrix** funcL1;
        Matrix** funcL2;
        Matrix** funcR1;
        Matrix** funcR2;

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

Rasklad* calcIndepPinceRasklad(WindPatternsProject* gfd, Form* F);

void SaveRasklad2(WindPatternsProject* gfd, Rasklad* rasklad);

void calcPatronPosNerv(WindPatternsProject* gfd, int noNerv1, bool sym1, int FaceDeb1, int FaceFin1, double Deb1, double Fin1,
        int noNerv2, bool sym2, int FaceDeb2, int FaceFin2, double Deb2, double Fin2,
        Matrix **Xdev1, Matrix **Ydev1, Matrix **Xdev2, Matrix **Ydev2,
        Matrix **X0, Matrix **Y0, Matrix **Z0, Matrix **P0,
        Matrix **X1, Matrix **Y1, Matrix **Z1, Matrix **P1);
void calcPatronKlapan(WindPatternsProject* gfd, int noNerv1, bool sym1, int noNerv2, bool sym2, double posKlapanInt, double posKlapanFin,
        Matrix **Xdev1, Matrix **Ydev1, Matrix **Xdev2, Matrix **Ydev2,
        Matrix **X0, Matrix **Y0, Matrix **Z0, Matrix **P0,
        Matrix **X1, Matrix **Y1, Matrix **Z1, Matrix **P1);

void calcPatron(WindPatternsProject* gfd, int noNerv1, bool sym1, int FaceDeb1, int FaceFin1, double Deb1, double Fin1,
        int noNerv2, bool sym2, int FaceDeb2, int FaceFin2, double Deb2, double Fin2,
        Matrix **Xdev1, Matrix **Ydev1, Matrix **Xdev2, Matrix **Ydev2,
        Matrix **X0, Matrix **Y0, Matrix **Z0, Matrix **P0,
        Matrix **X1, Matrix **Y1, Matrix **Z1, Matrix **P1);

double calcWidthNervs(WindPatternsProject* gfd, int noNerv1, int noNerv2, int face);

void GetMiddleProfile (WindPatternsProject* gfd, Form* F, int nerv1, int nerv2, int face, int realMashtab, Matrix** XProf, Matrix** YProf);

void GetMiddleProfileBal(WindPatternsProject* gfd, Form* F, int nerv1, int nerv2, int face,  int realMashtab, Matrix** XProf, Matrix** YProf);

ProfilGeom* getProfile(WindPatternsProject* gfd, Form* F, int nerv);

void goCalcIndepPinceNew(WindPatternsProject* gfd, int noNerv, int face, double *pLA, double *pLF, Matrix** fl, double *pRA, double *pRF, Matrix** fr, double *len);

Matrix* getReperPoints(WindPatternsProject* gfd, Matrix* Xd, Matrix* Yd, Matrix* P, 
						int nerv, float deb, int faceDeb, float fin, int faceFin, bool klapanShov);


Matrix* getReperPointsPince(WindPatternsProject* gfd,  Matrix* Xd0, Matrix* Yd0, double coeff, Matrix* Xd, Matrix* Yd, Matrix* P, 
								int nerv, float deb, int faceDeb, float fin, int faceFin, bool klapanShov);


void makePosProfile(Matrix* Ext, Matrix* Int, double percent, Matrix** ExtRes);

void goCalcNervureWithPince(WindPatternsProject* gfd, int noNerv1, int face1, int noNerv2, int face2, double lp1, double lp2, double *rcoeff);

void goCalcIndepPince(int noNerv, int face, double *pLA, double *pLF, double *pRA, double *pRF, double *len);

double calcPolyLength(Matrix* X, Matrix* Y, double xrp, double yrp);

int calcRepPointByLength(Matrix* Xd, Matrix* Yd, double l, double* x, double* y);

void printM(Matrix* m);