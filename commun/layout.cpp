#pragma warning(disable:4514)

#include "afx.h"		//class CString
#include "afxdlgs.h"	//class CFileDialog

#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>

#include "layout.h"
#include "dxf.h"
#include "geom.h"
#include "fichier.h"
#include "plot.h"
#include "pince.h"
#include "profil.h"
#include "matrice.h"
#include "logger.h"
#include "patron.h"

#define sqr(f1) ((f1)*(f1))
#define pi	3.141592675f

#ifndef DEBUG
    #define DEBUG false
#endif

void calculateLayout(WindPatternsProject* gfd, Layout* layout);

LayoutElement::~LayoutElement() {
}

LayoutElement::LayoutElement() {
    isVentHole = 0;
}

LayoutElementExport::LayoutElementExport(){
}
LayoutElementExport::~LayoutElementExport(){
}

KlapanLayoutElement::KlapanLayoutElement(){
    isVentHole = 0;
}
KlapanLayoutElement::~KlapanLayoutElement(){
}

ProfLayoutElement::ProfLayoutElement(){
    isVentHole = 0;
}

ProfLayoutElement::~ProfLayoutElement(){
}

PanelLayoutElement::PanelLayoutElement(){
    isVentHole = 0;
}
PanelLayoutElement::~PanelLayoutElement(){
}

DiagNervLayoutElement::DiagNervLayoutElement(){
    isVentHole = 0;
}
DiagNervLayoutElement::~DiagNervLayoutElement(){
}



Quad::Quad(double _pd1, double _pf1, double _pd2, double _pf2){
	pd1 = _pd1;
	pf1 = _pf1;
	pd2 = _pd2;
	pf2 = _pf2;
}

void KlapanLayoutElement::calculateExport(WindPatternsProject* gfd){
	bool debug = false;
    if (debug) printf ("\n KlapanLayoutElement::calculateExport");

    Matrix * Xd[2], *Yd[2];//,*newXd[2], *newYd[2], *Xdp[2], *Ydp[2], *rXdp[2], *rYdp[2], *rXd[2], *rYd[2];
    Matrix *X[2], *Y[2], *Z[2], *P[2];//, *rP[2];//, *newP[2];
    char* charname = new char[1];

    calcPatronKlapan(gfd, n1, s1, n2, s2, posKlapanIntDeb, posKlapanFin,
        &Xd[0], &Yd[0], &Xd[1], &Yd[1],
        &X[0], &Y[0], &Z[0], &P[0],
        &X[1], &Y[1], &Z[1], &P[1]);

    charname="K";
    double marge1 = gfd->Marge[0];
    double marge2 = gfd->Marge[1];
    double margeDeb = gfd->MargeDeb;
    double margeFin = gfd->MargeFin;
    int n = gfd->Form->m_nbProfils;
    char text[100];
    sprintf(text, "%s%dF%dt%dF%d", charname, n1, ff1, n2, ff2);
    leexport = new LayoutElementExport();
    GenerateCourbe(gfd, Xd[0], Yd[0], P[0], n1, posDeb1, fd1, posFin1, ff1,
        Xd[1], Yd[1], P[1], n2, posDeb2, fd2, posFin2, ff2, text,
        &(leexport->AxeP), &(leexport->AxePD), &(leexport->AxePTD), &(leexport->AxeMD), &(leexport->AxeCD), &(leexport->AxeRepD), vent, marge1, marge2, margeDeb, margeFin, false, debug);
    calcMaxWH(Xd[0], Yd[0], Xd[1], Yd[1], &(leexport->W), &(leexport->H));
    delete (X[0]);
    delete (Y[0]);
    delete (Z[0]);
    delete (P[0]);
    delete (X[1]);
    delete (Y[1]);
    delete (Z[1]);
    delete (P[1]);
    delete (Xd[0]);
    delete (Yd[0]);
    delete (Xd[1]);
    delete (Yd[1]);
}

void ProfLayoutElement::calculateExport(WindPatternsProject* gfd){
	bool debug = false;
    if (debug) printf ("\n ProfLayoutElement::calculateExport");

    Matrix * Xd[2], *Yd[2],*newXd[2], *newYd[2];//, *Xdp[2], *Ydp[2];//, *rXdp[2], *rYdp[2], *rXd[2], *rYd[2];
    Matrix *X[2], *Y[2], *Z[2], *P[2];//, *rP[2];//, *newP[2];
    char* charname = new char[1];

    calcPatron(gfd, n1, s1, fd1, ff1, posDeb1, posFin1,
        n2, s2, fd2, ff2, posDeb2, posFin2,
        &Xd[0], &Yd[0], &Xd[1], &Yd[1],
        &X[0], &Y[0], &Z[0], &P[0],
        &X[1], &Y[1], &Z[1], &P[1]);

    getLayoutLogger()->logprintf(" c[%f]", coeff);
    getLayoutLogger()->logprintf(" %d %d ff1=%d ff2=%d", n1, n2,  ff1, ff2);

    if (coeff != 0) {
        calcPatronWithCoeff(Xd[0], Yd[0], Xd[1], Yd[1], coeff, &newXd[0], &newYd[0], &newXd[1], &newYd[1]);
        delete (Xd[0]);
        delete (Yd[0]);
        delete(Xd[1]);
        delete (Yd[1]);
        Xd[0] = newXd[0];
        Yd[0] = newYd[0];
        Xd[1] = newXd[1];
        Yd[1] = newYd[1];
    }
    double marge1 = gfd->Marge[0];
    double marge2 = gfd->Marge[1];
    double margeDeb = gfd->MargeDeb;
    double margeFin = gfd->MargeFin;
    // nervures
    margeFin = gfd->margeFinNerv;
    vent = gfd->VentilationLayout;
    charname="N";
    int n = gfd->Form->m_nbProfils;
    char text[100];
    sprintf(text, "%s%dF%dt%dF%d", charname, n1, ff1, n2, ff2);
    leexport = new LayoutElementExport();
    GenerateCourbe(gfd, Xd[0], Yd[0], P[0], n1, posDeb1, fd1, posFin1, ff1,
        Xd[1], Yd[1], P[1], n2, posDeb2, fd2, posFin2, ff2, text,
        &(leexport->AxeP), &(leexport->AxePD), &(leexport->AxePTD), &(leexport->AxeMD), &(leexport->AxeCD), &(leexport->AxeRepD), vent, marge1, marge2, margeDeb, margeFin, true, debug);

    calcMaxWH(Xd[0], Yd[0], Xd[1], Yd[1], &(leexport->W), &(leexport->H));

    delete (X[0]);
    delete (Y[0]);
    delete (Z[0]);
    delete (P[0]);
    delete (X[1]);
    delete (Y[1]);
    delete (Z[1]);
    delete (P[1]);
    delete (Xd[0]);
    delete (Yd[0]);
    delete (Xd[1]);
    delete (Yd[1]);
}
void PanelLayoutElement::calculateExport(WindPatternsProject* gfd) {
    bool debug = false;

	if (debug) printf ("\n PanelLayoutElement::calculateExport");
    Matrix * Xd[2], *Yd[2],*newXd[2], *newYd[2], *Xdp[2], *Ydp[2], *rXdp[2], *rYdp[2], *rXd[2], *rYd[2];
    Matrix *X[2], *Y[2], *Z[2], *P[2], *rP[2];//, *newP[2];
    double _pa0, _pa00, _pf0, _pa1, _pa01, _pf1;
    Matrix *_f0, *_f1;
    char* charname = new char[1];
	if (debug) printf ("\n call calcPatron()");
    calcPatron(gfd, n1, s1, fd1, ff1, posDeb1, posFin1,
        n2, s2, fd2, ff2, posDeb2, posFin2,
        &Xd[0], &Yd[0], &Xd[1], &Yd[1],
        &X[0], &Y[0], &Z[0], &P[0],
        &X[1], &Y[1], &Z[1], &P[1]);
	if (debug) printf ("\n ...call calcPatron()");
    if (debug) printf ("\n n1=%d fd1=%d ff1=%d posDeb1=%f posFin1=%f", n1, fd1, ff1, posDeb1, posFin1);
    if (debug) printf ("\n n2=%d fd2=%d ff2=%d posDeb2=%f posFin2=%f", n2, fd2, ff2, posDeb2, posFin2);
    if (isPince) {
        if (ff1 == 1) {
            _pa00 = p1a00;
            _pa0 = p1a0;
            _pf0 = p1f0;

            _pa1 = p1a1;
            _pa01 = p1a01;
            _pf1 = p1f1;

            _f0 = func1f0;
            _f1 = func1f1;
        } else {
            _pa0 = p2a0;
            _pa00 = p2a00;
            _pf0 = p2f0;

            _pa1 = p2a1;
            _pa01 = p2a01;
            _pf1 = p2f1;

            _f0 = func2f0;
            _f1 = func2f1;
        }
        delete (X[0]);
        delete (Y[0]);
        delete (Z[0]);
        delete (P[0]);
        delete (X[1]);
        delete (Y[1]);
        delete (Z[1]);
        delete (P[1]);
        delete (Xd[0]);
        delete (Yd[0]);
        delete (Xd[1]);
        delete (Yd[1]);
        double myDeb = posDeb1;
		if (debug) printf ("\n posDeb1=%f posDeb2=%f ---> myDeb=%f", posDeb1, posDeb2, myDeb);

		double myFin = posFin1;
		if (posFin1 < posFin2) myFin = posFin2; else myFin = posFin1;

		if (debug) printf ("\n posFin1=%f posFin2=%f ---> myFin=%f", posFin1, posFin2, myFin);
         if (fd1 == ff1) {
            if (posDeb1 < posDeb2) myDeb = posDeb1; else myDeb = posDeb2;
        } else {
            if (posDeb1 > posDeb2) myDeb = posDeb1; else myDeb = posDeb2;
        }
        if (debug) printf ("\n calcPatron myDeb myFin");
		if (debug) printf ("\n ce3");
		if (debug) printf ("calcPatron(%d: %f(%d) %f(%d))-(%d: %f(%d) %f(%d)) ", n1, myDeb, fd1, myFin, ff1, n2, myDeb, fd2, myFin, ff2);
        calcPatron(gfd, n1, s1, fd1, ff1, myDeb, myFin,
                n2, s2, fd2, ff2, myDeb, myFin,
                &Xd[0], &Yd[0], &Xd[1], &Yd[1],
                &X[0], &Y[0], &Z[0], &P[0],
                &X[1], &Y[1], &Z[1], &P[1]);

		if (debug) printf ("\n ce4");
		if (debug) printf ("\n Xd[0]->GetLignes(): %d, Yd[0]->GetLignes(): %d", Xd[0]->GetLignes(), Yd[0]->GetLignes());
		if (debug) printf ("\n Xd[1]->GetLignes(): %d, Yd[1]->GetLignes(): %d", Xd[1]->GetLignes(), Yd[1]->GetLignes());
		if (debug) printf ("\n P[0]->GetLignes(): %d,  P[1]->GetLignes(): %d", P[0]->GetLignes(),  P[1]->GetLignes());
		if (debug) printf ("\n ********************************************************** ");
        if (debug) printf ("\n...calcPatron myDeb myFin");
        getLayoutLogger()->logprintf("\n p[%4.1f %4.1f %4.1f %4.1f]", _pa0 * 1000, _pf0 * 1000, _pa1 * 1000, _pf1 * 1000);
        getLayoutLogger()->logprintf(" \n %d %d [%5.1f(%d) %5.1f(%d)] [%5.1f(%d) %5.1f(%d)]", n1, n2, posDeb1, fd1, posFin1,  ff1, posDeb2, fd2, posFin2,  ff2);

        Pince* pince = new Pince();

        pince -> debug = false;
		//pince -> debug = true;

		if (debug) printf ("\n getPincePlus.1 (%d, %d, %f, %f, %f, %f, %d, %d)", n1, n2, myDeb, myFin,  myDeb, myFin,  fd1, ff1);
	
        Pince* tmpPince = getPincePlus(gfd, n1, n2, myDeb, myFin,  myDeb, myFin,  fd1, ff1);
		if (debug) printf ("\n ********************************************************** ");
		if (debug) printf ("\n tmpPince->function1->GetLignes(): %d    tmpPince->function2->GetLignes(): %d  ", tmpPince->function1->GetLignes(), tmpPince->function2->GetLignes());
		if (debug) printf ("\n getPincePlus.2");
        pince -> SetFunction1(tmpPince->function1);
        pince -> SetFunction2(tmpPince->function2);
		if (debug) printf ("\n ... - ");
        pince -> ipa1 = tmpPince -> ipa1;
        pince -> ipf1 = tmpPince -> ipf1;
        pince -> ipa2 = tmpPince -> ipa2;
        pince -> ipf2 = tmpPince -> ipf2;
        pince -> SetPos(gfd->PosPinceBA[0], gfd->PosPinceBF[0]);
        pince -> SetDiffAmps(_pa0, _pf0, _pa1, _pf1);
        pince -> SetDiffAmps0(_pa00, _pa01);
        pince -> i01 = tmpPince -> i01;
        pince -> i02 = tmpPince -> i02;
        if (debug) printf("\n pince->i01=%d", pince->i01);
        if (debug) printf("\n pince->ipa1=%d", pince->ipa1);
        if (debug) printf("\n pince->ipf1=%d", pince->ipf1);
		if (debug) printf ("\n ce5");
		if (debug) printf ("\n Xd[0]->GetLignes(): %d, Yd[0]->GetLignes(): %d", Xd[0]->GetLignes(), Yd[0]->GetLignes());
		if (debug) printf ("\n Xd[1]->GetLignes(): %d, Yd[1]->GetLignes(): %d", Xd[1]->GetLignes(), Yd[1]->GetLignes());
        calcPincePlusNew(Xd[0], Yd[0], Xd[1], Yd[1], pince, &Xdp[0], &Ydp[0], &Xdp[1], &Ydp[1]);
		if (debug) printf ("\n ce6");
        if ((myDeb != posDeb1) || (myDeb != posDeb2) || (myFin != posFin1) || (myFin != posFin2)) {
            // rezem Xd, newXd
            // rXdp[0] = GetFunctionSrezDeb(P[0], Xdp[0], posDeb1);
            // rYdp[0] = GetFunctionSrezDeb(P[0], Ydp[0], posDeb1);
            rXdp[0] = GetFunctionSrezDebFin(P[0], Xdp[0], posDeb1, posFin1);
            rYdp[0] = GetFunctionSrezDebFin(P[0], Ydp[0], posDeb1, posFin1);
            delete(Xdp[0]);
            Xdp[0]=rXdp[0];
            delete(Ydp[0]);
            Ydp[0]=rYdp[0];


            // rXdp[1] = GetFunctionSrezDeb(P[1], Xdp[1], posDeb2);
            // rYdp[1] = GetFunctionSrezDeb(P[1], Ydp[1], posDeb2);
            rXdp[1] = GetFunctionSrezDebFin(P[1], Xdp[1], posDeb2, posFin2);
            rYdp[1] = GetFunctionSrezDebFin(P[1], Ydp[1], posDeb2, posFin2);
            delete(Xdp[1]);
            Xdp[1]=rXdp[1];
            delete(Ydp[1]);
            Ydp[1]=rYdp[1];


            //rXd[0] = GetFunctionSrezDeb(P[0], Xd[0], posDeb1);
            //rYd[0] = GetFunctionSrezDeb(P[0], Yd[0], posDeb1);
            rXd[0] = GetFunctionSrezDebFin(P[0], Xd[0], posDeb1, posFin1);
            rYd[0] = GetFunctionSrezDebFin(P[0], Yd[0], posDeb1, posFin1);
            delete(Xd[0]);
            Xd[0]=rXd[0];
            delete(Yd[0]);
            Yd[0]=rYd[0];


            // rXd[1] = GetFunctionSrezDeb(P[1], Xd[1], posDeb2);
            // rYd[1] = GetFunctionSrezDeb(P[1], Yd[1], posDeb2);
            rXd[1] = GetFunctionSrezDebFin(P[1], Xd[1], posDeb2, posFin2);
            rYd[1] = GetFunctionSrezDebFin(P[1], Yd[1], posDeb2, posFin2);

            delete(Xd[1]);
            Xd[1]=rXd[1];
            delete(Yd[1]);
            Yd[1]=rYd[1];

            //rP[0] = GetFunctionSrezDeb(P[0], P[0], posDeb1);
			rP[0] = GetFunctionSrezDebFin(P[0], P[0], posDeb1, posFin1);
            delete(P[0]);
            P[0]=rP[0];

            //rP[1] = GetFunctionSrezDeb(P[1], P[1], posDeb2);
			rP[1] = GetFunctionSrezDebFin(P[1], P[1], posDeb2, posFin2);
            delete(P[1]);
            P[1]=rP[1];
        }

        //delete(P[0]);
        //delete(P[1]);
        //P[0] = newP[0];
        //P[1] = newP[1];
    }
    if (debug) printf ("\n PanelLayoutElement::calculateExport.1");
	if (debug) printf ("\n ce7");
    double marge1 = gfd->Marge[0];
    double marge2 = gfd->Marge[1];
    double margeDeb = gfd->MargeDeb;
    double margeFin = gfd->MargeFin;

    if (side == EXT_SIDE) {
        margeFin = gfd->margeFinExt;
        charname="P";
    }
    if (side == INT_SIDE) {
        // panels interior
        margeFin = gfd->margeFinInt;
        charname="P";
        if (gfd->VentHoles)
            if ((n1==-1) && gfd->VentCentralNerv) {
                 if (posFin1!=100.0f) {
                    margeFin = gfd->MargeFin;
                    charname="V";
                }
            } else
                if (gfd->noNervVH[n1]) {
                   if (posFin1!=100.0f) {
                        margeFin = gfd->MargeFin;
                        charname="V";
                   }
                }
    }
    if (debug) printf ("\n PanelLayoutElement::calculateExport.2");
	//printf ("\n ce8");
    int n = gfd->Form->m_nbProfils;
    char text[100];
    sprintf(text, "%s%d %.1fF%d-%.1fF%d t %d %.1fF%d-%.1fF%d", charname,
		n1, posDeb1, fd1, posFin1, ff1,
		n2, posDeb2, fd2, posFin2, ff2);
    TAxe *AxeP, *AxePD,  *AxePTD,  *AxeMD,  *AxeCD,  *AxeRepD;
    if (debug) printf ("\n PanelLayoutElement::calculateExport.3");
    leexport = new LayoutElementExport();

    GenerateCourbe(gfd, Xdp[0], Ydp[0], P[0], n1, posDeb1, fd1, posFin1, ff1,
        Xdp[1], Ydp[1], P[1], n2, posDeb2, fd2, posFin2, ff2, text,
        &(leexport->AxeP), &(leexport->AxePD), &(leexport->AxePTD), &(leexport->AxeMD), &(leexport->AxeCD), &(leexport->AxeRepD),
		vent, marge1, marge2, margeDeb, margeFin, true, debug,
		true, Xd[0], Yd[0], coeff1, Xd[1], Yd[1], coeff2);

    if (debug) printf ("\n PanelLayoutElement::calculateExport.3.2");

    //printf ("\n PanelLayoutElement::calculateExport.4");
    calcMaxWH(Xdp[0], Ydp[0], Xdp[1], Ydp[1], &(leexport->W), &(leexport->H));
    if (debug) printf ("\n PanelLayoutElement::calculateExport.5");
    delete (Xdp[0]);
    delete (Ydp[0]);
    delete (Xdp[1]);
    delete (Ydp[1]);
    //printf ("\n... PanelLayoutElement::calculateExport");
}
void DiagNervLayoutElement::calculateExport(WindPatternsProject* gfd) {

    bool debug = false;

    if (debug) printf ("\n DiagNervLayoutElement::calculateExport");
    Matrix * Xd[2], *Yd[2],*newXd[2], *newYd[2];
    Matrix *X[2], *Y[2], *Z[2], *P[2];
    char* charname = new char[1];

    calcPatron(gfd, n1, s1, fd1, ff1, posDeb1, posFin1,
        n2, s2, fd2, ff2, posDeb2, posFin2,
        &Xd[0], &Yd[0], &Xd[1], &Yd[1],
        &X[0], &Y[0], &Z[0], &P[0],
        &X[1], &Y[1], &Z[1], &P[1]);

    getLayoutLogger()->logprintf(" c[%f]", coeff);
    getLayoutLogger()->logprintf(" %d %d ff1=%d ff2=%d", n1, n2,  ff1, ff2);

    if (coeff != 0) {
        calcPatronWithCoeff(Xd[0], Yd[0], Xd[1], Yd[1], coeff, &newXd[0], &newYd[0], &newXd[1], &newYd[1]);
        delete (Xd[0]);
        delete (Yd[0]);
        delete(Xd[1]);
        delete (Yd[1]);
        Xd[0] = newXd[0];
        Yd[0] = newYd[0];
        Xd[1] = newXd[1];
        Yd[1] = newYd[1];
    }

    double marge1 = gfd->Marge[0];
    double marge2 = gfd->Marge[1];
    double margeDeb = gfd->MargeDeb;
    double margeFin = gfd->MargeFin;

    // diagonal nervures
    margeFin = gfd->margeFinDiagNerv;
    vent = gfd->VentilationLayout;
    charname="D";
    int n = gfd->Form->m_nbProfils;
    char text[100];
    sprintf(text, "%s%dF%dt%dF%d", charname, n1, ff1, n2, ff2);
    leexport = new LayoutElementExport();
    GenerateCourbe(gfd, Xd[0], Yd[0], P[0], n1, posDeb1, fd1, posFin1, ff1,
        Xd[1], Yd[1], P[1], n2, posDeb2, fd2, posFin2, ff2, text,
        &(leexport->AxeP), &(leexport->AxePD), &(leexport->AxePTD), &(leexport->AxeMD), &(leexport->AxeCD), &(leexport->AxeRepD), vent, marge1, marge2, margeDeb, margeFin,
        true, debug);

    calcMaxWH(Xd[0], Yd[0], Xd[1], Yd[1], &(leexport->W), &(leexport->H));

    delete (X[0]);
    delete (Y[0]);
    delete (Z[0]);
    delete (P[0]);
    delete (X[1]);
    delete (Y[1]);
    delete (Z[1]);
    delete (P[1]);
    delete (Xd[0]);
    delete (Yd[0]);
    delete (Xd[1]);
    delete (Yd[1]);
}

Layout::~Layout() {
}

Layout::Layout(int n) {
    funcL1 = new Matrix*[n - 1];
    funcR1 = new Matrix*[n - 1];
    funcL2 = new Matrix*[n - 1];
    funcR2 = new Matrix*[n - 1];

    pinceLAAmp1 = new double[n - 1];
    pinceRAAmp1 = new double[n - 1];

    pinceLA0Amp1 = new double[n - 1];
    pinceRA0Amp1 = new double[n - 1];

    pinceLFAmp1 = new double[n - 1];
    pinceRFAmp1 = new double[n - 1];

    lenp1 = new double[n - 1];

    pinceLA0Amp2 = new double[n - 1];
    pinceRA0Amp2 = new double[n - 1];

    pinceLAAmp2 = new double[n - 1];
    pinceRAAmp2 = new double[n - 1];

    pinceLFAmp2 = new double[n - 1];
    pinceRFAmp2 = new double[n - 1];

    lenp2 = new double[n - 1];
    coeffn = new double[n - 1];

    coeffd = new double[n - 1];
    noNervD = new int[n - 1];
    noNervD1 = new int[n - 1];
    noNervD2 = new int[n - 1];

    isCenterPanel = false;
}

Layout* calcIndepPinceLayout(WindPatternsProject* gfd, Form* F) {
	bool debug = false;
    if (debug) printf("\n calcIndepPinceLayout()");

    char name[255];
    if (getLayoutLogger() == NULL) {
	char* getlfn=getLogFileName();
        strcpy(name, getlfn);
        strcpy(gfd->logFileName, name);
        LayoutLoggerInit(name);
    }
    getLayoutLogger()->logprintf("\nPROJECT %s", gfd->name);
    getLayoutLogger()->logprintf("\nFORME %s", gfd->fileNameForm);
    getLayoutLogger()->logprintf("\nREP_POINTS %s", gfd->fileNameRepPoints);
    getLayoutLogger()->logprintf("\nDIAG_NERVS %s", gfd->fileNameDiagNerv);
    getLayoutLogger()->logprintf("\nVENT_HOLES %s", gfd->fileNameVentHoles);
    int n = F->m_nbProfils;
    bool isCenterPanel = (1 & F->NbCaiss);
    Layout* layout = new Layout(n);
    layout->isDesign = (gfd->layoutWithDesign == 1);
    layout->isCenterPanel = isCenterPanel;
    int startNerv;

    if (isCenterPanel) {
        startNerv = 0;
    } else {
        startNerv = 1;

        Pince* p0f1 = getPince(gfd, 0, 1, 1);
        layout->pinceLAAmp1[0] = p0f1->AmpA;
        layout->pinceLFAmp1[0] = p0f1->AmpF;
        layout->pinceRAAmp1[0] = p0f1->AmpA;
        layout->pinceRFAmp1[0] = p0f1->AmpF;
        layout->funcL1[0] = p0f1->function1;
        layout->funcR1[0] = p0f1->function1;
        goCalcPinceLen(gfd, 0, 1, &(layout->lenp1[0]));

        Pince* p0f2 = getPince(gfd, 0, 1, 2);
        layout->pinceLAAmp2[0] = p0f2->AmpA;
        layout->pinceRAAmp2[0] = p0f2->AmpA;
        layout->pinceLFAmp2[0] = p0f2->AmpF;
        layout->pinceRFAmp2[0] = p0f2->AmpF;
        layout->funcL2[0] = p0f2->function1;
        layout->funcR2[0] = p0f2->function1;
        goCalcPinceLen(gfd, 0, 2, &(layout->lenp2[0]));
    }

    for (int inoNerv = startNerv; inoNerv < n - 1; inoNerv++) {
        goCalcIndepPinceNew(gfd, inoNerv, 1,
                &(layout->pinceLAAmp1[inoNerv]), &(layout->pinceLFAmp1[inoNerv]), &(layout->funcL1[inoNerv]),
                &(layout->pinceRAAmp1[inoNerv]), &(layout->pinceRFAmp1[inoNerv]), &(layout->funcR1[inoNerv]),
                &(layout->lenp1[inoNerv]));
    }

    getLayoutLogger()->logprintf("\n");
    for (int inoNerv = startNerv; inoNerv < n - 1; inoNerv++) {
        goCalcIndepPinceNew(gfd, inoNerv, 2,
                &(layout->pinceLAAmp2[inoNerv]), &(layout->pinceLFAmp2[inoNerv]), &(layout->funcL2[inoNerv]),
                &(layout->pinceRAAmp2[inoNerv]), &(layout->pinceRFAmp2[inoNerv]), &(layout->funcR2[inoNerv]),
                &(layout->lenp2[inoNerv]));
    }

    // calc Profile Nervures
    for (int inoNerv = 0; inoNerv < n - 1; inoNerv++) {
        goCalcNervureWithPince(gfd, inoNerv, 2, inoNerv, 1,
                layout->lenp2[inoNerv], layout->lenp1[inoNerv], &(layout->coeffn[inoNerv]));
    }

    // calc Diagonal Nervures
	if (gfd->DiagNervs) {
		int isave = 0;
		for (int i = 0; i < gfd->quantDiag; i++) {
			goCalcNervureWithPince(gfd, gfd->noNervD[i], 2, gfd->noNervD[i] - 1, 1, layout->lenp2[gfd->noNervD[i]], layout->lenp1[gfd->noNervD[i] - 1], &(layout->coeffd[isave]));
			isave++;
			goCalcNervureWithPince(gfd, gfd->noNervD[i], 2, gfd->noNervD[i] + 1, 1, layout->lenp2[gfd->noNervD[i]], layout->lenp1[gfd->noNervD[i] + 1], &(layout->coeffd[isave]));
			isave++;
		}
	}
    if (debug) printf("\n...calcIndepPinceLayout()");
    return layout;
}


void Layout::preparePanelLayoutElementCoeff(LayoutElement* le, WindPatternsProject* gfd) {
   	double coeff1 = coeffn[le->n1];
	double coeff2 = coeffn[le->n2];
    int n = gfd->Form->m_nbProfils;
	if (le->n1 == -1) coeff1 = coeffn[0];
	if (le->n2 == n-1) coeff2 = -1;
    le->coeff1=coeff1;
    le->coeff2=coeff2;
}


void Layout::prepareKlapan(WindPatternsProject* gfd, int i) {
    if ((gfd->VentHoles) && (gfd->noNervVH[i])) {
        if (gfd->LayoutKlapans) {
            if (gfd->VentHolesDouble){
                KlapanLayoutElement* le4 = new KlapanLayoutElement();
                le4->side=INT_SIDE;
                le4->n1 = i;
                le4->n2 = i+1;
                le4->s1 = false;
                le4->s2 = false;
                le4->posKlapanIntDeb=0.0f;
                le4->posKlapanFin=gfd->PosKlapanFin;
                le4->isPince = 0;
                le4->isKlapan = 1;
		        le4->fd1 = 2;
		        le4->fd2 = 2;
		        le4->ff1 = 2;
		        le4->ff2 = 2;
                preparePanelLayoutElementCoeff(le4, gfd);
                klapans.push_back(le4);
            }

            KlapanLayoutElement* le5 = new KlapanLayoutElement();
            le5->side=INT_SIDE;
            le5->n1 = i;
            le5->n2 = i+1;
            le5->s1 = false;
            le5->s2 = false;
            le5->posKlapanIntDeb=gfd->VentHolesDeb;
            le5->posKlapanFin=gfd->PosKlapanFin;
            le5->isPince = 0;
            le5->isKlapan = 1;
	        le5->fd1 = 2;
	        le5->fd2 = 2;
	        le5->ff1 = 2;
	        le5->ff2 = 2;
            preparePanelLayoutElementCoeff(le5, gfd);
            klapans.push_back(le5);

            KlapanLayoutElement* le6 = new KlapanLayoutElement();
            le6->side=INT_SIDE;
            le6->n1 = i;
            le6->n2 = i+1;
            le6->s1 = false;
            le6->s2 = false;
            le6->posKlapanIntDeb=gfd->VentHolesFin;
            le6->posKlapanFin=gfd->PosKlapanFin;
            le6->isPince = 0;
            le6->isKlapan = 1;
	        le6->fd1 = 2;
	        le6->fd2 = 2;
	        le6->ff1 = 2;
	        le6->ff2 = 2;
            preparePanelLayoutElementCoeff(le6, gfd);
            klapans.push_back(le6);

        }
    }
}

void Layout::preparePanelInt(WindPatternsProject* gfd, int i) {
    PanelLayoutElement* le1 = new PanelLayoutElement();
    le1->side=INT_SIDE;
    int n = gfd->Form->m_nbProfils;
    double tpF1, tpF2;
    if ((gfd->VentHoles) && (gfd->noNervVH[i])) {
        tpF1 = gfd->VentHolesFin;
        tpF2 = gfd->VentHolesFin;
    } else {
        tpF1 = 100.0f;
        tpF2 = 100.0f;
    }
    //tututu
    if ((gfd->VentHoles) && (gfd->VentHolesDouble) && (gfd->noNervVH[i])) {
        le1->n1 = i;
        le1->n2 = i + 1;
        le1->s1 = false;
        le1->s2 = false;
        le1->posDeb1 = 0.0f;
        le1->posFin1 = debBorder;
        le1->posDeb2 = 0.0f;
        le1->posFin2 = debBorder;

        le1->p2a00 = 0.0f;
        le1->p2a0 = pinceRAAmp2 [i];
        le1->p2f0 = pinceRFAmp2 [i];

        le1->p2a01 = 0.0f;
        le1->p2a1 = pinceLAAmp2 [i + 1];
        le1->p2f1 = pinceLFAmp2 [i + 1];

        le1->func1f0 = 0;
        le1->func1f1 = 0;

        le1->func2f0 = funcR2[i];
        le1->func2f1 = funcL2[i + 1];

        le1->fd1 = 2;
        le1->fd2 = 2;
        le1->ff1 = 2;
        le1->ff2 = 2;
        le1->isPince = 1;
        le1->isVentHole = 1;
        le1->coeff = 0.0f;
        preparePanelLayoutElementCoeff(le1, gfd);
        panelsInt.push_back(le1);
    }

    PanelLayoutElement* le2 = new PanelLayoutElement();
    le2->side=INT_SIDE;
    le2->n1 = i;
    le2->n2 = i + 1;
    le2->s1 = false;
    le2->s2 = false;
    le2->posDeb1 = debBorder;
    le2->posFin1 = tpF1;
    le2->posDeb2 = debBorder;
    le2->posFin2 = tpF2;

    le2->p2a00 = 0.0f;
    le2->p2a0 = pinceRAAmp2 [i];
    le2->p2f0 = pinceRFAmp2 [i];

    le2->p2a01 = 0.0f;
    le2->p2a1 = pinceLAAmp2 [i + 1];
    le2->p2f1 = pinceLFAmp2 [i + 1];

    le2->func1f0 = 0;
    le2->func1f1 = 0;

    le2->func2f0 = funcR2[i];
    le2->func2f1 = funcL2[i + 1];

    le2->fd1 = 2;
    le2->fd2 = 2;
    le2->ff1 = 2;
    le2->ff2 = 2;
    le2->isPince = 1;
    
    le2->isVentHole = 0;
    if ((gfd->VentHoles) && (gfd->noNervVH[i])) {
        le2->isVentHole = 1;
    }

    le2->coeff = 0.0f;

    if (i == (n - 2)) {
        le2->func2f1 = Zeros((le2->func2f0)->GetLignes(), 1);
        le2->p2a01 = 0.0f;
        le2->p2a1 = 0.0f;
        le2->p2f1 = 0.0f;

        if (!((gfd->VentHoles) && (gfd->noNervVH[i]))  ) {
            le2->posFin1 = 100.0f;
            le2->posDeb2 = 0.0f;
            le2->posFin2 = 100.0f;
        }
    }
    preparePanelLayoutElementCoeff(le2, gfd);
    panelsInt.push_back(le2);

    if ((gfd->VentHoles) && (gfd->noNervVH[i])) {
        PanelLayoutElement* le3 = new PanelLayoutElement();
        le3->side=INT_SIDE;
        le3->n1 = i;
        le3->n2 = i + 1;
        le3->s1 = false;
        le3->s2 = false;
        le3->posDeb1 = gfd->VentHolesFin;
        le3->posFin1 = 100.0f;
        le3->posDeb2 = gfd->VentHolesFin;
        le3->posFin2 = 100.0f;

        le3->p2a00 = 0.0f;
        le3->p2a0 = pinceRAAmp2 [i];
        le3->p2f0 = pinceRFAmp2 [i];

        le3->func1f0 = 0;
        le3->func1f1 = 0;

        le3->func2f0 = funcR2[i];
        le3->func2f1 = funcL2[i + 1];

        le3->p2a01 = 0.0f;
        le3->p2a1 = pinceLAAmp2 [i + 1];
        le3->p2f1 = pinceLFAmp2 [i + 1];

        le3->fd1 = 2;
        le3->fd2 = 2;
        le3->ff1 = 2;
        le3->ff2 = 2;

        le3->isPince = 1;
        le3->isVentHole = 0;
        le3->coeff = 0.0f;

        //edge
        if (i == (n - 2)) {
            le3->func2f1 = Zeros((le3->func2f0)->GetLignes(), 1);
            le3->p2a01 = 0.0f;
            le3->p2a1 = 0.0f;
            le3->p2f1 = 0.0f;
        }
        preparePanelLayoutElementCoeff(le3, gfd);
        panelsInt.push_back(le3);
    }
}

void Layout::preparePanelExt(WindPatternsProject* gfd, int i) {
    int n = gfd->Form->m_nbProfils;
    // face == 1
    PanelLayoutElement* le = new PanelLayoutElement();
    le->side=EXT_SIDE;
    le->n1 = i;
    le->n2 = i + 1;
    le->s1 = false;
    le->s2 = false;
    le->posDeb1 = debBorder;
    le->posFin1 = 100.0f;
    le->posDeb2 = debBorder;
    le->posFin2 = 100.0f;

    le->p1a00 = pinceRAAmp2 [i];
    le->p1a0 = pinceRAAmp1 [i];
    le->p1f0 = pinceRFAmp1 [i];

    le->p1a01 = pinceLAAmp2 [i + 1];
    le->p1a1 = pinceLAAmp1 [i + 1];
    le->p1f1 = pinceLFAmp1 [i + 1];

    le->func1f0 = funcR1[i];
    le->func1f1 = funcL1[i + 1];

    le->func2f0 = 0;
    le->func2f1 = 0;

    le->ff1 = 1;
    le->ff2 = 1;
    le->fd1 = faceDebBorder;
    le->fd2 = faceDebBorder;

    le->isPince = 1;
	le->isKlapan = 0;
    le->isVentHole = 0;
    le->coeff = 0.0f;
    // tututu
    if ((gfd->VentHoles) && (gfd->VentHolesDouble) && (gfd->noNervVH[i])) {
        le->posDeb1 = 0.0f;
        le->posDeb2 = 0.0f;
        le->fd1 = 1;
        le->fd2 = 1;

    }

    if (i == (n - 2)) {
        le->func1f1 = Zeros((le->func1f0)->GetLignes(), 1);
        le->p1a01 = 0.0f;
        le->p1a1 = 0.0f;
        le->p1f1 = 0.0f;
        if (!((gfd->VentHoles) && (gfd->noNervVH[i]))) {
            le->posDeb2 = 0.0f;
            le->fd2 = 2;
        }
    }
    preparePanelLayoutElementCoeff(le, gfd);
    panelsExt.push_back(le);
}


void Layout::prepareCenterPanelInt(WindPatternsProject* gfd) {
        double tpF1, tpF2;
        if ((gfd->VentHoles) && (gfd->VentCentralNerv)) {
            tpF1 = gfd->VentHolesFin;
            tpF2 = gfd->VentHolesFin;
        } else {
            tpF1 = 100.0f;
            tpF2 = 100.0f;
        }
        int i = 0;
        if ((gfd->VentHoles) && (gfd->VentHolesDouble) && (gfd->VentCentralNerv)) {
            PanelLayoutElement* le = new PanelLayoutElement();
            le->side=INT_SIDE;
            le->n1 = -1;
            le->n2 = 0;
            le->s1 = false;
            le->s2 = false;
            le->posDeb1 = 0.0f;
            le->posFin1 = debBorder;
            le->posDeb2 = 0.0f;
            le->posFin2 = debBorder;

            le->p2a00 = 0.0f;
            le->p2a0 = pinceRAAmp2 [i];
            le->p2f0 = pinceRFAmp2 [i];

            le->p2a01 = 0.0f;
            le->p2a1 = pinceLAAmp2 [i + 1];
            le->p2f1 = pinceLFAmp2 [i + 1];

		    le->func1f0 = 0;
		    le->func1f1 = 0;
            le->func2f0 = funcR2[i];
            le->func2f1 = funcL2[i + 1];

            le->fd1 = 2;
            le->fd2 = 2;
            le->ff1 = 2;
            le->ff2 = 2;
            le->isPince = 1;
            le->coeff = 0.0f;
            panelsInt.push_back(le);
        }
        PanelLayoutElement* le1 = new PanelLayoutElement();
        le1->side=INT_SIDE;
        le1->n1 = -1;
        le1->n2 = 0;
        le1->s1 = false;
        le1->s2 = false;
        le1->fd1 = 2;
        le1->fd2 = 2;
        le1->ff1 = 2;
        le1->ff2 = 2;

        le1->isPince = 1;
        le1->coeff = 0.0f;

        le1->posDeb1 = debBorder;
        le1->posFin1 = tpF1;
        le1->posDeb2 = debBorder;
        le1->posFin2 = tpF2;

        le1->p2a00 = 0.0f;
        le1->p2a0 = pinceLAAmp2[0];
        le1->p2f0 = pinceLFAmp2[0];

        le1->p2a01 = 0.0f;
        le1->p2a1 = pinceLAAmp2[0];
        le1->p2f1 = pinceLFAmp2[0];

		le1->func1f0 = 0;
		le1->func1f1 = 0;

        le1->func2f0 = funcL2[0];
        le1->func2f1 = funcL2[0];
        panelsInt.push_back(le1);

        if ((gfd->VentHoles) && (gfd->VentCentralNerv)) {
            PanelLayoutElement* le2 = new PanelLayoutElement();
            le2->side=INT_SIDE;
            le2->n1 = -1;
            le2->n2 = 0;
            le2->s1 = false;
            le2->s2 = false;
            le2->fd1 = 2;
            le2->fd2 = 2;
            le2->ff1 = 2;
            le2->ff2 = 2;

            le2->isPince = 1;
            le2->coeff = 0.0f;

            le2->posDeb1 = gfd->VentHolesFin;
            le2->posFin1 = 100.0f;
            le2->posDeb2 = gfd->VentHolesFin;
            le2->posFin2 = 100.0f;


            le2->p2a00 = 0.0f;
            le2->p2a0 = pinceLAAmp2[0];
            le2->p2f0 = pinceLFAmp2[0];

            le2->p2a01 = 0.0f;
            le2->p2a1 = pinceLAAmp2[0];
            le2->p2f1 = pinceLFAmp2[0];

		    le2->func1f0 = 0;
		    le2->func1f1 = 0;

            le2->func2f0 = funcL2[0];
            le2->func2f1 = funcL2[0];
            panelsInt.push_back(le2);

            if (gfd->LayoutKlapans) {
                if (gfd->VentHolesDouble){
                    KlapanLayoutElement* le3 = new KlapanLayoutElement();
                    le3->side=INT_SIDE;
                    le3->n1 = -1;
                    le3->n2 = 0;
                    le3->s1 = false;
                    le3->s2 = false;
                    le3->posKlapanIntDeb = 0.0f;
                    le3->posKlapanFin = gfd->PosKlapanFin;
                    le3->isKlapan = 1;
                    le3->isPince = 0;
                    //printf ("\n klapan=%d", isave);
                    klapans.push_back(le3);
                }
                    KlapanLayoutElement* le4 = new KlapanLayoutElement();
                    le4->side=INT_SIDE;
                    le4->n1 = -1;
                    le4->n2 = 0;
                    le4->s1 = false;
                    le4->s2 = false;
                    le4->posKlapanIntDeb = gfd->VentHolesDeb;
                    le4->posKlapanFin = gfd->PosKlapanFin;
                    le4->isKlapan = 1;
                    le4->isPince = 0;
                    klapans.push_back(le4);


                    KlapanLayoutElement* le5 = new KlapanLayoutElement();
                    le5->side=INT_SIDE;
                    le5->n1 = -1;
                    le5->n2 = 0;
                    le5->s1 = false;
                    le5->s2 = false;
                    le5->posKlapanIntDeb = gfd->VentHolesFin;
                    le5->posKlapanFin = gfd->PosKlapanFin;
                    le5->isKlapan = 1;
                    le5->isPince = 0;
                    klapans.push_back(le5);
            }
        }
}

void Layout::prepareCenterPanelExt(WindPatternsProject* gfd) {
    PanelLayoutElement* le = new PanelLayoutElement();
    le->side=EXT_SIDE;
    le->n1 = -1;
    le->n2 = 0;
    le->s1 = false;
    le->s2 = false;
    le->fd1 = faceDebBorder;
    le->fd2 = faceDebBorder;
    le->ff1 = 1;
    le->ff2 = 1;
    le->isPince = 1;
	le->isKlapan = 0;
    le->coeff = 0.0f;
    le->posDeb1 = debBorder;
    le->posFin1 = 100.0f;
    le->posDeb2 = debBorder;
    le->posFin2 = 100.0f;
    le->p1a00 = pinceLAAmp2[0];
    le->p1a0 = pinceLAAmp1[0];
    le->p1f0 = pinceLFAmp1[0];
    le->p1a01 = pinceLAAmp2[0];
    le->p1a1 = pinceLAAmp1[0];
    le->p1f1 = pinceLFAmp1[0];

    le->func1f0 = funcL1[0];
    le->func1f1 = funcL1[0];
    le->func2f0 = 0;
    le->func2f1 = 0;

    if ((gfd->VentHoles) && (gfd->VentHolesDouble) && (gfd->VentCentralNerv)) {
        le->posDeb1 = 0.0f;
        le->posDeb2 = 0.0f;
        le->fd1 = 1;
        le->fd2 = 1;
    }
    panelsExt.push_back(le);
}

void Layout::prepareProfile(WindPatternsProject* gfd, int i) {
    ProfLayoutElement* lep = new ProfLayoutElement();
    lep->n1 = i;
    lep->n2 = i;
    lep->s1 = false;
    lep->s2 = false;
    lep->posDeb1 = 0.0f;
    lep->posFin1 = 100.0f;
    lep->posDeb2 = 0.0f;
    lep->posFin2 = 100.0f;
    lep->isPince = 0;
    lep->fd1 = 2;
    lep->ff1 = 2;
    lep->fd2 = 1;
    lep->ff2 = 1;
    lep->coeff = coeffn[i];
    profs.push_back(lep);
}

void Layout::prepareDiagNerv(WindPatternsProject* gfd, int i) {
    DiagNervLayoutElement* led1 = new DiagNervLayoutElement();
    int k = 2*i - 1;
    led1->n1 = gfd->noNervD[i];
    led1->n2 = gfd->noNervD[i] - 1;
    led1->s1 = false;
    led1->s2 = false;
    led1->isPince = 0;
    led1->fd1 = 2;
    led1->ff1 = 2;
    led1->fd2 = 1;
    led1->ff2 = 1;
    led1->posDeb1 = gfd->PosDiagNerv2A;
    led1->posFin1 = gfd->PosDiagNerv2F;

    led1->posDeb2 = gfd->PosDiagNerv1A;
    led1->posFin2 = gfd->PosDiagNerv1F;
    led1->coeff = coeffd[k];
    k++;

    DiagNervLayoutElement* led2 = new DiagNervLayoutElement();
    led2->n1 = gfd->noNervD[i];
    led2->n2 = gfd->noNervD[i] + 1;
    led2->s1 = false;
    led2->s2 = false;
    led2->isPince = 0;
    led2->fd1 = 2;
    led2->ff1 = 2;
    led2->fd2 = 1;
    led2->ff2 = 1;

    led2->posDeb1 = gfd->PosDiagNerv2A;
    led2->posFin1 = gfd->PosDiagNerv2F;
    led2->posDeb2 = gfd->PosDiagNerv1A;
    led2->posFin2 = gfd->PosDiagNerv1F;
    led2->coeff = coeffd[k];
    diagNervs.push_back(led2);
}



PanelLayoutElement* getPanelLayoutElementCopy(PanelLayoutElement* p){
    PanelLayoutElement* ple = new PanelLayoutElement();
    ple->side=p->side;
    ple->n1 = p->n1;
    ple->n2 = p->n2;
    ple->s1 = p->s1;
    ple->s2 = p->s2;
    ple->posDeb1 = p->posDeb1;
    ple->posFin1 = p->posFin1;
    ple->posDeb2 = p->posDeb2;
    ple->posFin2 = p->posFin2;

    ple->p1a00 = p->p1a00;
    ple->p1a0 = p->p1a0;
    ple->p1f0 = p->p1f0;

    ple->p1a01 = p->p1a01;
    ple->p1a1 = p->p1a1;
    ple->p1f1 = p->p1f1;

    ple->func1f0 = CloneMat(p->func1f0);
    ple->func1f1 = CloneMat(p->func1f1);

    ple->ff1 = p->ff1;
    ple->ff2 = p->ff2;
    ple->fd1 = p->fd1;
    ple->fd2 = p->fd2;

    ple->isPince = p->isPince;
    ple->isKlapan = p->isKlapan;
    ple->coeff = p->coeff;
    ple->isVentHole = p->isVentHole;
    return ple;
}

PanelLayoutElement* getPanelLayoutElementCopyWithNewQuad(PanelLayoutElement* p, Quad* q) {
	bool debug = false;
	if (debug) printf ("\n getPanelLayoutElementCopyWithNewQuad()");
    PanelLayoutElement* ple = new PanelLayoutElement();
    ple->side=p->side;
    ple->n1 = p->n1;
    ple->n2 = p->n2;
    ple->s1 = p->s1;
    ple->s2 = p->s2;
	if (debug) printf ("\n 1");
	ple->posDeb1 = q->pd1;
    ple->posFin1 = q->pf1;
    ple->posDeb2 = q->pd2;
    ple->posFin2 = q->pf2;
	if (debug) printf ("\n 2");
    ple->p1a00 = p->p1a00;
    ple->p1a0 = p->p1a0;
    ple->p1f0 = p->p1f0;
	if (debug) printf ("\n 3");
    ple->p1a01 = p->p1a01;
    ple->p1a1 = p->p1a1;
    ple->p1f1 = p->p1f1;
	if (debug) printf ("\n 4");

	if (debug) printf ("\n p->func1f0");
	if (p->func1f0 != 0) {
		if (debug) printf ("\n --p->func1f0");
		ple->func1f0 = CloneMat(p->func1f0);
	}
	
	if (debug) printf ("\n p->func1f1");
	if (p->func1f1 != 0) {
		if (debug) printf ("\n --p->func1f1");
		ple->func1f1 = CloneMat(p->func1f1);
	}

	if (debug) printf ("\n p->func2f0");
	if (p->func2f0 != 0) {
		if (debug) printf ("\n --p->func2f0");
		ple->func2f0 = CloneMat(p->func2f0);
	}

	if (debug) printf ("\n p->func2f1");
	if (p->func2f1 != 0) {
		if (debug) printf ("\n --p->func2f1");
		ple->func2f1 = CloneMat(p->func2f1);
	}


	if (debug) printf ("\n 5");
    ple->ff1 = p->ff1;
    ple->ff2 = p->ff2;
    ple->fd1 = p->fd1;
    ple->fd2 = p->fd2;
	if (debug) printf ("\n 6");
    ple->isPince = p->isPince;
    ple->isKlapan = p->isKlapan;
	if (debug) printf ("\n 7");
    ple->isVentHole = p->isVentHole;
    ple->coeff = p->coeff;
	if (debug) printf ("\n ...getPanelLayoutElementCopyWithNewQuad()");
    return ple;
}


PanelLayoutElement* getPanelLayoutElementCopyWithNewBorders(PanelLayoutElement* p, double d1, double f1, double d2, double f2, int face) {
    PanelLayoutElement* ple = new PanelLayoutElement();
    ple->side=p->side;
    ple->n1 = p->n1;
    ple->n2 = p->n2;
    ple->s1 = p->s1;
    ple->s2 = p->s2;
    ple->posDeb1 = d1;
    ple->posFin1 = f1;
    ple->posDeb2 = d2;
    ple->posFin2 = f2;

    ple->p1a00 = p->p1a00;
    ple->p1a0 = p->p1a0;
    ple->p1f0 = p->p1f0;

    ple->p1a01 = p->p1a01;
    ple->p1a1 = p->p1a1;
    ple->p1f1 = p->p1f1;

    ple->func1f0 = CloneMat(p->func1f0);
    ple->func1f1 = CloneMat(p->func1f1);

    ple->ff1 = p->ff1;
    ple->ff2 = p->ff2;
    ple->fd1 = p->fd1;
    ple->fd2 = p->fd2;

    ple->isPince = p->isPince;
    ple->isKlapan = p->isKlapan;
    ple->isVentHole = p->isVentHole;
    ple->coeff = p->coeff;
    return ple;
}

Quad* intersectQuads(Quad* qp, Quad* qcs) {
	bool debug = false;
	Quad* result = 0;
	double pd1 = qp->pd1;
	double pf1 = qp->pf1;
	double pd2 = qp->pd2;
	double pf2 = qp->pf2;
	if (debug) printf ("\n intersectQuads qp: (%f, %f) (%f, %f)", pd1, pf1, pd2, pf2);

	double csd1 = qcs->pd1;
	double csf1 = qcs->pf1;
	double csd2 = qcs->pd2;
	double csf2 = qcs->pf2;
	if (debug) printf ("\n intersectQuads qcs: (%f, %f) (%f, %f)", csd1, csf1, csd2, csf2);

	if ((csd1 <= pd1) && (csd2 <= pd2)) {
        if ((pf1  <= csf1) && (pf2 <= csf2)) {
            // ?
            // pf1  pf2 
            // pd1  pd2
            if ((pd1 <= pf1) && (pd2 <= pf2)){
				if (debug) printf ("\n!!! csd1 <= pd1 < pf1 <= csf1");
				result = new Quad(pd1, pf1, pd2, pf2);
                //ple = getPanelLayoutElementCopyWithNewBorders(p, pd1, pf1, pd2, pf2);

            }
        }

        if ((csf1  <= pf1) && (csf2 <= pf2)) {
            // ?
            // csf1  csf2 
            // pd1   pd2
            if ((pd1 <= csf1) && (pd2 <= csf2)){
				if (debug) printf ("\n!!! csd1 <= pd1 < csf1 <= pf1");
                //ple = getPanelLayoutElementCopyWithNewBorders(p, pd1, csf1, pd2, csf2);
				result = new Quad(pd1, csf1, pd2, csf2);
            }
        }
    }

    if ((pd1 <= csd1) && (pd2 <= csd2)) {
        if ((pf1  <= csf1) && (pf2 <= csf2)) {
            // ?
            // pf1  pf2 
            // csd1  csd2
            if ((csd1 <= pf1) && (csd2 <= pf2)){
				if (debug) printf ("\n!!! pd1 <= csd1 < pf1 <= csf1");
                //ple = getPanelLayoutElementCopyWithNewBorders(p, csd1, pf1,csd2, pf2);
				result = new Quad(csd1, pf1, csd2, pf2);
            }
        }

        if ((csf1  <= pf1) && (csf2 <= pf2)) {
            // ?
            // csf1  csf2 
            // csd1   csd2
            if ((csd1 <= csf1) && (csd2 <= csf2)){
				if (debug) printf ("\n!!! pd1 <= csd1 < csf1 <= pf1");
                //ple = getPanelLayoutElementCopyWithNewBorders(p, csd1, csf1, csd2, csf2);
				result = new Quad(csd1, csf1, csd2, csf2);
            }
        }
    }

	return result;
}

void Layout::intersectPanelExtWithColorSegment(PanelLayoutElement* p, ColorSegment* cs) {
	bool debug = false;
	if (debug) printf ("\n Layout::intersectPanelExtWithColorSegment()");

    double pd1 = p->posDeb1;
    double pf1 = p->posFin1;
    double pd2 = p->posDeb2;
    double pf2 = p->posFin2;
	if (debug) printf ("\nP  (%f(%d), %f(%d)) (%f(%d), %f(%d))", pd1, p->fd1,  pf1, p->ff1,pd2, p->fd2, pf2, p->ff2);

	int fd1 = p->fd1;
	int fd2 = p->fd2;

	bool was = false;
	if (fd1 == 2) {
		pd1 = 0.0f; // make fd1 = 1
		was = true;
	}
	if (fd2 == 2) {
		pd2 = 0.0f; // make fd2 = 1
		was = true;
	}
	if (debug) printf ("\nnewP  (%f(%d), %f(%d)) (%f(%d), %f(%d))", pd1, p->fd1, pf1, p->ff1, pd2, p->fd2,  pf2, p->ff2);
	
	Quad* qp = new Quad(pd1, pf1, pd2, pf2);

    double csd1 = cs->p00;
    double csf1 = cs->p10;
    double csd2 = cs->p01;
    double csf2 = cs->p11;
	Quad* qcs = new Quad(csd1, csf1, csd2, csf2);

	if (debug) printf ("\nCS (%f, %f) (%f, %f)", csd1, csf1, csd2, csf2);

	if (debug) printf ("\n IPE 1");
	Quad* qr = intersectQuads(qp, qcs);
	if (debug) printf ("\n IPE 2");
	if (qr == 0) {
		printf ("\n ERROR! quad result==0");
		return;
	}
	if (debug) printf ("Quad result (%f, %f) (%f, %f)", qr->pd1, qr->pf1, qr->pd2, qr->pf2);
	PanelLayoutElement* ple = getPanelLayoutElementCopyWithNewQuad(p, qr);
	if (debug) printf ("\n IPE 3");
	if (debug) printf ("\nP  (%f(%d), %f(%d)) (%f(%d), %f(%d))", ple->posDeb1, ple->fd1,  ple->posFin1, ple->ff1,ple->posDeb2, ple->fd2, ple->posFin2, ple->ff2);

	if ((was) && (qr->pd1 == 0.0f) && (qr->pd2 == 0.0f)) {
		if (debug) printf ("\n was IPE 4");
		if (debug) printf ("\n if ((qr->pd1 == 0.0f) && (qr->pd2 == 0.0f))...");
		ple->posDeb1 = p->posDeb1;
		ple->posDeb2 = p->posDeb2;
		ple->fd1 = p->fd1;
		ple->fd2 = p->fd2;
	} else {
		ple->fd1 = 1;
		ple->fd2 = 1;
	}
	if (debug) printf ("\n IPE 5");
	if (debug) printf ("\nresP  (%f(%d), %f(%d)) (%f(%d), %f(%d))", ple->posDeb1, ple->fd1, ple->posFin1, ple->ff1, ple->posDeb2, ple->fd2, ple->posFin2, ple->ff2);

	if (ple != 0 ) {
		if (debug) printf ("\n panelsExtDesign->push_back");
		panelsExtDesign.push_back(ple);
	}
}

void Layout::intersectPanelIntWithColorSegment(PanelLayoutElement* p, ColorSegment* cs) {
	bool debug = false;
	if (debug) printf ("\n Layout::intersectPanelIntWithColorSegment()");
    double pd1 = p->posDeb1;
    double pf1 = p->posFin1;
    double pd2 = p->posDeb2;
    double pf2 = p->posFin2;
	if (debug) printf ("\nP  (%f(%d), %f(%d)) (%f(%d), %f(%d))", pd1, p->fd1,  pf1, p->ff1,pd2, p->fd2, pf2, p->ff2);
	Quad* qp = new Quad(pd1, pf1, pd2, pf2);
    double csd1 = cs->p00;
    double csf1 = cs->p10;
    double csd2 = cs->p01;
    double csf2 = cs->p11;
	Quad* qcs = new Quad(csd1, csf1, csd2, csf2);
	if (debug) printf ("\nCS (%f, %f) (%f, %f)", csd1, csf1, csd2, csf2);

	if (debug) printf ("\n IPI 1");
	Quad* qr = intersectQuads(qp, qcs);
	if (debug) printf ("\n IPI 2");
	if (qr == 0) {
		printf ("\n ERROR! quad result==0");
		return;
	}
	if (debug) printf ("Quad result (%f, %f) (%f, %f)", qr->pd1, qr->pf1, qr->pd2, qr->pf2);
	PanelLayoutElement* ple = getPanelLayoutElementCopyWithNewQuad(p, qr);
	if (debug) printf ("\n IPI 3");
	if (debug) printf ("\nP  (%f(%d), %f(%d)) (%f(%d), %f(%d))", ple->posDeb1, ple->fd1,  ple->posFin1, ple->ff1,ple->posDeb2, ple->fd2, ple->posFin2, ple->ff2);
	if (debug) printf ("\n IPI 5");
	if (debug) printf ("\nresP  (%f(%d), %f(%d)) (%f(%d), %f(%d))", ple->posDeb1, ple->fd1, ple->posFin1, ple->ff1, ple->posDeb2, ple->fd2, ple->posFin2, ple->ff2);

	if (ple != 0 ) {
		if (debug) printf ("\n panelsIntDesign->push_back");
		panelsIntDesign.push_back(ple);
	}
}



void Layout::prepareDesignLayoutElements(WindPatternsProject* gfd) {
	bool debug = false;
	if (debug) printf ("\nprepareDesignLayoutElements()");
    // calculate ColorSegmentsTable here, and pass it (i)
    // to preparePanelExtDesign(gfd, i, KiteDesignExt . CST[i] // row );
    // to preparePanelIntDesign(gfd, i, KiteDesignInt . CST[i] // row );
    int n = gfd->Form->m_nbProfils;
	if (debug) printf ("\n gfd->kiteDesignExt->getColorSegmentsTable(%d)", n);
    ColorSegmentsTable* cst = gfd->kiteDesignExt->getColorSegmentsTable( n );
	
    // ------------------ Ext ------------------------------
	if (debug) printf ("\n panelsExt.size(): %d", panelsExt.size());
    for (int ipe = 0; ipe < panelsExt.size(); ipe++) {
	//for (int ipe = 3; ipe < 4; ipe++) {
		if (debug) printf("\n ipe: %d", ipe);
        PanelLayoutElement* p = panelsExt[ipe];
        int nerv = p->n1;
		if (debug) printf ("\nPLE nerv: %d", nerv);
		if (debug) printf ("\nPLE-- %f, %f, %f, %f ", p->fd1, p->ff1, p->fd2, p->ff2);

		if (nerv == -1) {
			if (debug) printf ("\ncontinue...");
			continue;
		}

        vector<ColorSegment*> vcs = cst->table[nerv];            
		if (debug) printf ("\n vcs.size()=%d", vcs.size());
		for (int i = 0; i < vcs.size(); i++) {
			ColorSegment* cs = vcs[i];
            intersectPanelExtWithColorSegment(p, cs);
        }
    }

	if (debug) printf ("\n result: ");
	if (debug) printf ("\n **************************************************************");
	if (debug) printf ("\n **************************************************************");
	if (debug) printf ("\n **************************************************************");
	for (int _i = 0; _i < panelsExtDesign.size(); _i++) {
		PanelLayoutElement* pd = panelsExtDesign[_i];
		
		if (debug) printf ("\n %d pd: {%d(%f, %d)-(%f, %d)  %d(%f, %d)-(%f, %d)}  ", _i, pd->n1, pd->posDeb1, pd->fd1, pd->posFin1, pd->ff1, pd->n2, pd->posDeb2, pd->fd2, pd->posFin2, pd->ff2);
	}
	if (debug) printf ("\n ***************************** ");
    
    // -----------------------------------------------------
	if (debug) printf ("\n INT INT INT INT INT");

	ColorSegmentsTable* csti = gfd->kiteDesignInt->getColorSegmentsTable( n );
	if (debug) printf ("\npanelsInt.size(): %d", panelsInt.size());
    //for (int ipi = 0; ipi < panelsInt.size(); ipi++) {
	for (int ipi = 0; ipi < 3; ipi++) {
		if (debug) printf("\n ipi: %d", ipi);
        PanelLayoutElement* p = panelsInt[ipi];
        int nerv = p->n1;
		if (nerv == -1) {
			if (debug) printf ("\ncontinue...");
			continue;
		}

        if (!p->isVentHole) {
			if (debug) printf ("\n!p->isVentHole");
            vector<ColorSegment*> vcs = csti->table[nerv];
			if (debug) printf ("\n int vcs.size()=%d", vcs.size());
		    for (int i = 0; i < vcs.size(); i++) {
			    ColorSegment* cs = vcs[i];
                intersectPanelIntWithColorSegment(p, cs);
            }
        } else {
			if (debug) printf ("\np->isVentHole");
            PanelLayoutElement* ple = getPanelLayoutElementCopy(p);
            panelsIntDesign.push_back(ple);
        }
    }
	//printf ("\n3");
	//printf ("\n...prepareDesignLayoutElements()");
}

void Layout::prepareLayoutElements(WindPatternsProject* gfd) {
	bool debug = false;
    if (debug) printf("\n prepareLayoutElements()");
    int i = 0, face = 0;
    int n = gfd->Form->m_nbProfils;
    float debBorder = 0.0f;
    int faceDebBorder = 0;
    if (gfd->VentHoles) {
        debBorder = gfd->VentHolesDeb;
        faceDebBorder = 2;
    } else {
        debBorder = 0.0f;
        faceDebBorder = 1;
    }

	this->debBorder = debBorder;
	this->faceDebBorder = faceDebBorder;
	
    if (isCenterPanel) {
        prepareCenterPanelExt(gfd);
        prepareCenterPanelInt(gfd);
    }


    for (i = 0; i < n - 1; i++) {
      preparePanelExt(gfd, i);
      preparePanelInt(gfd, i);

      prepareKlapan(gfd, i);
      prepareProfile(gfd, i);
    }
	
    if (gfd->DiagNervs) {
        for (i = 0; i < gfd->quantDiag; i++) {
          prepareDiagNerv(gfd, i);
        }
    }

    if (isDesign) {
        prepareDesignLayoutElements(gfd);
    }
    if (debug) printf("\n...prepareLayoutElements()");
}

void Layout::SaveLayout2(WindPatternsProject* gfd) {
	bool debug = false;
    if (debug) printf("\n SaveLayout2()");
    prepareLayoutElements(gfd);
	/* if (isDesign) {
		if (debug) printf("\n prepareDesignLayoutElements()");
		prepareDesignLayoutElements(gfd);
	} */
    calculateExport(gfd);
    saveCalcLayoutToFile(gfd);
    if (debug) printf("\n...SaveLayout2()");
}

void Layout::calculateExport(WindPatternsProject* gfd){
	bool debug = false;
    if (debug) printf("\n calculateExport()");
	if (debug) printf("\n calculateExport.1Panels");
    for (int i = 0; i < panelsExt.size(); i++) {
        panelsExt[i]->calculateExport(gfd);
    }

    for (int i = 0; i < panelsInt.size(); i++) {
        panelsInt[i]->calculateExport(gfd);
    }
		
	if (debug) printf("\n calculateExport.2PanelsDesign");
    if (isDesign) {
		if (debug) printf ("\n --- Design Ext panelsExtDesign.size()=%d", panelsExtDesign.size());
        for (int i = 0; i < panelsExtDesign.size(); i++) {
			if (debug) printf ("\next i=%d(%d)", i, panelsExtDesign.size());
            panelsExtDesign[i]->calculateExport(gfd);
        }
		if (debug) printf ("\n --- Design Int panelsIntDesign.size()=%d", panelsIntDesign.size());
        for (int i = 0; i < panelsIntDesign.size(); i++) {
			if (debug) printf ("\nint i=%d(%d)", i, panelsIntDesign.size());
            panelsIntDesign[i]->calculateExport(gfd);
        }
    }
	if (debug) printf("\n calculateExport.3diagNervs");
    for (int i = 0; i < diagNervs.size(); i++) {
        diagNervs[i]->calculateExport(gfd);
    }
	if (debug) printf("\n calculateExport.4profs");
    for (int i = 0; i < profs.size(); i++) {
        profs[i]->calculateExport(gfd);
    }
    if (debug) printf("\n...calculateExport()");
}

void calculateLayout(WindPatternsProject* gfd, Layout* layout) {
    /*
    TAxe **AxeP, **AxePD, **AxePTD, **AxeMD, **AxeCD, **AxeRepD;
    double *H, *W;
    int n = gfd->Form->m_nbProfils;
    char text[100];
    double _pa0, _pa00, _pf0, _pa1, _pa01, _pf1;
    Matrix *_f0, *_f1;
    bool debug = false;
    char* charname = new char[1];
    */

}

//void saveCalcLayoutToFile(int n, int q, TAxe **AxePD, TAxe **AxeMD, TAxe **RepD, int gfdVentilationLayout, TAxe **AxeCD, TAxe  **AxePTD, double *H, double *W, int* numncol) {

void Layout::writeLayoutWithDesignToDXF(char *fileName, WindPatternsProject* gfd) {
	bool debug = false;
    FILE *fid;
	if (debug) printf("\nwrite MANY POLY fichier de DXF: '%s'",fileName);
	if( (fid = fopen( fileName, "wt" )) == NULL ) {
		printf( "\nError write file '%s'", fileName);
		exit(0);
	}
    fprintf(fid,"0\nSECTION\n2\nENTITIES\n");
    double maxW1 = -1000000.0f, maxW2 = -1000000.0f, maxW3 = -1000000.0f;
    //panelsExt
    double dx = 0.0f, dy = 0.0f, dz = 0.0f;
    for (int i = 0; i < panelsExtDesign.size(); i++) {
        LayoutElementExport* lee = panelsExtDesign[i]->leexport;
        writeFichierPolyDXFDelta(fid, lee->AxePD, lee->AxeMD,
            1, lee->AxeRepD, gfd->VentilationLayout, lee->AxeCD, 1, lee->AxePTD, dx, dy, dz, 0);
        dy = dy + lee->H * 2.0f;
        if (lee->W > maxW1) maxW1 = lee->W;
    }

    //profs
    dx = maxW1*2.0; dy = 0.0f; maxW2 = -1000000.0f; maxW3 = -1000000.0f;
    for (int i = 0; i < profs.size(); i++) {
        LayoutElementExport* lee = profs[i]->leexport;
        writeFichierPolyDXFDelta(fid, lee->AxePD, lee->AxeMD,
            1, lee->AxeRepD, gfd->VentilationLayout, lee->AxeCD, 1, lee->AxePTD, dx, dy, dz, 0);
        dy = dy + lee->H * 2.0f;
        if (lee->W > maxW2) maxW2 = lee->W;
    }

    //diagNervs
    dx = maxW1 * 2.0 + maxW2 * 4.0; dy = 0.0f;
    for (int i = 0; i < diagNervs.size(); i++) {
        LayoutElementExport* lee = diagNervs[i]->leexport;
        writeFichierPolyDXFDelta(fid, lee->AxePD, lee->AxeMD,
            1, lee->AxeRepD, gfd->VentilationLayout, lee->AxeCD, 1, lee->AxePTD, dx, dy, dz, 0);
        dy = dy + lee->H * 2.0f;
        if (lee->W > maxW3) maxW3 = lee->W;
    }

    //panelsInt
    dx = maxW1 * 2.0 + maxW2 * 4.0 + maxW3 * 4.0; dy = 0.0f;
    for (int i = 0; i < panelsIntDesign.size(); i++) {
        LayoutElementExport* lee = panelsIntDesign[i]->leexport;
        writeFichierPolyDXFDelta(fid, lee->AxePD, lee->AxeMD,
            1, lee->AxeRepD, gfd->VentilationLayout, lee->AxeCD, 1, lee->AxePTD, dx, dy, dz, 0);
        dy = dy + lee->H * 2.0f;
    }

    fprintf(fid,"0\nENDSEC\n0\nEOF\n");
    if(fclose(fid))	{
	    printf("\nError close file");
	    exit(0);
    }
}


void Layout::writeLayoutToDXF(char *fileName, WindPatternsProject* gfd) {
	bool debug = false;
    FILE *fid;
	if (debug) printf("\nwrite MANY POLY fichier de DXF: '%s'",fileName);
	if( (fid = fopen( fileName, "wt" )) == NULL ) {
		printf( "\nError write file '%s'", fileName);
		exit(0);
	}
    fprintf(fid,"0\nSECTION\n2\nENTITIES\n");
    double maxW1 = -1000000.0f, maxW2 = -1000000.0f, maxW3 = -1000000.0f;
    //panelsExt
    double dx = 0.0f, dy = 0.0f, dz = 0.0f;
    for (int i = 0; i < panelsExt.size(); i++) {
        LayoutElementExport* lee = panelsExt[i]->leexport;
        writeFichierPolyDXFDelta(fid, lee->AxePD, lee->AxeMD,
            1, lee->AxeRepD, gfd->VentilationLayout, lee->AxeCD, 1, lee->AxePTD, dx, dy, dz, 0);
        dy = dy + lee->H * 2.0f;
        if (lee->W > maxW1) maxW1 = lee->W;
    }

    //profs
    dx = maxW1*2.0; dy = 0.0f; maxW2 = -1000000.0f; maxW3 = -1000000.0f;
    for (int i = 0; i < profs.size(); i++) {
        LayoutElementExport* lee = profs[i]->leexport;
        writeFichierPolyDXFDelta(fid, lee->AxePD, lee->AxeMD,
            1, lee->AxeRepD, gfd->VentilationLayout, lee->AxeCD, 1, lee->AxePTD, dx, dy, dz, 0);
        dy = dy + lee->H * 2.0f;
        if (lee->W > maxW2) maxW2 = lee->W;
    }

    //diagNervs
    dx = maxW1 * 2.0 + maxW2 * 4.0; dy = 0.0f;
    for (int i = 0; i < diagNervs.size(); i++) {
        LayoutElementExport* lee = diagNervs[i]->leexport;
        writeFichierPolyDXFDelta(fid, lee->AxePD, lee->AxeMD,
            1, lee->AxeRepD, gfd->VentilationLayout, lee->AxeCD, 1, lee->AxePTD, dx, dy, dz, 0);
        dy = dy + lee->H * 2.0f;
        if (lee->W > maxW3) maxW3 = lee->W;
    }

    //panelsInt
    dx = maxW1 * 2.0 + maxW2 * 4.0 + maxW3 * 4.0; dy = 0.0f;
    for (int i = 0; i < panelsInt.size(); i++) {
        LayoutElementExport* lee = panelsInt[i]->leexport;
        writeFichierPolyDXFDelta(fid, lee->AxePD, lee->AxeMD,
            1, lee->AxeRepD, gfd->VentilationLayout, lee->AxeCD, 1, lee->AxePTD, dx, dy, dz, 0);
        dy = dy + lee->H * 2.0f;
    }

    fprintf(fid,"0\nENDSEC\n0\nEOF\n");
    if(fclose(fid))	{
	    printf("\nError close file");
	    exit(0);
    }
}


void Layout::saveCalcLayoutToFile(WindPatternsProject* gfd) {
    //TAxe **AxeP, **AxePD, **AxePTD, **AxeMD, **AxeCD, **AxeRepD;
    CString fileName;
    LPTSTR PtrfileName;
    char ext[10];
    strcpy(ext, "*.dxf");
    CFileDialog DlgOpen(FALSE, NULL, ext, OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        fileName = DlgOpen.GetPathName();
        PtrfileName = fileName.GetBuffer(1);
        writeLayoutToDXF(PtrfileName, gfd);
        if (isDesign) {
            char filenameLayoutDesign[255];
            strcpy(filenameLayoutDesign, PtrfileName);
            strcat (filenameLayoutDesign,"_design.dxf");
            writeLayoutWithDesignToDXF(filenameLayoutDesign, gfd);
        }

        //writeManyFichierPolyDXF2(PtrfileName, n, q, AxePD, AxeMD, 1, AxeRepD, gfd->VentilationLayout, AxeCD, 1, AxePTD, W, H, numncol);
        char filename[255];
        strcpy(filename, PtrfileName);
        strcat (filename,".log");
        //printf ("\n logfilename=[%s]", filename);
        //printf ("\n gfd->logfilename=[%s]", gfd->logFileName);
        //printf ("\n ptrnomfichier=[%s]", PtrfileName);
        CopyFile(gfd->logFileName, filename, false);
    }

    /*for (i = 0; i < q; i++) {
        clearCourbesAxe(AxeP[i]);
        clearCourbesAxe(AxePD[i]);
        clearCourbesAxe(AxePTD[i]);
        clearCourbesAxe(AxeMD[i]);
        clearCourbesAxe(AxeCD[i]);
        clearCourbesAxe(AxeRepD[i]);
    }*/

}


void printXY(Matrix* X, Matrix* Y) {
	//printf ("\n printXY");
	for (int i = 0; i < X->GetLignes(); i++) {
		printf ("\n$ %d (%f, %f)", i, X->Element(i, 0), Y->Element(i, 0));
	}
}

void printM(Matrix* m) {
	printf ("\n printM");
	for (int i = 0; i < m->GetLignes(); i++) {
		printf ("\n$ %d (%f)", i, m->Element(i, 0));
	}
}


