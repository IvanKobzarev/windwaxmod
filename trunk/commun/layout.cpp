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
}

LayoutElementExport::LayoutElementExport(){
}
LayoutElementExport::~LayoutElementExport(){
}

KlapanLayoutElement::KlapanLayoutElement(){
}
KlapanLayoutElement::~KlapanLayoutElement(){
}

ProfLayoutElement::ProfLayoutElement(){
}
ProfLayoutElement::~ProfLayoutElement(){
}

PanelLayoutElement::PanelLayoutElement(){
}
PanelLayoutElement::~PanelLayoutElement(){
}

DiagNervLayoutElement::DiagNervLayoutElement(){
}
DiagNervLayoutElement::~DiagNervLayoutElement(){
}


void KlapanLayoutElement::calculateExport(WindPatternsProject* gfd){
    printf ("\n KlapanLayoutElement::calculateExport");
    bool debug = 0;
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
    printf ("\n ProfLayoutElement::calculateExport");
    bool debug = 0;
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
    printf ("\n PanelLayoutElement::calculateExport");
    bool debug = 0;
    Matrix * Xd[2], *Yd[2],*newXd[2], *newYd[2], *Xdp[2], *Ydp[2], *rXdp[2], *rYdp[2], *rXd[2], *rYd[2];
    Matrix *X[2], *Y[2], *Z[2], *P[2], *rP[2];//, *newP[2];
    double _pa0, _pa00, _pf0, _pa1, _pa01, _pf1;
    Matrix *_f0, *_f1;
    char* charname = new char[1];

    calcPatron(gfd, n1, s1, fd1, ff1, posDeb1, posFin1,
        n2, s2, fd2, ff2, posDeb2, posFin2,
        &Xd[0], &Yd[0], &Xd[1], &Yd[1],
        &X[0], &Y[0], &Z[0], &P[0],
        &X[1], &Y[1], &Z[1], &P[1]);
    printf ("\n n1=%d fd1=%d ff1=%d posDeb1=%f posFin1=%f", n1, fd1, ff1, posDeb1, posFin1);
    printf ("\n n2=%d fd2=%d ff2=%d posDeb2=%f posFin2=%f", n2, fd2, ff2, posDeb2, posFin2);
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
        double myFin = posFin1;
         if (fd1 == ff1) {
            if (posDeb1 < posDeb2) myDeb = posDeb1; else myDeb = posDeb2;
        } else {
            if (posDeb1 > posDeb2) myDeb = posDeb1; else myDeb = posDeb2;
        }
        if (debug) printf ("calcPatron myDeb myFin");
        calcPatron(gfd, n1, s1, fd1, ff1, myDeb, myFin,
                n2, s2, fd2, ff2, myDeb, myFin,
                &Xd[0], &Yd[0], &Xd[1], &Yd[1],
                &X[0], &Y[0], &Z[0], &P[0],
                &X[1], &Y[1], &Z[1], &P[1]);
        if (debug) printf ("...calcPatron myDeb myFin");
        getLayoutLogger()->logprintf("p[%4.1f %4.1f %4.1f %4.1f]", _pa0 * 1000, _pf0 * 1000, _pa1 * 1000, _pf1 * 1000);
        getLayoutLogger()->logprintf(" %d %d [%5.1f(%d) %5.1f(%d)] [%5.1f(%d) %5.1f(%d)]", n1, n2, posDeb1, fd1, posFin1,  ff1, posDeb2, fd2, posFin2,  ff2);
        Pince* pince = new Pince();
        pince -> debug = false;
        Pince* tmpPince = getPincePlus(gfd, n1, n2, myDeb, myFin,  myDeb, myFin,  fd1, ff1);
        pince -> SetFunction1(tmpPince->function1);
        pince -> SetFunction2(tmpPince->function2);
        pince -> ipa1 = tmpPince -> ipa1;
        pince -> ipf1 = tmpPince -> ipf1;
        pince -> ipa2 = tmpPince -> ipa2;
        pince -> ipf2 = tmpPince -> ipf2;
        pince -> SetPos(gfd->PosPinceBA[0], gfd->PosPinceBF[0]);
        pince -> SetDiffAmps(_pa0, _pf0, _pa1, _pf1);
        pince -> SetDiffAmps0(_pa00, _pa01);
        pince -> i01 = tmpPince -> i01;
        pince -> i02 = tmpPince -> i02;
        //printf("\n pince->i01=%d", pince->i01);
        //printf("\n pince->ipa1=%d", pince->ipa1);
        //printf("\n pince->ipf1=%d", pince->ipf1);
        calcPincePlusNew(Xd[0], Yd[0], Xd[1], Yd[1], pince, &Xdp[0], &Ydp[0], &Xdp[1], &Ydp[1]);

            if ((myDeb != posDeb1)||(myDeb != posDeb2)) {
                // rezem Xd, newXd
                rXdp[0] = GetFunctionSrezDeb(P[0], Xdp[0], posDeb1);
                rYdp[0] = GetFunctionSrezDeb(P[0], Ydp[0], posDeb1);
                delete(Xdp[0]);
                Xdp[0]=rXdp[0];
                delete(Ydp[0]);
                Ydp[0]=rYdp[0];


                rXdp[1] = GetFunctionSrezDeb(P[1], Xdp[1], posDeb2);
                rYdp[1] = GetFunctionSrezDeb(P[1], Ydp[1], posDeb2);
                delete(Xdp[1]);
                Xdp[1]=rXdp[1];
                delete(Ydp[1]);
                Ydp[1]=rYdp[1];


                rXd[0] = GetFunctionSrezDeb(P[0], Xd[0], posDeb1);
                rYd[0] = GetFunctionSrezDeb(P[0], Yd[0], posDeb1);
                delete(Xd[0]);
                Xd[0]=rXd[0];
                delete(Yd[0]);
                Yd[0]=rYd[0];


                rXd[1] = GetFunctionSrezDeb(P[1], Xd[1], posDeb2);
                rYd[1] = GetFunctionSrezDeb(P[1], Yd[1], posDeb2);
                delete(Xd[1]);
                Xd[1]=rXd[1];
                delete(Yd[1]);
                Yd[1]=rYd[1];

                rP[0] = GetFunctionSrezDeb(P[0], P[0], posDeb1);
                delete(P[0]);
                P[0]=rP[0];

                rP[1] = GetFunctionSrezDeb(P[1], P[1], posDeb2);
                delete(P[1]);
                P[1]=rP[1];
            }

        //delete(P[0]);
        //delete(P[1]);
        //P[0] = newP[0];
        //P[1] = newP[1];
    }
    printf ("\n PanelLayoutElement::calculateExport.1");
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
    printf ("\n PanelLayoutElement::calculateExport.2");

    int n = gfd->Form->m_nbProfils;
    char text[100];
    sprintf(text, "%s%dF%dt%dF%d", charname, n1, ff1, n2, ff2);
    TAxe *AxeP, *AxePD,  *AxePTD,  *AxeMD,  *AxeCD,  *AxeRepD;
    printf ("\n PanelLayoutElement::calculateExport.3");
    leexport = new LayoutElementExport();
    GenerateCourbe(gfd, Xdp[0], Ydp[0], P[0], n1, posDeb1, fd1, posFin1, ff1,
        Xdp[1], Ydp[1], P[1], n2, posDeb2, fd2, posFin2, ff2, text,
        &(leexport->AxeP), &(leexport->AxePD), &(leexport->AxePTD), &(leexport->AxeMD), &(leexport->AxeCD), &(leexport->AxeRepD), vent, marge1, marge2, margeDeb, margeFin, true, debug,
		true, Xd[0], Yd[0], coeff1, Xd[1], Yd[1], coeff2);
    printf ("\n PanelLayoutElement::calculateExport.3.2");

    printf ("\n PanelLayoutElement::calculateExport.4");
    calcMaxWH(Xdp[0], Ydp[0], Xdp[1], Ydp[1], &(leexport->W), &(leexport->H));
    printf ("\n PanelLayoutElement::calculateExport.5");
    delete (Xdp[0]);
    delete (Ydp[0]);
    delete (Xdp[1]);
    delete (Ydp[1]);
    printf ("\n... PanelLayoutElement::calculateExport");
}
void DiagNervLayoutElement::calculateExport(WindPatternsProject* gfd) {
    printf ("\n DiagNervLayoutElement::calculateExport");

    bool debug = 0;
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
        true, true);//debug);

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
    if ( gfd->debug ) printf("\n calcIndepPinceLayout()");

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
    layout->isCenterPanel = isCenterPanel;
    int startNerv;
    //printf ("\n 1");
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
	//printf ("\n 2");
    for (int inoNerv = startNerv; inoNerv < n - 1; inoNerv++) {
        goCalcIndepPinceNew(gfd, inoNerv, 1,
                &(layout->pinceLAAmp1[inoNerv]), &(layout->pinceLFAmp1[inoNerv]), &(layout->funcL1[inoNerv]),
                &(layout->pinceRAAmp1[inoNerv]), &(layout->pinceRFAmp1[inoNerv]), &(layout->funcR1[inoNerv]),
                &(layout->lenp1[inoNerv]));
    }
	//printf ("\n 3");
    getLayoutLogger()->logprintf("\n");
    for (int inoNerv = startNerv; inoNerv < n - 1; inoNerv++) {
        goCalcIndepPinceNew(gfd, inoNerv, 2,
                &(layout->pinceLAAmp2[inoNerv]), &(layout->pinceLFAmp2[inoNerv]), &(layout->funcL2[inoNerv]),
                &(layout->pinceRAAmp2[inoNerv]), &(layout->pinceRFAmp2[inoNerv]), &(layout->funcR2[inoNerv]),
                &(layout->lenp2[inoNerv]));
    }
	//printf ("\n 4");
    //if (DEBUG) printf("\n");

    // calc Profile Nervures
    for (int inoNerv = 0; inoNerv < n - 1; inoNerv++) {
        goCalcNervureWithPince(gfd, inoNerv, 2, inoNerv, 1,
                layout->lenp2[inoNerv], layout->lenp1[inoNerv], &(layout->coeffn[inoNerv]));
    }
    //if (DEBUG) printf("\n");
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

void Layout::preparePanelIntDesign(WindPatternsProject* gfd, int i) {
    //TODO: me
}

void Layout::preparePanelExtDesign(WindPatternsProject* gfd, int i) {
    //TODO: me
}

void Layout::prepareCenterPanelIntDesign(WindPatternsProject* gfd) {
    //TODO: me
}

void Layout::prepareCenterPanelExtDesign(WindPatternsProject* gfd) {
    //TODO: me
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

        le1->func2f0 = funcR2[i];
        le1->func2f1 = funcL2[i + 1];

        le1->fd1 = 2;
        le1->fd2 = 2;
        le1->ff1 = 2;
        le1->ff2 = 2;
        le1->isPince = 1;
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

    le2->func2f0 = funcR2[i];
    le2->func2f1 = funcL2[i + 1];

    le2->fd1 = 2;
    le2->fd2 = 2;
    le2->ff1 = 2;
    le2->ff2 = 2;
    le2->isPince = 1;
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

    le->ff1 = 1;
    le->ff2 = 1;
    le->fd1 = faceDebBorder;
    le->fd2 = faceDebBorder;

    le->isPince = 1;
	le->isKlapan = 0;
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


PanelLayoutElement* getPanelLayoutElementCopyWithNewBorders(PanelLayoutElement* p, double d1, double f1, double d2, double f2) {
    PanelLayoutElement* ple = new PanelLayoutElement();
    ple->side=p->side;
    ple->n1 = p->n1;
    ple->n2 = p->n2;
    ple->s1 = false;
    ple->s2 = false;
    ple->posDeb1 = d1;
    ple->posFin1 = f1;
    ple->posDeb2 = d2;
    ple->posFin2 = f2;

    /*
    ple->p1a00 = pinceRAAmp2 [i];
    ple->p1a0 = pinceRAAmp1 [i];
    ple->p1f0 = pinceRFAmp1 [i];

    ple->p1a01 = pinceLAAmp2 [i + 1];
    ple->p1a1 = pinceLAAmp1 [i + 1];
    ple->p1f1 = pinceLFAmp1 [i + 1];
    */

    ple->func1f0 = funcR1[i];
    ple->func1f1 = funcL1[i + 1];

    ple->ff1 = p->ff1;
    ple->ff2 = p->ff2;
    ple->fd1 = p->fd1;
    ple->fd2 = p->fd2;

    ple->isPince = p->isPince;
    ple->isKlapan = p->isKlapan;
    ple->coeff = p->coeff;
    return ple;
}
void Layout::intersectPanelExtWithColorSegment(PanelLayoutElement* p, ColorSegment* cs) {
    double pd1 = p->posDeb1;
    double pf1 = p->posFin1;
    double pd2 = p->posDeb2;
    double pf2 = p->posFin2;

    double csd1 = cs->p00;
    double csf1 = cs->p01;
    double csd2 = cs->p10;
    double csf2 = cs->p11;

    PanelLayoutElement* ple = 0;

    if ((csd1 <= pd1) && (csd2 <= pd2)) {
        if ((pf1  <= csf1) && (pf2 <= csf2)) {
            // ?
            // pf1  pf2 
            // pd1  pd2
            if ((pd1 < pf1) && (pd2 < pf2)){
                ple = getPanelLayoutElementCopyWithNewBorders(p, pd1, pf1,pd2, pf2);
            }
        }

        if ((csf1  <= pf1) && (csf2 <= pf2)) {
            // ?
            // csf1  csf2 
            // pd1   pd2
            if ((pd1 < csf1) && (pd2 < csf2)){
                ple = getPanelLayoutElementCopyWithNewBorders(p, pd1, csf1,pd2, csf2);
            }
        }
    }

    if ((pd1 <= csd1) && (pd2 <= csd2)) {
        if ((pf1  <= csf1) && (pf2 <= csf2)) {
            // ?
            // pf1  pf2 
            // csd1  csd2
            if ((csd1 < pf1) && (csd2 < pf2)){
                ple = getPanelLayoutElementCopyWithNewBorders(p, csd1, pf1,csd2, pf2);
            }
        }

        if ((csf1  <= pf1) && (csf2 <= pf2)) {
            // ?
            // csf1  csf2 
            // csd1   csd2
            if ((csd1 < csf1) && (csd2 < csf2)){
                ple = getPanelLayoutElementCopyWithNewBorders(p, csd1, csf1, csd2, csf2);
            }
        }
    }

    panelsExtDesign.push_back(ple);
}

void Layout::prepareDesignLayoutElements(WindPatternsProject* gfd) {
    // calculate ColorSegmentsTable here, and pass it (i)
    // to preparePanelExtDesign(gfd, i, KiteDesignExt . CST[i] // row );
    // to preparePanelIntDesign(gfd, i, KiteDesignInt . CST[i] // row );
    int n = gfd->Form->m_nbProfils;
    
    ColorSegmentsTable* cst = gfd->kiteDesignExt->getColorSegmentsTable( n );

    // ------------------ Ext ------------------------------

    for (int ipe = 0; ipe < panelsExt.size(); ipe++) {
        PanelLayoutElement* p = panelsExt[ipe];

        int nerv = p->n1;
        vector<ColorSegment*> vcs = cst->table[nerv];            
		for (int i = 0; i < vcs.size(); i++) {
			ColorSegment* cs = vcs[i];
            intersectPanelExtWithColorSegment(p, cs);

        }
    }
    
    // -----------------------------------------------------



    for (int i = 0; i < panelsInt.size(); i++) {
        panelsInt[i]->calculateExport(gfd);
    }

    // preparePanelExtDesign(gfd, i);
    // preparePanelIntDesign(gfd, i);

    // prepareCenterPanelExtDesign(gfd);
    // prepareCenterPanelIntDesign(gfd);
}

void Layout::prepareLayoutElements(WindPatternsProject* gfd) {
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

	debBorder = debBorder;
	faceDebBorder = faceDebBorder;

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

}

void Layout::SaveLayout2(WindPatternsProject* gfd) {
    if (gfd->debug) printf("\n SaveLayout2()");
    prepareLayoutElements(gfd);
    if (isDesign) prepareDesignLayoutElements(gfd);

    calculateExport(gfd);

    saveCalcLayoutToFile(gfd);
}

void Layout::calculateExport(WindPatternsProject* gfd){
    for (int i = 0; i < panelsExt.size(); i++) {
        panelsExt[i]->calculateExport(gfd);
    }
    for (int i = 0; i < panelsInt.size(); i++) {
        panelsInt[i]->calculateExport(gfd);
    }

    if (isDesign) {
        for (int i = 0; i < panelsExt.size(); i++) {
            panelsExtDesign[i]->calculateExport(gfd);
        }
        for (int i = 0; i < panelsInt.size(); i++) {
            panelsIntDesign[i]->calculateExport(gfd);
        }
    }

    for (int i = 0; i < diagNervs.size(); i++) {
        diagNervs[i]->calculateExport(gfd);
    }
    for (int i = 0; i < profs.size(); i++) {
        profs[i]->calculateExport(gfd);
    }
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
    FILE *fid;
	printf("\nwrite MANY POLY fichier de DXF: '%s'",fileName);
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


void Layout::writeLayoutToDXF(char *fileName, WindPatternsProject* gfd) {
    FILE *fid;
	printf("\nwrite MANY POLY fichier de DXF: '%s'",fileName);
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


