#pragma warning(disable:4514)

#include "afx.h"		//class CString
#include "afxdlgs.h"	//class CFileDialog

#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>

#include "layout.h"
#include "geom.h"
#include "fichier.h"
#include "plot.h"
#include "pince.h"
#include "profil.h"
#include "matrice.h"
#include "logger.h"


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

void KlapanLayoutElement::calculateExport(WindPatternsProject* gfd){
    bool debug = 0;
    Matrix * Xd[2], *Yd[2];//,*newXd[2], *newYd[2], *Xdp[2], *Ydp[2], *rXdp[2], *rYdp[2], *rXd[2], *rYd[2];
    Matrix *X[2], *Y[2], *Z[2], *P[2];//, *rP[2];//, *newP[2];

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
        &(leexport->AxeP), &&(leexport->AxePD), &(leexport->AxePTD), &(leexport->AxeMD), &(leexport->AxeCD), &(leexport->AxeRepD), vent, marge1, marge2, margeDeb, margeFin, false, debug);
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
    bool debug = 0;
    Matrix * Xd[2], *Yd[2],*newXd[2], *newYd[2];//, *Xdp[2], *Ydp[2];//, *rXdp[2], *rYdp[2], *rXd[2], *rYd[2];
    Matrix *X[2], *Y[2], *Z[2], *P[2];//, *rP[2];//, *newP[2];

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
    bool debug = 0;
    Matrix * Xd[2], *Yd[2],*newXd[2], *newYd[2], *Xdp[2], *Ydp[2], *rXdp[2], *rYdp[2], *rXd[2], *rYd[2];
    Matrix *X[2], *Y[2], *Z[2], *P[2], *rP[2];//, *newP[2];

    calcPatron(gfd, n1, s1, fd1, ff1, posDeb1, posFin1,
        n2, s2, fd2, ff2, posDeb2, posFin2,
        &Xd[0], &Yd[0], &Xd[1], &Yd[1],
        &X[0], &Y[0], &Z[0], &P[0],
        &X[1], &Y[1], &Z[1], &P[1]);

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

        double marge1 = gfd->Marge[0];
        double marge2 = gfd->Marge[1];
        double margeDeb = gfd->MargeDeb;
        double margeFin = gfd->MargeFin;

        /*
        if (t < numncol[1]) {
            // panel exterior
            margeFin = gfd->margeFinExt;
            charname="P";
        } else {
            if (t < numncol[2]) {
            } else {
                if (t < numncol[3]) {
                } else {
                    // panels interior
                    margeFin = gfd->margeFinInt;
                    charname="P";
                    if (gfd->VentHoles) 
                        if ((n1[t]==-1) && gfd->VentCentralNerv) {
                             if (posFin1[t]!=100.0f) {
                                margeFin = gfd->MargeFin;
                                charname="V";
                            } 
                        } else
                            if (gfd->noNervVH[n1[t]]) {
                               if (posFin1[t]!=100.0f) {
                                    margeFin = gfd->MargeFin;
                                    charname="V";
                               }
                            }
                    
                }
            }
        }
*/
    int n = gfd->Form->m_nbProfils;
    char text[100];
    sprintf(text, "%s%dF%dt%dF%d", charname, n1, ff1, n2, ff2);

    GenerateCourbe(gfd, Xdp[0], Ydp[0], P[0], n1, posDeb1, fd1, posFin1, ff1,
        Xdp[1], Ydp[1], P[1], n2, posDeb2, fd2, posFin2, ff2, text,
        &(leexport->AxeP), &(leexport->AxePD), &(leexport->AxePTD), &(leexport->AxeMD), &(leexport->AxeCD), &(leexport->AxeRepD), vent, marge1, marge2, margeDeb, margeFin, true, debug,
		true, Xd[0], Yd[0], coeff1, Xd[1], Yd[1], coeff2);

    calcMaxWH(Xdp[0], Ydp[0], Xdp[1], Ydp[1], &(leexport->W), &(leexport->H));

    if (debug) printf("\ngo delete isPince");
    delete (Xdp[0]);
    delete (Ydp[0]);
    delete (Xdp[1]);
    delete (Ydp[1]);
    if (debug) printf("\n... go delete isPince");
}
void DiagNervLayoutElement::calculateExport(WindPatternsProject* gfd) {
    bool debug = 0;
    Matrix * Xd[2], *Yd[2],*newXd[2], *newYd[2];
    Matrix *X[2], *Y[2], *Z[2], *P[2];

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


void preparePanelLayoutElementCoeff(LayoutElement* le, WindPatternsProject* gfd, Layout* layout) {
   	double coeff1 = layout->coeffn[le->n1];
	double coeff2 = layout->coeffn[le->n2];
    int n = gfd->Form->m_nbProfils;
	if (le->n1 == -1) coeff1 = layout->coeffn[0];
	if (le->n2 == n-1) coeff2 = -1;
    le->coeff1=coeff1;
    le->coeff2=coeff2;

}

void preparePanelInt(WindPatternsProject* gfd, Layout* layout, int i) {
    LayoutElement* le1 = new LayoutElement();
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
        le1->p2a0 = layout->pinceRAAmp2 [i];
        le1->p2f0 = layout->pinceRFAmp2 [i];

        le1->p2a01 = 0.0f;
        le1->p2a1 = layout->pinceLAAmp2 [i + 1];
        le1->p2f1 = layout->pinceLFAmp2 [i + 1];

        le1->func2f0 = layout->funcR2[i];
        le1->func2f1 = layout->funcL2[i + 1];

        le1->fd1 = 2;
        le1->fd2 = 2;
        le1->ff1 = 2;
        le1->ff2 = 2;
        le1->isPince = 1;
        le1->coeff = 0.0f;
        preparePanelLayoutElementCoeff(le1, gfd, layout);
        layout->panelsInt.push_back(le1);
    }

    LayoutElement* le2 = new LayoutElement();
    le2->n1 = i;
    le2->n2 = i + 1;
    le2->s1 = false;
    le2->s2 = false;
    le2->posDeb1 = debBorder;
    le2->posFin1 = tpF1;
    le2->posDeb2 = debBorder;
    le2->posFin2 = tpF2;

    le2->p2a00 = 0.0f;
    le2->p2a0 = layout->pinceRAAmp2 [i];
    le2->p2f0 = layout->pinceRFAmp2 [i];

    le2->p2a01 = 0.0f;
    le2->p2a1 = layout->pinceLAAmp2 [i + 1];
    le2->p2f1 = layout->pinceLFAmp2 [i + 1];

    le2->func2f0 = layout->funcR2[i];
    le2->func2f1 = layout->funcL2[i + 1];

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
    preparePanelLayoutElementCoeff(le2, gfd, layout);
    layout->panelsInt.push_back(le2);

    if ((gfd->VentHoles) && (gfd->noNervVH[i])) {
        LayoutElement* le3 = new LayoutElement();
        le3->n1 = i;
        le3->n2 = i + 1;
        le3->s1 = false;
        le3->s2 = false;
        le3->posDeb1 = gfd->VentHolesFin;
        le3->posFin1 = 100.0f;
        le3->posDeb2 = gfd->VentHolesFin;
        le3->posFin2 = 100.0f;

        le3->p2a00 = 0.0f;
        le3->p2a0 = layout->pinceRAAmp2 [i];
        le3->p2f0 = layout->pinceRFAmp2 [i];

        le3->func2f0 = layout->funcR2[i];
        le3->func2f1 = layout->funcL2[i + 1];

        le3->p2a01 = 0.0f;
        le3->p2a1 = layout->pinceLAAmp2 [i + 1];
        le3->p2f1 = layout->pinceLFAmp2 [i + 1];

        le3->fd1 = 2;
        le3->fd2 = 2;
        le3->ff1 = 2;
        le3->ff2 = 2;

        le3->isPince = 1;
        le3->coeff = 0.0f;

        if (i == (n - 2)) {
            le3->func2f1 = Zeros((le3->func2f0)->GetLignes(), 1);
            le3->p2a01 = 0.0f;
            le3->p2a1 = 0.0f;
            le3->p2f1 = 0.0f;
        }
        preparePanelLayoutElementCoeff(le3, gfd, layout);
        layout->panelsInt.push_back(le3);

        if (gfd->LayoutKlapans) {
            if (gfd->VentHolesDouble){
                LayoutElement* le4 = new LayoutElement();
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
                preparePanelLayoutElementCoeff(le4, gfd, layout);
                layout->panelsInt.push_back(le4);
            }
                LayoutElement* le5 = new LayoutElement();
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
                preparePanelLayoutElementCoeff(le5, gfd, layout);
                layout->panelsInt.push_back(le5);

               LayoutElement* le6 = new LayoutElement();
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
                preparePanelLayoutElementCoeff(le6, gfd, layout);
                layout->panelsInt.push_back(le6);

        }

    }
}

void preparePanelExt(WindPatternsProject* gfd, Layout* layout, int i) {
    // face == 1
    LayoutElement* le = new LayoutElement();
    le->n1 = i;
    le->n2 = i + 1;
    le->s1 = false;
    le->s2 = false;
    le->posDeb1 = debBorder;
    le->posFin1 = 100.0f;
    le->posDeb2 = debBorder;
    le->posFin2 = 100.0f;

    le->p1a00 = layout->pinceRAAmp2 [i];
    le->p1a0 = layout->pinceRAAmp1 [i];
    le->p1f0 = layout->pinceRFAmp1 [i];

    le->p1a01 = layout->pinceLAAmp2 [i + 1];
    le->p1a1 = layout->pinceLAAmp1 [i + 1];
    le->p1f1 = layout->pinceLFAmp1 [i + 1];

    le->func1f0 = layout->funcR1[i];
    le->func1f1 = layout->funcL1[i + 1];

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
    preparePanelLayoutElementCoeff(le, gfd, layout);
    layout->panelsExt.push_back(le);
}


void prepareCenterPanelInt(WindPatternsProject* gfd, Layout* layout) {
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
            LayoutElement* le = new LayoutElement();
            le->n1 = -1;
            le->n2 = 0;
            le->s1 = false;
            le->s2 = false;
            le->posDeb1 = 0.0f;
            le->posFin1 = debBorder;
            le->posDeb2 = 0.0f;
            le->posFin2 = debBorder;

            le->p2a00 = 0.0f;
            le->p2a0 = layout->pinceRAAmp2 [i];
            le->p2f0 = layout->pinceRFAmp2 [i];

            le->p2a01 = 0.0f;
            le->p2a1 = layout->pinceLAAmp2 [i + 1];
            le->p2f1 = layout->pinceLFAmp2 [i + 1];

            le->func2f0 = layout->funcR2[i];
            le->func2f1 = layout->funcL2[i + 1];

            le->fd1 = 2;
            le->fd2 = 2;
            le->ff1 = 2;
            le->ff2 = 2;
            le->isPince = 1;
            le->coeff = 0.0f;
            layout->panelsInt.push_back(le);
        }
        LayoutElement* le1 = new LayoutElement();
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
        le1->p2a0 = layout->pinceLAAmp2[0];
        le1->p2f0 = layout->pinceLFAmp2[0];

        le1->p2a01 = 0.0f;
        le1->p2a1 = layout->pinceLAAmp2[0];
        le1->p2f1 = layout->pinceLFAmp2[0];

        le1->func2f0 = layout->funcL2[0];
        le1->func2f1 = layout->funcL2[0];
        layout->panelsInt.push_back(le1);

        if ((gfd->VentHoles) && (gfd->VentCentralNerv)) {
            LayoutElement* le2 = new LayoutElement();
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
            le2->p2a0 = layout->pinceLAAmp2[0];
            le2->p2f0 = layout->pinceLFAmp2[0];

            le2->p2a01 = 0.0f;
            le2->p2a1 = layout->pinceLAAmp2[0];
            le2->p2f1 = layout->pinceLFAmp2[0];

            le2->func2f0 = layout->funcL2[0];
            le2->func2f1 = layout->funcL2[0];
            layout->panelsInt.push_back(le2);

            if (gfd->LayoutKlapans) {
                if (gfd->VentHolesDouble){
                    LayoutElement* le3 = new LayoutElement();
                    le3->n1 = -1;
                    le3->n2 = 0;
                    le3->s1 = false;
                    le3->s2 = false;
                    le3->posKlapanIntDeb = 0.0f;
                    le3->posKlapanFin = gfd->PosKlapanFin;
                    le3->isKlapan = 1;
                    le3->isPince = 0;
                    //printf ("\n klapan=%d", isave);
                    layout->panelsInt.push_back(le3);
                }
                    LayoutElement* le4 = new LayoutElement();
                    le4->n1 = -1;
                    le4->n2 = 0;
                    le4->s1 = false;
                    le4->s2 = false;
                    le4->posKlapanIntDeb = gfd->VentHolesDeb;
                    le4->posKlapanFin = gfd->PosKlapanFin;
                    le4->isKlapan = 1;
                    le4->isPince = 0;
                    layout->panelsInt.push_back(le4);


                    LayoutElement* le5 = new LayoutElement();
                    le5->n1 = -1;
                    le5->n2 = 0;
                    le5->s1 = false;
                    le5->s2 = false;
                    le5->posKlapanIntDeb = gfd->VentHolesFin;
                    le5->posKlapanFin = gfd->PosKlapanFin;
                    le5->isKlapan = 1;
                    le5->isPince = 0;
                    layout->panelsInt.push_back(le5);

            }

        }
}

void prepareCenterPanelExt(WindPatternsProject* gfd, Layout* layout) {
    LayoutElement* le = new LayoutElement();
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
    le->p1a00 = layout->pinceLAAmp2[0];
    le->p1a0 = layout->pinceLAAmp1[0];
    le->p1f0 = layout->pinceLFAmp1[0];
    le->p1a01 = layout->pinceLAAmp2[0];
    le->p1a1 = layout->pinceLAAmp1[0];
    le->p1f1 = layout->pinceLFAmp1[0];

    le->func1f0 = layout->funcL1[0];
    le->func1f1 = layout->funcL1[0];
    if ((gfd->VentHoles) && (gfd->VentHolesDouble) && (gfd->VentCentralNerv)) {
        le->posDeb1 = 0.0f;
        le->posDeb2 = 0.0f;
        le->fd1 = 1;
        le->fd2 = 1;
    }
    layout->panelsExt.push_back(le);
}

void prepareProfile(WindPatternsProject* gfd, Layout* layout, int i) {
    LayoutElement* lep = new LayoutElement();
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
    lep->coeff = layout->coeffn[i];
    layout->profs.push_back(lep);
}

void prepareDiagNerv(WindPatternsProject* gfd, Layout* layout, int i) {
    LayoutElement* led1 = new LayoutElement();
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
    led1->coeff = layout->coeffd[k];
    k++;

    LayoutElement* led2 = new LayoutElement();
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
    led2->coeff = layout->coeffd[k];
    layout->diagNervs.push_back(led2);
}


void prepareLayoutElements(WindPatternsProject* gfd, Layout* layout) {

    int startNoNerv = 0, noNerv = 0, face = 1;
    Matrix * Xd[2], *newXd[2], *Yd[2], *newYd[2], *Xdp[2], *Ydp[2], *rXdp[2], *rYdp[2], *rXd[2], *rYd[2];
    Matrix * X[2], *Y[2], *Z[2], *P[2], *rP[2];//, *newP[2];
    int n = gfd->Form->m_nbProfils;
    int size = 6 * n + 2 + 10;

    face = 1;
    int *n1, *n2, *fd1, *ff1, *fd2, *ff2, *isPince,*isKlapan, *vent;
    bool *s1, *s2;
    double *coeff, *p1a0, *p1a00, *p1a01, *p1f0, *p1a1, *p1f1, *p2a0, *p2a00, *p2f0, *p2a1, *p2a01, *p2f1, *posDeb1, *posDeb2, *posFin1, *posFin2, *posKlapanIntDeb, *posKlapanFin;
    Matrix **func1f0, **func1f1, **func2f0, **func2f1;

    fd1 = new int[size];
    fd2 = new int[size];
    ff1 = new int[size];
    ff2 = new int[size];
    n1 = new int[size];
    n2 = new int[size];
    s1 = new bool[size];
    s2 = new bool[size];

    posDeb1 = new double[size];
    posFin1 = new double[size];
    posDeb2 = new double[size];
    posFin2 = new double[size];

    coeff = new double[size];

    p1a0 = new double[size];
    p1a00 = new double[size];
    p1f0 = new double[size];
    p1a1 = new double[size];
    p1a01 = new double[size];
    p1f1 = new double[size];

    p2a0 = new double[size];
    p2a00 = new double[size];
    p2f0 = new double[size];
    p2a1 = new double[size];
    p2a01 = new double[size];
    p2f1 = new double[size];

    posKlapanIntDeb = new double[size];
    posKlapanFin = new double[size];


    func1f0 = new Matrix*[size];
    func1f1 = new Matrix*[size];
    func2f0 = new Matrix*[size];
    func2f1 = new Matrix*[size];

    isPince = new int[size];
    isKlapan = new int[size];
	memset(isPince, 0, sizeof(int)*size);
	memset(isKlapan, 0, sizeof(int)*size);

    vent = new int[size];
    int *numncol = new int[6];

    int isave = 0, i = 0;
    int col = 0;
    float debBorder = 0.0f;
    int faceDebBorder = 0;
    if (gfd->VentHoles) {
        debBorder = gfd->VentHolesDeb;
        faceDebBorder = 2;
    } else {
        debBorder = 0.0f;
        faceDebBorder = 1;
    }

	layout->debBorder = debBorder;
	layout->faceDebBorder = faceDebBorder;
    double tpF1, tpF2;

    for (face = 1; face <= 2; face++) {
        numncol[col] = isave;
        col++;
        /*
        if (gfd->LayoutSymetrique) {
	    //TOTEST not tested yet
            for (i = n - 1; i > 0; i--) {
                if (face == 1) {
					LayoutElement* le = new LayoutElement();
					layoutElementsList.push_back(le);

                    n1[isave] = i;
                    n2[isave] = i - 1;
                    s1[isave] = true;
                    s2[isave] = true;

                    ff1[isave] = face;
                    ff2[isave] = face;
                    fd1[isave] = faceDebBorder;
                    fd2[isave] = faceDebBorder;

                    posDeb1[isave] = debBorder;
                    posFin1[isave] = 100.0f;
                    posDeb2[isave] = debBorder;
                    posFin2[isave] = 100.0f;

                    p1a0[isave] = layout->pinceLAAmp1 [i];
                    p1a00[isave] = layout->pinceLAAmp2 [i];

                    p1f0[isave] = layout->pinceLFAmp1 [i];
                    func1f0[isave] = layout->funcL1[i];
                    func1f1[isave] = layout->funcR1[i - 1];
                    p1a1[isave] = layout->pinceRAAmp1 [i - 1];
                    p1a01[isave] = layout->pinceRAAmp2 [i - 1];

                    p1f1[isave] = layout->pinceRFAmp1 [i - 1];
                    isPince[isave] = 1;
                    coeff[isave] = 0.0f;
                    isave++;

                } else {

                    if ((gfd->VentHoles) && (gfd->noNervVH[i - 1])) {
                        tpF1 = gfd->VentHolesFin;
                        tpF2 = gfd->VentHolesFin;
                    } else {
                        tpF1 = 100.0f;
                        tpF2 = 100.0f;
                    }
                    n1[isave] = i;
                    n2[isave] = i - 1;
                    s1[isave] = true;
                    s2[isave] = true;
                    fd1[isave] = face;
                    fd2[isave] = face;
                    ff1[isave] = face;
                    ff2[isave] = face;
                    posDeb1[isave] = debBorder;
                    posFin1[isave] = tpF1;
                    posDeb2[isave] = debBorder;
                    posFin2[isave] = tpF2;
                    p2a0[isave] = layout->pinceLAAmp2 [i];
                    p2a00[isave] = 0.0f;

                    p2f0[isave] = layout->pinceLFAmp2 [i];
                    func2f0[isave] = layout->funcL2[i];
                    func2f1[isave] = layout->funcR2[i - 1];
                    p2a1[isave] = layout->pinceRAAmp2 [i - 1];
                    p2a01[isave] = 0.0f;

                    p2f1[isave] = layout->pinceRFAmp2 [i - 1];
                    isPince[isave] = 1;
                    coeff[isave] = 0.0f;
                    isave++;

                    if ((gfd->VentHoles) && (gfd->noNervVH[i - 1])) {
                        n1[isave] = i;
                        n2[isave] = i - 1;
                        s1[isave] = true;
                        s2[isave] = true;
                        fd1[isave] = face;
                        fd2[isave] = face;
                        ff1[isave] = face;
                        ff2[isave] = face;
                        posDeb1[isave] = gfd->VentHolesFin;
                        posFin1[isave] = 100.0f;
                        posDeb2[isave] = gfd->VentHolesFin;
                        posFin2[isave] = 100.0f;
                        p2a0[isave] = layout->pinceLAAmp2 [i];
                        p2a00[isave] = 0.0f;


                        p2f0[isave] = layout->pinceLFAmp2 [i];
                        func2f0[isave] = layout->funcL2[i];
                        func2f1[isave] = layout->funcR2[i - 1];
                        p2a1[isave] = layout->pinceRAAmp2 [i - 1];
                        p2a01[isave] = 0.0f;

                        p2f1[isave] = layout->pinceRFAmp2 [i - 1];
                        isPince[isave] = 1;
                        coeff[isave] = 0.0f;
                        isave++;
                    }

                }
            }
        }
        if layout symetric
        */

        if (layout->isCenterPanel) {
            if (face == 1) {
                prepareCenterPanelExt(gfd, layout);
                /*
                n1[isave] = -1;
                n2[isave] = 0;
                s1[isave] = false;
                s2[isave] = false;
                fd1[isave] = faceDebBorder;
                fd2[isave] = faceDebBorder;
                ff1[isave] = 1;
                ff2[isave] = 1;
                isPince[isave] = 1;
				isKlapan[isave] = 0;
                coeff[isave] = 0.0f;
                posDeb1[isave] = debBorder;
                posFin1[isave] = 100.0f;
                posDeb2[isave] = debBorder;
                posFin2[isave] = 100.0f;

                p1a00[isave] = layout->pinceLAAmp2[0];
                p1a0[isave] = layout->pinceLAAmp1[0];
                p1f0[isave] = layout->pinceLFAmp1[0];

                p1a01[isave] = layout->pinceLAAmp2[0];
                p1a1[isave] = layout->pinceLAAmp1[0];
                p1f1[isave] = layout->pinceLFAmp1[0];

                func1f0[isave] = layout->funcL1[0];
                func1f1[isave] = layout->funcL1[0];
                if ((gfd->VentHoles) && (gfd->VentHolesDouble) && (gfd->VentCentralNerv)) {
                    posDeb1[isave] = 0.0f;
                    posDeb2[isave] = 0.0f;
                    fd1[isave] = 1;
                    fd2[isave] = 1;

                }

                isave++;*/
            } else {
                prepareCenterPanelInt(gfd, layout);
                /*if ((gfd->VentHoles) && (gfd->VentCentralNerv)) {
                    tpF1 = gfd->VentHolesFin;
                    tpF2 = gfd->VentHolesFin;
                } else {
                    tpF1 = 100.0f;
                    tpF2 = 100.0f;
                }

                if ((gfd->VentHoles) && (gfd->VentHolesDouble) && (gfd->VentCentralNerv)) {
                    n1[isave] = -1;
                    n2[isave] = 0;
                    s1[isave] = false;
                    s2[isave] = false;
                    posDeb1[isave] = 0.0f;
                    posFin1[isave] = debBorder;
                    posDeb2[isave] = 0.0f;
                    posFin2[isave] = debBorder;

                    p2a00[isave] = 0.0f;
                    p2a0[isave] = layout->pinceRAAmp2 [i];
                    p2f0[isave] = layout->pinceRFAmp2 [i];

                    p2a01[isave] = 0.0f;
                    p2a1[isave] = layout->pinceLAAmp2 [i + 1];
                    p2f1[isave] = layout->pinceLFAmp2 [i + 1];

                    func2f0[isave] = layout->funcR2[i];
                    func2f1[isave] = layout->funcL2[i + 1];

                    fd1[isave] = 2;
                    fd2[isave] = 2;
                    ff1[isave] = 2;
                    ff2[isave] = 2;
                    isPince[isave] = 1;
                    coeff[isave] = 0.0f;
                    isave++;

                }
                
                n1[isave] = -1;
                n2[isave] = 0;
                s1[isave] = false;
                s2[isave] = false;
                fd1[isave] = 2;
                fd2[isave] = 2;
                ff1[isave] = 2;
                ff2[isave] = 2;

                isPince[isave] = 1;
                coeff[isave] = 0.0f;

                posDeb1[isave] = debBorder;
                posFin1[isave] = tpF1;
                posDeb2[isave] = debBorder;
                posFin2[isave] = tpF2;

                p2a00[isave] = 0.0f;
                p2a0[isave] = layout->pinceLAAmp2[0];
                p2f0[isave] = layout->pinceLFAmp2[0];

                p2a01[isave] = 0.0f;
                p2a1[isave] = layout->pinceLAAmp2[0];
                p2f1[isave] = layout->pinceLFAmp2[0];

                func2f0[isave] = layout->funcL2[0];
                func2f1[isave] = layout->funcL2[0];
                isave++;

                if ((gfd->VentHoles) && (gfd->VentCentralNerv)) {
                    n1[isave] = -1;
                    n2[isave] = 0;
                    s1[isave] = false;
                    s2[isave] = false;
                    fd1[isave] = 2;
                    fd2[isave] = 2;
                    ff1[isave] = 2;
                    ff2[isave] = 2;

                    isPince[isave] = 1;
                    coeff[isave] = 0.0f;

                    posDeb1[isave] = gfd->VentHolesFin;
                    posFin1[isave] = 100.0f;
                    posDeb2[isave] = gfd->VentHolesFin;
                    posFin2[isave] = 100.0f;


                    p2a00[isave] = 0.0f;
                    p2a0[isave] = layout->pinceLAAmp2[0];
                    p2f0[isave] = layout->pinceLFAmp2[0];

                    p2a01[isave] = 0.0f;
                    p2a1[isave] = layout->pinceLAAmp2[0];
                    p2f1[isave] = layout->pinceLFAmp2[0];

                    func2f0[isave] = layout->funcL2[0];
                    func2f1[isave] = layout->funcL2[0];
                    isave++;

                    if (gfd->LayoutKlapans) {
                        if (gfd->VentHolesDouble){
                            n1[isave] = -1;
                            n2[isave] = 0;
                            s1[isave] = false;
                            s2[isave] = false;
                            posKlapanIntDeb[isave]=0.0f;
                            posKlapanFin[isave]=gfd->PosKlapanFin;
                            isKlapan[isave] = 1;
                            isPince[isave] = 0;
                            //printf ("\n klapan=%d", isave);
                            isave++;
                        }
                            n1[isave] = -1;
                            n2[isave] = 0;
                            s1[isave] = false;
                            s2[isave] = false;
                            posKlapanIntDeb[isave]=gfd->VentHolesDeb;
                            posKlapanFin[isave]=gfd->PosKlapanFin;
                            isKlapan[isave] = 1;
                            isPince[isave] = 0;
                            //printf ("\n klapan=%d", isave);
                            isave++;

                            n1[isave] = -1;
                            n2[isave] = 0;
                            s1[isave] = false;
                            s2[isave] = false;
                            posKlapanIntDeb[isave]=gfd->VentHolesFin;
                            posKlapanFin[isave]=gfd->PosKlapanFin;
                            isKlapan[isave] = 1;
                            isPince[isave] = 0;
                            //printf ("\n klapan=%d", isave);
                            isave++;

                    }

                }*/
            }
        }

        for (i = 0; i < n - 1; i++) {
            if (face == 1) {
                preparePanelExt(gfd, layout, i);
                /*n1[isave] = i;
                n2[isave] = i + 1;
                s1[isave] = false;
                s2[isave] = false;
                posDeb1[isave] = debBorder;
                posFin1[isave] = 100.0f;
                posDeb2[isave] = debBorder;
                posFin2[isave] = 100.0f;

                p1a00[isave] = layout->pinceRAAmp2 [i];
                p1a0[isave] = layout->pinceRAAmp1 [i];
                p1f0[isave] = layout->pinceRFAmp1 [i];

                p1a01[isave] = layout->pinceLAAmp2 [i + 1];
                p1a1[isave] = layout->pinceLAAmp1 [i + 1];
                p1f1[isave] = layout->pinceLFAmp1 [i + 1];

                func1f0[isave] = layout->funcR1[i];
                func1f1[isave] = layout->funcL1[i + 1];

                ff1[isave] = 1;
                ff2[isave] = 1;
                fd1[isave] = faceDebBorder;
                fd2[isave] = faceDebBorder;

                isPince[isave] = 1;
				isKlapan[isave] = 0;
                coeff[isave] = 0.0f;
                // tututu
                if ((gfd->VentHoles) && (gfd->VentHolesDouble) && (gfd->noNervVH[i])) {
                    posDeb1[isave] = 0.0f;
                    posDeb2[isave] = 0.0f;
                    fd1[isave] = 1;
                    fd2[isave] = 1;

                }

                if (i == (n - 2)) {
                    func1f1[isave] = Zeros(func1f0[isave]->GetLignes(), 1);
                    p1a01[isave] = 0.0f;
                    p1a1[isave] = 0.0f;
                    p1f1[isave] = 0.0f;
                    if (!((gfd->VentHoles) && (gfd->noNervVH[i]))) {
                        posDeb2[isave] = 0.0f;
                        fd2[isave] = 2;
                    }

                }
                isave++;
                */
            } else {
                preparePanelInt(gfd, layout, i);
/*
                if ((gfd->VentHoles) && (gfd->noNervVH[i])) {
                    tpF1 = gfd->VentHolesFin;
                    tpF2 = gfd->VentHolesFin;
                } else {
                    tpF1 = 100.0f;
                    tpF2 = 100.0f;
                }
//tututu
                if ((gfd->VentHoles) && (gfd->VentHolesDouble) && (gfd->noNervVH[i])) {
                    n1[isave] = i;
                    n2[isave] = i + 1;
                    s1[isave] = false;
                    s2[isave] = false;
                    posDeb1[isave] = 0.0f;
                    posFin1[isave] = debBorder;
                    posDeb2[isave] = 0.0f;
                    posFin2[isave] = debBorder;

                    p2a00[isave] = 0.0f;
                    p2a0[isave] = layout->pinceRAAmp2 [i];
                    p2f0[isave] = layout->pinceRFAmp2 [i];

                    p2a01[isave] = 0.0f;
                    p2a1[isave] = layout->pinceLAAmp2 [i + 1];
                    p2f1[isave] = layout->pinceLFAmp2 [i + 1];

                    func2f0[isave] = layout->funcR2[i];
                    func2f1[isave] = layout->funcL2[i + 1];

                    fd1[isave] = 2;
                    fd2[isave] = 2;
                    ff1[isave] = 2;
                    ff2[isave] = 2;
                    isPince[isave] = 1;
                    coeff[isave] = 0.0f;
                    isave++;
                }

                n1[isave] = i;
                n2[isave] = i + 1;
                s1[isave] = false;
                s2[isave] = false;
                posDeb1[isave] = debBorder;
                posFin1[isave] = tpF1;
                posDeb2[isave] = debBorder;
                posFin2[isave] = tpF2;

                p2a00[isave] = 0.0f;
                p2a0[isave] = layout->pinceRAAmp2 [i];
                p2f0[isave] = layout->pinceRFAmp2 [i];

                p2a01[isave] = 0.0f;
                p2a1[isave] = layout->pinceLAAmp2 [i + 1];
                p2f1[isave] = layout->pinceLFAmp2 [i + 1];

                func2f0[isave] = layout->funcR2[i];
                func2f1[isave] = layout->funcL2[i + 1];

                fd1[isave] = 2;
                fd2[isave] = 2;
                ff1[isave] = 2;
                ff2[isave] = 2;
                isPince[isave] = 1;
                coeff[isave] = 0.0f;

                if (i == (n - 2)) {
                    func2f1[isave] = Zeros(func2f0[isave]->GetLignes(), 1);
                    p2a01[isave] = 0.0f;
                    p2a1[isave] = 0.0f;
                    p2f1[isave] = 0.0f;

                    if (!((gfd->VentHoles) && (gfd->noNervVH[i]))  ) {
                        posFin1[isave] = 100.0f;
                        posDeb2[isave] = 0.0f;
                        posFin2[isave] = 100.0f;
                    }
                }

                isave++;

                if ((gfd->VentHoles) && (gfd->noNervVH[i])) {
                    n1[isave] = i;
                    n2[isave] = i + 1;
                    s1[isave] = false;
                    s2[isave] = false;
                    posDeb1[isave] = gfd->VentHolesFin;
                    posFin1[isave] = 100.0f;
                    posDeb2[isave] = gfd->VentHolesFin;
                    posFin2[isave] = 100.0f;

                    p2a00[isave] = 0.0f;
                    p2a0[isave] = layout->pinceRAAmp2 [i];
                    p2f0[isave] = layout->pinceRFAmp2 [i];

                    func2f0[isave] = layout->funcR2[i];
                    func2f1[isave] = layout->funcL2[i + 1];

                    p2a01[isave] = 0.0f;
                    p2a1[isave] = layout->pinceLAAmp2 [i + 1];
                    p2f1[isave] = layout->pinceLFAmp2 [i + 1];

                    fd1[isave] = 2;
                    fd2[isave] = 2;
                    ff1[isave] = 2;
                    ff2[isave] = 2;

                    isPince[isave] = 1;
                    coeff[isave] = 0.0f;

                    if (i == (n - 2)) {
                        func2f1[isave] = Zeros(func2f0[isave]->GetLignes(), 1);
                        p2a01[isave] = 0.0f;
                        p2a1[isave] = 0.0f;
                        p2f1[isave] = 0.0f;
                    }
                    isave++;

                    if (gfd->LayoutKlapans) {
                        if (gfd->VentHolesDouble){
                            n1[isave] = i;
                            n2[isave] = i+1;
                            s1[isave] = false;
                            s2[isave] = false;
                            posKlapanIntDeb[isave]=0.0f;
                            posKlapanFin[isave]=gfd->PosKlapanFin;
                            isPince[isave] = 0;
                            isKlapan[isave] = 1;
                            //printf ("\n klapan=%d", isave);
							fd1[isave] = 2;
							fd2[isave] = 2;
							ff1[isave] = 2;
							ff2[isave] = 2;
                            isave++;
                        }
                            n1[isave] = i;
                            n2[isave] = i+1;
                            s1[isave] = false;
                            s2[isave] = false;
                            posKlapanIntDeb[isave]=gfd->VentHolesDeb;
                            posKlapanFin[isave]=gfd->PosKlapanFin;
                            isPince[isave] = 0;
                            isKlapan[isave] = 1;
                            //printf ("\n klapan=%d", isave);
							fd1[isave] = 2;
							fd2[isave] = 2;
							ff1[isave] = 2;
							ff2[isave] = 2;
                           isave++;

                            n1[isave] = i;
                            n2[isave] = i+1;
                            s1[isave] = false;
                            s2[isave] = false;
                            posKlapanIntDeb[isave]=gfd->VentHolesFin;
                            posKlapanFin[isave]=gfd->PosKlapanFin;
                            isPince[isave] = 0;
                            isKlapan[isave] = 1;
                            //printf ("\n klapan=%d", isave);
							fd1[isave] = 2;
							fd2[isave] = 2;
							ff1[isave] = 2;
							ff2[isave] = 2;
                            isave++;

                    }

                }
*/
            }
        }

        //profiles
        if (face == 1) {
            numncol[col] = isave;
            col++;

            // profile
            /*
            if (gfd->LayoutSymetrique) {
                for (i = n - 2; i > 0; i--) {
                    n1[isave] = i;
                    n2[isave] = i;
                    s1[isave] = true;
                    s2[isave] = true;
                    fd1[isave] = 2;
                    ff1[isave] = 2;

                    fd2[isave] = 1;
                    ff2[isave] = 1;

                    posDeb1[isave] = 0.0f;
                    posFin1[isave] = 100.0f;
                    posDeb2[isave] = 0.0f;
                    posFin2[isave] = 100.0f;

                    isPince[isave] = 0;
                    isave++;
                }
            }
            if layout sym
            */
            for (i = 0; i < n - 1; i++) {
                prepareProfile(gfd, layout, i);

                /*n1[isave] = i;
                n2[isave] = i;
                s1[isave] = false;
                s2[isave] = false;
                posDeb1[isave] = 0.0f;
                posFin1[isave] = 100.0f;
                posDeb2[isave] = 0.0f;
                posFin2[isave] = 100.0f;
                isPince[isave] = 0;
                fd1[isave] = 2;
                ff1[isave] = 2;

                fd2[isave] = 1;
                ff2[isave] = 1;
                coeff[isave] = layout->coeffn[i];
                isave++;*/

            }
            numncol[col] = isave;
            col++;

			if (gfd->DiagNervs) {
				int k = 2 * gfd->quantDiag - 1;
/*
				if (gfd->LayoutSymetrique) {
					for (int i = gfd->quantDiag - 1; i >= 0; i--) {
						n1[isave] = gfd->noNervD[i];
						n2[isave] = gfd->noNervD[i] + 1;
						s1[isave] = true;
						s2[isave] = true;
						isPince[isave] = 0;
						fd1[isave] = 2;
						ff1[isave] = 2;
						fd2[isave] = 1;
						ff2[isave] = 1;
						posDeb1[isave] = gfd->PosDiagNerv2A;
						posFin1[isave] = gfd->PosDiagNerv2F;
						posDeb2[isave] = gfd->PosDiagNerv1A;
						posFin2[isave] = gfd->PosDiagNerv1F;
						coeff[isave] = layout->coeffd[k];
						isave++;
						k--;
						n1[isave] = gfd->noNervD[i];
						n2[isave] = gfd->noNervD[i] - 1;
						s1[isave] = true;
						s2[isave] = true;
						isPince[isave] = 0;
						fd1[isave] = 2;
						ff1[isave] = 2;

						fd2[isave] = 1;
						ff2[isave] = 1;

						posDeb1[isave] = gfd->PosDiagNerv2A;
						posFin1[isave] = gfd->PosDiagNerv2F;
						posDeb2[isave] = gfd->PosDiagNerv1A;
						posFin2[isave] = gfd->PosDiagNerv1F;
						coeff[isave] = layout->coeffd[k];
						isave++;
						k--;
					}
				}
                if layout symetric
                */
            k = 0;
            for (i = 0; i < gfd->quantDiag; i++) {
                prepareDiagNerv(gfd, layout, i);

/*                n1[isave] = gfd->noNervD[i];
                n2[isave] = gfd->noNervD[i] - 1;
                s1[isave] = false;
                s2[isave] = false;
                isPince[isave] = 0;
                fd1[isave] = 2;
                ff1[isave] = 2;
                fd2[isave] = 1;
                ff2[isave] = 1;
                posDeb1[isave] = gfd->PosDiagNerv2A;
                posFin1[isave] = gfd->PosDiagNerv2F;

                posDeb2[isave] = gfd->PosDiagNerv1A;
                posFin2[isave] = gfd->PosDiagNerv1F;
                coeff[isave] = layout->coeffd[k];
                isave++;
                k++;

                n1[isave] = gfd->noNervD[i];
                n2[isave] = gfd->noNervD[i] + 1;
                s1[isave] = false;
                s2[isave] = false;
                isPince[isave] = 0;
                fd1[isave] = 2;
                ff1[isave] = 2;
                fd2[isave] = 1;
                ff2[isave] = 1;

                posDeb1[isave] = gfd->PosDiagNerv2A;
                posFin1[isave] = gfd->PosDiagNerv2F;
                posDeb2[isave] = gfd->PosDiagNerv1A;
                posFin2[isave] = gfd->PosDiagNerv1F;
                coeff[isave] = layout->coeffd[k];
                isave++;
                k++;*/
            }
            }
        }
    }


}

void SaveLayout2(WindPatternsProject* gfd, Layout* layout) {
    if (gfd->debug) printf("\n SaveLayout2()");
    prepareLayoutElements(gfd, layout);
    layout->calculateExport(gfd);
    saveCalcLayoutToFile(gfd, layout);
}

void Layout::calculateExport(WindPatternsProject* gfd){
    for (int i = 0; i < panelsExt.size(); i++) {
        panelsExt[i]->calculateExport(gfd);
    }
    for (int i = 0; i < panelsInt.size(); i++) {
        panelsInt[i]->calculateExport(gfd);
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
void saveCalcLayoutToFile(WindPatternsProject* gfd, Layout* layout) {
    //TAxe **AxeP, **AxePD, **AxePTD, **AxeMD, **AxeCD, **AxeRepD;
    CString fileName;
    LPTSTR PtrfileName;
    char ext[10];
    strcpy(ext, "*.dxf");
    CFileDialog DlgOpen(FALSE, NULL, ext, OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        fileName = DlgOpen.GetPathName();
        PtrfileName = fileName.GetBuffer(1);
        writeLayoutToDXF(PtrfileName, gfd, layout); 
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



void makeEqualSize(Matrix** X0,Matrix** Y0, Matrix** Z0,Matrix** P0,
				   Matrix** X1,Matrix** Y1, Matrix** Z1,Matrix** P1) {
    Matrix *iAv, *iAp, *Res;
    if ((*X0)->GetLignes() > (*X1)->GetLignes()) {
        iAv = LinSpace(0.0f, (double) (*X1)->GetLignes() - 1, (*X1)->GetLignes());
        iAp = LinSpace(0.0f, (double) (*X1)->GetLignes() - 1, (*X0)->GetLignes());
        Res = new Matrix((*X0)->GetLignes(), 1);
        InterpLinMat(iAv, (*X1), iAp, Res);
        delete(*X1);
        *X1 = Res;
        Res = new Matrix((*X0)->GetLignes(), 1);
        InterpLinMat(iAv, *Y1, iAp, Res);
        delete(*Y1);
        *Y1 = Res;
        Res = new Matrix((*X0)->GetLignes(), 1);
        InterpLinMat(iAv, *Z1, iAp, Res);
        delete(*Z1);
        *Z1 = Res;
        Res = new Matrix((*X0)->GetLignes(), 1);
        InterpLinMat(iAv, *P1, iAp, Res);
        delete(*P1);
        *P1 = Res;
        delete(iAv);
        delete(iAp);
    } else if ((*X1)->GetLignes() > (*X0)->GetLignes()) {
        iAv = LinSpace(0.0f, (double) (*X0)->GetLignes() - 1, (*X0)->GetLignes());
        iAp = LinSpace(0.0f, (double) (*X0)->GetLignes() - 1, (*X1)->GetLignes());
        Res = new Matrix((*X1)->GetLignes(), 1);
        InterpLinMat(iAv, *X0, iAp, Res);
        delete(*X0);
        *X0 = Res;
        Res = new Matrix((*X1)->GetLignes(), 1);
        InterpLinMat(iAv, *Y0, iAp, Res);
        delete(*Y0);
        *Y0 = Res;
        Res = new Matrix((*X1)->GetLignes(), 1);
        InterpLinMat(iAv, *Z0, iAp, Res);
        delete(*Z0);
        *Z0 = Res;
        Res = new Matrix((*X1)->GetLignes(), 1);
        InterpLinMat(iAv, *P0, iAp, Res);
        delete(*P0);
        *P0 = Res;
        delete(iAv);
        delete(iAp);
    }
}

void fillXYZP(int* nerv, double* DebT, int* FaceDebT, double* FinT,int* FaceFinT,
			  Matrix* xExtProf, Matrix* XExt,Matrix* YExt,Matrix* ZExt,
			  Matrix* xIntProf, Matrix* XInt,Matrix* YInt,Matrix* ZInt,
                   Matrix** X0,Matrix** Y0, Matrix** Z0,Matrix** P0,
		   Matrix** X1,Matrix** Y1, Matrix** Z1,Matrix** P1) {
    Matrix * X[2], *Y[2], *Z[2], *P[2];
	//printf ("\n DebT[%f, %f]", DebT[0], DebT[1]);
	//printf ("\n FinT[%f, %f]", FinT[0], FinT[1]);
    int iDeb, iFin, iCol;
	int i, j;
    for (i = 0; i < 2; i++) {
        if ((FaceDebT[i] == 1) && (FaceFinT[i] == 1)) //debut et fin en extrados
        {
            if (DebT[i] > FinT[i]) //inversion Deb<->Fin
            {
                Ind(FinT[i], xExtProf, &iDeb, &iCol);
                Ind(DebT[i], xExtProf, &iFin, &iCol);
            } else {
                Ind(DebT[i], xExtProf, &iDeb, &iCol);
                Ind(FinT[i], xExtProf, &iFin, &iCol);
            }
            X[i] = new Matrix(iFin - iDeb + 1, 1);
            Y[i] = new Matrix(iFin - iDeb + 1, 1);
            Z[i] = new Matrix(iFin - iDeb + 1, 1);
            P[i] = new Matrix(iFin - iDeb + 1, 1);
            for (j = 0; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Y[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Z[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < P[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xExtProf->Element(iDeb + j, 0));
        } else if (FaceDebT[i] == 1) //debut en extrados et fin en intrados
        {
            Ind(DebT[i], xExtProf, &iDeb, &iCol);
            Ind(FinT[i], xIntProf, &iFin, &iCol);
            X[i] = new Matrix(iFin + iDeb - 1, 1);
            Y[i] = new Matrix(iFin + iDeb - 1, 1);
            Z[i] = new Matrix(iFin + iDeb - 1, 1);
            P[i] = new Matrix(iFin + iDeb - 1, 1);
            for (j = 0; j <= iDeb; j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) P[i]->SetElement(j, 0, xExtProf->Element(iDeb - j, 0));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xIntProf->Element(j - iDeb, 0));
        } else if (FaceFinT[i] == 1) //debut en intrados et fin en extrados
        {
            Ind(DebT[i], xIntProf, &iDeb, &iCol);
            Ind(FinT[i], xExtProf, &iFin, &iCol);
			//printf ("\n iDeb=%d iFin=%d", iDeb, iFin);
            X[i] = new Matrix(iFin + iDeb + 1, 1);
            Y[i] = new Matrix(iFin + iDeb + 1, 1);
            Z[i] = new Matrix(iFin + iDeb + 1, 1);
            P[i] = new Matrix(iFin + iDeb + 1, 1);
            for (j = 0; j <= iDeb; j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) P[i]->SetElement(j, 0, xIntProf->Element(iDeb - j, 0));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xExtProf->Element(j - iDeb, 0));
        } else //debut et fin en intrados
        {
            if (DebT[i] > FinT[i]) //inversion Deb<->Fin
            {
                Ind(FinT[i], xIntProf, &iDeb, &iCol);
                Ind(DebT[i], xIntProf, &iFin, &iCol);
            } else {
                Ind(DebT[i], xIntProf, &iDeb, &iCol);
                Ind(FinT[i], xIntProf, &iFin, &iCol);
            }
            X[i] = new Matrix(iFin - iDeb + 1, 1);
            Y[i] = new Matrix(iFin - iDeb + 1, 1);
            Z[i] = new Matrix(iFin - iDeb + 1, 1);
            P[i] = new Matrix(iFin - iDeb + 1, 1);
            for (j = 0; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Y[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Z[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < P[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xIntProf->Element(iDeb + j, 0));
        }
    }
    (*X0) = X[0]; (*Y0) = Y[0]; (*Z0) = Z[0]; (*P0) = P[0];
    (*X1) = X[1]; (*Y1) = Y[1]; (*Z1) = Z[1]; (*P1) = P[1];
}

void makePosProfile(Matrix* Ext, Matrix* Int, double percent, Matrix** ExtRes) {
    //printf ("\n makePosProfile()");
    //IntRes = CloneMat(Int);
    (*ExtRes) = CloneMat(Ext);
    double perc = percent*0.01f;
    double yInt, yExt;
    for (int i = 0; i < Ext->GetLignes(); i++) {
        yExt = Ext->Element(i, 1);
        yInt = InterpLinX(Int, Ext->Element(i,0));
        (*ExtRes)->SetElement(i, 1, yInt + perc*(yExt-yInt));
        //printf ("\n%d -> (%f, %f) => (%f,%f)", i, Ext->Element(i,0), Ext->Element(i,1), (*ExtRes)->Element(i,0), (*ExtRes)->Element(i,1));
    }
    //printf ("\n ...makePosProfile()");
}
/*void makeChangePoints(Matrix* XExt,Matrix*  YExt,Matrix*  ZExt,
                      Matrix*  XInt,Matrix*  YInt,Matrix*  ZInt,
        Matrix**  XExt1,Matrix** YExt1, Matrix** ZExt1,
        Matrix** XInt1, Matrix** YInt1, Matrix** ZInt1, int nerv, float pos)
{
    XExt1=CloneMat(XExt); YExt1=CloneMat(YExt); ZExt1=CloneMat(ZExt);
    XInt1=CloneMat(XInt); YInt1=CloneMat(YInt); ZInt1=CloneMat(ZInt);
    int n = XInt->GetColonnes();
    Matrix* interpXInt = new Matrix(n, 2);
    Matrix* interpYInt = new Matrix(n, 2);
    Matrix* interpZInt = new Matrix(n, 2);
    for (int i = 0; i<n;i++) {
        interpXInt->SetElement(i,0,XInt->Element(nerv, i));
        interpXInt->SetElement(i,1,XInt->Element(nerv, i));
    }
    for (int i = 0; i < XExt->GetColonnes(); i++) {
        XExt1 -> SetElement(nerv, i, XInt->Element);
    }
}*/

void calcPatronPosNerv(WindPatternsProject* gfd, int noNerv1, bool sym1, int FaceDeb1, int FaceFin1, double Deb1, double Fin1,
        int noNerv2, bool sym2, int FaceDeb2, int FaceFin2, double Deb2, double Fin2,
        Matrix **Xdev1, Matrix **Ydev1, Matrix **Xdev2, Matrix **Ydev2,
        Matrix **X0, Matrix **Y0, Matrix **Z0, Matrix **P0,
        Matrix **X1, Matrix **Y1, Matrix **Z1, Matrix **P1)
{
    //printf ("\n calcPatronPosNerv()");
    Matrix *XExt, *YExt, *ZExt;
    Matrix *XInt, *YInt, *ZInt;
    Matrix *XExt1, *YExt1, *ZExt1;
    Matrix *XInt1, *YInt1, *ZInt1;
    Matrix *XExt2, *YExt2, *ZExt2;
    Matrix *XInt2, *YInt2, *ZInt2;

    int i, j;
    bool symetrique[2] = {sym1, sym2};
    int nerv[2];
    Matrix *X[2], *Y[2], *Z[2], *P[2];
	Matrix *Xt[2], *Yt[2], *Zt[2], *Pt[2];
    Matrix *xExtProf, *xIntProf;
    int NoNervT[2] = {noNerv1, noNerv2};
    double DebT[2] = {Deb1, Deb2}, FinT[2] = {Fin1, Fin2};
    int FaceDebT[2], FaceFinT[2];
    FaceDebT[0] = 1;
    FaceDebT[1] = 1;
    FaceFinT[0] = 1;
    FaceFinT[1] = 1;

    //Matrix *ExtProfCent1, *IntProfCent1, *ExtProfCent2, *IntProfCent2;
    //Matrix *ExtProfBout1, *IntProfBout1, *ExtProfBout2, *IntProfBout2;
/*
    printf ("\n before makePosProfile()");
    makePosProfile(gfd->ExtProfCent, gfd->IntProfCent, gfd->PosNerv[0], &ExtProfCent1);
    makePosProfile(gfd->ExtProfBout, gfd->IntProfCent, gfd->PosNerv[0], &ExtProfBout1);
    ExtProfCent1->print(0);
    ExtProfCent1->print(1);
    makePosProfile(gfd->ExtProfCent, gfd->IntProfCent, gfd->PosNerv[1], &ExtProfCent2);
    makePosProfile(gfd->ExtProfBout, gfd->IntProfCent, gfd->PosNerv[1], &ExtProfBout2);
    printf ("\n after makePosProfile()");*/
    //printf ("\n before Form0");
    calcForm3D(gfd->Form,0, 0.0f,
            gfd->ExtProfCent, gfd->IntProfCent, gfd->ExtProfBout, gfd->IntProfBout,
            &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);
    //printf ("\n before Form1");
    calcForm3D(gfd->Form, 1, gfd->PosNerv[0],
            gfd->ExtProfCent, gfd->IntProfCent, gfd->ExtProfBout, gfd->IntProfBout,
            &XExt1, &YExt1, &ZExt1, &XInt1, &YInt1, &ZInt1);
    //printf ("\n before Form1");
    calcForm3D(gfd->Form, 1, gfd->PosNerv[1],
            gfd->ExtProfCent, gfd->IntProfCent, gfd->ExtProfBout, gfd->IntProfBout,
            &XExt2, &YExt2, &ZExt2, &XInt2, &YInt2, &ZInt2);
    //printf ("\n after calcForm3D");
    xExtProf = new Matrix(gfd->ExtProfCent->GetLignes(), 1);
    xIntProf = new Matrix(gfd->IntProfCent->GetLignes(), 1);

    for (i = 0; i < gfd->ExtProfCent->GetLignes(); i++)   xExtProf->SetElement(i, 0, gfd->ExtProfCent->Element(i, 0));
    for (i = 0; i < gfd->IntProfCent->GetLignes(); i++)   xIntProf->SetElement(i, 0, gfd->IntProfCent->Element(i, 0));

    for (i = 0; i < 2; i++) {
        nerv[i] = NoNervT[i];
        if (NoNervT[i] == -1) {
            nerv[i] = 0;
            symetrique[i] = true;
        }
    }

    fillXYZP(nerv, DebT, FaceDebT, FinT, FaceFinT,
                    xExtProf, XExt1, YExt1, ZExt1,
                    xIntProf, XInt1, YInt1, ZInt1,
                   &X[0], &Y[0], &Z[0], &P[0],
		   &Xt[1], &Yt[1], &Zt[1], &Pt[1]) ;
    fillXYZP(nerv, DebT, FaceDebT, FinT, FaceFinT,
                    xExtProf, XExt2, YExt2, ZExt2,
                    xIntProf, XInt2, YInt2, ZInt2,
                   &Xt[0], &Yt[0], &Zt[0], &Pt[0],
		   &X[1], &Y[1], &Z[1], &P[1]);
    //printf ("\n after FillXYZP");
    makeEqualSize(&X[0], &Y[0], &Z[0], &P[0],&X[1], &Y[1], &Z[1], &P[1]);
    //printf ("\n after makeEqualSize");
    for (i = 0; i < 2; i++) {
        if (symetrique[i] == true) {
            for (j = 0; j < X[i]->GetLignes(); j++)
                X[i]->MultiplyElement(j, 0, -1.0f);
        }
    }
    calcDeveloppe(
            X[0], Y[0], Z[0], X[1], Y[1], Z[1],
            Xdev1, Ydev1, Xdev2, Ydev2);
    
    *X0 = CloneMat(X[0]);
    *Y0 = CloneMat(Y[0]);
    *Z0 = CloneMat(Z[0]);
    *P0 = CloneMat(P[0]);
    *X1 = CloneMat(X[1]);
    *Y1 = CloneMat(Y[1]);
    *Z1 = CloneMat(Z[1]);
    *P1 = CloneMat(P[1]);
    
    delete(X[0]);
    delete(Y[0]);
    delete(Z[0]);
    delete(P[0]);
    delete(XExt);
    delete(YExt);
    delete(ZExt);
    delete(XInt);
    delete(YInt);
    delete(ZInt);

    delete(XExt1);
    delete(YExt1);
    delete(ZExt1);
    delete(XInt1);
    delete(YInt1);
    delete(ZInt1);

    delete(XExt2);
    delete(YExt2);
    delete(ZExt2);
    delete(XInt2);
    delete(YInt2);
    delete(ZInt2);


    delete(xExtProf);
    delete(xIntProf);
    delete(X[1]);
    delete(Y[1]);
    delete(Z[1]);
    delete(P[1]);

}
/*        calcPatronKlapan(getWindPatternsProject(), nerv[0], symetrique[0], nerv[1], symetrique[1], PosKlapanInt, PosKlapanFin,
                    &Xd[0], &Yd[0], &Xd[1], &Yd[1],
                    &X[0], &Y[0], &Z[0], &P[0],
                    &X[1], &Y[1], &Z[1], &P[1]);*/
void calcPatronKlapan(WindPatternsProject* gfd, int noNerv1, bool sym1, int noNerv2, bool sym2, double posKlapanInt, double posKlapanFin,
        Matrix **Xdev1, Matrix **Ydev1, Matrix **Xdev2, Matrix **Ydev2,
        Matrix **X0, Matrix **Y0, Matrix **Z0, Matrix **P0,
        Matrix **X1, Matrix **Y1, Matrix **Z1, Matrix **P1)
{

    Matrix *XExt, *YExt, *ZExt,*XExt0, *YExt0, *ZExt0;
    Matrix *XInt, *YInt, *ZInt,*XInt0, *YInt0, *ZInt0;
    int i, j;
    bool symetrique[2] = {sym1, sym2};
    int nerv[2], iDebInt, iFinExt, iCol;
    int NoNervT[2] = {noNerv1, noNerv2};

    Matrix * X[2], *Y[2], *Z[2], *P[2];
    Matrix *xExtProf, *xIntProf;

    calcForm3D(gfd->Form,0, 0.0f,
            gfd->ExtProfCent, gfd->IntProfCent, gfd->ExtProfBout, gfd->IntProfBout,
            &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);

    xExtProf = new Matrix(gfd->ExtProfCent->GetLignes(), 1);
    xIntProf = new Matrix(gfd->IntProfCent->GetLignes(), 1);
    for (i = 0; i < gfd->ExtProfCent->GetLignes(); i++)   xExtProf->SetElement(i, 0, gfd->ExtProfCent->Element(i, 0));
    for (i = 0; i < gfd->IntProfCent->GetLignes(); i++)   xIntProf->SetElement(i, 0, gfd->IntProfCent->Element(i, 0));

    Ind(posKlapanInt, xIntProf, &iDebInt, &iCol);

    Matrix* zeroExtProfCent, *zeroExtProfBout;

    zeroExtProfCent = CloneMat(gfd-> ExtProfCent);
    for ( i = 0; i < gfd->ExtProfCent->GetLignes(); i++) zeroExtProfCent->SetElement(i, 1, 0.0f);
    zeroExtProfBout = CloneMat(gfd-> ExtProfBout);
    for ( i = 0; i < gfd->ExtProfBout->GetLignes(); i++) zeroExtProfBout->SetElement(i, 1, 0.0f);

    calcForm3D(gfd->Form,0, 0.0f,
            zeroExtProfCent, gfd->IntProfCent, zeroExtProfBout, gfd->IntProfBout,
            &XExt0, &YExt0, &ZExt0, &XInt0, &YInt0, &ZInt0);
    Ind(posKlapanFin, xExtProf, &iFinExt, &iCol);

    for (i = 0; i < 2; i++) {
        nerv[i] = NoNervT[i];
        if (NoNervT[i] == -1) {
            nerv[i] = 0;
            symetrique[i] = true;
        }
    }

    X[0] = new Matrix(2, 1);
    X[1] = new Matrix(2, 1);
    Y[0] = new Matrix(2, 1);
    Y[1] = new Matrix(2, 1);
    Z[0] = new Matrix(2, 1);
    Z[1] = new Matrix(2, 1);
    P[0] = new Matrix(2, 1);
    P[1] = new Matrix(2, 1);

    for ( i = 0; i < 2; i++) {
        X[i]->SetElement(0,0,XInt->Element(nerv[i], iDebInt));
        Y[i]->SetElement(0,0,YInt->Element(nerv[i], iDebInt));
        Z[i]->SetElement(0,0,ZInt->Element(nerv[i], iDebInt));
        P[i]->SetElement(0,0,posKlapanInt);

        X[i]->SetElement(1,0,XExt0->Element(nerv[i], iFinExt));
        Y[i]->SetElement(1,0,YExt0->Element(nerv[i], iFinExt));
        Z[i]->SetElement(1,0,ZExt0->Element(nerv[i], iFinExt));
        P[i]->SetElement(1,0,posKlapanFin);
    }

    for (i = 0; i < 2; i++) {
        if (symetrique[i] == true) {
            for (j = 0; j < X[i]->GetLignes(); j++)
                X[i]->MultiplyElement(j, 0, -1.0f);
        }
    }

    calcDeveloppe(
            X[0], Y[0], Z[0], X[1], Y[1], Z[1],
            Xdev1, Ydev1, Xdev2, Ydev2);

    *X0 = CloneMat(X[0]);
    *Y0 = CloneMat(Y[0]);
    *Z0 = CloneMat(Z[0]);
    *P0 = CloneMat(P[0]);
    *X1 = CloneMat(X[1]);
    *Y1 = CloneMat(Y[1]);
    *Z1 = CloneMat(Z[1]);
    *P1 = CloneMat(P[1]);
    delete(X[0]);
    delete(Y[0]);
    delete(Z[0]);
    delete(P[0]);
    delete(XExt);
    delete(YExt);
    delete(ZExt);
    delete(XInt);
    delete(YInt);
    delete(ZInt);
    delete(xExtProf);
    delete(xIntProf);
    delete(X[1]);
    delete(Y[1]);
    delete(Z[1]);
    delete(P[1]);
    
}


void calcPatron(WindPatternsProject* gfd, int noNerv1, bool sym1, int FaceDeb1, int FaceFin1, double Deb1, double Fin1,
        int noNerv2, bool sym2, int FaceDeb2, int FaceFin2, double Deb2, double Fin2,
        Matrix **Xdev1, Matrix **Ydev1, Matrix **Xdev2, Matrix **Ydev2,
        Matrix **X0, Matrix **Y0, Matrix **Z0, Matrix **P0,
        Matrix **X1, Matrix **Y1, Matrix **Z1, Matrix **P1)
{

    //printf ("\n calcPatron()");
	//printf ("\n FaceDeb1=%d FaceFin1=%d FaceDeb2=%d FaceFin2=%d", FaceDeb1, FaceFin1, FaceDeb2, FaceFin2);
	//printf ("\n[%f, %f] [%f, %f]", Deb1, Fin1, Deb2, Fin2);

    Matrix *XExt, *YExt, *ZExt;
    Matrix *XInt, *YInt, *ZInt;

    int i, j;
    bool symetrique[2] = {sym1, sym2};
    int nerv[2];

    Matrix * X[2], *Y[2], *Z[2], *P[2];
    //int iDeb, iFin, iCol;

    Matrix *xExtProf, *xIntProf;
    //Matrix *iAv, *iAp, *Res;

    int NoNervT[2] = {noNerv1, noNerv2};
    double DebT[2] = {Deb1, Deb2}, FinT[2] = {Fin1, Fin2};
    int FaceDebT[2], FaceFinT[2];
    FaceDebT[0] = FaceDeb1;
    FaceDebT[1] = FaceDeb2;
    FaceFinT[0] = FaceFin1;
    FaceFinT[1] = FaceFin2;

    calcForm3D(gfd->Form,0, 0.0f,
            gfd->ExtProfCent, gfd->IntProfCent, gfd->ExtProfBout, gfd->IntProfBout,
            &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);
    /*adde points de suspentage*/
    //	if(ViewPtsSuspentes){
    //		addPtsSuspentage(Axe3d, F, IntProfCent, XExt, YExt, ZExt, XInt, YInt, ZInt,
    //			ViewSymetrique, &PtsSuspentes);}

    //------ BEGIN fill X, Y, Z, P [i] ----------
    xExtProf = new Matrix(gfd->ExtProfCent->GetLignes(), 1);
    xIntProf = new Matrix(gfd->IntProfCent->GetLignes(), 1);
    for (i = 0; i < gfd->ExtProfCent->GetLignes(); i++)   xExtProf->SetElement(i, 0, gfd->ExtProfCent->Element(i, 0));
    for (i = 0; i < gfd->IntProfCent->GetLignes(); i++)   xIntProf->SetElement(i, 0, gfd->IntProfCent->Element(i, 0));


    for (i = 0; i < 2; i++) {
        nerv[i] = NoNervT[i];
        if (NoNervT[i] == -1) {
            nerv[i] = 0;
            symetrique[i] = true;
        }
    }

	//printf ("\n FaceDebT=%d FaceFinT=%d", FaceDebT, FaceFinT);
    fillXYZP(nerv, DebT, FaceDebT, FinT, FaceFinT,
				xExtProf, XExt, YExt, ZExt, xIntProf,
				XInt, YInt, ZInt, 
	           &X[0], &Y[0], &Z[0], &P[0],
			   &X[1], &Y[1], &Z[1], &P[1]) ;


/*
    for (i = 0; i < 2; i++) {
        nerv[i] = NoNervT[i];
        if (NoNervT[i] == -1) {
            nerv[i] = 0;
            symetrique[i] = true;
        }

        if (FaceDebT[i] == 1 && FaceFinT[i] == 1) //debut et fin en extrados
        {
            if (DebT[i] > FinT[i]) //inversion Deb<->Fin
            {
                Ind(FinT[i], xExtProf, &iDeb, &iCol);
                Ind(DebT[i], xExtProf, &iFin, &iCol);
            } else {
                Ind(DebT[i], xExtProf, &iDeb, &iCol);
                Ind(FinT[i], xExtProf, &iFin, &iCol);
            }
            X[i] = new Matrix(iFin - iDeb + 1, 1);
            Y[i] = new Matrix(iFin - iDeb + 1, 1);
            Z[i] = new Matrix(iFin - iDeb + 1, 1);
            P[i] = new Matrix(iFin - iDeb + 1, 1);
            for (j = 0; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Y[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Z[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < P[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xExtProf->Element(iDeb + j, 0));
        } else if (FaceDebT[i] == 1) //debut en extrados et fin en intrados
        {
            Ind(DebT[i], xExtProf, &iDeb, &iCol);
            Ind(FinT[i], xIntProf, &iFin, &iCol);
            X[i] = new Matrix(iFin + iDeb - 1, 1);
            Y[i] = new Matrix(iFin + iDeb - 1, 1);
            Z[i] = new Matrix(iFin + iDeb - 1, 1);
            P[i] = new Matrix(iFin + iDeb - 1, 1);
            for (j = 0; j <= iDeb; j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) P[i]->SetElement(j, 0, xExtProf->Element(iDeb - j, 0));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xIntProf->Element(j - iDeb, 0));
        } else if (FaceFinT[i] == 1) //debut en intrados et fin en extrados
        {
            Ind(DebT[i], xIntProf, &iDeb, &iCol);
            Ind(FinT[i], xExtProf, &iFin, &iCol);
            X[i] = new Matrix(iFin + iDeb + 1, 1);
            Y[i] = new Matrix(iFin + iDeb + 1, 1);
            Z[i] = new Matrix(iFin + iDeb + 1, 1);
            P[i] = new Matrix(iFin + iDeb + 1, 1);
            for (j = 0; j <= iDeb; j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) P[i]->SetElement(j, 0, xIntProf->Element(iDeb - j, 0));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xExtProf->Element(j - iDeb, 0));
        } else //debut et fin en intrados
        {
            if (DebT[i] > FinT[i]) //inversion Deb<->Fin
            {
                Ind(FinT[i], xIntProf, &iDeb, &iCol);
                Ind(DebT[i], xIntProf, &iFin, &iCol);
            } else {
                Ind(DebT[i], xIntProf, &iDeb, &iCol);
                Ind(FinT[i], xIntProf, &iFin, &iCol);
            }
            X[i] = new Matrix(iFin - iDeb + 1, 1);
            Y[i] = new Matrix(iFin - iDeb + 1, 1);
            Z[i] = new Matrix(iFin - iDeb + 1, 1);
            P[i] = new Matrix(iFin - iDeb + 1, 1);
            for (j = 0; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Y[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Z[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < P[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xIntProf->Element(iDeb + j, 0));
        }
    }
 */

    //---- END fill X, Y, Z, P [i] ----------
    //	for (int _i=0;_i < P[0]->GetLignes();_i++) {
    //		printf ("\nP[0]->Element(%d, 0)=%f", _i, P[0]->Element(_i, 0));
    //	}
    //---- BEGIN interpolation equal quantity points X[0] <-> X[1]; Y[0] <-> Y[1]
    makeEqualSize(&X[0], &Y[0], &Z[0], &P[0],&X[1], &Y[1], &Z[1], &P[1]);

    /*if (X[0]->GetLignes() > X[1]->GetLignes()) {
        iAv = LinSpace(0.0f, (double) X[1]->GetLignes() - 1, X[1]->GetLignes());
        iAp = LinSpace(0.0f, (double) X[1]->GetLignes() - 1, X[0]->GetLignes());
        Res = new Matrix(X[0]->GetLignes(), 1);
        InterpLinMat(iAv, X[1], iAp, Res);
        delete(X[1]);
        X[1] = Res;
        Res = new Matrix(X[0]->GetLignes(), 1);
        InterpLinMat(iAv, Y[1], iAp, Res);
        delete(Y[1]);
        Y[1] = Res;
        Res = new Matrix(X[0]->GetLignes(), 1);
        InterpLinMat(iAv, Z[1], iAp, Res);
        delete(Z[1]);
        Z[1] = Res;
        Res = new Matrix(X[0]->GetLignes(), 1);
        InterpLinMat(iAv, P[1], iAp, Res);
        delete(P[1]);
        P[1] = Res;
        delete(iAv);
        delete(iAp);
    } else if (X[1]->GetLignes() > X[0]->GetLignes()) {
        iAv = LinSpace(0.0f, (double) X[0]->GetLignes() - 1, X[0]->GetLignes());
        iAp = LinSpace(0.0f, (double) X[0]->GetLignes() - 1, X[1]->GetLignes());
        Res = new Matrix(X[1]->GetLignes(), 1);
        InterpLinMat(iAv, X[0], iAp, Res);
        delete(X[0]);
        X[0] = Res;
        Res = new Matrix(X[1]->GetLignes(), 1);
        InterpLinMat(iAv, Y[0], iAp, Res);
        delete(Y[0]);
        Y[0] = Res;
        Res = new Matrix(X[1]->GetLignes(), 1);
        InterpLinMat(iAv, Z[0], iAp, Res);
        delete(Z[0]);
        Z[0] = Res;
        Res = new Matrix(X[1]->GetLignes(), 1);
        InterpLinMat(iAv, P[0], iAp, Res);
        delete(P[0]);
        P[0] = Res;
        delete(iAv);
        delete(iAp);
    }*/

    //---- END OF interpolation equal quantity points X[0] <-> X[1]; Y[0] <-> Y[1]
    for (i = 0; i < 2; i++) {
        if (symetrique[i] == true) {
            for (j = 0; j < X[i]->GetLignes(); j++)
                X[i]->MultiplyElement(j, 0, -1.0f);
        }
    }

    //--- BEGIN calcate Developpe in Xd[0], Yd[0]   Xd[1], Yd[1]
    //        printf ("\n calcp |X0|=%d |Y0|=%d |Z0|=%d |P0|=%d ", X[0]->GetLignes(), Y[0]->GetLignes(), Z[0]->GetLignes(), P[0]->GetLignes());
    //        printf ("\n calcp |X1|=%d |Y1|=%d |Z1|=%d |P1|=%d ", X[1]->GetLignes(), Y[1]->GetLignes(), Z[1]->GetLignes(), P[1]->GetLignes());

    calcDeveloppe(
            X[0], Y[0], Z[0], X[1], Y[1], Z[1],
            Xdev1, Ydev1, Xdev2, Ydev2);
    //--- END calcate Developpe in Xd[0], Yd[0]   Xd[1], Yd[1]
    *X0 = CloneMat(X[0]);
    *Y0 = CloneMat(Y[0]);
    *Z0 = CloneMat(Z[0]);
    *P0 = CloneMat(P[0]);
    *X1 = CloneMat(X[1]);
    *Y1 = CloneMat(Y[1]);
    *Z1 = CloneMat(Z[1]);
    *P1 = CloneMat(P[1]);
    delete(X[0]);
    delete(Y[0]);
    delete(Z[0]);
    delete(P[0]);
    delete(XExt);
    delete(YExt);
    delete(ZExt);
    delete(XInt);
    delete(YInt);
    delete(ZInt);
    delete(xExtProf);
    delete(xIntProf);
    delete(X[1]);
    delete(Y[1]);
    delete(Z[1]);
    delete(P[1]);
    //printf ("\n ...calcPatron()");
}

double calcWidthNervs(WindPatternsProject* gfd, int noNerv1, int noNerv2, int face) {
    Matrix * Xd1[2], *Yd1[2]; //, *Xd1p[2], *Yd1p[2];
    Matrix * X1[2], *Y1[2], *Z1[2], *P1[2];
    calcPatron(gfd, noNerv1, false, face, face, 0.0f, 100.0f,
            noNerv2, false, face, face, 0.0f, 100.0f,
            &Xd1[0], &Yd1[0], &Xd1[1], &Yd1[1],
            &X1[0], &Y1[0], &Z1[0], &P1[0],
            &X1[1], &Y1[1], &Z1[1], &P1[1]);
    double width;
    calcWidth(Xd1[0], Yd1[0], Xd1[1], Yd1[1], &width);
    delete (Xd1[0]);
    delete (Yd1[0]);
    delete (Xd1[1]);
    delete (Yd1[1]);
    delete(X1[0]);
    delete(Y1[0]);
    delete(Z1[0]);
    delete(P1[0]);
    delete(X1[1]);
    delete(Y1[1]);
    delete(Z1[1]);
    delete(P1[1]);
    return width;
}

void GetProfileXY(WindPatternsProject* gfd, Form* F, int nerv, int face, Matrix** XProf, Matrix** YProf) 
{
	//printf ("\n getProfileXY()");
    Matrix *XExt, *YExt, *ZExt, *XInt, *YInt, *ZInt;
	//printf ("\n 1");
    calcForm3D(F, 0, 0.0f,
		gfd->ExtProfCent, gfd->IntProfCent, gfd->ExtProfBout, gfd->IntProfBout, &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);
	//printf ("\n 1.1");
    double LongNerv = 100.0f;
	//LongNerv = F->m_pProfils[nerv]->m_fLength;
	LongNerv = 100.0f;
	double coeffx = LongNerv/100.0f;
	//coeffx = 1.0;
	//printf ("\n 2");
    double EpaiRel = 0.0f;
    double m = -1.0f;
    if (nerv != -1) {
        EpaiRel = F->m_pProfils[nerv]->m_fWidth;
        m = F->m_pProfils[nerv]->m_fMorph;
    } else {
        EpaiRel = F->m_pProfils[0]->m_fWidth;
        m = 1.0f;
    }
	//printf ("\n 3");
    double EpaiRelProfCent = EpaisseurRelative(gfd->ExtProfCent, gfd->IntProfCent);
    double EpaiRelProfBout = EpaisseurRelative(gfd->ExtProfBout, gfd->IntProfBout);
	//printf ("\n inGetProfileXY: EpaiRelProfCent=%f", EpaiRelProfCent);
	//printf ("\n inGetProfileXY: EpaiRelProfBout=%f", EpaiRelProfBout);
    double coeffyCent = LongNerv * EpaiRel / (EpaiRelProfCent * 100.0f);
    double coeffyBout = LongNerv * EpaiRel / (EpaiRelProfBout * 100.0f);
	//printf ("\n coeffx=%f, coeffyCent=%f, coeffyBout=%f", coeffx, coeffyCent, coeffyBout);
    double xp, yp;
    int n = 0, j = 0;
    if (face == 1) n = gfd->ExtProfCent->GetLignes(); else n = gfd->IntProfCent->GetLignes();
    *XProf = new Matrix(n, 1);
    *YProf = new Matrix(n, 1);
	//printf ("\n 4");
    if (face == 1) {
        for (j = 0; j < gfd->ExtProfCent->GetLignes(); j++) {
            xp = gfd->ExtProfCent->Element(j, 0) * coeffx;
            yp = gfd->ExtProfCent->Element(j, 1) * coeffyCent * m + gfd->ExtProfBout->Element(j, 1) * coeffyBout * (1.0f - m);
            (*XProf)->SetElement(j, 0, xp);
            (*YProf)->SetElement(j, 0, yp);
        }
    } else {
        for (j = 0; j < gfd->IntProfCent->GetLignes(); j++) {
            xp = gfd->IntProfCent->Element(j, 0) * coeffx;
            yp = gfd->IntProfCent->Element(j, 1) * coeffyCent * m + gfd->IntProfBout->Element(j, 1) * coeffyBout * (1.0f - m);
            (*XProf)->SetElement(j, 0, xp);
            (*YProf)->SetElement(j, 0, yp);
        }
    }
	//printf ("\n ..getProfileXY()");
}

ProfilGeom* getProfile(WindPatternsProject* gfd, Form* F, int nerv) {
	//printf ("\n getProfileGeom()");
	ProfilGeom* pg = new ProfilGeom();

	Matrix* XExt, *YExt, *XInt, *YInt;
	GetProfileXY(gfd, F, nerv, 1, &XExt, &YExt);
	int n = XExt->GetLignes();
	pg->ExtProf = new Matrix (n,2);
	for (int i = 0; i < n; i++ ) {
		pg->ExtProf->SetElement(i, 0, XExt->Element(i, 0));
		pg->ExtProf->SetElement(i, 1, YExt->Element(i, 0));
	}

	GetProfileXY(gfd, F, nerv, 2, &XInt, &YInt);
	n = XInt->GetLignes();
	pg->IntProf = new Matrix (n,2);
	for (int i = 0; i < n; i++ ) {
		pg->IntProf->SetElement(i, 0, XInt->Element(i, 0));
		pg->IntProf->SetElement(i, 1, YInt->Element(i, 0));
	}
	//printf ("\n ...getProfileGeom()");
	return pg;
}

void GetMiddleProfile(WindPatternsProject* gfd, Form* F, int nerv1, int nerv2, int face, int realMashtab, Matrix** XProf, Matrix** YProf) {
    /*Matrix *XExt, *YExt, *ZExt, *XInt, *YInt, *ZInt;
    calcForm3D(F, 0, 0.0f,
		gfd->ExtProfCent, gfd->IntProfCent, gfd->ExtProfBout, gfd->IntProfBout, &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);*/
    int i1 = nerv1;
    int i2 = nerv2;
	if (i1 == -1) i1 = 0;
	if (i2 == -1) i2 = 0;

	double LongNerv = 100.0f;
	if (realMashtab == 1) {
		LongNerv = 0.5f*(F->m_pProfils[i1]->m_fLength + F->m_pProfils[i2]->m_fLength);
	} 
    double EpaiRel = 0.0f;
    double m = -1.0f;
    if (i1 != -1) {
        EpaiRel = 0.5f * (F->m_pProfils[i1]->m_fWidth + F->m_pProfils[i2]->m_fWidth);
        m = 0.5f * (F->m_pProfils[i1]->m_fMorph + F->m_pProfils[i2]->m_fMorph);
    } else {
        EpaiRel = F->m_pProfils[0]->m_fWidth;
        m = 1.0f;
    }
    double EpaiRelProfCent = EpaisseurRelative(gfd->ExtProfCent, gfd->IntProfCent);
    double EpaiRelProfBout = EpaisseurRelative(gfd->ExtProfBout, gfd->IntProfBout);
    double coeffyCent = LongNerv * EpaiRel / (EpaiRelProfCent * 100.0f);
    double coeffyBout = LongNerv * EpaiRel / (EpaiRelProfBout * 100.0f);
	double coeffx = LongNerv/100.0f;
	//printf ("\n in getMiddleProfile() LongNerv=%f", LongNerv);
    double xp, yp;
    int n = 0, j = 0;
    if (face == 1) n = gfd->ExtProfCent->GetLignes(); else n = gfd->IntProfCent->GetLignes();
    *XProf = new Matrix(n, 1);
    *YProf = new Matrix(n, 1);
    if (face == 1) {
        for (j = 0; j < gfd->ExtProfCent->GetLignes(); j++) {
            xp = gfd->ExtProfCent->Element(j, 0) * coeffx;
            yp = gfd->ExtProfCent->Element(j, 1) * coeffyCent * m + gfd->ExtProfBout->Element(j, 1) * coeffyBout * (1.0f - m);
            (*XProf)->SetElement(j, 0, xp);
            (*YProf)->SetElement(j, 0, yp);
        }
    } else {
        for (j = 0; j < gfd->IntProfCent->GetLignes(); j++) {
            xp = gfd->IntProfCent->Element(j, 0) * coeffx;
            yp = gfd->IntProfCent->Element(j, 1) * coeffyCent * m + gfd->IntProfBout->Element(j, 1) * coeffyBout * (1.0f - m);
            (*XProf)->SetElement(j, 0, xp);
            (*YProf)->SetElement(j, 0, yp);
        }
    }

}

void GetMiddleProfileBal(WindPatternsProject* gfd, Form* F, int nerv1, int nerv2, int face,  int realMashtab, Matrix** XProf, Matrix** YProf) {
	double EpaiRel, xp, yp;
	double EpaiRelProfCent, EpaiRelProfBout;
	double coeffx, coeffy, coeffyCent, coeffyBout;
	int i1, i2, j;
	i1 = nerv1;
	if (i1 == -1) i1=0;
	i2 = nerv2;
	if (i2 == -1) i2=0;
    //bool isCenterPanel = (1 & forme->NbCaiss);
	double LongNerv = 100.0f;
	if (realMashtab == 1) {
		LongNerv = (F->m_pProfils[i1]->m_fLength + F->m_pProfils[i2]->m_fLength) * 0.5f;
	} 
	EpaiRel = (F->m_pProfils[i1]->m_fWidth + F->m_pProfils[i2]->m_fWidth) * 0.5f;
    if (i1 != -1) {
        EpaiRel = 0.5f * (F->m_pProfils[i1]->m_fWidth + F->m_pProfils[i2]->m_fWidth);
    } else {
        EpaiRel = F->m_pProfils[0]->m_fWidth;
    }
	ProfilGeom* pgCur = getProfile(gfd, F, i1);
	ProfilGeom* pgCurBal = getBalloneProfilGeom(pgCur, gfd->ballonement->kChord->Element(i1, 0), gfd->ballonement->kMf->Element(i1, 0), EpaiRel, gfd->ballonement->wN->Element(i1, 0), gfd->ballonement->dyw->Element(i1, 0));
	double l = abs (pgCur->ExtProf->Element(pgCur->ExtProf->GetLignes() - 1, 0) - pgCur->ExtProf->Element(0, 0));
	double xv = l * (100.0f - (gfd->PosPinceBF[0])) * 0.01f;
	ProfilGeom* pg = getProfilGeomTailDown(pgCurBal, pgCur, xv, gfd->ballonement->powerTail->Element(i1, 0));
	coeffx = LongNerv/100.0f;
	double EpaiRelCur = EpaisseurRelative(pgCur->ExtProf, pgCur->IntProf);
	coeffy = LongNerv*EpaiRel/(EpaiRelCur*100.0f);

    int n = 0;
    if (face == 1) n = pg->ExtProf->GetLignes(); else n = pg->IntProf->GetLignes();
    *XProf = new Matrix(n, 1);
    *YProf = new Matrix(n, 1);
	if (face == 1) {
		for (j=0; j < pg->ExtProf->GetLignes(); j++)
		{
			xp = pg->ExtProf->Element(j, 0) * coeffx;
			yp = pg->ExtProf->Element(j, 1) * coeffy;
			(*XProf)->SetElement(j, 0, xp);
			(*YProf)->SetElement(j, 0, yp);
		}
	} else {
		for (j=0; j < pg->IntProf->GetLignes(); j++)
		{
			xp = pg->IntProf->Element(j, 0) * coeffx;
			yp = pg->IntProf->Element(j, 1) * coeffy;
			(*XProf)->SetElement(j, 0, xp);
			(*YProf)->SetElement(j, 0, yp);
		}
	}
}

void goCalcIndepPinceNew(WindPatternsProject* gfd, int noNerv, int face, double *pLA, double *pLF, Matrix** fl, double *pRA, double *pRF, Matrix** fr, double *len) {
    Matrix * Xd1[2], *Yd1[2], *Xd1p[2], *Yd1p[2];
    Matrix * Xd2[2], *Yd2[2], *Xd2p[2], *Yd2p[2];
    Matrix * X1[2], *Y1[2], *Z1[2], *P1[2], *newP1[2];
    Matrix * X2[2], *Y2[2], *Z2[2], *P2[2], *newP2[2];
    bool showPoints = false;
    int i = 0;
    calcPatron(gfd, noNerv - 1, false, face, face, 0.0f, 100.0f,
            noNerv, false, face, face, 0.0f, 100.0f,
            &Xd1[0], &Yd1[0], &Xd1[1], &Yd1[1],
            &X1[0], &Y1[0], &Z1[0], &P1[0],
            &X1[1], &Y1[1], &Z1[1], &P1[1]);
    Matrix *lenm1 = Longueur(Xd1[1], Yd1[1]);
    double len1 = lenm1->Element(lenm1->GetLignes() - 1, 0);
    double w1, w2;
    calcWidth(Xd1[0], Yd1[0], Xd1[1], Yd1[1], &w1);
    calcPatron(gfd, noNerv, false, face, face, 0.0f, 100.0f,
            noNerv + 1, false, face, face, 0.0f, 100.0f,
            &Xd2[0], &Yd2[0], &Xd2[1], &Yd2[1],
            &X2[0], &Y2[0], &Z2[0], &P2[0],
            &X2[1], &Y2[1], &Z2[1], &P2[1]);
    Matrix *lenm2 = Longueur(Xd2[0], Yd2[0]);
    double len2 = lenm2->Element(lenm2->GetLignes() - 1, 0);
    calcWidth(Xd2[0], Yd2[0], Xd2[1], Yd2[1], &w2);

    double ppA = gfd->PosPinceBA[0];
    double apAPercInit = gfd->AmpPinceBA[0];
    double ppF = gfd->PosPinceBF[0];
    double apFPercInit = gfd->AmpPinceBF[0];

    Pince* pince1 = getPince(gfd, noNerv - 1, noNerv, face);
    (*fl) = pince1->function2;
    Pince* pince2 = getPince(gfd, noNerv, noNerv + 1, face);
    (*fr) = pince2->function1;

    double apASumAbs = pince1 ->AmpA + pince2->AmpA;
    double apFSumAbs = pince1 ->AmpF + pince2->AmpF;

    double middleInit = w1 / (w1 + w2);
    double apAAbs1, apFAbs1, apAAbs2, apFAbs2;
    double middle = 0.0f;
    Matrix *long1p1, *long2p0;
    double minDelta = 1000000.0f, minMiddle = -1.0f, minDeltaMiddleInit = 1000000.0f, minlp1 = 1000000.0f, minlp2 = 1000000.0f;
    double l1p1 = -1.0f, l2p0 = -1.0f, delta = 0.0f, deltaMiddleInit = 0.0f, mina1, mina2;

    for (middle = 0.0f; middle <= 1.0f; middle += HACCUR) {
        apAAbs1 = apASumAbs*middle;
        apAAbs2 = apASumAbs - apAAbs1;

        apFAbs1 = apFSumAbs*middle;
        apFAbs2 = apFSumAbs - apFAbs1;

        calcPinceAloneAbsFunc(Xd1[1], Yd1[1], Xd1[0], Yd1[0], X1[1], Y1[1], Z1[1], P1[1], X1[0], Y1[0], Z1[0], P1[0],
                pince1, apAAbs1, apFAbs1, 1,
                &Xd1p[1], &Yd1p[1], &newP1[1]);

        long1p1 = Longueur(Xd1p[1], Yd1p[1]);
        l1p1 = long1p1->Element(long1p1->GetLignes() - 1, 0);

        calcPinceAloneAbsFunc(Xd2[0], Yd2[0], Xd2[1], Yd2[1], X2[0], Y2[0], Z2[0], P2[0], X2[1], Y2[1], Z2[1], P2[1],
                pince2, apAAbs2, apFAbs2, 0,
                &Xd2p[0], &Yd2p[0], &newP2[0]);
        long2p0 = Longueur(Xd2p[0], Yd2p[0]);

        l2p0 = long2p0->Element(long2p0->GetLignes() - 1, 0);
        delta = fabs(l1p1 - l2p0);
        deltaMiddleInit = fabs(middle - middleInit);

        if (minDelta < (gfd->tochnostLayout2 * 0.001f)) {
            if (deltaMiddleInit < minDeltaMiddleInit) {
                minDelta = delta;
                minlp1 = l1p1;
                minMiddle = middle;
                minDeltaMiddleInit = deltaMiddleInit;
                mina1 = apAAbs1;
                mina2 = apAAbs2;
                minlp2 = l2p0;
            }
        } else {
            if (delta < minDelta) {
                minDelta = delta;
                minlp1 = l1p1;
                minMiddle = middle;
                minDeltaMiddleInit = deltaMiddleInit;
                mina1 = apAAbs1;
                mina2 = apAAbs2;
                minlp2 = l2p0;
            }
        }
        delete (long1p1);
        delete (Xd1p[1]);
        delete (Yd1p[1]);
        delete (long2p0);
        delete (Xd2p[0]);
        delete (Yd2p[0]);
        delete (newP1[1]);
        delete (newP2[0]);
    }
    getLayoutLogger()->logprintf("\n%2df%d D=%6.4f(M=%5.3f)[%5.3f] %8.5f/%8.5f l%8.5f a%6.4f a%6.4f",
            noNerv, face, (minDelta * 1000), minMiddle, minMiddle - middleInit, minlp1, minlp2, len1, mina1, mina2);
    if (minDelta < (gfd->tochnostLayout2 * 0.001f)) getLayoutLogger()->logprintf(" OK:)"); else getLayoutLogger()->logprintf(" BAD!");
    *len = (minlp1 + minlp2) * 0.5f;

    *pLA = apASumAbs*minMiddle;
    *pRA = apASumAbs * (1.0f - minMiddle);
    *pLF = apFSumAbs*minMiddle;
    *pRF = apFSumAbs * (1.0f - minMiddle);
    delete (X1[0]);
    delete (Y1[0]);
    delete (Z1[0]);
    delete (P1[0]);
    delete (X1[1]);
    delete (Y1[1]);
    delete (Z1[1]);
    delete (P1[1]);
    delete (X2[0]);
    delete (Y2[0]);
    delete (Z2[0]);
    delete (P2[0]);
    delete (X2[1]);
    delete (Y2[1]);
    delete (Z2[1]);
    delete (P2[1]);
    delete (Xd1[0]);
    delete (Yd1[0]);
    delete (Xd2[0]);
    delete (Yd2[0]);
    delete (Xd1[1]);
    delete (Yd1[1]);
    delete (Xd2[1]);
    delete (Yd2[1]);
    delete (lenm1);
    delete (lenm2);
}


void goCalcNervureWithPince(WindPatternsProject* gfd, int noNerv1, int face1, int noNerv2, int face2, double lp1, double lp2, double *rcoeff) {
    Matrix * Xd[2], *Yd[2], *newXd[2], *newYd[2];
    Matrix * X[2], *Y[2], *Z[2], *P[2];
    calcPatron(gfd, noNerv1, false, face1, face1, 0.0f, 100.0f,
            noNerv2, false, face2, face2, 0.0f, 100.0f,
            &Xd[0], &Yd[0], &Xd[1], &Yd[1],
            &X[0], &Y[0], &Z[0], &P[0],
            &X[1], &Y[1], &Z[1], &P[1]);
    Matrix *long1 = Longueur(Xd[0], Yd[0]);
    double len1 = long1->Element(long1->GetLignes() - 1, 0);
    Matrix *long2 = Longueur(Xd[1], Yd[1]);
    double len2 = long2->Element(long2->GetLignes() - 1, 0);

    //if (DEBUG) printf("\n in (%f %f)", len1, len2);
    //if (DEBUG) printf(" out (%f %f)", lp1, lp2);
    getLayoutLogger()->logprintf("\n(%d, %d) _d1=%f _d2=%f", noNerv1, noNerv2, 1000 * (lp1 - len1), 1000 * (lp2 - len2));

    double coeff1 = lp1 / len1;
    double coeff2 = lp2 / len2;
    double coeff = (coeff1 + coeff2) / 2.0f;
    newXd[0] = MultReelMat(Xd[0], coeff);
    newYd[0] = MultReelMat(Yd[0], coeff);
    newXd[1] = MultReelMat(Xd[1], coeff);
    newYd[1] = MultReelMat(Yd[1], coeff);
    delete (long1);
    delete (long2);
    long1 = Longueur(newXd[0], newYd[0]);
    len1 = long1->Element(long1->GetLignes() - 1, 0);
    long2 = Longueur(newXd[1], newYd[1]);
    len2 = long2->Element(long2->GetLignes() - 1, 0);
    getLayoutLogger()->logprintf("... d1=%f d2=%f", 1000 * (lp1 - len1), 1000 * (lp2 - len2));
    if (1000 * (lp1 - len1) < gfd->tochnostLayout2) getLayoutLogger()->logprintf(" OK");
    else getLayoutLogger()->logprintf(" BAD");
    if (1000 * (lp2 - len2) < gfd->tochnostLayout2) getLayoutLogger()->logprintf(" OK");
    else getLayoutLogger()->logprintf(" BAD");
    delete (long1);
    delete (long2);
    *rcoeff = coeff;
}



void goCalcIndepPince(WindPatternsProject* gfd, int noNerv, int face, double *pLA, double *pLF, double *pRA, double *pRF, double *len) {
    double mode1 = gfd->modePinces;
    double mode2 = gfd->modePinces;
    Matrix * Xd1[2], *Yd1[2], *Xd1p[2], *Yd1p[2];
    Matrix * Xd2[2], *Yd2[2], *Xd2p[2], *Yd2p[2];
    Matrix * X1[2], *Y1[2], *Z1[2], *P1[2], *newP1[2];
    Matrix * X2[2], *Y2[2], *Z2[2], *P2[2], *newP2[2];
    bool showPoints = false;
    int i = 0;
    calcPatron(gfd, noNerv - 1, false, face, face, 0.0f, 100.0f,
            noNerv, false, face, face, 0.0f, 100.0f,
            &Xd1[0], &Yd1[0], &Xd1[1], &Yd1[1],
            &X1[0], &Y1[0], &Z1[0], &P1[0],
            &X1[1], &Y1[1], &Z1[1], &P1[1]);
    Matrix *lenm1 = Longueur(Xd1[1], Yd1[1]);
    double len1 = lenm1->Element(lenm1->GetLignes() - 1, 0);
    double w1, w2;
    calcWidth(Xd1[0], Yd1[0], Xd1[1], Yd1[1], &w1);
    calcPatron(gfd, noNerv, false, face, face, 0.0f, 100.0f,
            noNerv + 1, false, face, face, 0.0f, 100.0f,
            &Xd2[0], &Yd2[0], &Xd2[1], &Yd2[1],
            &X2[0], &Y2[0], &Z2[0], &P2[0],
            &X2[1], &Y2[1], &Z2[1], &P2[1]);
    Matrix *lenm2 = Longueur(Xd2[0], Yd2[0]);
    double len2 = lenm2->Element(lenm2->GetLignes() - 1, 0);
    calcWidth(Xd2[0], Yd2[0], Xd2[1], Yd2[1], &w2);

    double ppA = gfd->PosPinceBA[0];
    double apAPercInit = gfd->AmpPinceBA[0];
    double ppF = gfd->PosPinceBF[0];
    double apFPercInit = gfd->AmpPinceBF[0];

    double apASumAbs = (w1 + w2) * apAPercInit / 100.0f;
    double apFSumAbs = (w1 + w2) * apFPercInit / 100.0f;

    double middleInit = w1 / (w1 + w2);
    double apAAbs1, apFAbs1, apAAbs2, apFAbs2;
    double middle = 0.0f;
    Matrix *long1p1, *long2p0;
    double minDelta = 1000000.0f, minMiddle = -1.0f, minDeltaMiddleInit = 1000000.0f, minlp1 = 1000000.0f, minlp2 = 1000000.0f;
    double l1p1 = -1.0f, l2p0 = -1.0f, delta = 0.0f, deltaMiddleInit = 0.0f, mina1, mina2;
    for (middle = 0.0f; middle <= 1.0f; middle += 0.001f) {
        apAAbs1 = apASumAbs*middle;
        apAAbs2 = apASumAbs - apAAbs1;

        apFAbs1 = apFSumAbs*middle;
        apFAbs2 = apFSumAbs - apFAbs1;

        calcPinceAloneAbs(gfd, Xd1[1], Yd1[1], Xd1[0], Yd1[0],
                X1[1], Y1[1], Z1[1], P1[1], X1[0], Y1[0], Z1[0], P1[0],
                ppA, apAAbs1, ppF, apFAbs1, mode1, &Xd1p[1], &Yd1p[1], &newP1[1]);

        long1p1 = Longueur(Xd1p[1], Yd1p[1]);
        l1p1 = long1p1->Element(long1p1->GetLignes() - 1, 0);
        calcPinceAloneAbs(gfd, Xd2[0], Yd2[0], Xd2[1], Yd2[1], X2[0], Y2[0], Z2[0], P2[0], X2[1], Y2[1], Z2[1], P2[1],
                ppA, apAAbs2, ppF, apFAbs2, mode2, &Xd2p[0], &Yd2p[0], &newP2[0]);
        long2p0 = Longueur(Xd2p[0], Yd2p[0]);
        l2p0 = long2p0->Element(long2p0->GetLignes() - 1, 0);
        delta = fabs(l1p1 - l2p0);
        deltaMiddleInit = fabs(middle - middleInit);


        if (minDelta < (gfd->tochnostLayout2 / 1000.0f)) {
            if (deltaMiddleInit < minDeltaMiddleInit) {
                minDelta = delta;
                minlp1 = l1p1;
                minMiddle = middle;
                minDeltaMiddleInit = deltaMiddleInit;
                mina1 = apAAbs1;
                mina2 = apAAbs2;
                minlp2 = l2p0;
            }
        } else {
            if (delta < minDelta) {
                minDelta = delta;
                minlp1 = l1p1;
                minMiddle = middle;
                minDeltaMiddleInit = deltaMiddleInit;
                mina1 = apAAbs1;
                mina2 = apAAbs2;
                minlp2 = l2p0;
            }
        }

        //printf ("\n [%f] delta=%f  %f/%f", middle, delta,l1p1,l2p0);
        delete (long1p1);
        delete (Xd1p[1]);
        delete (Yd1p[1]);
        delete (long2p0);
        delete (Xd2p[0]);
        delete (Yd2p[0]);
        delete (newP1[1]);
        delete (newP2[0]);
    }
    getLayoutLogger()->logprintf("\n%2df%d %6.4f(%5.3f)[%5.3f] %9.6f/%9.6f l%9.6f a%6.4f a%6.4f",
            noNerv, face, (minDelta * 1000), minMiddle, minMiddle - middleInit, minlp1, minlp2, len1, mina1, mina2);
    if (minDelta < (0.5 * 0.001)) {
        getLayoutLogger()->logprintf(" OK:)");
    } else {
        getLayoutLogger()->logprintf(" BAD!");
    }
    *len = (minlp1 + minlp2) / 2.0f;
    *pLA = apASumAbs*minMiddle;
    *pRA = apASumAbs - *pLA;
    *pLF = apFSumAbs*minMiddle;
    *pRF = apFSumAbs - *pLF;

    delete (X1[0]);
    delete (Y1[0]);
    delete (Z1[0]);
    delete (P1[0]);
    delete (X1[1]);
    delete (Y1[1]);
    delete (Z1[1]);
    delete (P1[1]);
    delete (X2[0]);
    delete (Y2[0]);
    delete (Z2[0]);
    delete (P2[0]);
    delete (X2[1]);
    delete (Y2[1]);
    delete (Z2[1]);
    delete (P2[1]);
    delete (Xd1[0]);
    delete (Yd1[0]);
    delete (Xd2[0]);
    delete (Yd2[0]);
    delete (Xd1[1]);
    delete (Yd1[1]);
    delete (Xd2[1]);
    delete (Yd2[1]);
    delete (lenm1);
    delete (lenm2);
}


int calcRepPointByLength(Matrix* Xd, Matrix* Yd, double l, double* x, double* y) {

	if (l < 0) {
		getLayoutLogger()->logprintf("\n!!! calcRepPointByLength l=%f return -1", l);
		return -1;
	}
	if (l == 0) {
		*x = Xd -> Element(0, 0);
		*y = Yd -> Element(0, 0);
		return 0;
	}
	int i = 0;
	double sumLength = 0;

	while ( (i < (Xd->GetLignes() - 1)) && (sumLength < l) ) {
		sumLength += dist2d (Xd->Element(i, 0), Yd->Element(i, 0), Xd->Element(i + 1, 0), Yd->Element(i + 1, 0));
		i++;
	}

	if ((i == (Xd->GetLignes() - 1)) && (sumLength < l)) {
		//getLayoutLogger()->logprintf("\n calcRepPointByLength l=%f sumLength=%f return 1 (put point to the end:))", l, sumLength);
		*x = Xd -> Element(i, 0);
		*y = Yd -> Element(i, 0);
		return 1;
	}

	double x1 = Xd->Element(i - 1, 0);
	double y1 = Yd->Element(i - 1, 0);
	double x2 = Xd->Element(i, 0);
	double y2 = Yd->Element(i, 0);

	double l2 = sumLength;
	double l1 = sumLength - dist2d (x1, y1, x2, y2);

	*x = x1 + (x2 - x1)*(l - l1) / (l2 - l1);
	*y = y1 + (y2 - y1)*(l - l1) / (l2 - l1);
	return 0;
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


Matrix* getReperPoints(WindPatternsProject* gfd, Matrix* Xd, Matrix* Yd, Matrix* P, 
						int nerv, float deb, int faceDeb, float fin, int faceFin, bool klapanShov)
{
    //printf ("\nget Reper Points()");
	//printf ("\n P:");
	//printM (P);
	
	//printf ("\nP");
	//printM(P);
	//printf ("\n...");
    Matrix *interpXSuspente, *interpYSuspente, *interpXSuspente0, *interpYSuspente0;
    bool cerezNos = false;
    int i0 = 0;
    if ((faceDeb != faceFin) && (deb != 0.0f)) {
        cerezNos = true;
      //  printf ("\n cerezNos");
        while (P->Element(i0, 0) != 0.0f) {
        //    printf ("\n i0=%d", i0);
            i0++;
         }
     }
    //printf ("...1 i0=%d", i0);
    interpXSuspente = Zeros(P -> GetLignes() - i0, 2);
    interpYSuspente = Zeros(P -> GetLignes() - i0, 2);
    for (int j = i0; j < P->GetLignes(); j++) {
        interpXSuspente->SetElement(j - i0, 0, P->Element(j, 0));
        interpYSuspente->SetElement(j - i0, 0, P->Element(j, 0));
        interpXSuspente->SetElement(j - i0, 1, Xd->Element(j, 0));
        interpYSuspente->SetElement(j - i0, 1, Yd->Element(j, 0));
    }
    //printf ("...2");
    int isave = 0;
    if (cerezNos) {
		//printf ("\n cerezNos");
        interpXSuspente0 = Zeros(i0 + 1, 2);
        interpYSuspente0 = Zeros(i0 + 1, 2);
        for (int j = i0; j >= 0; j--) {
            interpXSuspente0->SetElement(isave, 0, P->Element(j, 0));
            interpYSuspente0->SetElement(isave, 0, P->Element(j, 0));
            interpXSuspente0->SetElement(isave, 1, Xd->Element(j, 0));
            interpYSuspente0->SetElement(isave, 1, Yd->Element(j, 0));
            isave++;
        }
    }
    //printf ("...3");
    //printf ("\n interpYSuspente0->GetLignes()=%d", interpYSuspente0->GetLignes());
    int ndn = 0;
    double posda = 0.0f, posdf = 100.0f;
    for (int _d = 0; _d < gfd->quantDiag; _d++) {
        int dnerv = gfd->noNervD[_d] - 1;
        int f = 1;
        if ((nerv == dnerv) && (faceDeb == f) && (faceFin == f)) {
            ndn = 2;
            posda = gfd->PosDiagNerv1A;
            posdf = gfd->PosDiagNerv1F;
            break;
        }
        dnerv = gfd->noNervD[_d];
        f = 2;
        if ((nerv == dnerv) && (faceDeb == f) && (faceFin == f)) {
            ndn = 2;
            posda = gfd->PosDiagNerv2A;
            posdf = gfd->PosDiagNerv2F;
            break;

        }
        dnerv = gfd->noNervD[_d] + 1;
        f = 1;
        if ((nerv == dnerv) && (faceDeb == f) && (faceFin == f)) {
            ndn = 2;
            posda = gfd->PosDiagNerv1A;
            posdf = gfd->PosDiagNerv1F;
            break;

        }
    }
    //printf ("...4");
    int nrp;
    if (gfd->ReperPointsFromFile) nrp = gfd->ReperPoints[1] -> GetLignes();
    else nrp = P -> GetLignes();
    int n = 20 * nrp + 12;
    Matrix* res = new Matrix(n, 5);
    isave = 0;
    double posSuspente, x, y;
    int mark, MARK_SUSPENTE=REP_TRIANGLE, MARK_REP=REP_CROSS, MARK_DIAG=REP_LINE;

    if (gfd->ReperPointsPlotterFormat) {
        MARK_SUSPENTE = REP_V;
        MARK_REP = REP_MIDDLE_LINE;
        MARK_DIAG = REP_LINE;
    }
    
    for (int j = 0; j < 11; j++) {
        switch (j) {
            case 0: posSuspente = gfd->Form->m_pProfils[nerv]->m_fPosA; break;
            case 1: posSuspente = gfd->Form->m_pProfils[nerv]->m_fPosB; break;
            case 2: posSuspente = gfd->Form->m_pProfils[nerv]->m_fPosC; break;
            case 3: posSuspente = gfd->Form->m_pProfils[nerv]->m_fPosD; break;
            case 4: posSuspente = gfd->Form->m_pProfils[nerv]->m_fPosE; break;
            case 5: posSuspente = posda; break;
            case 6: posSuspente = posdf; break;
            case 7: posSuspente = 0.0f; break;
            case 8: posSuspente = gfd->VentHolesDeb; break;
            case 9: posSuspente = gfd->VentHolesFin; break;
            case 10: posSuspente = gfd->PosKlapanFin; break;
        }
        //printf ("...5.1");
        mark = MARK_SUSPENTE;
        if ((j == 5) || (j == 6)) {
            if (ndn == 0) continue;
            mark = MARK_DIAG;
        }
        //printf ("\n-- j=%d", j);
        if (j == 7) {
            if (!klapanShov) break;
            mark=REP_0;
            //printf ("\nposSuspente=%f isave=%d", posSuspente, isave);
        }
        if (j == 8) {
            if (!klapanShov) continue;
            mark=REP_VENTDEB;
            //printf ("\nposSuspente=%f isave=%d", posSuspente, isave);
        }
        if (j == 9) {
            if (!klapanShov) continue;
            mark=REP_VENTFIN;
            //printf ("\nposSuspente=%f isave=%d", posSuspente, isave);
        }
        if (j == 10) {
            if (!klapanShov) continue;
            mark=REP_KLAPAN_FIN;
            //printf ("posSuspente=%f isave=%d", posSuspente, isave);
        }
        
        //printf ("...5.2");

        if ((!cerezNos && (deb <= posSuspente) && (posSuspente <= fin))
                || (cerezNos && (posSuspente <= fin))) {
            x = InterpLinX(interpXSuspente, posSuspente);
            y = InterpLinX(interpYSuspente, posSuspente);
            res->SetElement(isave, 0, x);
            res->SetElement(isave, 1, y);
            res->SetElement(isave, 2, mark);
            if (!cerezNos && (deb <= posSuspente) && (posSuspente <= fin)) {
                res->SetElement(isave, 3, j+1);
            } else {
                res->SetElement(isave, 3, 0);
            }
			res->SetElement(isave, 4, posSuspente);
            isave++;
        }
        //printf ("...5.3 posSuspente=%f deb=%f", posSuspente, deb);
        //if (cerezNos) printf ("\ncerezNos=1"); else printf ("\ncerezNos=2");

        if ((cerezNos) && (posSuspente <= deb)) {
            x = InterpLinX(interpXSuspente0, posSuspente);
            y = InterpLinX(interpYSuspente0, posSuspente);
            res->SetElement(isave, 0, x);
            res->SetElement(isave, 1, y);
            res->SetElement(isave, 2, mark);
            res->SetElement(isave, 3, 0);
			res->SetElement(isave, 4, posSuspente);
            isave++;
        }
    }
    
    //printf ("...6");
    //if (cerezNos) printf ("\n cerezNos"); else printf ("\n NOT cerezNos");
    //printf ("\nbefore Reper Points from file");
    if (gfd->ReperPointsFromFile) {
        //printf ("...6.1");
        //printf ("\nIN Reper Points from file");
        if ((cerezNos) && (interpXSuspente0->GetLignes()>1)) {
            for (int j = 0; j < gfd->ReperPoints[faceDeb]->GetLignes(); j++) {
                posSuspente = gfd->ReperPoints[faceDeb]->Element(j, 0);
                //printf ("\ncZ: posSuspente=%f deb=%f", posSuspente, deb);
                if (posSuspente <= deb) {
                    //printf ("\nbefore InterpLinX");
                    x = InterpLinX(interpXSuspente0, posSuspente);
                    y = InterpLinX(interpYSuspente0, posSuspente);
                    //printf ("...after InterpLinX");
                    res->SetElement(isave, 0, x);
                    res->SetElement(isave, 1, y);
                    res->SetElement(isave, 2, MARK_REP);
                    res->SetElement(isave, 3, 0);
					res->SetElement(isave, 4, posSuspente);
                    isave++;
                }
            }
        }
        //printf ("\nnot cerez nos");
        for (int j = 0; j < gfd->ReperPoints[faceFin]->GetLignes(); j++) {
            posSuspente = gfd->ReperPoints[faceFin]->Element(j, 0);
            if ((!cerezNos && (deb <= posSuspente) && (posSuspente <= fin))
                    || (cerezNos && (posSuspente <= fin))) {
                x = InterpLinX(interpXSuspente, posSuspente);
                y = InterpLinX(interpYSuspente, posSuspente);
                res->SetElement(isave, 0, x);
                res->SetElement(isave, 1, y);
                res->SetElement(isave, 2, MARK_REP);
                res->SetElement(isave, 3, 0);
				res->SetElement(isave, 4, posSuspente);
                isave++;
            }
        }
    } else {
        //printf ("...6.2");
      //  printf ("\nReper Points NOT FILE");
        if (cerezNos) {
            for (int j = i0; j > 0; j--) {
                posSuspente = P->Element(j, 0);
                x = InterpLinX(interpXSuspente0, posSuspente);
                y = InterpLinX(interpYSuspente0, posSuspente);
                res->SetElement(isave, 0, x);
                res->SetElement(isave, 1, y);
                res->SetElement(isave, 2, MARK_REP);
                res->SetElement(isave, 3, 0);
				res->SetElement(isave, 4, posSuspente);
                isave++;
            }
        }
        for (int j = i0; j < P->GetLignes(); j++) {
            posSuspente = P->Element(j, 0);
            x = InterpLinX(interpXSuspente, posSuspente);
            y = InterpLinX(interpYSuspente, posSuspente);
            res->SetElement(isave, 0, x);
            res->SetElement(isave, 1, y);
            res->SetElement(isave, 2, MARK_REP);
            res->SetElement(isave, 3, 0);
			res->SetElement(isave, 4, posSuspente);
            isave++;
        }
    }
    //printf ("...7");
    //printf ("\n isave=%d", isave);
    Matrix* resexit = new Matrix(isave, 5);
    for (int i = 0; i < isave; i++)
        for (int j = 0; j < 5; j++) resexit->SetElement(i, j, res->Element(i, j));
    //printf ("\n dels");
    delete (res);
    delete (interpXSuspente);
    delete (interpYSuspente);

    if (cerezNos) {
        delete (interpXSuspente0);
        delete (interpYSuspente0);
    }
   //printf ("\n...END get Reper Points()");
   return resexit;
}


double calcPolyLength(Matrix* X, Matrix* Y, double xrp, double yrp) {
	int i  = 0, res = 0;
	int n = X->GetLignes();
	double sumLength = 0;
	if ( (xrp==X->Element(0, 0)) && (yrp==Y->Element(0, 0)) ) return 0;

	if (X->Element(0,0) > X->Element(X->GetLignes() - 1,0)) getLayoutLogger()->logprintf ("\n calcPolyLength()... Achtung! X->Element(0,0) > X->Element(%d,0)", X->GetLignes() - 1);

	for (i = 0; i < n-1;i++ ) {
		res = pointAtSegment(xrp, yrp, X->Element(i,0), Y->Element(i,0), X->Element(i+1,0), Y->Element(i+1,0));
		if (1 == res ) break;
		sumLength += dist2d (X->Element(i,0), Y->Element(i,0), X->Element(i+1,0), Y->Element(i+1,0));
	}
	
	if (1 == res) 
		return sumLength + dist2d(xrp, yrp, X->Element(i,0), Y->Element(i,0));

	return -1;
}

Matrix* getReperPointsPince(WindPatternsProject* gfd,  Matrix* Xd0, Matrix* Yd0, double coeff, Matrix* Xd, Matrix* Yd, Matrix* P, 
								int nerv, float deb, int faceDeb, float fin, int faceFin, bool klapanShov)
{
	// printf ("\n\n\n Xd0 Yd0");
	// printXY (Xd0, Yd0);

	// printf ("\n\n\n Xd Yd");
	// printXY (Xd, Yd);

	if (Xd0->GetLignes() != P->GetLignes()) {
		printf ("\n!!!\n!!!\n!!!ACHTUNG! getRepPointsPince() Xd0->GetL[%d] != P->GetL[%d]\n!!!", Xd0->GetLignes(), P->GetLignes());
	}

	Matrix *reperPoints0 = getReperPoints(gfd, Xd0, Yd0, P, nerv, deb, faceDeb, fin, faceFin, klapanShov);

	//printf ("\n=1===============================================");
	Matrix *reperPointsP = getReperPoints(gfd, Xd, Yd, P, nerv, deb, faceDeb, fin, faceFin, klapanShov);
	//printf ("\n=2===============================================");
	
	//printf ("\n\n Xd0 Yd0");
    //printXY(Xd0, Yd0);

	printf ("\n(Xd0, Yd0).length=%f coeff=%f", 1000*calcCourbeLength(Xd0, Yd0), coeff);

	for (int i = 0; i < reperPoints0->GetLignes(); i++) {
		// printf ("\n ");
		double xrp = reperPoints0 -> Element(i, 0);
		double yrp = reperPoints0 -> Element(i, 1);

		double xrpP = reperPointsP -> Element(i, 0);
		double yrpP = reperPointsP -> Element(i, 1);

		double l = calcPolyLength(Xd0, Yd0, xrp, yrp);
		// printf ("\n l=%f l*coeff=%f", l, l*coeff);
		// printf ("\n xrp, yrp (%f, %f)", xrp, yrp);
		double lP0 = calcPolyLength(Xd, Yd, xrpP, yrpP);

		double newx = 0, newy = 0;
		
		int res = calcRepPointByLength(Xd, Yd, l*coeff, &newx, &newy);

		// printf ("\n res=%d newx, newy (%f, %f)", res, newx, newy);

		if (res) {
			//getLayoutLogger()->logprintf ("\n Ahtung! calcRepPointByLength() == %d ", res);
			//getLayoutLogger()->logprintf (" l=%f   coeff=%f", l, coeff);
			double XY0length =  calcCourbeLength(Xd0, Yd0);
			double XYlength =  calcCourbeLength(Xd, Yd);
			double realCoeff = XYlength / XY0length;

			if (res == -1) {

			}

			if (res == 1) {
				//getLayoutLogger()->logprintf ("\n! XY0length=%f  XYlength=%f coeff=%f", XY0length, XYlength, realCoeff);
				//getLayoutLogger()->logprintf ("\n delta=%f", abs(XYlength-XY0length));
				/*if (abs(XYlength - XY0length) < 0.0005)  
					getLayoutLogger()->logprintf(" OK:)!");
				else 
					getLayoutLogger()->logprintf(" BAD:(!");*/
			}
		}
		
		//getLayoutLogger()->logprintf ("\n%2d R res=%d pS=%6.2f l=%6.2f (%6.2f, %6.2f)", i, res, reperPoints0->Element(i,4), 1000*l, 1000*xrp, 1000*yrp); 

		/*getLayoutLogger()->logprintf ("\n%d R res=%d pS=%f cf=%f newl=%f  dxyrp=%f  dl=%f",
										  i,     res, reperPoints0->Element(i,4), coeff, 1000*l*coeff, 1000*dist2d(xrpP,yrpP,newx,newy), 1000*(lP0-l*coeff));
		getLayoutLogger()->logprintf ("\n l==%f l*coeff=%f (%f, %f)", 1000*l, 1000*l*coeff, xrp, yrp);*/
		reperPoints0->SetElement(i, 0, newx);
		reperPoints0->SetElement(i, 1, newy);
	}

	return reperPoints0;
}