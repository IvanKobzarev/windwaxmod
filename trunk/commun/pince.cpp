

#pragma warning(disable:4514)
#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>

#include "matrice.h"
#include "pince.h"
#include "geom.h"
#include "rasklad.h"
#define sqr(f1) ((f1)*(f1))
#define pi	3.141592675f

#ifndef DEBUG
#define DEBUG false
#endif


Pince :: Pince( int n1, int n2 ) : nerv1(n1), nerv2(n2)
{
    diffAmps=false;
    debug=false;
    glue = false;
    i01 = 0;
	i02 = 0;
}

Pince::~Pince()
{
}

void Pince::SetFunction1( Matrice *func )
{
    this->function1 = func;
}

void Pince::SetFunction2( Matrice *func )
{
    this->function2 = func;
}

void Pince::SetP1( Matrice *p )
{
    this->P1 = p;
}
void Pince::SetP2( Matrice *p )
{
    this->P2 = p;
}

void Pince::print()
{
    for (int i = 0; i < function1 -> GetLignes(); i++) {
        printf ("\n %d -> %f", i , function1->Element(i, 0));
    }
}

void Pince::SetAmps(double aa, int _ia,  double af, int _if) 
{
    this->AmpA = aa;
    this->ipa1=_ia;
    this->ipa2=_ia;
    this->AmpF = af;
    this->ipf1 =_if;
    this->ipf2 = _if;
}

void Pince::SetDiffAmps (double a1, double f1, double a2, double f2 )
{
    this->AmpA1 = a1;
    this->AmpF1 = f1;

    this->AmpA2 = a2;
    this->AmpF2 = f2;
    this->diffAmps=true;
}

void Pince::SetDiffAmps0 (double a1,  double a2 )
{
    this->AmpA01 = a1;
    this->AmpA02 = a2;
}

void Pince::SetPos (double pa, double pf) {
    this->PosA = pa;
    this->PosF = pf;
}

double CalculW0byHW1(double h, double w1) {
	if (h > (w1*0.5)) printf ("\n!!!Error: in Radius alg h < w1/2 => CalculW0byHW1 not work correct!");
    if (h == 0.0f) return w1;
    double r = (w1*w1)/(8*h) + h/2;
    //if (DEBUG) printf ("\nCalculW0byHW1 r=%f", r)
    double sinv=0.5f*w1/r;
    //printf (" sinv=%f", sinv);
    if (fabs(sinv) > 1.0f) printf ("\n!!! ASIN from > 1 (%f)!?!?!??!", sinv);
    double alpha = asin (sinv);
    //printf (" alpha=%f", alpha);
    double w0 = 2 * alpha * r;
    return w0;
}

double CalculHbyw0w1(double w0, double w1) {
    double minDelta = 1000000000.0f, delta = 0.0f;
    double rside = w1/w0;
    double alpha = 0.0f;
    double ih=0.00001;
    for (double ialpha = ih; ialpha < 3.14159265354f; ialpha += ih) {
        delta = (double) fabs (sin(ialpha)/ialpha - rside);
        if (delta < minDelta) {
            minDelta = delta;
            alpha = ialpha;
        }
    }
    double r = w0/alpha;
    return (r*(1-cos(alpha)));
}

double functionPincePower(double x, double power) {
    return pow(x, power);
}

double functionPinceArctanSym(double x, double k1, double k2) {
    return (
            (atan(k1*(2*k2*x-k2))
            /
            atan(k1*k2)
            +1.0f)/2.0f);
}

double functionPinceArctan(double x, double k1, double k2, double k3) {
    double x0 = -k3 * (2*k2);
    double x1 = (1-k3) * (2*k2);
    return (atan(k1*(2*k2*x+x0))-atan(k1*x0))
            /
            (atan(k1*x1)-atan(k1*x0));
}

void makePointPince(Matrice *Xd, Matrice *Yd, Matrice *P,
        double PosPinceBA, double PosPinceBF, Matrice **Xr, Matrice **Yr, Matrice **Pr, int *ipair, int *ipfir)
{
    Matrice *interpXSuspente, *interpYSuspente;
    interpXSuspente = Zeros(P->GetLignes(), 2);
    interpYSuspente = Zeros(P->GetLignes(), 2);
    for (int j = 0; j < P->GetLignes(); j++) {
        interpXSuspente->SetElement(j, 0, P->Element(j, 0));
        interpYSuspente->SetElement(j, 0, P->Element(j, 0));
        interpXSuspente->SetElement(j, 1, Xd->Element(j, 0));
        interpYSuspente->SetElement(j, 1, Yd->Element(j, 0));
    }
    double xdPBA = InterpLinX(interpXSuspente, PosPinceBA);
    double ydPBA = InterpLinX(interpYSuspente, PosPinceBA);
    double xdPBF = InterpLinX(interpXSuspente, 100.0f - PosPinceBF);
    double ydPBF = InterpLinX(interpYSuspente, 100.0f - PosPinceBF);
    int ipai = 0;
    int ipfi = 0;
    bool pxba = false, pxbf = false;
    int k = 2, isave = 0;
    int i = 0;
    for (i = 0; i < Xd->GetLignes(); i++) {
        if (xdPBA == Xd->Element(i, 0)) {
            pxba = true;
            ipai=i;
            k--;
        }
        if (xdPBF == Xd->Element(i, 0)) {
            pxbf = true;
            ipfi=i;
            k--;
        }
    }

    (*Xr) = new Matrice(Xd->GetLignes() + k, 1);
    (*Yr) = new Matrice(Xd->GetLignes() + k, 1);
    (*Pr) = new Matrice(Xd->GetLignes() + k, 1);

    i = 0;
    while (i < Xd->GetLignes()) {
        if (!pxba)
            if (Xd->Element(i, 0) > xdPBA) {
                (*Xr)->SetElement(isave, 0, xdPBA);
                (*Yr)->SetElement(isave, 0, ydPBA);
                (*Pr)->SetElement(isave, 0, PosPinceBA);
                ipai = isave;
                isave++;
                pxba = true;
            }
        if (!pxbf)
            if (Xd->Element(i, 0) > xdPBF) {
                (*Xr)->SetElement(isave, 0, xdPBF);
                (*Yr)->SetElement(isave, 0, ydPBF);
                (*Pr)->SetElement(isave, 0, 100.0f - PosPinceBF);
                ipfi = isave;
                isave++;
                pxbf = true;
            }
        (*Xr)->SetElement(isave, 0, Xd->Element(i, 0));
        (*Yr)->SetElement(isave, 0, Yd->Element(i, 0));
        (*Pr)->SetElement(isave, 0, P->Element(i, 0));
        i++;
        isave++;
    }
    (*ipair) = ipai;
    (*ipfir) = ipfi;
    delete (interpXSuspente);
    delete (interpYSuspente);
}


double pinceFunctionA(WindPatternsProject* gfd, double x) {
    double sum = 0.0f;
    double k=0.0f;
    if (gfd->PincePowerA) {
        sum=sum + functionPincePower(x, (1.0f/gfd->PincePowerValueA));
        k = k + 1.0f;
    }
    if (gfd->PinceArctanA) {
        sum=sum + functionPinceArctan(x, gfd->PinceArctanK1ValueA, gfd->PinceArctanK2ValueA, gfd->PinceArctanK3ValueA);
        k = k + 1.0f;
    }
    if (k == 0.0f) return 0.0f;
    return (sum / k);
}

double pinceFunctionF(WindPatternsProject* gfd, double x) {
    double sum = 0.0f;
    double k=0.0f;
    if (gfd->PincePowerF) {
        sum = sum + functionPincePower(x, (1.0f/gfd->PincePowerValueF));
        k = k + 1.0f;
    }
    if (gfd->PinceArctanF) {
        sum = sum + functionPinceArctan(x, gfd->PinceArctanK1ValueF, gfd->PinceArctanK2ValueF, gfd->PinceArctanK3ValueF);
        k = k + 1.0f;
    }
    if (k==0.0f) return 0.0f;
    return (sum/k);
}


void CalculPinceAloneAbs(WindPatternsProject* gfd, Matrice *Xd, Matrice *Yd,
        Matrice *Xd0, Matrice *Yd0,
        Matrice *X, Matrice *Y, Matrice *Z, Matrice *P,
        Matrice *X0, Matrice *Y0, Matrice *Z0, Matrice *P0,
        double PosPinceBA, double AmpPinceBA, double PosPinceBF, double AmpPinceBF, double modePincesT,
        Matrice **Xdp, Matrice **Ydp, Matrice **Pp) {

    double amp, larg;
    Matrice *newXd, *newYd, *newP;
    Matrice *newXd0, *newYd0, *newP0;
    int ipai = 0;
    int ipfi = 0;
    int i = 0;
    int ipai0 = 0;
    int ipfi0 = 0;

    makePointPince(Xd, Yd, P, PosPinceBA, PosPinceBF, &newXd, &newYd, &newP, &ipai, &ipfi);
    makePointPince(Xd0, Yd0, P0, PosPinceBA, PosPinceBF, &newXd0, &newYd0, &newP0, &ipai0, &ipfi0);
    double ampPinceAabs = AmpPinceBA;
    double ampPinceFabs = AmpPinceBF;
    Matrice* long1 = Longueur(newXd, newYd);
    double l1 = long1->Element(long1->GetLignes() - 1, 0);
    delete (long1);
    double x1, y1, x2, y2;
	double argument=-10000.0f;
    bool showAmp = 0;
    for (i = 0; i < newXd->GetLignes(); i++) {
        if (i < ipai) {
	    argument = (newXd->Element(i, 0) - newXd->Element(0, 0))  /  (newXd->Element(ipai, 0) - newXd->Element(0, 0));
            amp = ampPinceAabs * pinceFunctionA(gfd, argument);
            if (showAmp) printf("\n_1_%d -> %8.5f", i, amp);
        } else {
            if (i > ipfi) {
                argument = (newXd->Element(newXd->GetLignes() - 1, 0) - newXd->Element(i, 0))  /   (newXd->Element(newXd->GetLignes() - 1, 0) - newXd->Element(ipfi, 0));
                amp = ampPinceFabs * pinceFunctionF(gfd, argument);
                if (showAmp) printf("\n_3_%d -> %8.5f", i, amp);
            } else {
                x1 = ipai;
                x2 = ipfi;
                y1 = ampPinceAabs;
                y2 = ampPinceFabs;
                amp = (i * (y2 - y1) + (y1 * (x2 - x1) - x1 * (y2 - y1))) / (x2 - x1);
                if (showAmp) printf("\n_2_%d -> %8.5f", i, amp);
                /*if (ampPinceAabs <= ampPinceFabs) {

                        amp = ampPinceAabs + (ampPinceFabs - ampPinceAabs) * pow(
                        (Xd->Element(i,0)-Xd->Element(ipai,0))
                        /
                        (Xd->Element(ipfi,0)-Xd->Element(ipai,0)), modePincesT );

                } else {
                amp = ampPinceFabs + (ampPinceAabs - ampPinceFabs) * pow(
                        (Xd->Element(ipfi,0)-Xd->Element(i,0))
                        /
                        (Xd->Element(ipfi,0)-Xd->Element(ipai,0)), modePincesT );

                }*/
            }
        }
        larg = sqrt(sqr(newXd->Element(i, 0) - newXd0->Element(i, 0)) + sqr(newYd->Element(i, 0) - newYd0->Element(i, 0)));
        double dx = (newXd->Element(i, 0) - newXd0->Element(i, 0)) * amp / larg;
        double dy = (newYd->Element(i, 0) - newYd0->Element(i, 0)) * amp / larg;
        if (showAmp) printf(" dx=%f dy=%f", dx, dy);
        newXd->SetElement(i, 0, newXd->Element(i, 0) + (newXd->Element(i, 0) - newXd0->Element(i, 0)) * amp / larg);
        newYd->SetElement(i, 0, newYd->Element(i, 0) + (newYd->Element(i, 0) - newYd0->Element(i, 0)) * amp / larg);
    }
    *Xdp = CloneMat(newXd);
    *Ydp = CloneMat(newYd);
    *Pp = CloneMat(newP);
    delete(newXd);
    delete(newYd);
    delete(newP);
    delete(newXd0);
    delete(newYd0);
    delete(newP0);

}
void CalculWidth(Matrice *Xd, Matrice *Yd, Matrice *Xd0, Matrice *Yd0, double *width) {
    double max = -1000000.0f;
    double val = 0.0f;
    for (int i = 0; i < Xd->GetLignes(); i++) {
        val = sqrt(sqr(Xd->Element(i, 0) - Xd0->Element(i, 0))
                + sqr(Yd->Element(i, 0) - Yd0->Element(i, 0)));
        if (val > max) {
            max = val;
        }
    }
    *width = max;
}

void CalculPinceAlone(WindPatternsProject* gfd, Matrice *Xd, Matrice *Yd, Matrice *Xd0, Matrice *Yd0, Matrice *X, Matrice *Y, Matrice *Z, Matrice *P, Matrice *X0, Matrice *Y0, Matrice *Z0, Matrice *P0,
        double PosPinceBA, double AmpPinceBA, double PosPinceBF, double AmpPinceBF, double modePincesT,
        Matrice **Xdp, Matrice **Ydp, Matrice **Pp) {
    double width;
    CalculWidth(Xd, Yd, Xd0, Yd0, &width);
    double ampPinceAabs = AmpPinceBA * width / 100.0f;
    double ampPinceFabs = AmpPinceBF * width / 100.0f;
    CalculPinceAloneAbs(gfd, Xd, Yd, Xd0, Yd0,
            X, Y, Z, P, X0, Y0, Z0, P0, PosPinceBA, ampPinceAabs, PosPinceBF, ampPinceFabs, modePincesT,
            Xdp, Ydp, Pp);
}


void CalculPinceAloneAbsFunc(Matrice *Xd, Matrice *Yd,
                             Matrice *Xd0, Matrice *Yd0,
                            Matrice *X, Matrice *Y, Matrice *Z, Matrice *P,
        Matrice *X0, Matrice *Y0, Matrice *Z0, Matrice *P0, Pince* p, double ampa, double ampf, int no,
        Matrice **Xdp, Matrice **Ydp, Matrice **Pp)
{
    if (p->debug) printf ("\n CalculPinceAloneAbsFunc()...");
    double amp, larg;
    Matrice *newXd, *newYd, *newP;
    Matrice *newXd0, *newYd0, *newP0;
    Matrice* func;
    int ipai;// = p->ipa;
    int ipfi;// = p->ipf;

    if (no==0) { 
		func = p->function1; 
		ipai = p->ipa1;
		ipfi = p->ipf1;
	}
	else {
		func=p->function2;
		ipai = p->ipa2;
		ipfi = p->ipf2;
	}
    int i = 0;
	int ipai0 = 0;
	int ipfi0 = 0;
    double ampPinceAabs = ampa;
    double ampPinceFabs = ampf;
    makePointPince(Xd, Yd, P, p->PosA, p->PosF, &newXd, &newYd, &newP, &ipai, &ipfi);
    makePointPince(Xd0, Yd0, P0, p->PosA, p->PosF, &newXd0, &newYd0, &newP0, &ipai0, &ipfi0);
    double x1, y1, x2, y2;
    bool showAmp = 0;
    for (i = 0; i < newXd->GetLignes(); i++) {
        if (p->debug) printf ("\n %d ", i);
        if (i < ipai) {
            amp = ampPinceAabs * func->Element(i,0);
            if (p->debug)  printf("_1_%d -> %8.5f", i, amp);
        } else {
            if (i > ipfi) {
                amp = ampPinceFabs * func->Element(i,0);
                if (p->debug)  printf("_3_%d -> %8.5f", i, amp);
            } else {
                x1 = ipai;
                x2 = ipfi;
                y1 = ampPinceAabs;
                y2 = ampPinceFabs;
                amp = (i * (y2 - y1) + (y1 * (x2 - x1) - x1 * (y2 - y1))) / (x2 - x1);
                if (p->debug)  printf("_2_%d -> %8.5f", i, amp);
            }
        }
        larg = sqrt(sqr(newXd->Element(i, 0) - newXd0->Element(i, 0)) + sqr(newYd->Element(i, 0) - newYd0->Element(i, 0)));
        double dx = (newXd->Element(i, 0) - newXd0->Element(i, 0)) * amp / larg;
        double dy = (newYd->Element(i, 0) - newYd0->Element(i, 0)) * amp / larg;
        if (showAmp) printf(" dx=%f dy=%f", dx, dy);
        newXd->SetElement(i, 0, newXd->Element(i, 0) + (newXd->Element(i, 0) - newXd0->Element(i, 0)) * amp / larg);
        newYd->SetElement(i, 0, newYd->Element(i, 0) + (newYd->Element(i, 0) - newYd0->Element(i, 0)) * amp / larg);
    }
    *Xdp = CloneMat(newXd);
    *Ydp = CloneMat(newYd);
    *Pp = CloneMat(newP);
    delete(newXd);
    delete(newYd);
    delete(newP);
    delete(newXd0);
    delete(newYd0);
    delete(newP0);
}

void CalculPincePlusAloneAbsFunc(Matrice *Xd, Matrice *Yd,
                             Matrice *Xd0, Matrice *Yd0, Pince* p, int no,
        Matrice **Xdp, Matrice **Ydp)
{
    if (p->debug)  printf ("\n CalculPinceAloneAbsFunc()...p->ipa=%d no=%d", p->ipa1, no);
    double amp, larg;
    Matrice *newXd, *newYd;
    Matrice* func;
    double ampPinceA0abs = 0.0f;
    double ampPinceAabs = 0.0f;
    double ampPinceFabs = 0.0f;
    int i0=0, ipai=0, ipfi=0;
    if (no==0) 
    {
        func = p -> function1;
        ampPinceA0abs = p -> AmpA01;
        ampPinceAabs = p -> AmpA1;
        ampPinceFabs = p -> AmpF1;
        i0  = p-> i01;
        ipai = p->ipa1;
        ipfi = p->ipf1;

    } else {
        func=p->function2;
        ampPinceA0abs = p -> AmpA02;
        ampPinceAabs = p -> AmpA2;
        ampPinceFabs = p -> AmpF2;
        i0  = p-> i02;
        ipai = p->ipa2;
        ipfi = p->ipf2;

    }
    int i = 0;


    double x1, y1, x2, y2;
    bool showAmp = 0;

    
    printf ("\n Xd->GetLignes()=%d", Xd->GetLignes());
    printf (" func->GetLignes()=%d", func->GetLignes());
    

    newXd = new Matrice (Xd->GetLignes(), 1);
    newYd = new Matrice (Yd->GetLignes(), 1);
    for (i = 0; i < Xd->GetLignes(); i++) {
        if (p->debug) printf ("\n %d ", i);
        if (i < i0) {
            amp = ampPinceA0abs * func->Element(i,0);
            if (p->debug)  printf("_0_%d -> %8.5f", i, amp);
        }
        if ((i >= i0) && (i < ipai) ) {
            amp = ampPinceAabs * func->Element(i,0);
            if (p->debug)  printf("_1_%d -> %8.5f", i, amp);
        }
        if ((i >= ipai) && (i <= ipfi) ) {
               x1 = ipai;
               x2 = ipfi;
               y1 = ampPinceAabs;
               y2 = ampPinceFabs;
               amp = (i * (y2 - y1) + (y1 * (x2 - x1) - x1 * (y2 - y1))) / (x2 - x1);
               if (p->debug)  printf("_2_%d -> %8.5f", i, amp);
        }
        if (i > ipfi) {
                    amp = ampPinceFabs * func->Element(i,0);
                    if (p->debug)  printf("_3_%d -> %8.5f", i, amp);
        } 
            
        larg = sqrt(sqr(Xd->Element(i, 0) - Xd0->Element(i, 0)) + sqr(Yd->Element(i, 0) - Yd0->Element(i, 0)));
        double dx = (Xd->Element(i, 0) - Xd0->Element(i, 0)) * amp / larg;
        double dy = (Yd->Element(i, 0) - Yd0->Element(i, 0)) * amp / larg;
        if (showAmp) printf(" dx=%f dy=%f", dx, dy);
        newXd->SetElement(i, 0, Xd->Element(i, 0) + (Xd->Element(i, 0) - Xd0->Element(i, 0)) * amp / larg);
        newYd->SetElement(i, 0, Yd->Element(i, 0) + (Yd->Element(i, 0) - Yd0->Element(i, 0)) * amp / larg);
    }
    *Xdp = CloneMat(newXd);
    *Ydp = CloneMat(newYd);
    if (p->debug) printf ("\n ... end calculPincePlusAloneAbsFunc()...");
}

void CalculPinceNew(Matrice *Xd1, Matrice *Yd1, Matrice *Xd2, Matrice *Yd2,
        Matrice *X1, Matrice *Y1, Matrice *Z1, Matrice *P1,
        Matrice *X2, Matrice *Y2, Matrice *Z2, Matrice *P2, Pince* p,
        Matrice **Xdp1, Matrice **Ydp1, Matrice **Pp1,
        Matrice **Xdp2, Matrice **Ydp2, Matrice **Pp2)
{
    if (p->debug) printf ("\n CalculPinceNew...");
    CalculPinceAloneAbsFunc(Xd1, Yd1, Xd2, Yd2,
            X1, Y1, Z1, P1,
            X2, Y2, Z2, P2,p, p->AmpA1, p->AmpF1, 0,
            Xdp1, Ydp1, Pp1);
    CalculPinceAloneAbsFunc( Xd2, Yd2, Xd1, Yd1,
            X2, Y2, Z2, P2,
            X1, Y1, Z1, P1, p, p->AmpA2, p->AmpF2, 1,
            Xdp2, Ydp2, Pp2);
    if (p->debug) printf ("\n end calculPinceNew...");

}

void CalculPincePlusNew(Matrice *Xd1, Matrice *Yd1, Matrice *Xd2, Matrice *Yd2,
        Pince* p,
        Matrice **Xdp1, Matrice **Ydp1,
        Matrice **Xdp2, Matrice **Ydp2)
{
    if (p->debug) printf ("\n CalculPincePlusNew p->ipa1=%d", p->ipa1);
    CalculPincePlusAloneAbsFunc( Xd1, Yd1, Xd2, Yd2, p, 0, Xdp1, Ydp1);
    CalculPincePlusAloneAbsFunc( Xd2, Yd2, Xd1, Yd1, p, 1, Xdp2, Ydp2);
    if (p->debug) printf ("\n ...calculPincePlusNew");
}

void CalculPince(WindPatternsProject* gfd,Matrice *Xd1, Matrice *Yd1, Matrice *Xd2, Matrice *Yd2,
        Matrice *X1, Matrice *Y1, Matrice *Z1, Matrice *P1,
        Matrice *X2, Matrice *Y2, Matrice *Z2, Matrice *P2,
        double PosPinceBA, double AmpPinceBA, double PosPinceBF, double AmpPinceBF, double modePincesT,
        Matrice **Xdp1, Matrice **Ydp1, Matrice **Pp1, Matrice **Xdp2, Matrice **Ydp2, Matrice **Pp2) {

    CalculPinceAlone(gfd, Xd1, Yd1, Xd2, Yd2, X1, Y1, Z1, P1, X2, Y2, Z2, P2,
            PosPinceBA, AmpPinceBA, PosPinceBF, AmpPinceBF, modePincesT,
            Xdp1, Ydp1, Pp1);
    CalculPinceAlone(gfd, Xd2, Yd2, Xd1, Yd1, X2, Y2, Z2, P2, X1, Y1, Z1, P1,
            PosPinceBA, AmpPinceBA, PosPinceBF, AmpPinceBF, modePincesT,
            Xdp2, Ydp2, Pp2);
}


void CalculPinceDiff(WindPatternsProject* gfd,Matrice *Xd1, Matrice *Yd1, Matrice *Xd2, Matrice *Yd2,
        Matrice *X1, Matrice *Y1, Matrice *Z1, Matrice *P1,
        Matrice *X2, Matrice *Y2, Matrice *Z2, Matrice *P2,
        double PosPinceBA1, double AmpPinceBA1, double PosPinceBF1, double AmpPinceBF1,
        double PosPinceBA2, double AmpPinceBA2, double PosPinceBF2, double AmpPinceBF2,
        double mode,
        Matrice **Xdp1, Matrice **Ydp1, Matrice **Pp1, Matrice **Xdp2, Matrice **Ydp2, Matrice **Pp2) {

    CalculPinceAlone(gfd, Xd1, Yd1, Xd2, Yd2, X1, Y1, Z1, P1, X2, Y2, Z2, P2,
            PosPinceBA1, AmpPinceBA1, PosPinceBF1, AmpPinceBF1, mode,
            Xdp1, Ydp1, Pp1);
    CalculPinceAlone(gfd, Xd2, Yd2, Xd1, Yd1, X2, Y2, Z2, P2, X1, Y1, Z1, P1,
            PosPinceBA2, AmpPinceBA2, PosPinceBF2, AmpPinceBF2, mode,
            Xdp2, Ydp2, Pp2);
}

Matrice* GetFunctionSrez(Matrice*X, Matrice *Y, double pa, double pf, int *ia) {
    int i1 = 0, i2 = 0;
    for (int i = 0; i < X -> GetLignes(); i++) {
        if (X->Element(i, 0) == pa) i1 = i;
        if (X->Element(i, 0) == pf) i2 = i;
    }
    *ia = i1;
    Matrice* res = new Matrice(i2-i1+1, 1);
    for (int i = i1; i <= i2; i++) res -> SetElement(i-i1, 0, Y->Element(i, 0));
    return res;
}

Matrice* GetFunctionSrezDeb(Matrice *X, Matrice *Y, double pa) {
    //printf ("\n GetFunctionSrezDeb() ");
    int i1 = 0, i = 0;
    //printf (" XL=%d YL=%d pa=%f", X->GetLignes(), Y->GetLignes(), pa);
    while (X->Element(i1, 0) != pa)  i1++; 
    //printf ("\n  result i1=%d resL=%d", i1, (Y->GetLignes()-i1));
    Matrice* res = new Matrice(Y->GetLignes()-i1, 1);
    for (i = i1; i < Y->GetLignes(); i++) res -> SetElement(i-i1, 0, Y->Element(i, 0));
    //printf ("\n ...GetFunctionSrezDeb()");
    return res;
}

Pince* getGluePinceOnlyFuncs (Pince* pince1, double pos1, Pince* pince2, double pos2)
{
    Pince* res = new Pince();
    int i = 0;
    while (pince1 -> P1 -> Element(i, 0) != pos1) i++;
    int i1 = i;
    i = pince2->P1->GetLignes()-1;
    while (pince2 -> P1->Element(i, 0) != pos2) i--;
   
    int i2 = i;
    
    int n = i2 + i1 + 1;
    Matrice* func1 = new Matrice (n, 1);
    Matrice* func2 = new Matrice (n, 1);
    Matrice* P1 = new Matrice (n, 1);
    Matrice* P2 = new Matrice (n, 1);
    int isave = 0;

    for (i = i1; i > 0; i--) {
        func1 -> SetElement(isave,0, pince1->function1->Element(i, 0));
        func2 -> SetElement(isave,0, pince1->function2->Element(i, 0));
        P1 -> SetElement(isave, 0, pince1->P1->Element(i, 0));
        P2 -> SetElement(isave, 0, pince1->P1->Element(i, 0));
        isave++;
    }
    for (i = 0; i <= i2; i++) {
        func1 -> SetElement(isave,0, pince2->function1->Element(i, 0));
        func2 -> SetElement(isave,0, pince2->function2->Element(i, 0));
        P1 -> SetElement(isave, 0, pince2->P1->Element(i, 0));
        P2 -> SetElement(isave, 0, pince2->P1->Element(i, 0));
        isave++;
    }
    res -> function1 = func1;
    res -> function2 = func2;
    res -> P1 = P1;
    res -> P2 = P2;
    res -> i01 = i1;
    res -> i02 = i1;
    res -> ipa1 = pince2 -> ipa1 + i1;
    res -> ipa2 = pince2 -> ipa2 + i1;

    res -> ipf1 = pince2 -> ipf1 + i1;
    res -> ipf2 = pince2 -> ipf2 + i1;
    return res;
}

Pince* getPincePlus(WindPatternsProject* gfd,
                        int nerv1, int nerv2,
                        float deb1, float fin1, float deb2, float fin2, int faceDeb, int faceFin )
{
    //printf ("\ngetPincePlus()");
    Pince *res, *res1, *res2;
    if (faceDeb == 1 && faceFin == 1) 
        {
            // Ext -> Ext
            res = getPince(gfd, nerv1, nerv2, 1);
            res -> makeSrez(deb1, fin1, deb2, fin2);
        } else if (faceDeb == 1)
        {
            // Ext -> Int
            Pince* pext = getPince(gfd, nerv1, nerv2, 1);
            Pince* pint = getPince(gfd, nerv1, nerv2, 2);

            res1 = getGluePinceOnlyFuncs (pext, deb1, pint, fin1);
            res2 = getGluePinceOnlyFuncs (pext, deb2, pint, fin2);
            res = new Pince();
            res -> nerv1 = nerv1;
            res -> nerv2 = nerv2;
            res -> glue = true;

            res->function1 = res1->function1;
            res->function2 = res2->function2;
            res->P1 = res1->P1;
            res->P2 = res2->P2;
            res->ipa1 = res1->ipa1;
            res->ipa2 = res2->ipa2;
            res->ipf1 = res1->ipf1;
            res->ipf2 = res2->ipf2;
            res->i01 = res1->i01;
            res->i02 = res2->i02;
            res -> AmpA = pint -> AmpA;
            res -> AmpF = pint -> AmpF;
            res -> AmpA0 = pext -> AmpA;
       
        } else if (faceFin == 1)
        {
            // Int -> Ext
            Pince* pext = getPince(gfd, nerv1, nerv2, 1);
            Pince* pint = getPince(gfd, nerv1, nerv2, 2);
            //res = getGluePinceOnlyFuncs (pint, deb, pext, fin);

            res1 = getGluePinceOnlyFuncs (pint, deb1, pext, fin1);
            res2 = getGluePinceOnlyFuncs (pint, deb2, pext, fin2);

        	res = new Pince();

            res -> nerv1 = nerv1;
            res -> nerv2 = nerv2;
            res -> glue = true;

            res->function1 = res1->function1;
            res->function2 = res2->function2;
            res->P1 = res1->P1;
            res->P2 = res2->P2;
            res->ipa1 = res1->ipa1;
            res->ipa2 = res2->ipa2;
            res->ipf1 = res1->ipf1;
            res->ipf2 = res2->ipf2;
            res->i01 = res1->i01;
            res->i02 = res2->i02;
            res -> AmpA = pint -> AmpA;
            res -> AmpF = pint -> AmpF;
            res -> AmpA0 = pext -> AmpA;
        } else 
        {
            // Int -> Int
            res = getPince(gfd, nerv1, nerv2, 2);
            res->makeSrez(deb1, fin1, deb2, fin2);
        }

        //printf ("\n...getPincePlusNew()");
        return res;
}

Pince* getPince(WindPatternsProject* gfd, int nerv1, int nerv2, int face)
{
    //printf ("\ngetPince(%d, %d, %d)", nerv1, nerv2, face);
    Matrice *ResFunction1, *ResFunction2, *RadiusFunction, *FuncFunction1, *FuncFunction2;
    double ampFfA=0.0f, ampFfF=0.0f, ampRfA=0.0f, ampRfF=0.0f;
    int ipa = 0, ipf = 0;

    Matrice *P;
    getPinceFunctions(gfd, nerv1, nerv2, face,
                        &ipa, &ipf, &P, &FuncFunction1,&FuncFunction2, &ampFfA, &ampFfF,
                        &RadiusFunction, &ampRfA, &ampRfF);
    printf ("\n   in getPince(): P->GetLignes()=%d", P->GetLignes());
	printf ("\n   in getPince(): RadiusFunction->GetLignes()=%d", RadiusFunction->GetLignes());
    int n = FuncFunction1->GetLignes();
    int _ipai=ipa, _ipfi=ipf;

    ResFunction1 = new Matrice(n, 1);
    ResFunction2 = new Matrice(n, 1);

    double AmpA, AmpF;
    bool goHvost = true;
    if (gfd->PinceNosRadio == 0) {
        AmpA = ampRfA;
        if (gfd->PinceNosEqualAmp) {
            AmpF = AmpA;
            goHvost = false;
        }
    } else AmpA = ampFfA;

    if (goHvost){
        if (gfd->PinceHvostRadio == 0) {
            AmpF = ampRfF;
            if (gfd->PinceHvostEqualAmp)
                AmpA = AmpF;
        } else AmpF = ampFfF;
    }
 

    for (int i = 0; i < n; i++) {
        if (i <= _ipai) {
            // i in [0, ipai]
            if (gfd->PinceNosRadio==0) {
                ResFunction1->SetElement(i,0, RadiusFunction->Element(i,0));
                ResFunction2->SetElement(i,0, RadiusFunction->Element(i,0));
            } else {
                ResFunction1->SetElement(i,0, FuncFunction1->Element(i,0));
                ResFunction2->SetElement(i,0, FuncFunction2->Element(i,0));
            }
        } else {
            if (i < _ipfi) {
                // i in [ipai+1, ipfi-1]
                ResFunction1->SetElement(i,0, 0.0f);//mpipai + (i - _ipai)*(mpipfi-mpipai)/(_ipfi - _ipai));
                ResFunction2->SetElement(i,0, 0.0f);//mpipai + (i - _ipai)*(mpipfi-mpipai)/(_ipfi - _ipai));
            } else {
                // i in [ipfi,n-1]
                if (gfd->PinceHvostRadio==0) {
                    ResFunction1->SetElement(i,0,RadiusFunction->Element(i,0));
                    ResFunction2->SetElement(i,0,RadiusFunction->Element(i,0));
                } else {
                    ResFunction1->SetElement(i,0,FuncFunction1->Element(i,0));
                    ResFunction2->SetElement(i,0,FuncFunction2->Element(i,0));
                }

            }
        }

    }
    double PosA = gfd->PosPinceBA[0];
    double PosF = gfd->PosPinceBF[0];

    Pince* res = new Pince(nerv1, nerv2);
    res -> SetFunction1 (ResFunction1);
    res -> SetFunction2 (ResFunction2);
    //delete ResFunction1;
    //delete ResFunction2;
    res -> SetP1 (CloneMat(P));
    res -> SetP2 (CloneMat(P));
    delete (P);
    res -> SetAmps(AmpA, _ipai, AmpF, _ipfi);
    res -> SetPos (PosA, PosF);
    res -> debug = false;
    return res;
}


void Pince::makeSrez(double deb1, double fin1, double deb2, double fin2)
{
    if (debug) printf ("\n makeSrez()");
    Matrice* func1 = this->function1;
    Matrice* func2 = this->function2;
    Matrice* Pn1 = this -> P1;
    Matrice* Pn2 = this -> P2;
    int ia1 = 0;
    int ia2 = 0;
    //printf ("\nbefore P1->GetLignes");
    //printf ("\nP1->GetLignes()=%d", P1->GetLignes());
    //printf ("\nP2->GetLignes()=%d", P2->GetLignes());
    //printf ("\nafter P1->GetLignes");

    function1 = GetFunctionSrez(P1, func1, deb1, fin1, &ia1);
    //printf ("\nF1 ipa1=%d ipf1=%d ia1=%d function1->GetLignes()=%d",ipa1, ipf1, ia1, function1->GetLignes());
    function2 = GetFunctionSrez(P2, func2, deb2, fin2, &ia2);
    //printf ("\nF2 ipa2=%d ipf2=%d ia2=%d",ipa2, ipf2, ia2);
    P1 = GetFunctionSrez(Pn1, Pn1, deb1, fin1, &ia1);
    P2 = GetFunctionSrez(Pn2, Pn2, deb2, fin2, &ia2);
    ipa1 = ipa1 - ia1;
    ipf1 = ipf1 - ia1;
    ipa2 = ipa2 - ia2;
    ipf2 = ipf2 - ia2;
    delete (func1);
    delete (func2);
    delete (Pn1);
    delete (Pn2);
    if (debug) printf ("\n .... makeSrez()");
}

void getPinceFunctions(WindPatternsProject* gfd, int nerv1, int nerv2, int face, int* ipa, int* ipf, Matrice** Pout,
                        Matrice** FuncFunction1,Matrice** FuncFunction2, double* ampFfA, double* ampFfF,
                        Matrice** RadiusFunction, double* ampRfA, double* ampRfF)
{
    printf ("\n getPinceFunctions()");
    int nerv[2] = {nerv1, nerv2};
    int i = 0;
    Matrice *X[2], *Y[2], *Z[2], *P[2], *Xd[2], *Yd[2];
    CalculPatron(gfd, nerv1, false, face, face, 0.0f, 100.0f,
             nerv2, false, face, face, 0.0f, 100.0f,
             &Xd[0], &Yd[0], &Xd[1], &Yd[1],
             &X[0], &Y[0], &Z[0], &P[0],
             &X[1], &Y[1], &Z[1], &P[1]);
    //printf ( " P[0]->GetLignes()=%d", P[0]->GetLignes());
/*
    printf ("\ngetPF Xd[0]");
    Xd[0]->print(0);
    printf ("\ngetPF Yd[0]");
    Yd[0]->print(0);
    printf ("\ngetPF Xd[1]");
    Xd[1]->print(0);
    printf ("\ngetPF Yd[1]");
    Yd[1]->print(0);
  */
    Matrice *newXd[2], *newYd[2], *newP[2];
    int ipai[2] = {0, 0};
    int ipfi[2] = {0, 0};
    double PosA = gfd->PosPinceBA[0];
    double PosF = gfd->PosPinceBF[0];

    makePointPince(Xd[0], Yd[0], P[0], PosA, PosF, &newXd[0], &newYd[0], &newP[0], &ipai[0], &ipfi[0]);
    makePointPince(Xd[1], Yd[1], P[1], PosA, PosF, &newXd[1], &newYd[1], &newP[1], &ipai[1], &ipfi[1]);

    (*ipa)=ipai[0];
    (*ipf)=ipfi[0];

    int n = newXd[0]->GetLignes();

    (*FuncFunction1) = new Matrice (n, 1);
    (*FuncFunction2) = new Matrice (n, 1);
    (*RadiusFunction) = new Matrice (n, 1);

    // ==============FuncFunction=============
    double width=0.0f, argument=0.0f;
    CalculWidth(Xd[0], Yd[0], Xd[1], Yd[1], &width);
    (*ampFfA)=gfd->AmpPinceBA[0] * width / 100.0f;
    (*ampFfF)=gfd->AmpPinceBF[0] * width / 100.0f;

    for (i = 0; i < n; i++) {
        if (i <= ipai[0]) {
	     argument = (newXd[0]->Element(i, 0) - newXd[0]->Element(0, 0))  /  (newXd[0]->Element(ipai[0], 0) - newXd[0]->Element(0, 0));
            (*FuncFunction1) -> SetElement (i, 0, pinceFunctionA(gfd,argument));
        } else {
            if (i >= ipfi[0]) {
                argument = (newXd[0]->Element(newXd[0]->GetLignes() - 1, 0) - newXd[0]->Element(i, 0))  /   (newXd[0]->Element(newXd[0]->GetLignes() - 1, 0) - newXd[0]->Element(ipfi[0], 0));
                (*FuncFunction1) -> SetElement (i, 0, pinceFunctionF(gfd,argument));
            } else {
                (*FuncFunction1) -> SetElement (i, 0, 0.0f);
            }
        }
    }

    for (i = 0; i < n; i++) {
        if (i <= ipai[1]) {
	     argument = (newXd[1]->Element(i, 0) - newXd[1]->Element(0, 0))
                            / (newXd[1]->Element(ipai[1], 0) - newXd[1]->Element(0, 0));
            (*FuncFunction2) -> SetElement (i, 0, pinceFunctionA(gfd,argument));
        } else {
            if (i >= ipfi[1]) {
                argument = (newXd[1]->Element(newXd[1]->GetLignes() - 1, 0) - newXd[1]->Element(i, 0))
                                /   (newXd[1]->Element(newXd[1]->GetLignes() - 1, 0) - newXd[1]->Element(ipfi[1], 0));
                (*FuncFunction2) -> SetElement (i, 0, pinceFunctionF(gfd,argument));
            } else {
                (*FuncFunction2) -> SetElement (i, 0, 0.0f);
            }
        }

    }

    // ==============RadiusFunction=============
	    printf ("\n RadiusFunction()");
    //double kn = gfd->PinceRadiusAlgKNos;
    //double kh = gfd->PinceRadiusAlgKHvost;

    double w10a = dist2d(newXd[0] -> Element(ipai[0], 0), newYd[0] -> Element (ipai[0], 0),
                        newXd[1] -> Element(ipai[1], 0), newYd[1] -> Element(ipai[1], 0));

    double w10f = dist2d(newXd[0] -> Element(ipfi[0], 0), newYd[0] -> Element (ipfi[0], 0),
                        newXd[1] -> Element(ipfi[1], 0), newYd[1] -> Element(ipfi[1], 0));
	// oldSchool Radius Alg
    //(*ampRfA) = (w10a/kn - w10a)*0.5f;
    //(*ampRfF) = (w10f/kh - w10f)*0.5f;

    //double h0a = CalculHbyw0w1 (w10a/kn, w10a);
    //double h0f = CalculHbyw0w1 (w10f/kn, w10f);

    Matrice *XProfil0, *YProfil0;
    printf ("\n before getMiddleProfile()");
    GetMiddleProfile(gfd, gfd->Forme, nerv[0], nerv[1], face, 1,  &XProfil0, &YProfil0);
    Matrice *XProfil, *YProfil, *newPtmp;
    int i1 = 0, i2 = 0;
    //makePointPince(XProfil0, YProfil0, XProfil0, PosA, PosF, &XProfil, &YProfil, &newPtmp, &i1, &i2);
    XProfil = CloneMat(XProfil0);
    for (i = 0; i < XProfil->GetLignes();i++) {
        if (XProfil->Element(i, 0) == PosA) i1 = i;
        if (XProfil->Element(i, 0) == (100.0f-PosF)) i2 = i;
    }
    //printf ("\n i1=%d i2=%d", i1, i2);
    YProfil = CloneMat(YProfil0);


    //double ypa = YProfil->Element(i1,0);
    //double kynos = (ypa+h0a)/ypa;
    //double ypf = YProfil->Element(i2,0);
    //double kyhvost = (ypf+h0f)/ypf;

    Matrice *XProfil2n, *YProfil2n;
    //XProfil2n = CloneMat(XProfil);
    //YProfil2n = MultReelMat( CloneMat(YProfil), kynos );
	//printf ("\n$$1");
	GetMiddleProfileBal(gfd, gfd->Forme, nerv[0], nerv[1], face, 1,  &XProfil2n, &YProfil2n);
	//printf ("\n$$2");
    (*Pout)=XProfil2n;
    Matrice *XProfil2h, *YProfil2h;
    //  XProfil2h = CloneMat(XProfil);
    //  YProfil2h = MultReelMat( CloneMat(YProfil), kyhvost );
	//  printf ("\n$$$1");
	GetMiddleProfileBal(gfd, gfd->Forme, nerv[0], nerv[1], face, 1, &XProfil2h, &YProfil2h);
	//  printf ("\n$$$2");

	//	printf ("\n CalculNormalRasst");
    Matrice* Hn = CalculNormalRasst(XProfil, YProfil, XProfil2n, YProfil2n);
	//	printf ("\n Hn:");
	//  Hn->print(0);
	//  printf ("\n ...CalculNormalRasst");
    Matrice* Hh = CalculNormalRasst(XProfil, YProfil, XProfil2h, YProfil2h);
	//	printf ("\n Hh:");
	//  Hh->print(0);

	/*if (nerv1 == -1) {
	XProfil->print(0);
		printf ("\n...");
		XProfil2n->print(0);
	}*/
	printf ("\n XProfil->GL=%d XProfil2n->GL=%d", XProfil->GetLignes(), XProfil2n->GetLignes());
	printf ("\n newXd[0]->GetLignes()=%d", newXd[0]->GetLignes());

    double w = dist2d(newXd[0]->Element(i1,0),newYd[0]->Element(i1,0),  newXd[1]->Element(i1,0), newYd[1]->Element(i1,0));
    double dkn = 1.0f/(CalculW0byHW1(Hn->Element(i1, 0), w)-w);
    w = dist2d(newXd[0]->Element(i2,0),newYd[0]->Element(i2,0),  newXd[1]->Element(i2,0), newYd[1]->Element(i2,0));
    double dkh = 1.0f/(CalculW0byHW1(Hh->Element(i2, 0), w)-w);

    double wi = 0.0f;

    for (i=0; i < newXd[0]->GetLignes(); i++) {
        wi = dist2d(newXd[0]->Element(i,0), newYd[0]->Element(i,0),  newXd[1]->Element(i,0), newYd[1]->Element(i,0));
        if ( i <= i1) {
            (*RadiusFunction) -> SetElement(i,0,
                    (CalculW0byHW1(Hn->Element(i, 0), wi)-wi) * dkn);
        } else {
            if ( i >= i2){
            //    printf (" hh=%f", Hh->Element(i, 0));
                (*RadiusFunction) -> SetElement(i,0,
                        (CalculW0byHW1(Hh->Element(i, 0), wi)-wi) * dkh);
            } else (*RadiusFunction) -> SetElement(i,0, 0.0f);
        }
       // printf ("\n RF=%f", (*RadiusFunction) -> Element(i, 0));
    }
	// new School!
	printf ("\n Hn->Element(%d, 0)=%f w10a=%f", ipai[0], Hn->Element(ipai[0], 0), w10a);
	printf ("\n Hn->Element(%d, 0)=%f w10f=%f", ipfi[0], Hn->Element(ipfi[0], 0), w10f);
	printf ("\n ampRfA=%f ampRfF=%f", (CalculW0byHW1(Hn->Element(ipai[0], 0), w10a) - w10a)*0.5f, (CalculW0byHW1(Hn->Element(ipfi[0], 0), w10f) - w10f)*0.5f);

	printf ("\n RF->GetLignes()=%d", (*RadiusFunction)->GetLignes());

	(*ampRfA) = (CalculW0byHW1(Hn->Element(ipai[0], 0), w10a) - w10a)*0.5f;
    (*ampRfF) = (CalculW0byHW1(Hn->Element(ipfi[0], 0), w10f) - w10f)*0.5f;

	 printf ("\n ...RadiusFunction()");
     printf ("\n ...getPinceFunctions()");
}

void getPinceRadiusFunction(WindPatternsProject* gfd, int nerv1, int nerv2, int face) {

    //printf ("\n goCalcPinceRadiusFunction(%d, %d, %d)", nerv1, nerv2, face);
    double k = gfd->PinceRadiusAlgKNos;
    double k2 = gfd->PinceRadiusAlgKHvost;
    //Matrice* PinceFunction;
    int nerv[2] = {nerv1, nerv2};
    double *funcPince;
    Matrice *X[2], *Y[2], *Z[2], *P[2], *Xd[2], *Yd[2];

    CalculPatron(gfd,nerv[0], false, face, face, 0.0f, 100.0f,
             nerv[1], false, face, face, 0.0f, 100.0f,
             &Xd[0], &Yd[0], &Xd[1], &Yd[1],
             &X[0], &Y[0], &Z[0], &P[0],
             &X[1], &Y[1], &Z[1], &P[1]);
    Matrice *newXd[2], *newYd[2], *newP[2];

    int ipai[2] = {0, 0};
    int ipfi[2] = {0, 0};
    double PosA = gfd->PosPinceBA[0];
    double PosF = gfd->PosPinceBF[0];
    makePointPince(Xd[0], Yd[0], P[0], PosA, PosF, &newXd[0], &newYd[0], &newP[0], &ipai[0], &ipfi[0]);
    makePointPince(Xd[1], Yd[1], P[1], PosA, PosF, &newXd[1], &newYd[1], &newP[1], &ipai[1], &ipfi[1]);

    double ampPinceAabs = gfd->AmpPinceBA[0];
    double ampPinceFabs = gfd->AmpPinceBF[0];

    funcPince = new double[Xd[0]->GetLignes()];

    double Xr[2], Yr[2];
    for (int i = 0; i < 2; i++) {
            getPointByPos(Xd[i], Yd[i], P[i], PosA, &Xr[i], &Yr[i]);
    }

    /*for (i = 0; i < P[0]->GetLignes(); i++) {
        printf ("\nXYP %d (%f, %f) %f", i, Xd[0]->Element(i,0),Yd[0]->Element(i,0), P[0]->Element(i,0));
    }*/

    //printf ("\n\n Xr[0]=%f Yr[0]=%f Xr[1]=%f Yr[1]=%f ", Xr[0], Yr[0], Xr[1], Yr[1]);
    double w10 = dist2d(Xr[0], Yr[0], Xr[1], Yr[1]);
    //printf ("\n\nW10=%f", w10);

    double funcPosPince = 1/k;
    double funcPosPinceAmp = w10/k - w10;
    //printf ("\nBEG -> %f (%f)", funcPosPince, (funcPosPince-1)*w10);

    double h0 = CalculHbyw0w1 (w10/k, w10);

    Matrice *XProfil1, *YProfil1;
    GetMiddleProfile(gfd, gfd->Forme,nerv[0], nerv[1], face, 1, &XProfil1, &YProfil1);
    /*for (i = 0; i < XProfil1->GetLignes() ; i++) {
        printf ("\n MiddleProfil: %d (%f, %f)", i, XProfil1->Element(i,0), YProfil1->Element(i,0));
    }*/


    double xp, yp;
    getPointByPos(XProfil1, YProfil1, XProfil1, PosA, &xp, &yp);

    double kymorda=(yp+h0)/yp;
    //printf ("\nxp=%f yp %f kymorda=%f", xp, yp, kymorda);
    Matrice *XProfil2, *YProfil2;
    //XProfil2 = CloneMat(XProfil1);
    //YProfil2 = MultReelMat( CloneMat(YProfil1), kymorda );
	printf ("\n$1");
	GetMiddleProfileBal(gfd, gfd->Forme,nerv[0], nerv[1], face, 1, &XProfil2, &YProfil2);
	printf ("\n$2");
	int _ipai = 0;

    for (int i = 0; i < XProfil1->GetLignes(); i++) {
		if (XProfil1->Element(i,0) > xp) {
			_ipai = i - 1;
			break;
		}
    }

    Matrice* H = CalculNormalRasst(XProfil1, YProfil1, XProfil2, YProfil2);
    /*for (i = 0; i < XProfil1->GetLignes() ; i++) {
        printf ("\n H(Profil): %d (%f)", i, H->Element(i,0));
    }*/

    //printf ("\nipai = %d", ipai);
    for (int i = _ipai; i >= 0; i--) {
        double hi = H->Element(i, 0);
        // need to calculate on normal!!!!
        double wi = dist2d(Xd[0]->Element(i,0),Yd[0]->Element(i,0), Xd[1]->Element(i,0), Yd[1]->Element(i,0));
        double funcPosPincei = CalculW0byHW1(hi, wi);
    }


}


void goCalcPinceLen(WindPatternsProject* gfd, int noNerv, int face, double *len) {
    //double mode1 = modePinces;
    Matrice * Xd1[2], *Yd1[2], *Xd1p[2], *Yd1p[2], *long1;
    //Matrice * Xd2[2], *Yd2[2];
    Matrice * X1[2], *Y1[2], *Z1[2], *P1[2], *newP1[2];
    CalculPatron(gfd, noNerv, false, face, face, 0.0f, 100.0f,
            noNerv + 1, false, face, face, 0.0f, 100.0f,
            &Xd1[0], &Yd1[0], &Xd1[1], &Yd1[1],
            &X1[0], &Y1[0], &Z1[0], &P1[0],
            &X1[1], &Y1[1], &Z1[1], &P1[1]);
    Pince* pince = getPince(gfd, noNerv, noNerv+1, face);
    /*CalculPinceAlone(getWindPatternsProject(),Xd1[0], Yd1[0], Xd1[1], Yd1[1],
            X1[0], Y1[0], Z1[0], P1[0],
            X1[1], Y1[1], Z1[1], P1[1],
            PosPinceBA[0], AmpPinceBA[0], PosPinceBF[0], AmpPinceBF[0], mode1,
            &Xd1p[0], &Yd1p[0], &newP1[1]);*/
    CalculPinceAloneAbsFunc(Xd1[0], Yd1[0], Xd1[1], Yd1[1],
            X1[0], Y1[0], Z1[0], P1[0],
            X1[1], Y1[1], Z1[1], P1[1], pince, pince->AmpA, pince->AmpF, 0,
            &Xd1p[0], &Yd1p[0], &newP1[0]);
    long1 = Longueur(Xd1p[0], Yd1p[0]);
    *len = long1->Element(long1->GetLignes() - 1, 0);
    delete (X1[0]);
    delete (Y1[0]);
    delete (Z1[0]);
    delete (P1[0]);

    delete (X1[1]);
    delete (Y1[1]);
    delete (Z1[1]);
    delete (P1[1]);

    delete (Xd1[0]);
    delete (Yd1[0]);
    delete (Xd1[1]);
    delete (Yd1[1]);

    delete (long1);
}
