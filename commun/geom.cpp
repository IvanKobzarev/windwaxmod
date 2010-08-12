/*************************************/
/* geom.c                            */
/*************************************/

#pragma warning(disable:4514)

#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>

#include "plot.h"
#include "rasklad.h"
#include "matrice.h"
#include "geom.h"
#include "fichier.h"
#include "profil.h"
#include "patternsproject.h"

#define sqr(l1) ((l1)*(l1))
#define pi	3.141592675f

#ifndef DEG2RAD
#define DEG2RAD	(3.141592675f/180.0f)
#endif

#ifndef DEBUG
#define DEBUG false
#endif


/********************************************************
* res=MonBezier(tab,NbrPts,tabi)
* calcul Bezier de NbrPts selon la technique geometrique
* tab contient les 4 points ABCD selon le format
* [XA YA
*  XB YB
*  XC YC
*  XD YD]
********************************************************/

Matrice* MonBezier(Matrice *tab, int NbrPts)
{
	Matrice *res = NULL;
	double XA,YA,XB,YB,XC,YC,XD,YD,XE,YE,XF,YF,XG,YG,XH,YH,XI,YI,XJ,YJ;
	double t,pas;
	int i;
	/*allocation matrice resultat*/
	res=new Matrice(NbrPts,2);
	/*pour simpifier l'ecriture*/
	XA=tab->Element(0,0); YA=tab->Element(0,1);
	XB=tab->Element(1,0); YB=tab->Element(1,1);
	XC=tab->Element(2,0); YC=tab->Element(2,1);
	XD=tab->Element(3,0); YD=tab->Element(3,1);
	/*boucle remplissage*/
	pas=1.0F/(NbrPts-1.0F); t=0.0F;
	for (i=0; i<NbrPts; i++)
	{ 
		XE=XA+(XB-XA)*t;
		YE=YA+(YB-YA)*t;
		XF=XB+(XC-XB)*t;
		YF=YB+(YC-YB)*t;
		XG=XC+(XD-XC)*t;
		YG=YC+(YD-YC)*t;
		XH=XE+(XF-XE)*t;
		YH=YE+(YF-YE)*t;
		XI=XF+(XG-XF)*t;
		YI=YF+(YG-YF)*t;
		XJ=XH+(XI-XH)*t;
		YJ=YH+(YI-YH)*t;
		res->SetElement(i,0,XJ);
		res->SetElement(i,1,YJ);
		t+=pas;
	}
	return res;
}



/*****************************************************/
/* CalculVecteurNormal                               */
/* calcul vecteur BC normal au vecteur AB            */
/* direct (cotï¿½=+1) ou indirect (cotï¿½=-1)            */
/* de longeur l.                                     */
/*****************************************************/

void CalculVecteurNormal(double xa,double ya,double xb,double yb,
						 double *xc,double *yc,double l,int cote)
{
	double xab, yab, a11, a12, a21, a22, b1, b2; /*var intermï¿½daires*/
	double delta, delta1, delta2; /*determinant principal et secondaires*/
	/*pour simplifier ecriture*/
	xab=xb-xa; yab=yb-ya;
	/*preparation matrices A et B pour resolution AX=B*/
	a11=-yab; a12=xab;
	a21=xab;  a22=yab;
	b1=l*(double)sqrt(xab*xab+yab*yab)*cote +yb*xab -xb*yab;
	b2=xb*xab +yb*yab;
	/*calcul determinants*/
	delta=a11*a22-a21*a12;
	if (delta !=0)
	{
		delta1=b1*a22-b2*a12;
		delta2=a11*b2-a21*b1;
        	/*calcul resultat*/
		*xc=delta1/delta;
		*yc=delta2/delta;
	}
	else
	{
		printf("\nIMPSBL CalculVecteurNormal (%f, %f)(%f, %f)", xa, ya, xb, yb);
		*xc=0.0; *yc=0.0;
	}

}



/*****************************************************/
/* CalculContour                                     */
/* xy vecteurs colonnes ...						     */
/* res: matrice contenant le contour                 */
/* d : distance a la courbe en chaque point          */
/* cote: +1:contour gauche, -1:contour droite        */
/*****************************************************/

Matrice* CalculContour(Matrice* xy, Matrice *d, int cote)
{
	int i; 
	/*test parametres ok*/
	/* a faire ...*/
	/*creation matrices*/
	Matrice *res = new Matrice(xy->GetLignes(),2);
	/*calcul matrices contour xc,yc a partir du 2eme point*/
	double xc, yc;
	for(i=0; i<xy->GetLignes()-1; i++)
	{
		xc = res->Element(i+1,0);
		yc = res->Element(i+1,1);
                //printf ("\n CalculContour1()");
		CalculVecteurNormal(xy->Element(i,0), xy->Element(i,1), xy->Element(i+1,0), xy->Element(i+1,1),	&xc, &yc, d->Element(i,0), cote);
		res->SetElement(i+1,0,xc);
		res->SetElement(i+1,1,yc);
	}

	/*calcul premier point xc,yc*/
	xc = res->Element(0,0);
	yc = res->Element(0,1);
        //printf ("\n CalculContour2()");
	CalculVecteurNormal(xy->Element(1,0), xy->Element(1,1), xy->Element(0,0), xy->Element(0,1), &xc, &yc, d->Element(0,0), -cote);
	res->SetElement(0,0,xc);
	res->SetElement(0,1,yc);
	return res;

}

void AddPointsToCourb(Matrice* courb, Matrice* pts, Matrice** newCourb)
{
    //printf ("\n add Points to courb()");
    int n = courb->GetLignes();
    int np = pts->GetLignes();
    //printf ("\n n=%d np=%d",n,np);
    int q = np;
    int i, j ,s;
    for (i = 0; i < n; i++ ) 
        for (j = 0; j < np; j++) 
            if (courb->Element(i,0) == pts->Element(j, 0)) q--;
      //  printf ("\nq=%d",q);
    (*newCourb) = new Matrice(n + q, 2);
    i = 0;
    j = 0;
    s = 0;
    double xc, yc, xp;
    while ((i + j - s) < (n+q)) {

        if (i < n) {
            xc = courb -> Element(i, 0);
            yc = courb -> Element(i, 1);
        } else {
            xc = courb -> Element(n-1, 0);
            yc = courb -> Element(n-1, 1);
        }

        if (j < np) xp = pts -> Element(j, 0); else xp = pts -> Element(np-1, 0);
        //printf ("\n (%d, %d, %d)    (%f) (%f)", i, j, s, xc, xp);

        if (  (i < n) && ((xc < xp) || ((j == np)))) {
          //  printf ("\n < i");
            (*newCourb) -> SetElement(i + j - s, 0, xc);
            (*newCourb) -> SetElement(i + j - s, 1, yc);
            i++;
        } 
        if (  (j < np)  && ((xc > xp) || (i == n)) ) {
            //printf ("\n > j");
            (*newCourb) -> SetElement(i + j - s, 0, xp);
            (*newCourb) -> SetElement(i + j - s, 1, InterpLinX(courb, xp));
            j++;
         }
         if (xc==xp) {
              //  printf ("\n =");
                 (*newCourb) -> SetElement(i + j - s, 0, xc);
                 (*newCourb) -> SetElement(i + j - s, 1, yc);
                 i++;
                 j++;
                 s++;
         }
       

        
    }
    //printf ("\n ...add Points to courb()");
}

Matrice* CalculNormalRasst(Matrice* X1, Matrice* Y1,Matrice* X2, Matrice* Y2) {
    int n = X1->GetLignes();
    Matrice *res = new Matrice(n, 1);
    int i = 0;
    res->SetElement(0, 0,dist2d (X1->Element(0,0), Y1->Element(0,0), X2->Element(0,0), Y2->Element(0,0)));
    res->SetElement(n-1, 0,dist2d (X1->Element(n-1,0), Y1->Element(n-1,0), X2->Element(n-1,0), Y2->Element(n-1,0)));
    double xn, yn, x , y;
    for (i=1; i<n-1;i++) {
        double x0 = X1->Element(i, 0);
        double y0 = Y1->Element(i, 0);
        CalculVecteurBissec(X1->Element(i-1, 0), Y1->Element(i-1, 0),
                        x0, y0,
                        X1->Element(i+1, 0), Y1->Element(i+1, 0),
		         &xn, &yn, 10.0f, +1);
        
        for (int j = 0; j < X2->GetLignes()-1; j++) {
            double x1 = X2->Element(j,0);
            double y1 = Y2->Element(j,0);
            double x2 = X2->Element(j+1,0);
            double y2 = Y2->Element(j+1,0);
            Inter2Vecteurs(x0,y0,xn,yn, x1,y1,x2,y2,&x, &y);
            if ((x >= x1) && (x <= x2)) {
                res->SetElement(i, 0, dist2d (x,y,x0,y0));
                break;
            }
        }
    }
    return res;
}

/*****************************************************/
/* CalculVecteurBissec                               */
/*****************************************************/

void CalculVecteurBissec(double x1,double y1,double x2,double y2, double x3,double y3,
						 double *x,double *y,double l,int cote)
{
	double x4,y4,x5,y5;
	double x12,y12,x23,y23;
	double a11,a12,a21,a22,b1,b2;
	double delta, delta1, delta2;
	//pour simplifier calcul
	x12=x2-x1; y12=y1-y2;
	x23=x3-x2; y23=y3-y2;
	//test pts 1,2,3 alignï¿½s
	if (x12*y23-x23*y12 == 0.0f)
	{
//printf ("\nCalculVecteurBissec cote1");
		CalculVecteurNormal(x1,y1,x2,y2,x,y,l,cote);
	}
	else
	{
//printf ("\nCalculVecteurBissec cote2");
		CalculVecteurNormal(x1,y1,x2,y2,&x4,&y4,l,cote);
//printf ("\nCalculVecteurBissec -cote");
		CalculVecteurNormal(x3,y3,x2,y2,&x5,&y5,l,-cote);
		/*preparation matrices A et B pour resolution AX=B*/
		a11=y12; a12=-x12;
		a21=y23;  a22=-x23;
		b1=x4*y12-y4*x12;
		b2=x5*y23-y5*x23;
		/*calcul determinants*/
		delta=a11*a22-a21*a12;
		if (delta !=0)
		{
			delta1=b1*a22-b2*a12;
			delta2=a11*b2-a21*b1;
        		/*calcul resultat*/
			*x=delta1/delta;
			*y=delta2/delta;
		}
		else
		{
			printf("\npb: impossibilitï¿½ de calcul dans la procï¿½dure 'CalculVecteurBissec'");
			*x=0.0; *y=0.0;
		}
	}
}

/*********************/
/* EpaisseurRelative */
/*********************/

double EpaisseurRelative(Matrice* extrados, Matrice* intrados)

{
	int i;
	double MaxExt=-1000.0, MinInt=+1000.0;
	/*calcul epaisseur relative -> a refaire car calcul approchï¿½*/
	for(i=0; i<extrados->GetLignes(); i++)
		if(MaxExt < extrados->Element(i,1)) MaxExt=extrados->Element(i,1);
	for(i=0; i<intrados->GetLignes(); i++)
		if(MinInt > intrados->Element(i,1)) MinInt=intrados->Element(i,1);
	return MaxExt-MinInt;
}



/*****************/
/* CalculForme3D */
/*****************/
void CalculForme3D(Forme *forme, int isPercent, double percent,
				   Matrice *ExtProfCent, Matrice *IntProfCent,
				   Matrice *ExtProfBout, Matrice *IntProfBout,
				   Matrice **XExt, Matrice **YExt, Matrice **ZExt,
				   Matrice **XInt, Matrice **YInt, Matrice **ZInt)

{
    //printf ("\n CalculForme3D");
    Matrice *ExtProfCentN, *ExtProfBoutN;
	double LongNerv, EpaiRel, xp,yp, xo,yo,zo, a,v,m;
	double EpaiRelProfCent, EpaiRelProfBout;
	double coeffx, coeffyCent, coeffyBout;
	int i,j;
    bool isCenterPanel = (1 & forme->NbCaiss);

	/*calcul epaisseur relative profil central et bout*/
	EpaiRelProfCent = EpaisseurRelative(ExtProfCent, IntProfCent);
	//printf ("\nEpaiRelProfCent=%f", EpaiRelProfCent);
	EpaiRelProfBout = EpaisseurRelative(ExtProfBout, IntProfBout);
	//printf ("\nEpaiRelProfBout=%f", EpaiRelProfBout);
	/*init matrices extrados*/
	//delete(*XExt); delete(*YExt); delete(*ZExt); 
	//*XExt = Zeros(5, 5);
	*XExt = Zeros(forme->m_nbProfils, ExtProfCent->GetLignes());
	*YExt = Zeros(forme->m_nbProfils, ExtProfCent->GetLignes());
	*ZExt = Zeros(forme->m_nbProfils, ExtProfCent->GetLignes());

	/*init matrices intrados*/
	//delete(*XInt); delete(*YInt); delete(*ZInt); 
	*XInt = Zeros(forme->m_nbProfils, IntProfCent->GetLignes());
	*YInt = Zeros(forme->m_nbProfils, IntProfCent->GetLignes());
	*ZInt = Zeros(forme->m_nbProfils, IntProfCent->GetLignes());
	/*boucle ï¿½ partir de la 1ere nervure du centre vers l'extrï¿½mitï¿½*/


    if (isPercent) {
        //printf ("\n isPercent");
        makePosProfile(ExtProfCent, IntProfCent, percent, &ExtProfCentN);
        makePosProfile(ExtProfBout, IntProfBout, percent, &ExtProfBoutN);
    }
        
	for (i=0; i<forme->m_nbProfils; i++)
	{
		//longueur nervure courante
		LongNerv = forme->m_pProfils[i]->m_fLength;
		//epaisseur relative
		EpaiRel = forme->m_pProfils[i]->m_fWidth;
		//position xo,yo,zo du nez
		xo = forme->m_pProfils[i]->m_fNezX;
		yo = forme->m_pProfils[i]->m_fNezY;
		zo = forme->m_pProfils[i]->m_fNezZ;
        //printf ("\n %d long=%f width=%f (%f, %f, %f)", i, forme->m_pProfils[i]->m_fLength, forme->m_pProfils[i]->m_fWidth, xo, yo, zo);
		//inclinaison de la nervure par rapport a l'horizontale
		a = forme->m_pProfils[i]->m_fInclin;
        if ((isCenterPanel) && (i == 0)) a = forme->m_pProfils[0]->m_fInclin - 0.33333333f * fabs(forme->m_pProfils[0]->m_fInclin - forme->m_pProfils[1]->m_fInclin);
                
		//angle vrillage
		v = forme->m_pProfils[i]->m_fWash;
		//coeff morphing
		m = forme->m_pProfils[i]->m_fMorph;
               
		//printf ("\n%3d -> LNerv=%f EpRel=%f a=%f v=%f m=%f",
		//				i,  LongNerv, EpaiRel, a * 180.0f/pi, v, m);

		//calcul coeffx et coeffy des points du profil en
		//fonction de l'ï¿½paisseur relative et de la longueur de nervure
		coeffx = LongNerv/100.0f;
		coeffyCent = LongNerv*EpaiRel/(EpaiRelProfCent*100.0f);
		coeffyBout = LongNerv*EpaiRel/(EpaiRelProfBout*100.0f);

        /*boucle sur les points du profil en intrados*/
		for (j=0; j<IntProfCent->GetLignes(); j++)
		{
			xp = IntProfCent->Element(j,0)*coeffx;
			yp = IntProfCent->Element(j,1)*coeffyCent*m
				+ IntProfBout->Element(j,1)*coeffyBout*(1.0f-m);
			(*XInt)->SetElement(i,j, xo+(yp-xp*(double)sin(v))*(double)cos(a));
			(*YInt)->SetElement(i,j, yo+(yp-xp*(double)sin(v))*(double)sin(a));
			(*ZInt)->SetElement(i,j, zo-(xp*(double)cos(v)));

		}

		/*boucle sur les points du profil en extrados*/
		for (j=0; j<ExtProfCent->GetLignes(); j++)
		{
            if (isPercent) {
                xp = ExtProfCentN->Element(j,0)*coeffx;
                yp = ExtProfCentN->Element(j,1)*coeffyCent*m
                        + ExtProfBoutN->Element(j,1)*coeffyBout*(1.0f-m);
            } else {
                xp = ExtProfCent->Element(j,0)*coeffx;
                yp = ExtProfCent->Element(j,1)*coeffyCent*m
                        + ExtProfBout->Element(j,1)*coeffyBout*(1.0f-m);
            }
			(*XExt)->SetElement(i,j, xo+(yp-xp*(double)sin(v))*(double)cos(a));
			(*YExt)->SetElement(i,j, yo+(yp-xp*(double)sin(v))*(double)sin(a));
			(*ZExt)->SetElement(i,j, zo-(xp*(double)cos(v)));
          /*  if (i == 0) {
                printf ("\n%d DELTA(%f, %f, %f)",j, (yp-xp*(double)sin(v))*(double)cos(a), (yp-xp*(double)sin(v))*(double)sin(a), -(xp*(double)cos(v)));
            } */
		}

	}
        //printf ("\n ...CalculForme3D");

}



void CalculForme3DBallonement
				(WindPatternsProject* gfd, Forme *forme, Ballonement* bal, int isPercent, double percent,
				   Matrice *ExtProfCent, Matrice *IntProfCent,
				   Matrice *ExtProfBout, Matrice *IntProfBout,
				   Matrice **XExt, Matrice **YExt, Matrice **ZExt,
				   Matrice **XInt, Matrice **YInt, Matrice **ZInt)

{
	// !!! support of isPercent, percent not implemented yet!!!

    //printf ("\n CalculForme3DBallonement");
    Matrice *ExtProfCentN, *ExtProfBoutN;
	double LongNerv, EpaiRel, xp, yp, xo, yo, zo, a, a0, v, m;
	double EpaiRelProfCent, EpaiRelProfBout;
	double coeffx, coeffy, coeffyCent, coeffyBout;
	int i,j;
    bool isCenterPanel = (1 & forme->NbCaiss);

	/*calcul epaisseur relative profil central et bout*/
	EpaiRelProfCent = EpaisseurRelative(ExtProfCent, IntProfCent);
	EpaiRelProfBout = EpaisseurRelative(ExtProfBout, IntProfBout);

	*XExt = Zeros(forme->m_nbProfils*2 - 1, ExtProfCent->GetLignes());
	*YExt = Zeros(forme->m_nbProfils*2 - 1, ExtProfCent->GetLignes());
	*ZExt = Zeros(forme->m_nbProfils*2 - 1, ExtProfCent->GetLignes());

	*XInt = Zeros(forme->m_nbProfils*2 - 1, IntProfCent->GetLignes());
	*YInt = Zeros(forme->m_nbProfils*2 - 1, IntProfCent->GetLignes());
	*ZInt = Zeros(forme->m_nbProfils*2 - 1, IntProfCent->GetLignes());

    if (isPercent) {
        makePosProfile(ExtProfCent, IntProfCent, percent, &ExtProfCentN);
        makePosProfile(ExtProfBout, IntProfBout, percent, &ExtProfBoutN);
    }
	//printf ("\n CalculForme3DBallonement go in FOR");        
	for (i=0; i<forme->m_nbProfils; i++)
	{
		// normal nervure profile
		//printf ("\n\n\n %d nervure normal", i);        
		LongNerv = forme->m_pProfils[i]->m_fLength;
		EpaiRel = forme->m_pProfils[i]->m_fWidth;
		xo = forme->m_pProfils[i]->m_fNezX;
		yo = forme->m_pProfils[i]->m_fNezY;
		zo = forme->m_pProfils[i]->m_fNezZ;
		a = forme->m_pProfils[i]->m_fInclin;
        if ((isCenterPanel) && (i == 0)) a = forme->m_pProfils[0]->m_fInclin - 0.33333333f * fabs(forme->m_pProfils[0]->m_fInclin - forme->m_pProfils[1]->m_fInclin);
              
		//angle vrillage
		v = forme->m_pProfils[i]->m_fWash;
		//coeff morphing
		m = forme->m_pProfils[i]->m_fMorph;
		coeffx = LongNerv/100.0f;
		coeffyCent = LongNerv*EpaiRel/(EpaiRelProfCent*100.0f);
		coeffyBout = LongNerv*EpaiRel/(EpaiRelProfBout*100.0f);
		//printf ("\n %d N longNerv=%f", i, LongNerv);
		//printf ("\n %d N coeffx=%f coeffy=(%f, %f)", i, coeffx, coeffyCent, coeffyBout);
		//printf ("\n %d N(init) EpaiRelProfCent=%f", i, EpaiRelProfCent);
		//printf ("\n %d N(init) EpaiRelProfBout=%f", i, EpaiRelProfBout);
		for (j=0; j<IntProfCent->GetLignes(); j++)
		{
			xp = IntProfCent->Element(j,0)*coeffx;
			yp = IntProfCent->Element(j,1)*coeffyCent*m
				+ IntProfBout->Element(j,1)*coeffyBout*(1.0f-m);
			(*XInt)->SetElement(2*i,j, xo+(yp-xp*(double)sin(v))*(double)cos(a));
			(*YInt)->SetElement(2*i,j, yo+(yp-xp*(double)sin(v))*(double)sin(a));
			(*ZInt)->SetElement(2*i,j, zo-(xp*(double)cos(v)));

		}
		for (j=0; j<ExtProfCent->GetLignes(); j++)
		{
            if (isPercent) {
                xp = ExtProfCentN->Element(j,0)*coeffx;
                yp = ExtProfCentN->Element(j,1)*coeffyCent*m
                        + ExtProfBoutN->Element(j,1)*coeffyBout*(1.0f-m);
            } else {
                xp = ExtProfCent->Element(j,0)*coeffx;
                yp = ExtProfCent->Element(j,1)*coeffyCent*m
                        + ExtProfBout->Element(j,1)*coeffyBout*(1.0f-m);
            }
			(*XExt)->SetElement(2*i,j, xo+(yp-xp*(double)sin(v))*(double)cos(a));
			(*YExt)->SetElement(2*i,j, yo+(yp-xp*(double)sin(v))*(double)sin(a));
			(*ZExt)->SetElement(2*i,j, zo-(xp*(double)cos(v)));
		}
		//printf ("\n len (%d) ext(%d) int(%d)",i, ExtProfCent->GetLignes(), IntProfCent->GetLignes());
		// --  normal nervure profile
		if ( i < (forme->m_nbProfils - 1)) {
			// nervure ballonement
			LongNerv = (forme->m_pProfils[i]->m_fLength + forme->m_pProfils[i+1]->m_fLength) * 0.5f;
			EpaiRel = (forme->m_pProfils[i]->m_fWidth + forme->m_pProfils[i+1]->m_fWidth) * 0.5f;
			xo = (forme->m_pProfils[i]->m_fNezX + forme->m_pProfils[i+1]->m_fNezX) * 0.5f;
			yo = (forme->m_pProfils[i]->m_fNezY + forme->m_pProfils[i+1]->m_fNezY) * 0.5f;
			zo = (forme->m_pProfils[i]->m_fNezZ + forme->m_pProfils[i+1]->m_fNezZ) * 0.5f;
			a0 = (forme->m_pProfils[i]->m_fInclin + forme->m_pProfils[i+1]->m_fInclin) * 0.5f;
			if ((isCenterPanel) && (i == 0)) a0 = forme->m_pProfils[0]->m_fInclin - 0.33333333f * fabs(forme->m_pProfils[0]->m_fInclin - forme->m_pProfils[1]->m_fInclin);
			a = (a0 + forme->m_pProfils[i+1]->m_fInclin) * 0.5f;
			//angle vrillage
			v = (forme->m_pProfils[i]->m_fWash + forme->m_pProfils[i]->m_fWash) * 0.5f;
			//coeff morphing
			m = (forme->m_pProfils[i]->m_fMorph + forme->m_pProfils[i]->m_fMorph) * 0.5f;
			// time to calculate profile ballone
			ProfilGeom* pgCur = getProfile(gfd, forme, i);
			ProfilGeom* pgCurBal = getBalloneProfilGeom(pgCur, bal->kChord->Element(i, 0), bal->kMf->Element(i, 0), EpaiRel, bal->wN->Element(i, 0), bal->dyw->Element(i, 0));
			double l = abs (pgCur->ExtProf->Element(pgCur->ExtProf->GetLignes() - 1, 0) - pgCur->ExtProf->Element(0, 0));
			double xv = l * (100.0f - (gfd->PosPinceBF[0])) * 0.01f;
			ProfilGeom* pg = getProfilGeomTailDown(pgCurBal, pgCur, xv, bal->powerTail->Element(i, 0));
			//  -- time to calculate profile ballone
			coeffx = LongNerv/100.0f;
			double EpaiRelCur = EpaisseurRelative(pgCur->ExtProf, pgCur->IntProf);
			coeffy = LongNerv*EpaiRel/(EpaiRelCur*100.0f);
			//printf ("\n B longNerv=%f", LongNerv);
			//printf ("\n B coeffx=%f coeffy=(%f)", coeffx, coeffy);
			//printf ("\n B EpaiRelCur=%f", EpaiRelCur);
			double _x=0, _y=0, _z=0;
			for (j=0; j < pg->IntProf->GetLignes(); j++)
			{
				xp = pg->IntProf->Element(j, 0) * coeffx;
				yp = pg->IntProf->Element(j, 1) * coeffy;
				_x = xo+(yp-xp*(double)sin(v))*(double)cos(a);
				_y = yo+(yp-xp*(double)sin(v))*(double)sin(a);
				_z = zo-(xp*(double)cos(v));
				(*XInt)->SetElement(2*i+1, j, _x);
				(*YInt)->SetElement(2*i+1, j, _y);
				(*ZInt)->SetElement(2*i+1, j, _z);
			//	printf ("\n int %d (%f, %f, %f)", j, (*XInt)->Element(2*i+1, j), (*YInt)->Element(2*i+1, j), (*ZInt)->Element(2*i+1, j));
			}
			if (pg->IntProf->GetLignes() < IntProfCent->GetLignes()){
				for ( j = pg->IntProf->GetLignes(); j < IntProfCent->GetLignes(); j++) {
					(*XInt)->SetElement(2*i+1, j, _x);
					(*YInt)->SetElement(2*i+1, j, _y);
					(*ZInt)->SetElement(2*i+1, j, _z);
			//		printf ("\n +int %d (%f, %f, %f)", j, (*XInt)->Element(2*i+1, j), (*YInt)->Element(2*i+1, j), (*ZInt)->Element(2*i+1, j));
				}
			}

			//printf ("\n\n ext(%d) int(%d)", pg->ExtProf->GetLignes(), pg->IntProf->GetLignes());
			//printf ("\n\n BECOME ext(%d) int(%d)", ExtProfCent->GetLignes(), IntProfCent->GetLignes());
			for (j=0; j < pg->ExtProf->GetLignes(); j++)
			{
/*				if (isPercent) {
					xp = ExtProfCentN->Element(j,0)*coeffx;
					yp = ExtProfCentN->Element(j,1)*coeffyCent*m
							+ ExtProfBoutN->Element(j,1)*coeffyBout*(1.0f-m);
				} else { */
				xp = pg->ExtProf->Element(j, 0) * coeffx;
				yp = pg->ExtProf->Element(j, 1) * coeffy;
				//}
				_x = xo+(yp-xp*(double)sin(v))*(double)cos(a);
				_y = yo+(yp-xp*(double)sin(v))*(double)sin(a);
				_z = zo-(xp*(double)cos(v));
				(*XExt)->SetElement(2*i+1, j,_x);
				(*YExt)->SetElement(2*i+1, j, _y);
				(*ZExt)->SetElement(2*i+1, j, _z);
			//	printf ("\n ext %d (%f, %f, %f)",j, (*XExt)->Element(2*i+1, j), (*YExt)->Element(2*i+1, j), (*ZExt)->Element(2*i+1, j));
			}
			if (pg->ExtProf->GetLignes() < ExtProfCent->GetLignes()){
				for ( j = pg->ExtProf->GetLignes(); j < ExtProfCent->GetLignes(); j++) {
					(*XExt)->SetElement(2*i+1, j, _x);
					(*YExt)->SetElement(2*i+1, j, _y);
					(*ZExt)->SetElement(2*i+1, j, _z);
			//		printf ("\n +ext %d (%f, %f, %f)",j, (*XExt)->Element(2*i+1, j), (*YExt)->Element(2*i+1, j), (*ZExt)->Element(2*i+1, j));
				}
			}

			//printf ("\n --%d nervure ballone", i);        
		}
	}
}

/***********************/
/* InterpoleProfilBout */
/***********************/
void InterpoleProfilBout(Matrice** XYBout, Matrice* XYCent)

{
	/*interpolation des Y du profil du bout avec les X du profil central*/

	Matrice *Xi=Zeros(XYCent->GetLignes(),1); 
	Matrice *Yi=Zeros(XYCent->GetLignes(),1);

	Matrice *X=Zeros((*XYBout)->GetLignes(),1); 
	Matrice *Y=Zeros((*XYBout)->GetLignes(),1);

	for(int i=0; i<(*XYBout)->GetLignes(); i++)
	{
		X->SetElement(i,0,(*XYBout)->Element(i,0)); 
		Y->SetElement(i,0,(*XYBout)->Element(i,1));
	}

	for(int i=0; i<XYCent->GetLignes(); i++)
	{
		Xi->SetElement(i,0,XYCent->Element(i,0));
	}

	InterpLinMat(X, Y, Xi, Yi); 


	/*affectation au profil du bout*/

	if ( *XYBout != NULL )
		delete((*XYBout));

	(*XYBout)=Zeros(XYCent->GetLignes(),2);

	for(int i=0; i<XYCent->GetLignes(); i++)
	{
		(*XYBout)->SetElement(i,0,Xi->Element(i,0)); 
		(*XYBout)->SetElement(i,1,Yi->Element(i,0));
	}

	/*liberation matrice intermediaire*/

	if ( X != NULL )
		delete(X); 
	if ( Y != NULL )
		delete(Y);
	if ( Xi != NULL )
		delete(Xi); 
	if ( Yi != NULL )
		delete(Yi);

}



/**********/
/* dist3d */
/**********/
double dist3d(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return (double) sqrt(sqr(x1-x2)+sqr(y1-y2)+sqr(z1-z2));
}


/**********/
/* dist2d */
/**********/
double dist2d(double x1, double y1, double x2, double y2)
{
    return (double) sqrt(sqr(x1-x2)+sqr(y1-y2));
}

/*********************************************************************
* int_cer
* role:
*	retourne les points d'intersection de deux cercles
* parametres en entree:
*	coordonnees x,y du centre o1 et rayon r1 du premier cercle
*	coordonnees x,y du centre o2 et rayon r2 du second cercle
* resultat:
*	nombre de points solutions, coordonnees x,y de deux points p1 et p2
***********************************************************************/

void int_cer(
		double xo1, double yo1, double r1,
		double xo2, double yo2, double r2,
		int *nbr_sol, double *xp1, double *yp1, double *xp2, double *yp2)
{
	double normeO1O2, normeO1P, normeO2P;
	double a,b,c,A,B,C,delta;
	double x[4],y[4],distR[4];
	double distmin1, distmin2;
	int i, imin1, imin2;
	//init par defaut
	*xp1=0.0f; *yp1=0.0f; *xp2=0.0f; *yp2=0.0f; 
	normeO1O2 = (double)sqrt(sqr(xo1-xo2)+sqr(yo1-yo2));
	////////test cercles non secants ?
	if ( (normeO1O2>(r1+r2)) || (normeO1O2<fabs(r1-r2)) )
	{
		*nbr_sol=0;
	}
	////////test cercles tangents exterieurs ?
	else if (normeO1O2==(r1+r2))
	{
		*nbr_sol=1;
		*xp1=xo1+r1/(r1+r2)*(xo2-xo1);
		*yp1=yo1+r1/(r1+r2)*(yo2-yo1);
	}
	////////test cercles tangents interieurs ?
	else if (normeO1O2==fabs(r1-r2))
	{
		////////test cercles confondus ?
		if (normeO1O2==0.0f)
		{
			*nbr_sol=0;
		}
		////////cercles tangents interieurs en un point
		else 
		{
			*nbr_sol=1;
			*xp1=xo1+r1/(r1-r2)*(xo2-xo1);
			*yp1=yo1+r1/(r1-r2)*(yo2-yo1);
		}
	}
	////////cercles secants en 2 points ...
	else
	{
		*nbr_sol=2;
		////////resolution du systeme dans le repere centre sur o1
		A=xo1-xo2;
		B=yo1-yo2;
		C=r1*r1-r2*r2+A*A+B*B;
		a=4*A*A+4*B*B;
		b=2*A*C;
		c=C*C-4*B*B*r1*r1;
		////////resolution du systeme d'ordre 2: a*x^2+2*b*x+c=0
		delta=b*b-a*c;
		x[0]=(-b+(double)sqrt(delta))/a;
		x[1]=(-b-(double)sqrt(delta))/a;
		////////les quatres solutions resultantes
		x[2]=x[0];
		x[3]=x[1];
		y[0]=(double)sqrt(r1*r1-x[0]*x[0]);
		y[1]=(double)sqrt(r1*r1-x[1]*x[1]);
		y[2]=-y[0];
		y[3]=-y[1];
		////////elimination de 2 solutions
		////////remarque: repere centre sur O1
		//remarques: P point de coordonnees x[i],y[i]
		for(i=0; i<4; i++)
		{
			normeO1P = (double)sqrt( sqr(x[i]) + sqr(y[i]) );
			normeO2P = (double)sqrt( sqr(x[i]+xo1-xo2) + sqr(y[i]+yo1-yo2) );
			distR[i] = (double)sqrt( sqr(normeO1P-r1) + sqr(normeO2P-r2) );
		}
		// choix de la solution au sens des moindres carrï¿½s ...
		//test normeO1P=r1 et normeO2P=r2
		//choix des deux valeurs mini ...
		distmin1=1000000000000.0f; imin1=0;
		for(i=0; i<4; i++){
			if (distmin1 > distR[i]){
				distmin1 = distR[i]; imin1 = i;}}
		distmin2=1000000000000.0f; imin2=0;
		for(i=0; i<4; i++){
			if ( (distmin2 > distR[i]) && (imin1 != i) ){
				distmin2 = distR[i]; imin2 = i;}}
		*xp1 = x[imin1]; *yp1 = y[imin1]; 
		*xp2 = x[imin2]; *yp2 = y[imin2]; 
		////////changement de repere
		*xp1=*xp1+xo1; *yp1=*yp1+yo1;
		*xp2=*xp2+xo1; *yp2=*yp2+yo1;
	}
	//resultat=[nbr_sol xp1 yp1 xp2 yp2];
}


/*********************************************************************
* int_cer_bis
* role:
*	retourne les points d'intersection de deux cercles
* parametres en entree:
*	coordonnees x,y du centre o1 et rayon r1 du premier cercle
*	coordonnees x,y du centre o2 et rayon r2 du second cercle
* resultat:
*	nombre de points solutions, coordonnees x,y de deux points p1 et p2
***********************************************************************/

void int_cer_bis(double xo1, double yo1, double r1,
		double xo2, double yo2, double r2,
		int *nbr_sol, double *xp1, double *yp1, double *xp2, double *yp2)
{
    //printf ("\n int_cer_bis(%f, %f) r%f (%f, %f) r%f", xo1, yo1, r1, xo2, yo2, r2 );
	double normeO1O2;
	double M,N,a,b,c,delta;
	//init par defaut
	*xp1=0.0f; *yp1=0.0f; *xp2=0.0f; *yp2=0.0f; 
	normeO1O2 = (double)sqrt(sqr(xo1-xo2)+sqr(yo1-yo2));
    //printf ("\n norme0102=%f (r1+r2)=%f", normeO1O2, (r1+r2));
	////////test cercles non secants ?
	if ( (normeO1O2>(r1+r2)) || (normeO1O2<fabs(r1-r2)) )
	{
		*nbr_sol=0;
	}
	////////test cercles tangents exterieurs ?
	else if (normeO1O2==(r1+r2))
	{
		*nbr_sol=1;
		*xp1=xo1+r1/(r1+r2)*(xo2-xo1);
		*yp1=yo1+r1/(r1+r2)*(yo2-yo1);
	}
	////////test cercles tangents interieurs ?
	else if (normeO1O2==fabs(r1-r2))
	{
		////////test cercles confondus ?
		if (normeO1O2==0.0f)
		{
			*nbr_sol=0;
		}
		////////cercles tangents interieurs en un point
		else 
		{
			*nbr_sol=1;
			*xp1=xo1+r1/(r1-r2)*(xo2-xo1);
			*yp1=yo1+r1/(r1-r2)*(yo2-yo1);
		}
	}
	////////cercles secants en 2 points ...
	else
	{
		*nbr_sol=2;
		//on se ramene ï¿½ un systï¿½me a 2 equations: une d'ordre 1 et l'autre d'ordre 2
		//  y=Mx+N
		//	(x-xo1)^2+(y-yo1)^2-r1^2=0
		// avec pour M et N:
		M=(xo1-xo2)/(yo2-yo1);

		N=(sqr(xo2)-sqr(xo1)+sqr(yo2)-sqr(yo1)-sqr(r2)+sqr(r1))/(2*(yo2-yo1));
		a=1+sqr(M);
		b=2*M*(N-yo1)-2*xo1;
		c=sqr(xo1)+sqr(N-yo1)-sqr(r1);
		////////resolution du systeme d'ordre 2: a*x^2+2*b*x+c=0
		delta=b*b-4*a*c;
                if (delta < 0) printf ("\n ERROR! Delta<0 in int_cer_bis");
                //printf ("\nM=%f N=%f delta=%f", M, N, delta);
                //printf (" xp1=%f", (-b+(double)sqrt(delta))/(2*a));
		*xp1=(-b+(double)sqrt(delta))/(2*a);
		*xp2=(-b-(double)sqrt(delta))/(2*a);
                //printf (" *xp1=%f", *xp1);
		////////les quatres solutions resultantes
		*yp1=(double)(M*(*xp1)+N);
		*yp2=(double)(M*(*xp2)+N);
	}
	//resultat=[nbr_sol xp1 yp1 xp2 yp2];
}



/*******************/
/* CalculDeveloppe */
/*******************/

void CalculDeveloppe(
					 //coordonnï¿½es 3D des 2 cotes de la surface 
					 Matrice *Xs1, Matrice *Ys1, Matrice *Zs1,
					 Matrice *Xs2, Matrice *Ys2, Matrice *Zs2,
					 //coordonnï¿½es 2D du dï¿½veloppï¿½
					 Matrice **Xb1, Matrice **Yb1,
					 Matrice **Xb2, Matrice **Yb2)
{
	double diag, dX, dY, X1, Y1, X2, Y2, Scal1, Scal2;
	int i, N;
	//pour rotation
	Matrice *Thb1=NULL, *Thb2=NULL, *Rb1=NULL, *Rb2=NULL;
	double Th;
	//test nb pts identique
	if (Xs1->GetLignes() != Xs2->GetLignes())
	{
		printf("\n erreur procedure CalculDeveloppe");
		fflush(stdin); getch(); exit(0);
	}
	//init tableaux
	*Xb1=Zeros(Xs1->GetLignes(),1);	*Yb1=Zeros(Xs1->GetLignes(),1);
	*Xb2=Zeros(Xs1->GetLignes(),1);	*Yb2=Zeros(Xs1->GetLignes(),1);
	(*Yb2)->SetElement(0,0, dist3d(
		Xs2->Element(0,0),Ys2->Element(0,0),Zs2->Element(0,0),
        	Xs1->Element(0,0),Ys1->Element(0,0),Zs1->Element(0,0))
            );
        //printf ("\nCD (*Yb2)->(0,0)=%f", (*Yb2)->Element(0,0));
	//boucle 
	for(i=0; i<Xs1->GetLignes()-1; i++)
	{
		/////////////////////////////////////
		//calcul point suivant sur le bord 2
		/////////////////////////////////////
		diag=dist3d(
			Xs2->Element(i+1,0),Ys2->Element(i+1,0),Zs2->Element(i+1,0),
			Xs1->Element(i,0),Ys1->Element(i,0),Zs1->Element(i,0));     
		dX=dist3d(
			Xs2->Element(i+1,0),Ys2->Element(i+1,0),Zs2->Element(i+1,0),
			Xs2->Element(i,0),Ys2->Element(i,0),Zs2->Element(i,0));     
		//test b1(i) et b2(i) confondus ?
		//printf ("\n_1_ i=%d diag=%f dX=%f", i, diag, dX);
		if ((Xs2->Element(i,0)==Xs1->Element(i,0))
			&&(Ys2->Element(i,0)==Ys1->Element(i,0))
			&&(Zs2->Element(i,0)==Zs1->Element(i,0)))
		{
                    //printf ("\n Xs2->Element(i,0)==Xs1->Element(i,0)");
        		(*Xb2)->SetElement(i+1,0,(*Xb2)->Element(i,0));
			(*Yb2)->SetElement(i+1,0,(*Yb2)->Element(i,0)+dX);
                }
		else
		{
                    //printf ("\n Xs2->Element(i,0)!=Xs1->Element(i,0)");
			//test points b2->Element(i,0) et b2->Element(i+1,0) confondus ?
			if (dX!=0.0f)
			{
                            //printf ("\n dX!=0.0f");
				//int_cer(
				int_cer_bis(
					(*Xb1)->Element(i,0), (*Yb1)->Element(i,0), diag,
					(*Xb2)->Element(i,0), (*Yb2)->Element(i,0),dX,
					&N, &X1, &Y1, &X2, &Y2);
				//printf ("\nN=%d (%f, %f) (%f, %f)", N, X1, Y1, X2, Y2);
				//choix du point resultat ?
				//calcul de produits scalaires ...
				if (i>0)
				{
					Scal1 =
						(X1-(*Xb2)->Element(i,0))*((*Xb2)->Element(i,0)-(*Xb2)->Element(i-1,0))
						+(Y1-(*Yb2)->Element(i,0))*((*Yb2)->Element(i,0)-(*Yb2)->Element(i-1,0));
					Scal2 =
						(X2-(*Xb2)->Element(i,0))*((*Xb2)->Element(i,0)-(*Xb2)->Element(i-1,0))
						+(Y2-(*Yb2)->Element(i,0))*((*Yb2)->Element(i,0)-(*Yb2)->Element(i-1,0));
					if (Scal1>Scal2)
        				{
						(*Xb2)->SetElement(i+1,0,X1);
						(*Yb2)->SetElement(i+1,0,Y1);
					}
					else
					{
						(*Xb2)->SetElement(i+1,0,X2); 
						(*Yb2)->SetElement(i+1,0,Y2);
					}
				}
				else //(i=0) 1er point calcule du bord 2
				{
					if (X1>X2)
					{
						(*Xb2)->SetElement(i+1,0,X1);
						(*Yb2)->SetElement(i+1,0,Y1);
					}
					else
					{
						(*Xb2)->SetElement(i+1,0,X2);
						(*Yb2)->SetElement(i+1,0,Y2);
					}
				}
			}
			else //points b2->Element(i,0) et b2->Element(i+1,0) confondus
			{
                            //printf ("\n dX==0.0f");
				(*Xb2)->SetElement(i+1,0,(*Xb2)->Element(i,0));
				(*Yb2)->SetElement(i+1,0,(*Yb2)->Element(i,0));
			} //test points b2->Element(i,0) et b2->Element(i+1,0) confondus ?   
		} //test b1->Element(i,0) et b2->Element(i,0) confondus ?
		/////////////////////////////////////
		//calcul point suivant sur le bord 1
		/////////////////////////////////////
		dY=dist3d(
			Xs2->Element(i+1,0),Ys2->Element(i+1,0),Zs2->Element(i+1,0),
			Xs1->Element(i+1,0),Ys1->Element(i+1,0),Zs1->Element(i+1,0));     
		dX=dist3d(
			Xs1->Element(i+1,0),Ys1->Element(i+1,0),Zs1->Element(i+1,0),
			Xs1->Element(i,0),Ys1->Element(i,0),Zs1->Element(i,0));
		//test points b1->Element(i,0) et b1->Element(i+1,0) confondus ?
                //printf ("\n_i=%d dY=%f dX=%f", i, dY, dX);
		if (dX!=0)
		{
			//test points b1->Element(i+1,0) et b2->Element(i+1,0) confondus ?
			if (dY!=0)
			{
				//int_cer(
				int_cer_bis(
					(*Xb2)->Element(i+1,0),(*Yb2)->Element(i+1,0),dY,
					(*Xb1)->Element(i,0),(*Yb1)->Element(i,0),dX,
					&N, &X1, &Y1, &X2, &Y2);
                //printf ("\n_N=%d (%f, %f) (%f, %f)", N, X1, Y1, X2, Y2);
				//choix du point resultat ?
				//calcul de produits scalaires ...
				if (i>1)
        			{
                			Scal1 =
						(X1-(*Xb1)->Element(i,0))*((*Xb1)->Element(i,0)-(*Xb1)->Element(i-1,0))
						+(Y1-(*Yb1)->Element(i,0))*((*Yb1)->Element(i,0)-(*Yb1)->Element(i-1,0));
					Scal2 =
						(X2-(*Xb1)->Element(i,0))*((*Xb1)->Element(i,0)-(*Xb1)->Element(i-1,0))
						+(Y2-(*Yb1)->Element(i,0))*((*Yb1)->Element(i,0)-(*Yb1)->Element(i-1,0));
					if (Scal1>Scal2)
					{
						(*Xb1)->SetElement(i+1,0,X1); 
						(*Yb1)->SetElement(i+1,0,Y1);
					}
					else
					{
						(*Xb1)->SetElement(i+1,0,X2); 
						(*Yb1)->SetElement(i+1,0,Y2);
					}
				}
				else //1er point calcule bord 1
				{
					if (X1>X2)
					{
						(*Xb1)->SetElement(i+1,0,X1);
						(*Yb1)->SetElement(i+1,0,Y1);
					}
					else
					{
						(*Xb1)->SetElement(i+1,0,X2);
						(*Yb1)->SetElement(i+1,0,Y2);
					}
				}
			}
			else //b1->Element(i+1,0) et b2->Element(i+1,0) confondus
			{
				(*Xb1)->SetElement(i+1,0,(*Xb2)->Element(i+1,0));
				(*Yb1)->SetElement(i+1,0,(*Yb2)->Element(i+1,0));
			} //test points b1->Element(i+1,0) et b2->Element(i+1,0) confondus
		}
		else
		{
			(*Xb1)->SetElement(i+1,0,(*Xb1)->Element(i,0));
			(*Yb1)->SetElement(i+1,0,(*Yb1)->Element(i,0));
		}//test b1->Element(i,0) et b1->Element(i+1,0) points confondus   
	}	//boucle point
	//////////////////////////////////////////
	// aligne extremite du bord 1 sur l'axe X
	//////////////////////////////////////////
	//passage en coordonnees polaire
        for (i=0;i<(*Xb1)->GetLignes();i++){
            //printf ("\nCDBP%d 1(%f, %f)  2(%f, %f)",i, (*Xb1)->Element(i,0), (*Yb1)->Element(i,0), (*Xb2)->Element(i,0), (*Yb2)->Element(i,0));
        }
	Cart2Pol(*Xb1, *Yb1, &Thb1, &Rb1);
	Cart2Pol(*Xb2, *Yb2, &Thb2, &Rb2);
	//rotation
	for(i=0; i<Thb1->GetLignes(); i++)
	{
		Th = Thb1->Element(Thb1->GetLignes()-1,0);
		Thb1->RemoveElement(i,0,Th);
		Thb2->RemoveElement(i,0,Th);
	}
	//passage en coordonnees cartesienne
	if ( *Xb1 != NULL )
		delete(*Xb1);
	*Xb1 = NULL;
	if ( *Yb1 != NULL )
		delete(*Yb1);
	*Yb1 = NULL;
	if ( *Xb2 != NULL )
		delete(*Xb2);
	*Xb2 = NULL;
	if ( *Yb2 != NULL )
		delete(*Yb2);
	*Yb2 = NULL;
	Pol2Cart(Thb1, Rb1, Xb1, Yb1);
	Pol2Cart(Thb2, Rb2, Xb2, Yb2);
	if ( Thb1 != NULL )
		delete(Thb1); 
	if ( Rb1 != NULL )
		delete(Rb1);
	if ( Thb2 != NULL )
		delete(Thb2); 
	if ( Rb2 != NULL )
		delete(Rb2);
}



/*************************************************/
/* passage de coordonnï¿½es cartesiennes ï¿½ polaire */
/*************************************************/

void Cart2Pol(Matrice *X, Matrice *Y, Matrice **T, Matrice **R)

{

	int i;

	*R=new Matrice(X->GetLignes(),1); 
	*T=new Matrice(X->GetLignes(),1);

	for(i=0; i<X->GetLignes(); i++)

	{

		(*T)->SetElement(i,0, (double)atan2(Y->Element(i,0), X->Element(i,0)));

		(*R)->SetElement(i,0, (double)sqrt(sqr(Y->Element(i,0))+sqr(X->Element(i,0))));

	}

}

/*************************************************/
/* passage de coordonnï¿½es polaire ï¿½ cartesiennes */
/*************************************************/

void Pol2Cart(Matrice *T, Matrice *R, Matrice **X, Matrice **Y)

{

	int i;

	*X=new Matrice(R->GetLignes(),1); 
	*Y=new Matrice(R->GetLignes(),1);

	for(i=0; i<R->GetLignes(); i++)

	{

		(*X)->SetElement(i,0, (double)(R->Element(i,0)*cos(T->Element(i,0))));

		(*Y)->SetElement(i,0, (double)(R->Element(i,0)*sin(T->Element(i,0))));

	}

}

/*******************************/
/* intersectione de 2 vecteurs */
/*******************************/
void Inter2Vecteurs(double xa, double ya,
         		double xb, double yb,
			double xc, double yc,
			double xd, double yd,
			double *x, double *y)
{
	double xab, yab, xcd, ycd;
	double a11,a12,a21,a22,b1,b2,delta,delta1,delta2;
	//poursimplifier ecriture
	xab=xb-xa; yab=yb-ya;
	xcd=xd-xc; ycd=yd-yc;
	//preparation matrices A et B pour resolution AX=B
	a11=yab; a12=-xab; b1=xa*yab-ya*xab;
	a21=ycd; a22=-xcd; b2=xc*ycd-yc*xcd;
	//calcul determinants
	delta=a11*a22-a21*a12;
	if (delta !=0)
	{
		delta1=b1*a22-b2*a12;
		delta2=a11*b2-a21*b1;
		/*calcul resultat*/
		*x=delta1/delta;
		*y=delta2/delta;
	}
	else
	{
		printf("\npb: impossibilitï¿½ de calcul dans la procï¿½dure 'Inter2Vecteurs'");
		*x=0.0f; *y=0.0f;
	}
}

/*******************************************************************************
*
*     Linear vorticity surface panel method for airfoils.
*
*     Adapted from Kuethe and Chow 4th Edition
*     This version by I. Kroo  2-2-87
*
*     No attempt has been made to make this particularly fast or efficient.
*
*     Inputs:
*     -------
*     XB, YB     x and y coordinates of panel edges starting at trailing edge
*                proceeding forward on lower surface, wrapping around
*                leading edge and then running back to the trailingedge.
*     alphaD     Angle of attack in degrees
*
*     Outputs:
*     --------
*     X, Y       coordinates of panel centers
*     Cp         incompressible pressure coefficient at panel center
*     
********************************************************************************/

// version originale en fortran, traduit en C et Matlab par Thierry Pï¿½bayle, 2001
// function [X,Y,Cp]=LVFoil(XB,YB,alphaD);

void LVFoil(Matrice *XB, Matrice *YB, double alphaD, Matrice **X, Matrice **Y, Matrice **Cp )
{
	int n, np1, i, ip1, j;
	double alpha, SINA, COSA;
	Matrice *S=NULL, *THETA=NULL, *SINT=NULL, *COST=NULL;
	Matrice *CN1=NULL, *CN2=NULL, *CT1=NULL, *CT2=NULL, *AN=NULL, *AT=NULL, *RHS=NULL, *XX=NULL, *V=NULL;
	double A,B,C,D,E,F,G,P1,P2,P3,P4,P,Q;
	double ANmax;
	//
	//     Calculation of geometric data:
	//     ------------------------------

	n=XB->GetLignes()-1; // n:Number of panels
	np1 = n+1;

	

	//init matrices

	*X=Zeros(n,1); *Y=Zeros(n,1); *Cp=Zeros(n,1);

	S=Zeros(n,1); THETA=Zeros(n,1); SINT=Zeros(n,1); COST=Zeros(n,1);

	

	//calcul position centre et taille des panneaux 

	for(i=0;i<n; i++)

	{

		ip1=i+1;

		(*X)->SetElement(i,0, .5f * (XB->Element(i,0) + XB->Element(ip1,0)));

		(*Y)->SetElement(i,0, .5f * (YB->Element(i,0) + YB->Element(ip1,0)));

		S->SetElement(i,0, (double)sqrt( sqr(XB->Element(ip1,0)-XB->Element(i,0)) + sqr(YB->Element(ip1,0)-YB->Element(i,0)) ));

		THETA->SetElement(i,0, (double)atan2( (YB->Element(ip1,0)-YB->Element(i,0)), (XB->Element(ip1,0)-XB->Element(i,0)) ));

		SINT->SetElement(i,0, (double)sin(THETA->Element(i,0)));

		COST->SetElement(i,0, (double)cos(THETA->Element(i,0)));

	}	

	alpha = alphaD * DEG2RAD;

	SINA = (double)sin(alpha); 

	COSA = (double)cos(alpha);

	

	//

	//     Calculation of influence coefficients

	//     -------------------------------------

	//

	

	// init tableaux

	CN1=Zeros(n,n); CN2=Zeros(n,n); CT1=Zeros(n,n); CT2=Zeros(n,n);

	AN=Zeros(np1,np1); AT=Zeros(np1,np1);

	

	for (i=0; i<n; i++)

	{

		for (j=0; j<n; j++)

		{

			if (i==j) 

			{

				CN1->SetElement(i,j, -1.0);

				CN2->SetElement(i,j, 1.0);

				CT1->SetElement(i,j, 0.5*pi);

				CT2->SetElement(i,j, CT1->Element(i,j));

			}

			else

			{

				A = -((*X)->Element(i,0)-XB->Element(j,0))*COST->Element(j,0)

					- ((*Y)->Element(i,0)-YB->Element(j,0))*SINT->Element(j,0);

				B = sqr((*X)->Element(i,0)-XB->Element(j,0)) + sqr((*Y)->Element(i,0)-YB->Element(j,0));

				C = SINT->Element(i,0)*COST->Element(j,0)-COST->Element(i,0)*SINT->Element(j,0);

				D = COST->Element(i,0)*COST->Element(j,0)+SINT->Element(i,0)*SINT->Element(j,0);

				E = ((*X)->Element(i,0)-XB->Element(j,0))*SINT->Element(j,0)

					- ((*Y)->Element(i,0)-YB->Element(j,0))*COST->Element(j,0);

				F = (double)log( 1. + S->Element(j,0)*(S->Element(j,0)+2.*A)/B );

				G = (double)atan2(E*S->Element(j,0), B+A*S->Element(j,0));

				P1 = 1.0f-2.0f*SINT->Element(j,0)*SINT->Element(j,0);

				P2 = 2.0f*SINT->Element(j,0)*COST->Element(j,0);

				P3 = SINT->Element(i,0)*P1 - COST->Element(i,0)*P2;

				P4 = COST->Element(i,0)*P1 + SINT->Element(i,0)*P2;

				P = ((*X)->Element(i,0)-XB->Element(j,0)) * P3 + ((*Y)->Element(i,0)-YB->Element(j,0)) * P4;

				Q = ((*X)->Element(i,0)-XB->Element(j,0)) * P4 - ((*Y)->Element(i,0)-YB->Element(j,0)) * P3;

				CN2->SetElement(i,j, D + ( .5f*Q*F - (A*C+D*E)*G )/S->Element(j,0));

				CN1->SetElement(i,j, .5f*D*F + C*G - CN2->Element(i,j));

				CT2->SetElement(i,j, C + ( .5f*P*F + (A*D-C*E)*G )/S->Element(j,0));

				CT1->SetElement(i,j, .5f*C*F - D*G - CT2->Element(i,j));

			}

		}

	}

	

	ANmax = 0.0;

	for (i=0; i<n; i++)

	{

		AN->SetElement(i,0, CN1->Element(i,0));

		AN->SetElement(i,n, CN2->Element(i,n-1));

		AT->SetElement(i,0, CT1->Element(i,0));

		AT->SetElement(i,n, CT2->Element(i,n-1));

		for (j=1; j<n; j++)

		{

			AN->SetElement(i,j, CN1->Element(i,j) + CN2->Element(i,j-1));

			AT->SetElement(i,j, CT1->Element(i,j) + CT2->Element(i,j-1));

			if(fabs(AN->Element(i,j))>ANmax)

			{

				ANmax = (double)fabs(AN->Element(i,j));

			}

		}

	}

				

	//     The Kutta Condition is imposed by the relation:

	//     A*gamma(1) + A*gamma(n+1) = 0.  A would be 1 but for matrix

	//     conditioning problems.

			

	AN->SetElement(n,0, 1.0);
	AN->SetElement(n,n, 1.0);

	for(j=1; j<n; j++)
	{
		AN->SetElement(n,j, 0.0);
	}
				

	//

	//     Decompose and solve the system:

	//     -------------------------------

	//

	//     AN.XX = RHS

	//     avec  RHS(i) = sin( THETA(i)-alpha )

	//			

	RHS=Zeros(np1,1);

	for(i=0; i<n; i++)

	{

		RHS->SetElement(i,0, SINT->Element(i,0)*COSA - COST->Element(i,0)*SINA);

	}

	//XX=inv(AN)*RHS;

	ResolutionGauss(AN, RHS, &XX);

				

	//

	//     Compute derived quantities:

	//     ---------------------------

	//

	//     V(i) = cos( THETA(i) - alpha )

	V=Zeros(n,1);

	*Cp=Zeros(n,1);

	for (i=0; i<n; i++)

	{

		V->SetElement(i,0,  COST->Element(i,0)*COSA + SINT->Element(i,0)*SINA);

		for (j=0; j<np1; j++)

		{

			V->SetElement(i,0, V->Element(i,0) + AT->Element(i,j) * XX->Element(j,0));

		}

		(*Cp)->SetElement(i,0, 1-V->Element(i,0)*V->Element(i,0));

	}



	//liberation matrices de calcul

	delete(S); delete(THETA); delete(SINT); delete(COST);

	delete(CN1); delete(CN2); delete(CT1); delete(CT2);

	delete(AN); delete(AT); delete(RHS); delete(XX);

	delete(V);



}

/*********************************************/
/* calcul longueur dï¿½veloppï¿½ d'une courbe XY */
/*********************************************/

Matrice* Longueur(Matrice* x, Matrice *y)

{

	int i; 

	Matrice *res;

	

	res=new Matrice(x->GetLignes(),1);

	res->SetElement(0,0,0.0f);//premiï¿½re valeur =0.0 

	for(i=1; i<x->GetLignes(); i++)

		res->SetElement(i,0, res->Element(i-1,0)

		+ (double)sqrt(sqr(x->Element(i,0)-x->Element(i-1,0))+sqr(y->Element(i,0)-y->Element(i-1,0))));

	return res;

}

/*************************************/
/* calcul coordonï¿½es XY d'un cercle  */
/*************************************/

Matrice* Cercle(double xo, double yo, double rayon, int nbp)
{
	int i; 
	Matrice *res;
	res=new Matrice(nbp+1,2);
	for(i=0; i<nbp; i++)
	{

		res->SetElement(i,0, xo + (double)cos(2*pi*i/nbp)*rayon);

		res->SetElement(i,1, yo + (double)sin(2*pi*i/nbp)*rayon);

	}

	res->SetElement(nbp,0, res->Element(0,0)); 
	res->SetElement(nbp,1, res->Element(0,1));

	return res;

}



void CalculMaxWH(Matrice *Xd0, Matrice *Yd0, Matrice *Xd1, Matrice *Yd1, double *width, double *height) {
    // calculate width and height of (Xd[0],Yd[0]) (Xd[1],Yd[1])
    // calcul minX, maxX, minY, maxY
    //printf("\n CalculMaxWH()");
    double minX = 10000000.0f, minY = 10000000.0f, maxX = -10000000.0f, maxY = -10000000.0f, x = 0.0f, y = 0.0f;
    for (int i = 0; i < Xd0->GetLignes(); i++) {
        x = Xd0->Element(i, 0);
        if (x < minX) minX = x;
        if (x > maxX) maxX = x;
    }
    for (int i = 0; i < Xd1->GetLignes(); i++) {
        x = Xd1->Element(i, 0);
        if (x < minX) minX = x;
        if (x > maxX) maxX = x;
    }

    for (int i = 0; i < Yd0->GetLignes(); i++) {
        y = Yd0->Element(i, 0);
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
    }

    for (int i = 0; i < Yd1->GetLignes(); i++) {
        y = Yd1->Element(i, 0);
        if (y < minY) minY = y;
        if (y > maxY) maxY = y;
    }
    *width = fabs(maxX - minX);
    *height = fabs(maxY - minY);
}


void getPointByPos (Matrice *Xd, Matrice *Yd, Matrice *P, double Pos, double *xr, double *yr) {
    Matrice *interpXSuspente, *interpYSuspente;
    interpXSuspente = Zeros(P->GetLignes(), 2);
    interpYSuspente = Zeros(P->GetLignes(), 2);

    for (int j = 0; j < P->GetLignes(); j++) {
        interpXSuspente->SetElement(j, 0, P->Element(j, 0));
        interpYSuspente->SetElement(j, 0, P->Element(j, 0));
        interpXSuspente->SetElement(j, 1, Xd->Element(j, 0));
        interpYSuspente->SetElement(j, 1, Yd->Element(j, 0));
    }
    *xr = InterpLinX(interpXSuspente, Pos);
    *yr = InterpLinX(interpYSuspente, Pos);
}


void CalculPatronWithCoeff(Matrice *Xd0, Matrice *Yd0, Matrice *Xd1, Matrice *Yd1, double coeff, Matrice **newXd0, Matrice **newYd0, Matrice **newXd1, Matrice **newYd1) {
    /*    Matrice *Xd[2], *Yd[2], *newXd[2], *newYd[2];
        Matrice *X[2], *Y[2], *Z[2], *P[2];
        CalculPatron(noNerv, false, 2, 2, Deb[0], Fin[0],
                noNerv, false, 1, 1, Deb[1], Fin[1],
                &Xd[0], &Yd[0], &Xd[1], &Yd[1],
                &X[0], &Y[0], &Z[0], &P[0],
                &X[1], &Y[1], &Z[1], &P[1]);*/
    *newXd0 = MultReelMat(Xd0, coeff);
    *newYd0 = MultReelMat(Yd0, coeff);
    *newXd1 = MultReelMat(Xd1, coeff);
    *newYd1 = MultReelMat(Yd1, coeff);
    //delete (Xd[0]);delete (Yd[0]);delete (X[0]); delete (Y[0]);delete (Z[0]); delete(P[0]);
    //delete (Xd[1]);delete (Yd[1]);delete (X[1]); delete (Y[1]);delete (Z[1]); delete(P[1]);
}

double calculCourbeLength(Matrice* X, Matrice* Y) {
	double res = 0;

	for (int i = 0; i < X->GetLignes()-1; i++) {
		res += dist2d(X->Element(i,0), Y->Element(i,0), X->Element(i+1,0), Y->Element(i+1,0));
	}

	return res;
}

int pointAtSegment(double x, double y, double x1, double y1, double x2, double y2 ) {
	if ((x >= x1) && (x <= x2)) {
		if (  abs(y - ( y1 + (y2-y1)*(x-x1)/(x2-x1) )) < 0.000001  )
			return 1;
	} else
		return -1;

	return 0;
}