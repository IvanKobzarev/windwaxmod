#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>

#include "plot.h"
#include "const.h"
#include "layout.h"
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
* calc Bezier de NbrPts selon la technique geometrique
* tab contient les 4 points ABCD selon le format
* [XA YA
*  XB YB
*  XC YC
*  XD YD]
********************************************************/
Point2d::Point2d(){
}

Point2d::Point2d(double _x, double _y) {
	x = _x;
	y = _y;
}

Vector2d::Vector2d(Point2d p1, Point2d p2) {
	x = p1.x - p2.x;
	y = p1.y - p2.y;
}

double Vector2d::length() {
	return sqrt(x*x+y*y);
}

Segment2d::Segment2d(Point2d* _p1, Point2d* _p2) {
	p1 = _p1;
	p2 = _p2;
	A = (p1->y - p2->y);
	B = (p2->x - p1->x); 
	C = (p1->x*p2->y - p2->x*p1->y);
}

void Segment2d::print() {
	//printf ("\n segment2d::print");
}

bool Segment2d::contains(Point2d* pt) {
	double x0 = pt->x; double y0 = pt->y;
	double x1 = p1->x; double y1 = p1->y;
	double x2 = p2->x; double y2 = p2->y;
	double v = fabs((x0 - x1) * (y2 - y1) - (y0 - y1) * (x2 - x1));
	//printf ("\n v=%f", v);
	bool onLine = (v <= 0.0001f);
	//printf ("\n onLine=%d", onLine);
	double vx = (x0-x1)*(x0-x2);
	//printf ("\n vx=%f", vx);
	bool btwX = (vx <= EPSILON);
	//printf ("\n btwX=%d", btwX);

	double vy = (y0-y1)*(y0-y2);
	//printf ("\n vy=%f", vy);
	bool btwY = (vy <= EPSILON);
	//printf ("\n btwY=%d", btwY);
	return (onLine && btwX && btwY);
}

ResultIntersect2d::ResultIntersect2d() {
}

ResultIntersect2d* intersectSegments2d (Segment2d* s1, Segment2d* s2) {
	double d = s1->A * s2->B - s2->A * s1->B;
	double dx = (-s1->C) * s2->B - (-s2->C) * s1->B;
	double dy = s1->A * (-s2->C) - s2->A * (-s1->C);
	ResultIntersect2d* res = new ResultIntersect2d();
	//printf ("\nintersect d=%f dx=%f dy=%f", d, dx, dy);
	if (fabs(d) < EPSILON) {
		if ((fabs(dx) < EPSILON) && (fabs(dy) < EPSILON)) {
			//printf ("\nSEG");
			res->type = SEG;
		} else {
			//printf ("\nNO_");
			res->type = NO;
		}
	} else {
		double x = dx / d;
		double y = dy / d;
		//printf ("\nintersect (%f, %f)", x, y);
		Point2d* p = new Point2d(x, y);
		if (s1->contains(p) && s2->contains(p)) {
			//printf("ONE");
			res->type = ONE;
			res->p1 = p;
		} else {
			//printf ("\n s1->contains(p)=%d s2->contains(p)=%d",s1->contains(p),s2->contains(p));
			//printf ("\n p (%f, %f)", p->x, p->y );
			//printf ("\n s1 (%f, %f)_(%f, %f)", s1->p1->x, s1->p1->y, s1->p2->x, s1->p2->y );
			//printf ("\n s2 (%f, %f)_(%f, %f)", s2->p1->x, s2->p1->y, s2->p2->x, s2->p2->y );
			//printf(" NO");
			res->type = NO;
		}
	}
	return res;
}


Matrix* MonBezier(Matrix *tab, int NbrPts)
{
	Matrix *res = NULL;
	double XA,YA,XB,YB,XC,YC,XD,YD,XE,YE,XF,YF,XG,YG,XH,YH,XI,YI,XJ,YJ;
	double t,pas;
	int i;
	/*allocation matrice resultat*/
	res=new Matrix(NbrPts,2);
	/*pour simpifier l'write*/
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
/* calcVecteurNormal                               */
/* calc vecteur BC normal au vecteur AB            */
/* direct (cotï¿½=+1) ou indirect (cotï¿½=-1)            */
/* de longeur l.                                     */
/*****************************************************/

void calcVecteurNormal(double xa,double ya,double xb,double yb,
						 double *xc,double *yc,double l,int cote)
{
	double xab, yab, a11, a12, a21, a22, b1, b2; /*var intermï¿½daires*/
	double delta, delta1, delta2; /*determinant principal et secondaires*/
	/*pour simplifier write*/
	xab=xb-xa; yab=yb-ya;
	/*preparation matrices A et B pour resolution AX=B*/
	a11=-yab; a12=xab;
	a21=xab;  a22=yab;
	b1=l*(double)sqrt(xab*xab+yab*yab)*cote +yb*xab -xb*yab;
	b2=xb*xab +yb*yab;
	/*calc determinants*/
	delta=a11*a22-a21*a12;
	if (delta != 0.0f)
	{
		delta1=b1*a22-b2*a12;
		delta2=a11*b2-a21*b1;
        	/*calc resultat*/
		*xc=delta1/delta;
		*yc=delta2/delta;
	}
	else
	{
		printf ("\n * ");
		printf ("\n * ");
		printf ("\n * ");
		printf ("\n geom.cpp calcVecteurNormal delta == 0");
		printf ("\nA(%f, %f) B(%f, %f)", xa, ya, xb, yb);
		printf ("\n !!! calcVecteurNormal return 0.0f 0.0f");
		printf ("\nIMPSBL calcVecteurNormal (%f, %f)(%f, %f)", xa, ya, xb, yb);
		*xc=0.0; *yc=0.0;
	}

}



/*****************************************************/
/* calcContour                                     */
/* xy vecteurs colonnes ...						     */
/* res: matrice contenant le contour                 */
/* d : distance a la courbe en chaque point          */
/* cote: +1:contour gauche, -1:contour droite        */
/*****************************************************/

Matrix* calcContour(Matrix* xy, Matrix *d, int cote)
{
	int i; 
	int n = xy->GetLignes();
	Matrix *res = new Matrix(n, 2);
	double xc, yc;
	if (n > 1) {
		for(i = 0; i < xy->GetLignes()-1; i++)
		{
			xc = res->Element(i+1,0);
			yc = res->Element(i+1,1);
			calcVecteurNormal(xy->Element(i,0), xy->Element(i,1), xy->Element(i+1,0), xy->Element(i+1,1), &xc, &yc, d->Element(i,0), cote);
			res->SetElement(i+1,0,xc);
			res->SetElement(i+1,1,yc);
		}

		/*calc premier point xc,yc*/
		xc = res->Element(0,0);
		yc = res->Element(0,1);
		//printf ("\n calcContour2()");
		calcVecteurNormal(xy->Element(1,0), xy->Element(1,1), xy->Element(0,0), xy->Element(0,1), &xc, &yc, d->Element(0,0), -cote);
		res->SetElement(0,0,xc);
		res->SetElement(0,1,yc);
	} else {
		res->SetElement(0,0,xy->Element(0,0));
		res->SetElement(0,1,xy->Element(0,1));
	}
	return res;
}

Matrix* calcContour1(Matrix* xy, Matrix* xy1, Matrix *d, int cote)
{
	double x = xy->Element(0,0);
	double y = xy->Element(0,1);

	double x1 = xy1->Element(0,0);
	double y1 = xy1->Element(0,1);
	int n1 = xy1->GetLignes();
	double x2 = xy1->Element(n1-1,0);
	double y2 = xy1->Element(n1-1,1);

	double _vx = (x-x1) + (x-x2);
	double _vy = (y-y1) + (y-y2);
	
	double norme = (double) sqrt(sqr(_vx) + sqr(_vy));

	double vx = _vx/norme;
	double vy = _vy/norme;

	Matrix *res = new Matrix(1, 2);
	double dist = d->Element(0,0);
	res->SetElement(0,0, x + dist*vx);
	res->SetElement(0,1, y + dist*vy);

	return res;
}


void AddPointsToCourb(Matrix* courb, Matrix* pts, Matrix** newCourb)
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
    (*newCourb) = new Matrix(n + q, 2);
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

        if ( (i < n) && ((xc < xp) || ((j == np)))) {
          //  printf ("\n < i");
            (*newCourb) -> SetElement(i + j - s, 0, xc);
            (*newCourb) -> SetElement(i + j - s, 1, yc);
            i++;
        } 
        if ( (j < np)  && ((xc > xp) || (i == n)) ) {
            //printf ("\n > j");
            (*newCourb) -> SetElement(i + j - s, 0, xp);
            (*newCourb) -> SetElement(i + j - s, 1, InterpLinX(courb, xp));
            j++;
         }

         if (fabs(xc-xp) < EPSILON) {
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

Matrix* calcNormalRasst(Matrix* X1, Matrix* Y1,Matrix* X2, Matrix* Y2) {
	/*printf ("\nXY1");
	X1->print(0);
	Y1->print(0);
	printf ("\nXY2");
	X2->print(0);
	Y2->print(0);*/

    int n = X1->GetLignes();
	/*for (int i =0; i < n-5; i++) {
		printf ("\n %d (%f, %f)    (%f, %f)", i, X1->Element(i, 0), Y1->Element(i, 0), X2->Element(i, 0), Y2->Element(i, 0));
	}*/
    Matrix *res = new Matrix(n, 1);
    int i = 0;
    res->SetElement(0, 0,dist2d (X1->Element(0,0), Y1->Element(0,0), X2->Element(0,0), Y2->Element(0,0)));
	res->SetElement(n-1, 0,dist2d (X1->Element(n-1,0), Y1->Element(n-1,0), X2->Element(X2->GetLignes()-1,0), Y2->Element(Y2->GetLignes()-1,0)));
    double xn, yn, x , y;
    for (i=1; i < n-1; i++) {
        double x0 = X1->Element(i, 0);
        double y0 = Y1->Element(i, 0);
        calcVecteurBissec(X1->Element(i-1, 0), Y1->Element(i-1, 0),
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
				//printf ("\n set elem %f", dist2d (x,y,x0,y0));
                break;
            }
        }
    }
    return res;
}

/*****************************************************/
/* calcVecteurBissec                               */
/*****************************************************/

void calcVecteurBissec(double x1,double y1,double x2,double y2, double x3,double y3,
						 double *x,double *y,double l,int cote)
{
	double x4,y4,x5,y5;
	double x12,y12,x23,y23;
	double a11,a12,a21,a22,b1,b2;
	double delta, delta1, delta2;
	//pour simplifier calc
	x12=x2-x1; y12=y1-y2;
	x23=x3-x2; y23=y3-y2;
	//test pts 1,2,3 alignï¿½s
	if (fabs(x12*y23-x23*y12)  < EPSILON)
	{
		//printf ("\ncalcVecteurBissec cote1");
		calcVecteurNormal(x1,y1,x2,y2,x,y,l,cote);
	}
	else
	{
		//printf ("\ncalcVecteurBissec cote2");
		calcVecteurNormal(x1,y1,x2,y2,&x4,&y4,l,cote);
		//printf ("\ncalcVecteurBissec -cote");
		calcVecteurNormal(x3,y3,x2,y2,&x5,&y5,l,-cote);
		/*preparation matrices A et B pour resolution AX=B*/
		a11=y12; a12=-x12;
		a21=y23;  a22=-x23;
		b1=x4*y12-y4*x12;
		b2=x5*y23-y5*x23;
		/*calc determinants*/
		delta=a11*a22-a21*a12;
		if (fabs(delta) > EPSILON)
		{
			delta1=b1*a22-b2*a12;
			delta2=a11*b2-a21*b1;
        		/*calc resultat*/
			*x=delta1/delta;
			*y=delta2/delta;
		}
		else
		{
			printf ("\n * ");
			printf ("\n * ");
			printf ("\n * ");
			printf ("\n geom.cpp calcVecteurBissec delta == 0");
			printf ("\n1(%f, %f) 2(%f, %f) 3(%f, %f)", x1, y1, x2, y2, x3, y3);
			printf("\n !!! calcVecteurBissec return 0.0f 0.0f");
			*x=0.0f; *y=0.0f;
		}
	}
}

/*********************/
/* EpaisseurRelative */
/*********************/

double EpaisseurRelative(Matrix* extrados, Matrix* intrados)

{
	int i;
	double MaxExt=-1000.0, MinInt=+1000.0;
	/*calc epaisseur relative -> a refaire car calc approchï¿½*/
	for(i=0; i<extrados->GetLignes(); i++)
		if(MaxExt < extrados->Element(i,1)) MaxExt=extrados->Element(i,1);
	for(i=0; i<intrados->GetLignes(); i++)
		if(MinInt > intrados->Element(i,1)) MinInt=intrados->Element(i,1);
	return MaxExt-MinInt;
}

/***********************/
/* InterpoleProfilBout */
/***********************/
void InterpoleProfilBout(Matrix** XYBout, Matrix* XYCent)
{
	/*interpolation des Y du profil du bout avec les X du profil central*/
	Matrix *Xi=Zeros(XYCent->GetLignes(),1); 
	Matrix *Yi=Zeros(XYCent->GetLignes(),1);
	Matrix *X=Zeros((*XYBout)->GetLignes(),1); 
	Matrix *Y=Zeros((*XYBout)->GetLignes(),1);
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
	if ( X != NULL ) delete(X); 
	if ( Y != NULL ) delete(Y);
	if ( Xi != NULL ) delete(Xi); 
	if ( Yi != NULL ) delete(Yi);
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
		if (fabs(normeO1O2) < EPSILON)
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
		//choix des deux Value mini ...
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
/* calcDeveloppe */
/*******************/

void calcDeveloppe(
					 Matrix *Xs1, Matrix *Ys1, Matrix *Zs1,
					 Matrix *Xs2, Matrix *Ys2, Matrix *Zs2,
					 Matrix **Xb1, Matrix **Yb1,
					 Matrix **Xb2, Matrix **Yb2)
{
	double diag, dX, dY, X1, Y1, X2, Y2, Scal1, Scal2;
	int i, N;
	//pour rotation
	Matrix *Thb1=NULL, *Thb2=NULL, *Rb1=NULL, *Rb2=NULL;
	double Th;
	//test nb pts identique
	if (Xs1->GetLignes() != Xs2->GetLignes())
	{
		printf("\n erreur procedure calcDeveloppe");
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
		//calc point suivant sur le bord 2
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
				//calc de produits scalaires ...
				if (i>0)
				{
					Scal1 = (X1-(*Xb2)->Element(i,0))*((*Xb2)->Element(i,0)-(*Xb2)->Element(i-1,0))
								+ (Y1-(*Yb2)->Element(i,0))*((*Yb2)->Element(i,0)-(*Yb2)->Element(i-1,0));
					Scal2 = (X2-(*Xb2)->Element(i,0))*((*Xb2)->Element(i,0)-(*Xb2)->Element(i-1,0))
								+ (Y2-(*Yb2)->Element(i,0))*((*Yb2)->Element(i,0)-(*Yb2)->Element(i-1,0));
					if (Scal1 > Scal2) {
						(*Xb2)->SetElement(i+1,0,X1);
						(*Yb2)->SetElement(i+1,0,Y1);
					} else {
						(*Xb2)->SetElement(i+1,0,X2); 
						(*Yb2)->SetElement(i+1,0,Y2);
					}
				}
				else //(i=0) 1er point calce du bord 2
				{
					if (X1>X2) {
						(*Xb2)->SetElement(i+1,0,X1);
						(*Yb2)->SetElement(i+1,0,Y1);
					} else {
						(*Xb2)->SetElement(i+1,0,X2);
						(*Yb2)->SetElement(i+1,0,Y2);
					}
				}
			} else //points b2->Element(i,0) et b2->Element(i+1,0) confondus 
			{
                            //printf ("\n dX==0.0f");
				(*Xb2)->SetElement(i+1,0,(*Xb2)->Element(i,0));
				(*Yb2)->SetElement(i+1,0,(*Yb2)->Element(i,0));
			} //test points b2->Element(i,0) et b2->Element(i+1,0) confondus ?   
		} //test b1->Element(i,0) et b2->Element(i,0) confondus ?
		/////////////////////////////////////
		//calc point suivant sur le bord 1
		/////////////////////////////////////
		dY=dist3d(
			Xs2->Element(i+1,0),Ys2->Element(i+1,0),Zs2->Element(i+1,0),
			Xs1->Element(i+1,0),Ys1->Element(i+1,0),Zs1->Element(i+1,0));     
		dX=dist3d(
			Xs1->Element(i+1,0),Ys1->Element(i+1,0),Zs1->Element(i+1,0),
			Xs1->Element(i,0),Ys1->Element(i,0),Zs1->Element(i,0));
		//test points b1->Element(i,0) et b1->Element(i+1,0) confondus ?
                //printf ("\n_i=%d dY=%f dX=%f", i, dY, dX);
		if (dX != 0.0f)
		{
			//test points b1->Element(i+1,0) et b2->Element(i+1,0) confondus ?
			if (dY != 0.0f)
			{
				//int_cer(
				int_cer_bis(
					(*Xb2)->Element(i+1,0),(*Yb2)->Element(i+1,0),dY,
					(*Xb1)->Element(i,0),(*Yb1)->Element(i,0),dX,
					&N, &X1, &Y1, &X2, &Y2);
                //printf ("\n_N=%d (%f, %f) (%f, %f)", N, X1, Y1, X2, Y2);
				//choix du point resultat ?
				//calc de produits scalaires ...
				if (i>1) {
                	Scal1 = (X1-(*Xb1)->Element(i,0))*((*Xb1)->Element(i,0)-(*Xb1)->Element(i-1,0))
								+(Y1-(*Yb1)->Element(i,0))*((*Yb1)->Element(i,0)-(*Yb1)->Element(i-1,0));
					Scal2 =	(X2-(*Xb1)->Element(i,0))*((*Xb1)->Element(i,0)-(*Xb1)->Element(i-1,0))
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
				else //1er point calce bord 1
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

void Cart2Pol(Matrix *X, Matrix *Y, Matrix **T, Matrix **R)
{
	int i;
	*R=new Matrix(X->GetLignes(),1); 
	*T=new Matrix(X->GetLignes(),1);

	for(i=0; i<X->GetLignes(); i++)
	{
		(*T)->SetElement(i,0, (double)atan2(Y->Element(i,0), X->Element(i,0)));
		(*R)->SetElement(i,0, (double)sqrt(sqr(Y->Element(i,0))+sqr(X->Element(i,0))));
	}
}

/*************************************************/
/* passage de coordonnï¿½es polaire ï¿½ cartesiennes */
/*************************************************/

void Pol2Cart(Matrix *T, Matrix *R, Matrix **X, Matrix **Y)
{
	int i;
	*X=new Matrix(R->GetLignes(),1); 
	*Y=new Matrix(R->GetLignes(),1);
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
					double *x, double *y, bool debug)
{
	if (debug) printf ("\nInter2Vecteur A(%f, %f)->B(%f, %f)", xa, ya, xb, yb);
	if (debug) printf ("\nInter2Vecteur C(%f, %f)->D(%f, %f)", xc, yc, xd, yd);
	if (debug) printf ("\nInter2Vecteur AB(%f, %f)   CD(%f, %f)", (xb - xa), (yb - ya), (xd - xc), (yd - yc));
	double xab, yab, xcd, ycd;
	double a11, a12, a21, a22, b1, b2, delta, delta1, delta2;
	//poursimplifier write
	xab = xb - xa; yab = yb - ya;
	xcd = xd-xc; ycd=  yd - yc;
	//preparation matrices A et B pour resolution AX=B
	a11 = yab; a12 = -xab; b1 = xa*yab - ya*xab;
	a21 = ycd; a22 = -xcd; b2 = xc*ycd - yc*xcd;
	//calc determinants
	delta = a11*a22 - a21*a12;
	if (debug) printf ("\n delta=%f", delta);
	if (fabs(delta) > EPSILON )
	{
		delta1 = b1*a22 - b2*a12;
		delta2 = a11*b2 - a21*b1;
		/*calc resultat*/
		*x = delta1/delta;
		*y = delta2/delta;
	} else {
		printf ("\n * ");
		printf ("\n * ");
		printf ("\n * ");
		printf ("\n geom.cpp Inter2Vecteurs delta == 0, no ponints ");
		printf ("A(%f, %f)->B(%f, %f)   C(%f, %f)->D(%f, %f)", xa, ya, xb, yb, xc, yc, xd, yd);
		printf("\npb: impossibilitï¿½ de calc dans la procï¿½dure 'Inter2Vecteurs'");
		printf("\n !!! Inter2Vecteurs return 0.0f 0.0f");
		*x=0.0f; *y=0.0f;
	}
	if (debug) printf ("\n res x=%f y=%f", *x, *y);
}

Matrix* Longueur(Matrix* x, Matrix *y)
{
	int i; 
	Matrix *res;
	res=new Matrix(x->GetLignes(),1);
	res->SetElement(0,0,0.0f);//premiï¿½re valeur =0.0 
	for(i=1; i<x->GetLignes(); i++)
		res->SetElement(i,0, res->Element(i-1,0)
		+ (double)sqrt(sqr(x->Element(i,0)-x->Element(i-1,0))+sqr(y->Element(i,0)-y->Element(i-1,0))));
	return res;
}

/*************************************/
/* calc coordonï¿½es XY d'un cercle  */
/*************************************/

Matrix* Cercle(double xo, double yo, double rayon, int nbp)
{
	int i; 
	Matrix *res;
	res=new Matrix(nbp+1,2);
	for(i=0; i<nbp; i++)
	{
		res->SetElement(i,0, xo + (double)cos(2*pi*i/nbp)*rayon);
		res->SetElement(i,1, yo + (double)sin(2*pi*i/nbp)*rayon);
	}
	res->SetElement(nbp,0, res->Element(0,0)); 
	res->SetElement(nbp,1, res->Element(0,1));
	return res;
}

void calcMaxWH(Matrix *Xd0, Matrix *Yd0, Matrix *Xd1, Matrix *Yd1, double *width, double *height) {
    // calcate width and height of (Xd[0],Yd[0]) (Xd[1],Yd[1])
    // calc minX, maxX, minY, maxY
    //printf("\n calcMaxWH()");
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


void getPointByPos (Matrix *Xd, Matrix *Yd, Matrix *P, double Pos, double *xr, double *yr) {
    Matrix *interpXSuspente, *interpYSuspente;
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

void getPoint3dByPos (Matrix *X, Matrix *Y, Matrix *Z, Matrix *P, double pos, double *xr, double *yr, double *zr) {
	//printf ("\ngetPoint3dByPos ()");
    Matrix *interpX, *interpY, *interpZ;
    interpX = Zeros(P->GetLignes(), 2);
    interpY = Zeros(P->GetLignes(), 2);
    interpZ = Zeros(P->GetLignes(), 2);

    for (int i = 0; i < P->GetLignes(); i++) {
        interpX->SetElement(i, 0, P->Element(i, 0));
		interpX->SetElement(i, 1, X->Element(i, 0));

        interpY->SetElement(i, 0, P->Element(i, 0));
        interpY->SetElement(i, 1, Y->Element(i, 0));

        interpZ->SetElement(i, 0, P->Element(i, 0));
        interpZ->SetElement(i, 1, Z->Element(i, 0));

		//printf ("%f -> %f   ", P->Element(i,0), interpZ->Element(i,1));
    }
    *xr = InterpLinX(interpX, pos);
    *yr = InterpLinX(interpY, pos);
    *zr = InterpLinX(interpZ, pos);
}


void getPoint3dFormByPosDNerv(Form3D* f3d, double nerv, int side, double pos, double *xr, double *yr, double *zr){
	int nerv0 = floor(nerv);
	double nervPos = nerv - nerv0;
	double x0, y0, z0, x1,y1,z1;
	getPoint3dFormByPosNerv(f3d, nerv0, side, pos, &x0, &y0, &z0);
	getPoint3dFormByPosNerv(f3d, nerv0+1, side, pos, &x1, &y1, &z1);

	*xr = x0 + nervPos * (x1 - x0);
	*yr = y0 + nervPos * (y1 - y0); 
	*zr = z0 + nervPos * (z1 - z0);
}

void getPoint3dFormByPosNerv(Form3D* f3d, int nerv, int side, double pos, double *xr, double *yr, double *zr) 
{
	Matrix *X, *Y, *Z, *P;
	if (side == EXT_SIDE) {
		P = f3d->forme->getExtProf(nerv, false);
	} else {
		// INT_SIDE
		P = f3d->forme->getIntProf(nerv, false);
	}


	X = Zeros(P->GetLignes(), 2);
	Y = Zeros(P->GetLignes(), 2);
	Z = Zeros(P->GetLignes(), 2);

	for (int i = 0; i < P->GetLignes(); i++) {
		if (side == EXT_SIDE) {
			X->SetElement(i, 0, f3d->XExt->Element(nerv, i));
			Y->SetElement(i, 0, f3d->YExt->Element(nerv, i));
			Z->SetElement(i, 0, f3d->ZExt->Element(nerv, i));
		} else {
			//INT_SIDE
			X->SetElement(i, 0, f3d->XInt->Element(nerv, i));
			Y->SetElement(i, 0, f3d->YInt->Element(nerv, i));
			Z->SetElement(i, 0, f3d->ZInt->Element(nerv, i));
		}
	}
	double x, y, z;
	getPoint3dByPos(X, Y, Z, P, pos, xr, yr, zr);
}



void calcPatronWithCoeff(Matrix *Xd0, Matrix *Yd0, Matrix *Xd1, Matrix *Yd1, double coeff, Matrix **newXd0, Matrix **newYd0, Matrix **newXd1, Matrix **newYd1) {
    /*    Matrix *Xd[2], *Yd[2], *newXd[2], *newYd[2];
        Matrix *X[2], *Y[2], *Z[2], *P[2];
        calcPatron(noNerv, false, 2, 2, Deb[0], Fin[0],
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

double calcCourbeLength(Matrix* X, Matrix* Y) {
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
