#ifndef __GEOM_H__
#define __GEOM_H__

class Forme;
class Matrice;
class Ballonement;
class WindPatternsProject;

Matrice* MonBezier(Matrice *tab, int NbrPts);

Matrice* CalculNormalRasst(Matrice* X1, Matrice* Y1,Matrice* X2, Matrice* Y2) ;
void AddPointsToCourb(Matrice* courb, Matrice* pts, Matrice** newCourb);
void CalculVecteurNormal(double xa,double ya,double xb,double yb,
						 double *xc,double *yc,double l,int cote);

Matrice* CalculContour(Matrice* xy, Matrice *d, int cote);

void CalculVecteurBissec(double x1,double y1,double x2,double y2, double x3,double y3,
						 double *x,double *y,double l,int cote);

void CalculForme3D(Forme *forme, int isPercent, double percent,
				   Matrice *ExtProfCent, Matrice *IntProfCent,
				   Matrice *ExtProfBout, Matrice *IntProfBout,
				   Matrice **XExt, Matrice **YExt, Matrice **ZExt,
				   Matrice **XInt, Matrice **YInt, Matrice **ZInt);

void CalculForme3DBallonement
				(WindPatternsProject* gfd, Forme *forme, Ballonement* bal, int isPercent, double percent,
				   Matrice *ExtProfCent, Matrice *IntProfCent,
				   Matrice *ExtProfBout, Matrice *IntProfBout,
				   Matrice **XExt, Matrice **YExt, Matrice **ZExt,
				   Matrice **XInt, Matrice **YInt, Matrice **ZInt);

double EpaisseurRelative(Matrice* extrados, Matrice* intrados);

double dist2d(double x1, double y1, double x2, double y2);

void InterpoleProfilBout(Matrice** XYBout, Matrice* XYCent);

void CalculDeveloppe(
					 //coordonn�es 3D des 2 cotes de la surface 
					 Matrice *Xs1, Matrice *Ys1, Matrice *Zs1,
					 Matrice *Xs2, Matrice *Ys2, Matrice *Zs2,
					 //coordonn�es 2D du d�velopp�
					 Matrice **Xb1, Matrice **Yb1,
					 Matrice **Xb2, Matrice **Yb2);

void Cart2Pol(Matrice *X, Matrice *Y, Matrice **T, Matrice **R);

void Pol2Cart(Matrice *T, Matrice *R, Matrice **X, Matrice **Y);

void Inter2Vecteurs(
					double xa, double ya,
					double xb, double yb, 
					double xc, double yc,
					double xd, double yd,
					double *x, double *y);

void LVFoil(Matrice *XB, Matrice *YB, double alphaD, Matrice **X, Matrice **Y, Matrice **Cp );

Matrice* Longueur(Matrice* x, Matrice *y);

Matrice* Cercle(double xo, double yo, double rayon, int nbp);

void CalculMaxWH(Matrice *Xd0, Matrice *Yd0, Matrice *Xd1, Matrice *Yd1, double *width, double *height);

void getPointByPos (Matrice *Xd, Matrice *Yd, Matrice *P, double Pos, double *xr, double *yr);

void CalculPatronWithCoeff(Matrice *Xd0, Matrice *Yd0, Matrice *Xd1, Matrice *Yd1, double coeff, Matrice **newXd0, Matrice **newYd0, Matrice **newXd1, Matrice **newYd1);

double dist3d(double x1, double y1, double z1, double x2, double y2, double z2);

double calculCourbeLength(Matrice* X, Matrice* Y);

int pointAtSegment(double x, double y, double x1, double y1, double x2, double y2 );
#endif