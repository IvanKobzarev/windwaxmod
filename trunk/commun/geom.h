#pragma once

#define NO 0
#define ONE 1
#define SEG 2

class Form;
class Form3D;
class FormProjection;
class Matrix;
class Ballonement;
class WindPatternsProject;

class Point2d {
public:
	Point2d();
	double x, y;
	Point2d(double _x, double _y);
};

class Vector2d {
public:
	double x, y;
	Vector2d(Point2d p1, Point2d p2);
	double length();
};

class Segment2d {
public:
	Point2d* p;
	Point2d *p1, *p2;
	double A, B, C;

	bool contains(Point2d* pt);
	Segment2d(Point2d* _p1, Point2d* _p2);
	void print();
};

class ResultIntersect2d {
public:
	ResultIntersect2d();
	int type;
	Point2d *p1, *p2;

};

ResultIntersect2d* intersectSegments2d (Segment2d* s1, Segment2d* s2);

Matrix* MonBezier(Matrix *tab, int NbrPts);

Matrix* calcNormalRasst(Matrix* X1, Matrix* Y1,Matrix* X2, Matrix* Y2) ;
void AddPointsToCourb(Matrix* courb, Matrix* pts, Matrix** newCourb);
void calcVecteurNormal(double xa,double ya,double xb,double yb,
						 double *xc,double *yc,double l,int cote);

Matrix* calcContour(Matrix* xy, Matrix *d, int cote);

void calcVecteurBissec(double x1,double y1,double x2,double y2, double x3,double y3,
						 double *x,double *y,double l,int cote);

void calcForm3D(Form *forme, int isPercent, double percent,
				   Matrix *ExtProfCent, Matrix *IntProfCent,
				   Matrix *ExtProfBout, Matrix *IntProfBout,
				   Matrix **XExt, Matrix **YExt, Matrix **ZExt,
				   Matrix **XInt, Matrix **YInt, Matrix **ZInt);

void calcForm3DBallonement
				(WindPatternsProject* gfd, Form *forme, int isPercent, double percent,
				   Matrix *ExtProfCent, Matrix *IntProfCent,
				   Matrix *ExtProfBout, Matrix *IntProfBout,
				   Matrix **XExt, Matrix **YExt, Matrix **ZExt,
				   Matrix **XInt, Matrix **YInt, Matrix **ZInt);

Form3D* getForm3D(Form *forme, int isPercent, double percent);
				   //Matrix *ExtProfCent, Matrix *IntProfCent,
				   //Matrix *ExtProfBout, Matrix *IntProfBout);

FormProjection* getFormProjection(Form3D* f3d);

double EpaisseurRelative(Matrix* extrados, Matrix* intrados);

double dist2d(double x1, double y1, double x2, double y2);

void InterpoleProfilBout(Matrix** XYBout, Matrix* XYCent);

void calcDeveloppe(
					 //coordonn�es 3D des 2 cotes de la surface 
					 Matrix *Xs1, Matrix *Ys1, Matrix *Zs1,
					 Matrix *Xs2, Matrix *Ys2, Matrix *Zs2,
					 //coordonn�es 2D du d�velopp�
					 Matrix **Xb1, Matrix **Yb1,
					 Matrix **Xb2, Matrix **Yb2);

void Cart2Pol(Matrix *X, Matrix *Y, Matrix **T, Matrix **R);

void Pol2Cart(Matrix *T, Matrix *R, Matrix **X, Matrix **Y);

void Inter2Vecteurs(
					double xa, double ya,
					double xb, double yb, 
					double xc, double yc,
					double xd, double yd,
					double *x, double *y);

void LVFoil(Matrix *XB, Matrix *YB, double alphaD, Matrix **X, Matrix **Y, Matrix **Cp );

Matrix* Longueur(Matrix* x, Matrix *y);

Matrix* Cercle(double xo, double yo, double rayon, int nbp);

void calcMaxWH(Matrix *Xd0, Matrix *Yd0, Matrix *Xd1, Matrix *Yd1, double *width, double *height);

void getPointByPos (Matrix *Xd, Matrix *Yd, Matrix *P, double Pos, double *xr, double *yr);
void getPoint3dByPos (Matrix *X, Matrix *Y, Matrix *Z, Matrix *P, double pos, double *xr, double *yr, double *zr);
void getPoint3dFormByPosDNerv(Form3D* f3d, double nerv, int side, double pos, double *xr, double *yr, double *zr);
void getPoint3dFormByPosNerv(Form3D* f3d, int nerv, int side, double pos, double *xr, double *yr, double *zr);

void calcPatronWithCoeff(Matrix *Xd0, Matrix *Yd0, Matrix *Xd1, Matrix *Yd1, double coeff, Matrix **newXd0, Matrix **newYd0, Matrix **newXd1, Matrix **newYd1);

double dist3d(double x1, double y1, double z1, double x2, double y2, double z2);

double calcCourbeLength(Matrix* X, Matrix* Y);

int pointAtSegment(double x, double y, double x1, double y1, double x2, double y2 );

void makeEqualSize(Matrix** X0,Matrix** Y0, Matrix** Z0,Matrix** P0,
				   Matrix** X1,Matrix** Y1, Matrix** Z1,Matrix** P1);
