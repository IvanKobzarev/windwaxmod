//design.cpp
#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

#include "design.h"
#include "plot.h"
#include "fichier.h"
#include "profil.h"
#include "geom.h"


KiteDesignElement::KiteDesignElement(){
}

KiteDesignElement::~KiteDesignElement(){
}

Line::Line(){
	n_points = 0;
}

Line::Line(std::ifstream& in){
	n_points = 0;
	in >> n_points;
	cout << "Line n_points: "<< n_points << endl;
	pointsNervs = new int[n_points];
	pointsPercents = new double[n_points];
	for (int j = 0; j < n_points; j++) {
		int pn;
		double pp;
		in >> pointsNervs[j] >> pointsPercents[j];
		cout << "Line p "<< pointsNervs[j] << " "<< pointsPercents[j] <<endl;
	}
}

void Line::ajoutCourbesToAxe(TAxe* axe, FormeProjection* fp, int symetric, double dy, double ymult, double ymin) 
{
	cout << "Line::ajoutCourbesToAxe()" << endl;
	Courbe* courbe = new Courbe("CourbeLine");
	courbe->points = OFF;
	if (symetric) courbe->symX = ON;
	courbe->pts = new Matrice(n_points, 2);

	for (int i = 0; i < n_points; i++) {
		
		int nerv = pointsNervs[i];

		//printf ("\ni=%d nerv=%d", i, nerv);
		double percent = pointsPercents[i];
		
		double x0 = fp->X->Element(i, 0);
		double y0 = fp->Y->Element(i, 0);
		
		double x1 = fp->X->Element(i, 1);
		double y1 = fp->Y->Element(i, 1);
		//printf ("\n(%f,%f) (%f,%f)", x0, y0, x1, y1);
		double xp = x0 + percent * 0.01 * (x1 - x0);
		double yp = y0 + percent * 0.01 * (y1 - y0);
		//printf ("\n(%f,%f)", xp, yp);
		
		courbe->pts->SetElement(i, 0, xp);
		courbe->pts->SetElement(i, 1, dy + ymult*yp - ymin);
	}

	AjoutCourbe(axe, courbe);
}

void Line::ajoutCourbesToAxe3d(TAxe* axe, Forme3D* f3d, KiteDesign* kd, int side, int symetric) { 
	printf ("\nLine::ajoutCourbesToAxe3d %d", side);
}

Line::~Line(){

}

void Line::print(){
	printf ("\nline print()");
	printf ("\n np=%d", n_points);
	for (int i = 0; i < n_points; i++){
		printf ("\n pn:%d pp:%f", pointsNervs[i], pointsPercents[i]);
	}
}

KiteDesign::KiteDesign(){
	n_elements = 0;
}

KiteDesign::~KiteDesign(){
}