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

Color::Color(){
	r = 0.0f; 
	g = 0.0f; 
	b = 0.0f;
}

void Color::print() {
	printf ("\n(%f, %f, %f)", r, g, b);
}

Color::~Color(){
}

ColorTable::ColorTable(){
}

ColorTable::~ColorTable() {
}



KiteDesignElement::KiteDesignElement(){
}

KiteDesignElement::~KiteDesignElement(){
}

Line::Line(){
	n_points = 0;
}

Line::Line(std::ifstream& in){
	in >> colorR >> colorG >> colorB;
	//cout << "Line color: "<< colorR << colorG << colorB  ;
	n_points = 0;
	in >> n_points;
	//cout << "Line n_points: "<< n_points << endl;
	pointsNervs = new int[n_points];
	pointsPercents = new double[n_points];
	for (int j = 0; j < n_points; j++) {
		int pn;
		double pp;
		in >> pointsNervs[j] >> pointsPercents[j];
		//cout << "Line p "<< pointsNervs[j] << " "<< pointsPercents[j] <<endl;
	}
}

ColorTable::ColorTable(std::ifstream& in){
	int n_nervs = 0;
	in >> n_nervs;
	for (int i = 0; i < n_nervs; i++) {
		vector <Color*> row;
		int n_row = 0;
		int nerv = 0;
		in >> nerv >> n_row;
		for (int j = 0; j < n_row; j++) {
			Color* c = new Color(in);
			row.push_back(c);
		}
		table.push_back(row);
	}
}

Color::Color(std::ifstream& in){
	in >> r >> g >> b;
}

void Line::initColor(Courbe* courbe) {
	courbe->CouleurSegments[0] = colorR;
	courbe->CouleurSegments[1] = colorG;
	courbe->CouleurSegments[2] = colorB;
}

void Line::ajoutCourbesToAxe(TAxe* axe, FormeProjection* fp, int symetric, double dy, double ymult, double ymin) 
{
	Courbe* courbe = new Courbe("CourbeLine");
	initColor(courbe);
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

void Line::ajoutCourbesToAxe3d(TAxe* axe, Forme3D* f3d, int side, int symetric) { 
	Courbe* courbe = new Courbe("CourbeLine");
	initColor(courbe);
	courbe->points = OFF;
	if (symetric) courbe->symX = ON;
	courbe->pts = new Matrice(n_points, 2);
	double x=0, y=0, z=0;
	for (int i = 0; i < n_points; i++) {
		int nerv = pointsNervs[i];
		double perc = pointsPercents[i];
		getPoint3dFormeByPosNerv(f3d, nerv, side, perc, &x, &y, &z);
		courbe->pts->SetElement(i, 0, x);
		courbe->pts->SetElement(i, 1, y);
		courbe->pts->SetElement(i, 2, z);
	}
	AjoutCourbe(axe, courbe);
}

Line::~Line(){
}

void Line::print(){
	printf ("\nline print()");
	printf ("\n n_points=%d", n_points);
	for (int i = 0; i < n_points; i++){
		printf ("\n pointsNerv:%d pointPercent:%f", pointsNervs[i], pointsPercents[i]);
	}
}

void ColorTable::print() {
	printf ("\nColorTable::print");
	vector<vector<Color*>>::iterator iter_vc;
	vector<Color*>::iterator iter_c;

	for (int i = 0; i < table.size(); i++) {
		printf ("\n%d:", i);
		for (int j = 0; j < table[i].size(); j++) {
			table[i][j]->print();
		}
	}
}

KiteDesign::KiteDesign(){
	n_elements = 0;
}

KiteDesign::~KiteDesign(){
}

void KiteDesign::ajoutMeshesToAxe3d( TAxe *Axe3d, Forme3D* f3d, int side, int symetric){
	printf ("\nKiteDesign::ajoutMeshesToAxe3d");
	Axe3d->eclairage = ON;
	TMesh *Mesh;
	Mesh = CreerMesh(); 

	//Mesh->points = ON;
	//Mesh->segments = ON;  Mesh->faces = OFF;
	Mesh->segments=OFF; Mesh->faces=ON;
	if (side == INT_SIDE) Mesh->InvNormales=ON;
	
	int nerv = 1;
	int nPtsExt = f3d->XExt->GetLignes();

	Mesh->x=new Matrice(nPtsExt, 2);
	Mesh->y=new Matrice(nPtsExt, 2);
	Mesh->z=new Matrice(nPtsExt, 2);
	Mesh->CouleurFaces[0] = 1.0f;
	Mesh->CouleurFaces[1] = 0.0f;
	Mesh->CouleurFaces[2] = 0.0f;
	//Mesh->CouleurSegments[0] = 0.0f;
	//Mesh->CouleurSegments[1] = 0.0f;
	//Mesh->CouleurSegments[2] = 1.0f;
	//Mesh->CouleurPoints[0] = 0.0f;
	//Mesh->CouleurPoints[1] = 0.0f;
	//Mesh->CouleurPoints[2] = 1.0f;
	if(symetric==1)
	{
		Mesh->symX = ON;
	}
	for (int i=0; i<nPtsExt; i++)
	{
		for (int j=nerv; j<nerv+2; j++)
		{
			Mesh->x->SetElement(i,j-nerv, f3d->XExt->Element(j, i));
			Mesh->y->SetElement(i,j-nerv, f3d->YExt->Element(j, i));
			Mesh->z->SetElement(i,j-nerv, f3d->ZExt->Element(j, i));
		}
	}
	AjoutMesh(Axe3d, Mesh);
}