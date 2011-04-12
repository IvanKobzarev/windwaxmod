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

Color::Color(double _r, double _g, double _b) {
	r = _r;
	g = _g;
	b = _b;
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

ColorSegment::ColorSegment(){
	nerv = 0;
	p00 = 0.0f;
	p01 = 0.0f;
	p10 = 0.0f;
	p11 = 0.0f;
}

ColorSegment::~ColorSegment(){
}

int indexBeforePosProf(Matrice* Prof, double pos) {
	int i = 0;
	int n = Prof->GetLignes() - 1;
	while (i < n) {
		if (Prof->Element(i+1, 0) > pos) break;
		i++;
	}
	return i;
}

int indexAfterPosProf(Matrice* Prof, double pos) {
	int i = 0;
	int n = Prof->GetLignes()-1;
	while (Prof->Element(i, 0) < pos) {
		i++;
		if ( i == n ) break;
	}
	return i;
}


void ajoutColorSegmentToAxe3d(TAxe *Axe3d, Forme3D* f3d, ColorSegment* colorSegment, int side, int symetric) {
	printf ("\ndesign:ajoutColorSegmentToAxe3d()");
	Axe3d->eclairage = ON;
	TMesh *Mesh;
	Mesh = CreerMesh(); 

	Mesh->segments=OFF; Mesh->faces=ON;
	if (side == INT_SIDE) Mesh->InvNormales=ON;
	
	int nerv = colorSegment->nerv;
	
	double x00, y00, z00, x01, y01, z01, x10, y10, z10, x11, y11, z11;
	
	getPoint3dFormeByPosNerv(f3d, nerv, side, colorSegment->p00, &x00, &y00, &z00) ;
	getPoint3dFormeByPosNerv(f3d, nerv+1, side, colorSegment->p01, &x01, &y01, &z01) ;

	getPoint3dFormeByPosNerv(f3d, nerv, side, colorSegment->p10, &x10, &y10, &z10) ;
	getPoint3dFormeByPosNerv(f3d, nerv+1, side, colorSegment->p11, &x11, &y11, &z11) ;

	Matrice* mProf;
	if (side == INT_SIDE) mProf = f3d->forme->getIntProf(nerv,false); else mProf = f3d->forme->getExtProf(nerv,false);
	
	int i00 = indexAfterPosProf(mProf,  colorSegment->p00);
	int i01 = indexAfterPosProf(mProf,  colorSegment->p01);

	int i10 = indexBeforePosProf(mProf,  colorSegment->p10);
	int i11 = indexBeforePosProf(mProf,  colorSegment->p11);
	printf ("\ni00=%d", i00);
	printf ("\ni01=%d", i01);
	printf ("\ni10=%d", i10);
	printf ("\ni11=%d", i11);

	int max = 0;
	if ((i10-i00) > (i11-i01)) max= (i10-i00); else max = (i11-i01);

	int nPts = max + 3;
	printf ("\nnPts=%d", nPts);
	Mesh->x=new Matrice(nPts, 2);
	Mesh->y=new Matrice(nPts, 2);
	Mesh->z=new Matrice(nPts, 2);
	Mesh->CouleurFaces[0] = colorSegment->color->r;
	Mesh->CouleurFaces[1] = colorSegment->color->g;
	Mesh->CouleurFaces[2] = colorSegment->color->b;

	if(symetric==1)
	{
		Mesh->symX = ON;
	}

	Mesh->x->SetElement(0,0, x00);
	Mesh->y->SetElement(0,0, y00);
	Mesh->z->SetElement(0,0, z00);

	Mesh->x->SetElement(0,1, x01);
	Mesh->y->SetElement(0,1, y01);
	Mesh->z->SetElement(0,1, z01);

	int il = i00;
	int ir = i01;
	int _i = 1;
	while ((il <= i10) && (ir <= i11))
	{
		if (side == INT_SIDE) {
			Mesh->x->SetElement(_i, 0, f3d->XInt->Element(nerv, il));
			Mesh->y->SetElement(_i, 0, f3d->YInt->Element(nerv, il));
			Mesh->z->SetElement(_i, 0, f3d->ZInt->Element(nerv, il));

			Mesh->x->SetElement(_i, 1, f3d->XInt->Element(nerv + 1, ir));
			Mesh->y->SetElement(_i, 1, f3d->YInt->Element(nerv + 1, ir));
			Mesh->z->SetElement(_i, 1, f3d->ZInt->Element(nerv + 1, ir));
		} else {
			Mesh->x->SetElement(_i, 0, f3d->XExt->Element(nerv, il));
			Mesh->y->SetElement(_i, 0, f3d->YExt->Element(nerv, il));
			Mesh->z->SetElement(_i, 0, f3d->ZExt->Element(nerv, il));

			Mesh->x->SetElement(_i, 1, f3d->XExt->Element(nerv + 1, ir));
			Mesh->y->SetElement(_i, 1, f3d->YExt->Element(nerv + 1, ir));
			Mesh->z->SetElement(_i, 1, f3d->ZExt->Element(nerv + 1, ir));
		}
		printf ("\n%d %d %d -> (%f, %f, %f) (%f, %f, %f)", _i, il, ir, Mesh->x->Element(_i,0), Mesh->y->Element(_i,0), Mesh->z->Element(_i,0), Mesh->x->Element(_i,1), Mesh->y->Element(_i,1), Mesh->z->Element(_i,1));
		_i++;

		if ((il == i10) && (ir == i11)) break;
		if (il < i10) il++;
		if (ir < i11) ir++;
	}

	Mesh->x->SetElement(_i, 0, x10);
	Mesh->y->SetElement(_i, 0, y10);
	Mesh->z->SetElement(_i, 0, z10);

	Mesh->x->SetElement(_i, 1, x11);
	Mesh->y->SetElement(_i, 1, y11);
	Mesh->z->SetElement(_i, 1, z11);
	printf ("\n(%f, %f, %f) (%f, %f, %f)",x10, y10,z10,x11, y11,z11);

	AjoutMesh(Axe3d, Mesh);
}

void KiteDesign::ajoutMeshesToAxe3d( TAxe *Axe3d, Forme3D* f3d, int side, int symetric){
	printf ("\nKiteDesign::ajoutMeshesToAxe3d");
	ColorSegment* cs = new ColorSegment();
	cs->p00=0.0f;
	cs->p01=0.0f;

	cs->p10=40.0f;
	cs->p11=60.0f;

	cs->nerv=2;
	Color* c = new Color(1.0f, 0.0f, 0.0f);
	cs->color = c;


	ajoutColorSegmentToAxe3d(Axe3d, f3d, cs, side, symetric);
	/*Axe3d->eclairage = ON;
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
	AjoutMesh(Axe3d, Mesh);*/
}