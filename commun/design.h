//design.h
#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <algorithm>

#include "plot.h"
#include "pince.h"
#include "matrice.h"
#include "rasklad.h"
#include "patternsproject.h"

class Color;

class KiteDesignElement {
public:
	KiteDesignElement();
	virtual ~KiteDesignElement();

	virtual void ajoutCourbesToAxe(TAxe* axe, FormeProjection* fp, int symetric, double dy, double ymult, double ymin) { }
	virtual void ajoutCourbesToAxe3d(TAxe* axe, Forme3D* f3d, int side, int symetric) { }
};

class Line : public KiteDesignElement {
public:
	Line();
	Line(std::ifstream& in);
	virtual ~Line();
	int n_points;
	int* pointsNervs;
	double* pointsPercents;
	//double colorR, colorG, colorB;
	Color* c;
	void print();
	void ajoutCourbesToAxe(TAxe* axe, FormeProjection* fp, int symetric, double dy, double ymult, double ymin);
	void ajoutCourbesToAxe3d(TAxe* axe, Forme3D* f3d, int side, int symetric);
private:
	void initColor(Courbe* courbe);
};

class Color {
public:
	Color();
	Color(std::ifstream& in);
	Color(double r, double g, double b);
	virtual ~Color();

	double r, g, b, a;
	void print();
};

class ColorTable {
public:
	ColorTable();
	ColorTable(std::ifstream& in);
	virtual ~ColorTable();
	std::vector<std::vector <Color*> > table;

	void print();
};

class ColorSegment
{
public:
    ColorSegment();
    virtual ~ColorSegment();
	void print();

    int nerv;
    Color* color;
    double p00, p01, p10, p11;
};

class ColorSegmentsTable
{
public:
    std::vector< std::vector <ColorSegment*> > table;

	ColorSegmentsTable();
	ColorSegmentsTable(int n);
};


class PanelLine
{
public:

	PanelLine();
	PanelLine(int _nerv1, double _pos1, int _nerv2, double _pos2);
    virtual ~PanelLine();
	
	int nerv1;
	int nerv2;
    double pos1;
    double pos2;
};

class PanelLinesTable
{
public:
    std::vector< std::vector<PanelLine*>> table;

    PanelLinesTable();
	PanelLinesTable(int n_profils);
	void addPanelLine(int nerv, PanelLine* pl);
	void sortNervsVectors();
};

class KiteDesign
{

public:

	KiteDesign();
    KiteDesign(std::ifstream& in);

	virtual ~KiteDesign();

	int n_elements;
	KiteDesignElement** kiteDesignElements;
	ColorTable* colorTable;

	void ajoutMeshesToAxe3d( TAxe *Axe3d, Forme3D* f3d, float opac, int side, int symetric);
	void ajoutMeshesToAxeProjection( TAxe *AxeProjection, FormeProjection* fp, int symetric, double dy, double ymult, double ymin);
	
    ColorSegmentsTable* getColorSegmentsTable(int n_profils );
    PanelLinesTable* getPanelLinesTable();
};


void ajoutColorSegmentToAxe3d(TAxe *Axe3d, Forme3D* f3d, ColorSegment* colorSegment, int side, int symetric);

void ajoutColorSegmentToAxeProjection(TAxe *AxeProjection, FormeProjection* fp, ColorSegment* colorSegment, int symetric,  double dy, double ymult, double ymin);

void setMeshPoint(TMesh *Mesh, int i, int j, double x, double y, double z);

void getCourbeFromProfilGeom(ProfilGeom* pg, Courbe** courbeExt, Courbe** courbeInt);

void ajoutFormeProjectionCourbesToAxe(TAxe* axe, FormeProjection* fp, KiteDesign* kd, int symetrie, double dy, int dir);

void AjoutForme3DKiteDesign( TAxe *Axe3d, Forme3D* f3d, KiteDesign* kdExt, float transpExt, KiteDesign* kdInt, float transpInt, int mesh, int symetric);

void AjoutFormeProjectionKiteDesign(TAxe* axe, FormeProjection* fp, KiteDesign* kd, int symetric, double dy, int dir);
/* vector< vector<int> > vI2Matrix(3, vector<int>(2,0));
   vI2Matrix[0][0] = 0;
   vI2Matrix[0][1] = 1;
   vI2Matrix[1][0] = 10;
   vI2Matrix[1][1] = 11;
   vI2Matrix[2][0] = 20;
   vI2Matrix[2][1] = 21; */
