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

		virtual void addCourbesToAxe(TAxe* axe, FormProjection* fp, int symetric, double dy, double ymult, double ymin) { }
		virtual void addCourbesToAxe3d(TAxe* axe, Form3D* f3d, int side, int symetric) { }
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
		void addCourbesToAxe(TAxe* axe, FormProjection* fp, int symetric, double dy, double ymult, double ymin);
		void addCourbesToAxe3d(TAxe* axe, Form3D* f3d, int side, int symetric);
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

		void addMeshesToAxe3d2d( TAxe *Axe3d, TAxe *Axe2d, Form3D* f3d, float opac, int side, int symetric, double dy, double ymult, double ymin);
		void addMeshesToAxeProjection( TAxe *AxeProjection, FormProjection* fp, int symetric, double dy, double ymult, double ymin);
		
		ColorSegmentsTable* getColorSegmentsTable(int n_profils );
		PanelLinesTable* getPanelLinesTable();
};



void addColorSegmentToAxe3d2d(TAxe *Axe3d, TAxe *Axe2d, Form3D* f3d, ColorSegment* colorSegment, int side, int symetric);

void addColorSegmentToAxeProjection(TAxe *AxeProjection, FormProjection* fp, ColorSegment* colorSegment, int symetric,  double dy, double ymult, double ymin);

void setMeshPoint(TMesh *Mesh, int i, int j, double x, double y, double z);

void getCourbeFromProfilGeom(ProfilGeom* pg, Courbe** courbeExt, Courbe** courbeInt);

void addFormProjectionCourbesToAxe(TAxe* axe, FormProjection* fp, KiteDesign* kd, int symetrie, double dy, int dir);

void addForm3d2dKiteDesign( TAxe *Axe3d,TAxe *Axe2d, Form3D* f3d, FormProjection* fp, KiteDesign* kdExt, int dire, float opacExt, KiteDesign* kdInt, int diri, float opacInt, int mesh, int symetric);

void addFormProjectionKiteDesign(TAxe* axe, FormProjection* fp, KiteDesign* kd, int symetric, double dy, int dir);
