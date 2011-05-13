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
	a = 0.5f;
}

Color::Color(double _r, double _g, double _b) {
	r = _r;
	g = _g;
	b = _b;
	a = 0.5f;
}

void Color::print() {
	printf ("\n(%f, %f, %f)(a=%f)", r, g, b,a);
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
	//in >> colorR >> colorG >> colorB;
	c = new Color(in);
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
	a = 0.5f;
}

void Line::initColor(Courbe* courbe) {
	courbe->colorSegments[0] = c->r;
	courbe->colorSegments[1] = c->g;
	courbe->colorSegments[2] = c->b;
	courbe->colorSegments[3] = c->a;
}

void Line::addCourbesToAxe(TAxe* axe, FormProjection* fp, int symetric, double dy, double ymult, double ymin)
{
	Courbe* courbe = new Courbe("CourbeLine");
	initColor(courbe);
	courbe->points = OFF;
	if (symetric) courbe->symX = ON;
	courbe->pts = new Matrix(n_points, 2);

	for (int i = 0; i < n_points; i++) {
		int nerv = pointsNervs[i];
		//printf ("\ni=%d nerv=%d", i, nerv);
		double percent = pointsPercents[i];

		//double x0 = fp->X->Element(i, 0);
		//double y0 = fp->Y->Element(i, 0);
		double x0 = fp->X->Element(nerv, 0);
		double y0 = fp->Y->Element(nerv, 0);

		//double x1 = fp->X->Element(i, 1);
		//double y1 = fp->Y->Element(i, 1);
		double x1 = fp->X->Element(nerv, 1);
		double y1 = fp->Y->Element(nerv, 1);

		//printf ("\n(%f,%f) (%f,%f)", x0, y0, x1, y1);
		double xp = x0 + percent * 0.01 * (x1 - x0);
		double yp = y0 + percent * 0.01 * (y1 - y0);
		//printf ("\n(%f,%f)", xp, yp);

		courbe->pts->SetElement(i, 0, xp);
		courbe->pts->SetElement(i, 1, dy + ymult*yp - ymin);
	}

	addCourbe(axe, courbe);
}

void Line::addCourbesToAxe3d(TAxe* axe, Form3D* f3d, int side, int symetric) {
	Courbe* courbe = new Courbe("CourbeLine");
	initColor(courbe);
	courbe->points = OFF;
	if (symetric) courbe->symX = ON;
	courbe->pts = new Matrix(n_points, 2);
	double x=0, y=0, z=0;
	for (int i = 0; i < n_points; i++) {
		int nerv = pointsNervs[i];
		double perc = pointsPercents[i];
		getPoint3dFormByPosNerv(f3d, nerv, side, perc, &x, &y, &z);
		courbe->pts->SetElement(i, 0, x);
		courbe->pts->SetElement(i, 1, y);
		courbe->pts->SetElement(i, 2, z);
	}
	addCourbe(axe, courbe);
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

KiteDesign::KiteDesign(std::ifstream& in){
  //TODO:
}

KiteDesign::~KiteDesign(){
}

ColorSegmentsTable::ColorSegmentsTable(){
}

ColorSegmentsTable::ColorSegmentsTable(int n){
	table.resize(n);
}


ColorSegment::ColorSegment(){
	nerv = 0;
	p00 = 0.0f;
	p01 = 0.0f;
	p10 = 0.0f;
	p11 = 0.0f;
}

void ColorSegment::print() {
	printf ("\nColorSegment::print ");
	printf (" %d [%f, %f]-[%f, %f]", nerv, p00, p01, p10, p11);
}


ColorSegment::~ColorSegment(){
}

int indexBeforePosProf(Matrix* Prof, double pos) {
	int i = 0;
	int n = Prof->GetLignes() - 1;
	while (i < n) {
		if (Prof->Element(i+1, 0) > pos) break;
		i++;
	}
	return i;
}

int indexAfterPosProf(Matrix* Prof, double pos) {
	int i = 0;
	int n = Prof->GetLignes()-1;
	while (Prof->Element(i, 0) < pos) {
		i++;
		if ( i == n ) break;
	}
	return i;
}


PanelLine::PanelLine(){
}

PanelLine::PanelLine(int _nerv1, double _pos1, int _nerv2, double _pos2)
	: nerv1(_nerv1), pos1(_pos1), nerv2(_nerv2), pos2(_pos2)
{
}



PanelLine::~PanelLine(){
}

PanelLinesTable::PanelLinesTable(){
}

bool PanelLinePredicate(const PanelLine* pl1, const PanelLine* pl2)
{
	if (pl1->pos1 == pl2->pos1)
		return (pl1->pos2 < pl2->pos2);
	else 
		return (pl1->pos1 < pl2->pos1);
}

PanelLinesTable::PanelLinesTable(int n_profils) {
	printf ("\nPanelLinesTable::PanelLinesTable(%d)", n_profils);
	table.resize(n_profils);
}

void PanelLinesTable::addPanelLine(int nerv, PanelLine* pl){
	//printf ("\nPanelLinesTable::addPanelLine(%d, ...)", nerv);
	//printf ("\n v.size()=%d", table[nerv].size());
	vector<PanelLine*>::iterator it; 
	it=table[nerv].begin();
	while ((it < table[nerv].end()) && ((*it)->pos1 <= pl->pos1) && ((*it)->pos2 <= pl->pos2))
	{
	//	printf ("\nit++ (%f <= %f) (%f <= %f) ", (*it)->pos1,  pl->pos1, (*it)->pos2, pl->pos2);
		it++;
	}
	table[nerv].insert(it, pl);
	//printf ("\nPanelLinesTable::...addPanelLine(%d) v.size()=%d", nerv, table[nerv].size());
}

void PanelLinesTable::sortNervsVectors(){
	for (int i = 0; i < table.size(); i++) {
		std::vector<PanelLine*> vpl = table[i];
		//std::sort(vpl.begin(), vpl.end(), PanelLinePredicate);
	}
}
PanelLinesTable* KiteDesign::getPanelLinesTable() {
	printf ("\n KiteDesign::getPanelLinesTable()");
    PanelLinesTable* panelLinesTable = new PanelLinesTable(200);

    for (int i = 0; i < n_elements; i++) {
		//printf ("\n KiteDesignElement i=%d/%d", i, n_elements);
        KiteDesignElement* kde = kiteDesignElements[i];
        Line* line = (Line*) kde;

        for (int j = 0; j < line->n_points - 1; j++){
		  //printf ("\n nerv j=%d/%d", j, line->n_points);
          int nerv1 = line->pointsNervs[j];
          double perc1 = line->pointsPercents[j];
          int nerv2 = line->pointsNervs[j+1];
          double perc2 = line->pointsPercents[j+1];
		  //printf ("\n new PanelLine(%d, %f, %d, %f);",nerv1, perc1, nerv2, perc2);
          PanelLine* pl = new PanelLine(nerv1, perc1, nerv2, perc2);
          panelLinesTable->addPanelLine(nerv1, pl);
        }
    }
	//panelLinesTable->sortNervsVectors();
	printf ("\n...KiteDesign::getPanelLinesTable()");
    return panelLinesTable;
}

ColorSegmentsTable* KiteDesign::getColorSegmentsTable(int n_profils ){
    printf ("\nKiteDesign::getColorSegmentsTable()");
    PanelLinesTable* panelLinesTable = getPanelLinesTable();
    ColorSegmentsTable* colorSegmentsTable = new ColorSegmentsTable(n_profils);
	//printf ("\n nbprofils=%d colorTable->table.size()=%d", n_profils,colorTable->table.size());
	for (int nerv = 0; (nerv < n_profils) && (nerv < colorTable->table.size()); nerv++) {
		//printf ("\n===================================================");
		//printf ("\n nerv=%d", nerv);

		vector<PanelLine*> v = panelLinesTable->table[nerv];
		double prev_p0 = 0.0f;
		double prev_p1 = 0.0f;
		vector<Color*> vc = colorTable->table[nerv];
		//printf ("\n colorTable->table[%d] vc.size()=%d", nerv, vc.size());
		//printf ("\n panelLinesTable->table[%d] v.size()=%d", nerv, v.size());

		for (int i = 0; i < v.size(); i++) {
			if ( i >= vc.size()) {
				//printf ("\n colors for %d : vc.size(): %d", nerv, vc.size());
				break;
			}			
			ColorSegment* colorSegment = new ColorSegment();
			colorSegment->nerv=nerv;
			colorSegment->p00 = prev_p0;
			colorSegment->p01 = prev_p1;

			colorSegment->p10 = v[i]->pos1;
			colorSegment->p11 = v[i]->pos2;
			//printf ("\nv[%d]->pos1=%f pos2=%f", i, v[i]->pos1, v[i]->pos2);
			colorSegment->color = vc[i];

			//printf ("\n push cs (%f, %f)-(%f, %f)", colorSegment->p00, colorSegment->p01, colorSegment->p10, colorSegment->p11);
			colorSegmentsTable->table[nerv].push_back(colorSegment);
			prev_p0 = v[i]->pos1;
			prev_p1 = v[i]->pos2;
		
		}
		if ( vc.size() >= (v.size()+1))
		{
			if ((prev_p0 != 100.0f) || (prev_p1 != 100.0f)) {
				ColorSegment* colorSegment = new ColorSegment();
				colorSegment->p00 = prev_p0;
				colorSegment->p01 = prev_p1;
				colorSegment->nerv=nerv;
				colorSegment->p10 = 100.0f;
				colorSegment->p11 = 100.0f;
				colorSegment->color = vc[v.size()];
				colorSegmentsTable->table[nerv].push_back(colorSegment);
			}
		}
    }
	printf ("\n...KiteDesign::getColorSegmentsTable()");
    return colorSegmentsTable;

}


void setMeshPoint(TMesh *Mesh, int i, int j, double x, double y, double z) {
	Mesh->x->SetElement(i, j, x);
	Mesh->y->SetElement(i, j, y);
	Mesh->z->SetElement(i, j, z);
}

void setMeshPoint3d2d(TMesh *mesh3,TMesh *mesh2, int i, int j, double x, double y, double z, double dy, double ymult, double ymin) {
	mesh3->x->SetElement(i, j, x);
	mesh3->y->SetElement(i, j, y);
	mesh3->z->SetElement(i, j, z);

	mesh2->x->SetElement(i, j, x);
	mesh2->y->SetElement(i, j, dy + ymult*z - ymin);
	mesh2->z->SetElement(i, j, 0.0f);

}



void addColorSegmentToAxe3d2d(TAxe *Axe3d, TAxe *Axe2d, Form3D* f3d, ColorSegment* colorSegment, int side, int symetric, double dy, double ymult, double ymin) {
	//printf ("\ndesign:addColorSegmentToAxe3d()");
	Axe3d->eclairage = ON;
	TMesh *mesh3d;
	TMesh *mesh2d;

	mesh3d = createMesh();
	mesh2d = createMesh();

	mesh3d->segments=OFF; mesh3d->faces=ON;
	mesh2d->segments=OFF; mesh2d->faces=ON;

	if (side == INT_SIDE) mesh3d->InvNormales=ON;
	mesh3d->side = side;
	mesh2d->side = side;

	int nerv = colorSegment->nerv;
	//printf ("\ndesign: addColorSegmentToAxe3d() nerv=%d", nerv);
	double x00, y00, z00, x01, y01, z01, x10, y10, z10, x11, y11, z11;

	getPoint3dFormByPosNerv(f3d, nerv, side, colorSegment->p00, &x00, &y00, &z00) ;
	getPoint3dFormByPosNerv(f3d, nerv+1, side, colorSegment->p01, &x01, &y01, &z01) ;

	getPoint3dFormByPosNerv(f3d, nerv, side, colorSegment->p10, &x10, &y10, &z10) ;
	getPoint3dFormByPosNerv(f3d, nerv+1, side, colorSegment->p11, &x11, &y11, &z11) ;

	Matrix* mProf;
	if (side == INT_SIDE) mProf = f3d->forme->getIntProf(nerv,false); else mProf = f3d->forme->getExtProf(nerv,false);

	mesh3d->colorFaces[0] = colorSegment->color->r;
	mesh3d->colorFaces[1] = colorSegment->color->g;
	mesh3d->colorFaces[2] = colorSegment->color->b;
	mesh3d->colorFaces[3] = colorSegment->color->a;

	mesh2d->colorFaces[0] = colorSegment->color->r;
	mesh2d->colorFaces[1] = colorSegment->color->g;
	mesh2d->colorFaces[2] = colorSegment->color->b;
	mesh2d->colorFaces[3] = colorSegment->color->a;

	if(symetric == 1)
	{
		mesh3d->symX = ON;
		mesh2d->symX = ON;
	}
	//---------------------------------------------
	//printf ("\n colorSegment: ");
	//printf ("\n nerv=%d", nerv);

	//printf ("\n %f(%d) %f(%d)", colorSegment->p00, nerv, colorSegment->p01, nerv+1);
	Point2d* pnt00 = new Point2d(nerv, colorSegment->p00);
	Point2d* pnt01 = new Point2d(nerv+1, colorSegment->p01);
	//printf ("\n %f(%d) %f(%d)", colorSegment->p10, nerv, colorSegment->p11, nerv+1);
	Point2d* pnt10 = new Point2d(nerv, colorSegment->p10);
	Point2d* pnt11 = new Point2d(nerv+1, colorSegment->p11);

	Segment2d* sl0 = new Segment2d(pnt00, pnt01);
	Segment2d* sl1 = new Segment2d(pnt10, pnt11);
	Segment2d* se0 = new Segment2d(pnt00, pnt10);
	Segment2d* se1 = new Segment2d(pnt01, pnt11);

	std::vector<double> pcse;
	pcse.push_back(colorSegment->p00);
	pcse.push_back(colorSegment->p01);
	pcse.push_back(colorSegment->p10);
	pcse.push_back(colorSegment->p11);

	std::sort(pcse.begin(), pcse.end());

	double minp = pcse[0];
	//printf ("\npcse[0]=%f", pcse[0]);
	int iAfterMin = indexAfterPosProf(mProf,  minp);
	//printf ("\niAfterMin=%d -- %f", iAfterMin, mProf->Element(iAfterMin, 0));
	double maxp = pcse[3];
	//printf ("\npcse[3]=%f", pcse[3]);
	int iBeforeMax = indexBeforePosProf(mProf,  maxp);
	//printf ("\niBeforeMax=%d -- %f", iBeforeMax, mProf->Element(iBeforeMax, 0));

	std::vector<double> pcs;
	pcs.push_back(colorSegment->p00);
	pcs.push_back(colorSegment->p01);
	pcs.push_back(colorSegment->p10);
	pcs.push_back(colorSegment->p11);
	int i = 0;
	
	for (i = iAfterMin; i <= iBeforeMax; i++) {
		pcs.push_back(mProf->Element(i, 0));
	}
	std::sort(pcs.begin(), pcs.end());
	
	/*for (int ii = 0; ii < pcs.size(); ii++) {
		printf ("\nPPCCSS pcs %d - %f",ii, pcs[ii]);
	}*/

	int k = 0;
	double prev = -1000;

	for (i = 0; i < pcs.size(); i++) {
		if (pcs[i] != prev) {
		
			k++;
			prev = pcs[i];
		}
	}

	int nPts = k;
	//printf ("\nnPts=%d", nPts);

	mesh3d->x=new Matrix(nPts, 2);
	mesh3d->y=new Matrix(nPts, 2);
	mesh3d->z=new Matrix(nPts, 2);

	mesh2d->x=new Matrix(nPts, 2);
	mesh2d->y=new Matrix(nPts, 2);
	mesh2d->z=new Matrix(nPts, 2);

	prev = -1000;
	k = 0;
	for (i = 0; i < pcs.size(); i++) {
		//printf ("\n==%d=======================================================", i);
		//printf ("\npcs[%d]=%f", i, pcs[i]);
		if (pcs[i] != prev) {
			//printf ("\n!=prev go");
			double p = pcs[i];
			Point2d* p1 = new Point2d(nerv, p);
			Point2d* p2 = new Point2d(nerv+1, p);

			Segment2d* s = new Segment2d(p1, p2);
			//printf ("\n LINE_0");
			ResultIntersect2d* rsl0 = intersectSegments2d (s, sl0);
			//printf ("\n rsl0->type=%d", rsl0->type);

			//printf ("\n LINE_1");
			ResultIntersect2d* rsl1 = intersectSegments2d (s, sl1);
			//printf ("\n rsl1->type=%d", rsl1->type);

			//printf ("\n E_0");
			ResultIntersect2d* rse0 = intersectSegments2d (s, se0);
			//printf ("\n rse0->type=%d", rse0->type);

			//printf ("\n E_1");
			ResultIntersect2d* rse1 = intersectSegments2d (s, se1);
			//printf ("\n rse1->type=%d", rse1->type);

			if ((rsl0->type == SEG) || (rsl1->type == SEG))  {
				//printf ("\nSEG");
				if (rsl0->type == SEG) {
					//printf ("\nrsl0->type == SEG");
					//printf ("\nset %d (%f, %f, %f) - (%f, %f, %f) ",k, x00, y00, z00, x01, y01, z01);
					setMeshPoint3d2d(mesh3d, mesh2d, k, 0, x00, y00, z00, dy, ymult, ymin);
					setMeshPoint3d2d(mesh3d, mesh2d, k, 1, x01, y01, z01, dy, ymult, ymin);
				} else
				if (rsl1->type == SEG) {
					//printf ("\nrsl1->type == SEG");
					//printf ("\nset %d (%f, %f, %f) - (%f, %f, %f) ",k,  x10, y10, z10, x11, y11, z11);
					setMeshPoint3d2d(mesh3d, mesh2d, k, 0, x10, y10, z10, dy, ymult, ymin);
					setMeshPoint3d2d(mesh3d, mesh2d, k, 1, x11, y11, z11, dy, ymult, ymin);
				}
			} else {
				std::vector<double> res;
				vector<double>::iterator it;
				if (rsl0->type == ONE) {
					double xi = rsl0->p1->x;
					//printf ("\nrsl0->type == ONE : %f", xi);
					it = std::find(res.begin(),res.end(), xi);

					//printf ("\n[%d]", res.size());
					bool inres = false;
					for (int _ir = 0; _ir < res.size(); _ir++) {
						//printf ("%f ", res[_ir]);
						if (fabs(res[_ir]-xi)<0.00001f) inres = true;
					}

					if (!inres) {
						res.push_back(xi);
						//printf ("\n push!");	
					} else {
						//printf ("\n also in res");	
					}
				}
				if (rsl1->type == ONE) {
					double xi = rsl1->p1->x;
					//printf ("\nrsl1->type == ONE : %f", xi);
					it = std::find(res.begin(),res.end(), xi);

					//printf ("\n[%d]", res.size());
					bool inres = false;
					for (int _ir = 0; _ir < res.size(); _ir++) {
						//printf ("%f ", res[_ir]);
						if (fabs(res[_ir]-xi)<0.00001f) inres = true;
					}

					if (!inres) {
						res.push_back(xi);
						//printf ("\n push!");	
					} else {
						//printf ("\n also in res");	
					}

				}
				if (rse0->type == ONE) {
					double xi = rse0->p1->x;
					it = std::find(res.begin(),res.end(), xi);
					//printf ("\nrse0->type == ONE : %f", xi);

					//printf ("\n[%d]", res.size());
					bool inres = false;
					for (int _ir = 0; _ir < res.size(); _ir++) {
						//printf ("%f ", res[_ir]);
						if (fabs(res[_ir]-xi)<0.00001f) inres = true;
					}


					if (!inres) {
						res.push_back(xi);
						//printf ("\n push!");	
					} else {
						//printf ("\n also in res");	
					}
				}
				if (rse1->type == ONE) {
					double xi = rse1->p1->x;
					it = std::find(res.begin(),res.end(), xi);

					//printf ("\n[%d]", res.size());
					bool inres = false;
					for (int _ir = 0; _ir < res.size(); _ir++) {
						//printf ("%f ", res[_ir]);
						if (fabs(res[_ir]-xi)<0.00001f) inres = true;
					}


					//printf ("\nrse1->type == ONE : %f", rse1->p1->x);
					if (!inres) {
						res.push_back(xi);
						//printf ("\n push!");	
					} else {
						//printf ("\n also in res");	
					}

				}
				std::sort(res.begin(), res.end());
				int rsize = res.size();
				//printf ("\n*** res.size()=%d ", rsize);
				/*for (int _ir = 0; _ir < res.size(); _ir++) {
					printf ("%f ", res[_ir]);
				}*/

				if ((rsize == 2) || (rsize == 1)) {
					double x0=0.0f, y0=0.0f, z0=0.0f, x1=0.0f, y1=0.0f, z1=0.0f;
					getPoint3dFormByPosDNerv(f3d, res[0], side, pcs[i], &x0, &y0, &z0) ;
					setMeshPoint3d2d(mesh3d, mesh2d, k, 0, x0, y0, z0, dy, ymult, ymin);	
					//printf ("\nset %d (%f, %f, %f)", k, x0, y0, z0);
					if (rsize == 2) {
						getPoint3dFormByPosDNerv(f3d, res[1], side, pcs[i], &x1, &y1, &z1) ;
						setMeshPoint3d2d(mesh3d, mesh2d, k, 1, x1, y1, z1, dy, ymult, ymin);
						//printf(" - (%f, %f, %f) ", x1, y1, z1);
					} else {
						setMeshPoint3d2d(mesh3d, mesh2d, k, 1, x0, y0, z0, dy, ymult, ymin);
						//printf(" - (%f, %f, %f) ", x0, y0, z0);
					}
				} 
			}
			k++;
			prev = pcs[i];
		}
	}
	//=============================================
	addMesh(Axe3d, mesh3d);
	addMesh(Axe2d, mesh2d);
}
void KiteDesign::addMeshesToAxe3d2d( TAxe *Axe3d, TAxe *Axe2d, Form3D* f3d, float opac, int side, int symetric, double dy, double ymult, double ymin){
	printf ("\nKiteDesign::addMeshesToAxe3d");
	/*
	// test 
	ColorSegment* cs = new ColorSegment();
	cs->p00 = 0.0f;
	cs->p01 = 50.0f;
	cs->nerv=1;
	cs->p10 = 10.0f;
	cs->p11 = 60.0f;
	cs->color = new Color(1.0f, 0.0f, 0.0f);
	addColorSegmentToAxe3d(Axe3d, f3d, cs, side, symetric);
	*/
	ColorSegmentsTable* cst = getColorSegmentsTable( f3d->forme->m_nbProfils );
	for (int nerv = 0; nerv < f3d->forme->m_nbProfils; nerv++) {
	//for (int nerv = 3; nerv < 4; nerv++) {
		vector<ColorSegment*> vcs = cst->table[nerv];
		//printf ("\nKiteDesign::addMeshesToAxe3d nerv=%d vcs.size()=%d", nerv, vcs.size());
		for (int i = 0; i < vcs.size(); i++) {
		//for (int i = 1; i < 2; i++) {
		//	printf ("\nKiteDesign::addMeshesToAxe3d colorSegment i=%d", i);
			ColorSegment* cs = vcs[i];
			cs->color->a = opac;
			//printf ("\nColorSegmentKiteDesign::addMeshesToAxe3d colorSegment i=%d", i);
			//cs->print();
			addColorSegmentToAxe3d2d(Axe3d, Axe2d, f3d, cs, side, symetric, dy, ymult, ymin);
		}
    }
	printf ("\n...KiteDesign::addMeshesToAxe3d");
}
 
double calcYminFp(FormProjection* fp, double ymult) {
	int n = fp->X->GetLignes();
	int m = fp->X->GetColonnes();

	double ymin = 1000000;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			double y = fp->Y->Element(i, j);
			if (y*ymult < ymin) ymin = y*ymult;
		}
	}

	return ymin;
}


void addForm3d2dKiteDesign( TAxe *Axe3d,TAxe *Axe2d, Form3D* f3d, FormProjection* fp, KiteDesign* kdExt, int dire, float opacExt, KiteDesign* kdInt, int diri, float opacInt, int mesh, int symetric) 
{
	double ymulte = 1.0f, ymulti = 1.0f;
	
	if (dire == DIR_NOSE_DOWN) ymulte=-1;
	if (diri == DIR_NOSE_DOWN) ymulti=-1;
	
	double ymine = calcYminFp( fp, ymulte );
	double ymini = calcYminFp( fp, ymulti );
	//------------------------------------------------

	// KiteDesign Int Side
	for (int i = 0; i < kdInt->n_elements; i++) {
		KiteDesignElement* kde = kdInt->kiteDesignElements[i];
		kde -> addCourbesToAxe3d(Axe3d, f3d, INT_SIDE, symetric);
		kde -> addCourbesToAxe( Axe2d, fp, symetric, 2.0f, ymulti, ymini );
	}
	kdInt->addMeshesToAxe3d2d( Axe3d, Axe2d, f3d, opacInt, INT_SIDE, symetric, 0.0f, ymulti, ymini );

	// KiteDesign Ext Side
	for (int i = 0; i < kdExt->n_elements; i++) {
		KiteDesignElement* kde = kdExt->kiteDesignElements[i];
		kde -> addCourbesToAxe3d(Axe3d, f3d, EXT_SIDE, symetric);
		kde -> addCourbesToAxe( Axe2d, fp, symetric, 0.0f, ymulte, ymine );
	}
	ColorTable* ct = kdExt->colorTable;
	kdExt->addMeshesToAxe3d2d( Axe3d, Axe2d, f3d, opacExt, EXT_SIDE, symetric, 2.0f, ymulte, ymine );

}

void getCourbeFromProfilGeom(ProfilGeom* pg, Courbe** courbeExt, Courbe** courbeInt){
	*courbeExt = new Courbe("ProfileExt");
    (*courbeExt)->points = OFF;
    (*courbeExt)->symX = OFF;

	int nExt = pg->ExtProf->GetLignes();
	int nInt = pg->IntProf->GetLignes();
	(*courbeExt)->pts = new Matrix(nExt,2);
	for (int i = 0; i < nExt; i++) {
		(*courbeExt)->pts->SetElement(i, 0, pg->ExtProf->Element(i, 0));
		(*courbeExt)->pts->SetElement(i, 1, pg->ExtProf->Element(i, 1));
	}
	(*courbeInt) = new Courbe("ProfileInt");
    (*courbeInt)->points = OFF;
    (*courbeInt)->symX = OFF;
	(*courbeInt)->pts = new Matrix(nInt,2);
	for (int i = 0; i < nInt; i++) {
		(*courbeInt)->pts->SetElement( i, 0, pg->IntProf->Element(i, 0));
		(*courbeInt)->pts->SetElement( i, 1, pg->IntProf->Element(i, 1));
	}
}


void addFormProjectionCourbesToAxe(TAxe* axe, FormProjection* fp, KiteDesign* kd, int symetric, double dy, int dir) {
	//printf ("\n addFormProjectionCourbesToAxe()");
	double ymult = 1;
	if (dir == DIR_NOSE_DOWN) ymult=-1;
	int n = fp->X->GetLignes();
	int m = fp->X->GetColonnes();
	//printf ("\n n=%d, m=%d", n, m);
	
	Courbe* courbe0 = new Courbe("CourbeUp");
	courbe0->points = OFF;
	if (symetric) courbe0->symX = ON;
	courbe0->pts = new Matrix(n, 2);

	Courbe* courbe1 = new Courbe("CourbeDown");
	courbe1->points = OFF;
	if (symetric) courbe1->symX = ON;
	courbe1->pts = new Matrix(n, 2);
	

	double ymin = 1000000;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			double y = fp->Y->Element(i, j);
			if (y*ymult < ymin) ymin = y*ymult;
		}
	}

	for (int i = 0; i < n; i++) {
		//printf ("\n i=%d", i);
		Courbe* courbe = new Courbe("Courbe");
		courbe->points = OFF;
		courbe->symX = OFF;
		if (symetric) courbe->symX = ON;
		courbe0->pts->SetElement(i, 0, fp->X->Element(i, 0));
		courbe0->pts->SetElement(i, 1, dy + ymult*fp->Y->Element(i, 0)- ymin);
		courbe->pts = new Matrix(m, 2);
		for (int j = 0; j < m; j++) {
			double x = fp->X->Element(i, j);
			courbe->pts->SetElement(j, 0, x);
			double y = fp->Y->Element(i, j);
			courbe->pts->SetElement(j, 1, dy + (ymult*y - ymin));
			//printf ("\n (%d, %d)    %f, %f", i, j, fp->X->Element(i, j), fp->Y->Element(i, j));
		}
		courbe1->pts->SetElement(i, 0, fp->X->Element(i, m-1) );
		courbe1->pts->SetElement(i, 1, dy + ymult*fp->Y->Element(i, m-1)-ymin);

		//printf ("\n addCourbe()");
		addCourbe(axe, courbe);
	}
	
	addCourbe(axe, courbe0);
	addCourbe(axe, courbe1);
}
