//design.h
#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>

#include "plot.h"
#include "pince.h"
#include "matrice.h"
#include "rasklad.h"
#include "patternsproject.h"

class KiteDesignElement {
public:
	KiteDesignElement();
	virtual ~KiteDesignElement();

	virtual void ajoutCourbesToAxe(TAxe* axe, FormeProjection* fp, double ymin, int symetric, double dy, double ymult, double yoffset) { }
};

class Line : public KiteDesignElement {
public:
	Line();
	Line(std::ifstream& in);
	virtual ~Line();
	int n_points;
	int* pointsNervs;
	double* pointsPercents;
	void print();
	void ajoutCourbesToAxe(TAxe* axe, FormeProjection* fp, double ymin, int symetric, double dy, double ymult, double yoffset);
};

class KiteDesign
{

public:
	KiteDesign();
	virtual ~KiteDesign();

	int n_elements;
	KiteDesignElement** kiteDesignElements;

};