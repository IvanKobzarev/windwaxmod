//design.h
#pragma once

#include <string>

#include "plot.h"
#include "pince.h"
#include "matrice.h"
#include "rasklad.h"
#include "patternsproject.h"

int test();

class KiteDesignElement {
public:
	KiteDesignElement();
	virtual ~KiteDesignElement();

};

class Line : public KiteDesignElement {
public:
	Line();
	Line(ifstream& in);
	virtual ~Line();
	int n_points;
	int* pointsNervs;
	double* pointsPercents;
	void print();
};

class KiteDesign
{

public:
	KiteDesign();
	virtual ~KiteDesign();

	int n_elements;
	KiteDesignElement** kiteDesignElements;

};