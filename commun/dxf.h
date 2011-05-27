#pragma once

#include <string>

#include "plot.h"
#include "pince.h"
#include "matrice.h"
#include "layout.h"
#include "patternsproject.h"


void writeFichierDXF(char *fileName, TAxe *axe);

bool writeFichierWpa(char *fileName, Form *forme);

void writeFichierPolyDXF(char *fileName, TAxe *axe, TAxe *axe2, int rep, TAxe *axeR, int vent, TAxe *axeC, int num, TAxe *axeT);

void writeFichierPolyDXFDelta(FILE *fid, TAxe *axe, TAxe *axe2, int rep, TAxe *axeR, int vent, TAxe *axeC, int num, TAxe *axeT, double dx, double dy, double dz, int n );

void writeManyFichierPolyDXF(char *fileName, int np, int n, TAxe **axe, TAxe **axe2, int rep, TAxe **axeR, int vent, TAxe **axeC, int num, TAxe **axeT, double* W, double* H);

void writeManyFichierPolyDXF2(char *fileName, int np, int n, TAxe **axe, TAxe **axe2, int rep, TAxe **axeR, int vent, TAxe **axeC, int num, TAxe **axeT, double* W, double* H, int* numncon);

void writeLayoutToDXF(char *fileName, WindPatternsProject *gfd, Layout *layout);

