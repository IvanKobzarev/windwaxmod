

#pragma warning(disable:4514)
#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>

#include "matrice.h"
#include "patternsproject.h"
#include "geom.h"
#include "layout.h"
#define sqr(f1) ((f1)*(f1))
#define pi 3.141592675f

#ifndef DEBUG
#define DEBUG false
#endif


WindPatternsProject::~WindPatternsProject()
{

}

void WindPatternsProject::print()
{
    printf ("\nWindPatternsProject::print() Name=[%s]", name);
    printf ("\n[%s] [%s] [%s]", fileNameDiagNerv, fileNameRepPoints, fileNameVentHoles);
    printf ("\n XMastab=%f", XMashtab);
}
