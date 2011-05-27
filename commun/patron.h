#pragma once

#include <string>

#include "plot.h"
#include "pince.h"
#include "matrice.h"
#include "layout.h"
#include "patternsproject.h"


void GenerateCourbe(WindPatternsProject* gfd, Matrix *Xd1,
					Matrix *Yd1, Matrix *P1,
					int nerv1, double deb1, int faceDeb1, double fin1, int faceFin1,
					Matrix *Xd2, Matrix *Yd2, Matrix *P2,
					int nerv2, double deb2, int faceDeb2, double fin2, int faceFin2, char *text,
			        TAxe **AxePatronP, TAxe **AxePatronDXFP, TAxe **AxePatronTextDXFP, TAxe **AxeMarginDXFP, TAxe **AxeCercleDXFP, TAxe **AxeRepDXFP, int Ventilation,
                    double marge1, double marge2, double margeDeb, double margeFin,bool makeRep, bool debug,
					bool isPince=false, Matrix *Xd01=0, Matrix *Yd01=0,double coeff1=0, Matrix *Xd02=0, Matrix *Yd02=0,double coeff2=0);
