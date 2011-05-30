#pragma once

#include "matrice.h"
#include "fichier.h"
#include "patternsproject.h"

Matrix* getReperPoints(WindPatternsProject* gfd, Matrix* Xd, Matrix* Yd, Matrix* P, 
						int nerv, float deb, int faceDeb, float fin, int faceFin, bool klapanShov);

double calcPolyLength(Matrix* X, Matrix* Y, double xrp, double yrp);

Matrix* getReperPointsPince(WindPatternsProject* gfd,  Matrix* Xd0, Matrix* Yd0, double coeff, Matrix* Xd, Matrix* Yd, Matrix* P, 
								int nerv, float deb, int faceDeb, float fin, int faceFin, bool klapanShov);

int calcRepPointByLength(Matrix* Xd, Matrix* Yd, double l, double* x, double* y);

