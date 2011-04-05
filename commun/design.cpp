//design.cpp


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <conio.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <iostream>


#include "fichier.h"
#include "profil.h"
#include "geom.h"
#include "design.h"

int test() {
	printf ("\n test()");
	return 0;
}

KiteDesignElement::KiteDesignElement(){

}

KiteDesignElement::~KiteDesignElement(){

}

Line::Line(){
	n_points = 0;
}

Line::Line(ifstream& in){
	n_points = 0;
	in >> n_points;
	cout << "Line n_points: "<< n_points << endl;
	pointsNervs = new int[n_points];
	pointsPercents = new double[n_points];
	for (int j = 0; j < n_points; j++) {
		int pn;
		double pp;
		in >> pointsNervs[j] >> pointsPercents[j];
		cout << "Line p "<< pointsNervs[j] << " "<< pointsPercents[j] <<endl;
	}
}

Line::~Line(){

}

void Line::print(){
	printf ("\nline print()");
	printf ("\n np=%d", n_points);
	for (int i = 0; i < n_points; i++){
		printf ("\n pn:%d pp:%f", pointsNervs[i], pointsPercents[i]);
	}
}



KiteDesign::KiteDesign(){
	n_elements = 0;

}

KiteDesign::~KiteDesign(){

}