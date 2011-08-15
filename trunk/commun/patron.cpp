#include "patron.h"

#include "matrice.h"
#include "geom.h"
#include "plot.h"
#include "patternsproject.h"

#define sqr(l1) ((l1)*(l1))
#define pi	3.141592675f

#ifndef DEG2RAD
#define DEG2RAD	(3.141592675f/180.0f)
#endif


void GenerateCourbe(WindPatternsProject* gfd, Matrix *Xd1,
					Matrix *Yd1, Matrix *P1,
					int nerv1, double deb1, int faceDeb1, double fin1, int faceFin1,
					Matrix *Xd2, Matrix *Yd2, Matrix *P2,
					int nerv2, double deb2, int faceDeb2, double fin2, int faceFin2, char *text,
			        TAxe **AxePatronP, TAxe **AxePatronDXFP, TAxe **AxePatronTextDXFP, TAxe **AxeMarginDXFP, TAxe **AxeCercleDXFP, TAxe **AxeRepDXFP, int Ventilation,
                    double marge1, double marge2, double margeDeb, double margeFin,bool makeRep, bool debug,
					bool isPince, Matrix *Xd01, Matrix *Yd01,double coeff1, Matrix *Xd02, Matrix *Yd02,double coeff2)
{
    //printf ("\n GenerateCourbe()");
    double Marge[2]={marge1, marge2};
    *AxePatronP = createAxe(gfd->windowPatron);
    *AxePatronDXFP = createAxe(gfd->windowPatron);
    *AxeRepDXFP = createAxe(gfd->windowPatron);
    *AxePatronTextDXFP = createAxe(gfd->windowPatron);
    *AxeMarginDXFP = createAxe(gfd->windowPatron);
    *AxeCercleDXFP = createAxe(gfd->windowPatron);
    //printf ("\n GenerateCourbe def 1");
    int ajCoAr = 0, ajCoMaAr = 0, ajCoAv = 0, ajCoMaAv = 0;
    int ajCo1[2] = {0, 0};
    int ajCo2[2] = {0, 0};
    int nerv[2] = {nerv1, nerv2};
    if (nerv1 == -1) nerv[0]=0;
    int face[2] = {faceFin1, faceFin2};
    int FaceDeb[2]= {faceDeb1, faceDeb2};
    int FaceFin[2]= {faceFin1, faceFin2};
    double Deb[2] = {deb1, deb2};
    double Fin[2] = {fin1, fin2};
    //printf ("\n GenerateCourbe def 2");
    Matrix *distance;
    int n;
    double xText, yText;
    char texte[100], texteExtInt[3] = "EI";
    Matrix *interpXSuspente, *interpYSuspente;
    double posSuspente, xSuspente[2][205], ySuspente[2][205];
    //printf ("\n GenerateCourbe def 3");
    Matrix * Xd[2] = {Xd1, Xd2};

    //printf ("\n Xd1->GetLignes()=%d", Xd1->GetLignes());
    //printf ("\n Xd2->GetLignes()=%d", Xd2->GetLignes());

	Matrix * Xd0[2] = {Xd01, Xd02};
    Matrix * Yd[2] = {Yd1, Yd2};
	Matrix * Yd0[2] = {Yd01, Yd02};
    Matrix * P[2] = {P1, P2};
	double coeff[2]= {coeff1, coeff2};
    int i = 0, j = 0;
    //printf ("\n GenerateCourbe def 4");
    Courbe * CourbPatron[2], *CourbPatronDXF[2], *CourbPatronBack[2], *CourbMarge[2], *CourbMargeBack[2], *CourbMargeDXF[2];
    Courbe *CourbAv, *CourbAvBack, *CourbAr, *CourbArDXF, *CourbMargeAv, *CourbMargeAvBack, *CourbMargeAr, *CourbCoin, *CourbMargeArDXF, *CourbCoin1[2], *CourbCoin1Back[2], *CourbCoin2[2], *CourbCoin2Back[2], *CourbRep, *CourbRepDXF;
    Courbe *CourbCercle, *CourbCercleDXF;
    //printf ("\n GenerateCourbe def 5");
    for (i = 0; i < 2; i++) {
       // printf ("\n GenerateCourbe 1for: %d", i);
        CourbPatron[i] = new Courbe("Patron");
        CourbPatronDXF[i] = new Courbe("Patron");
        CourbPatronBack[i] = new Courbe("PatronBack");
       // printf ("\n G0.1");
        CourbPatron[i]->points = OFF;
        CourbPatronDXF[i]->points = OFF;
        CourbPatronBack[i]->points = OFF;
       // printf ("\n G0.2");
        CourbPatron[i]->symX = OFF;
        CourbPatronDXF[i]->symX = OFF;
        CourbPatronBack[i]->symX = OFF;
       // printf ("\n G0.3");
        CourbPatron[i]->pts = new Matrix(Xd[i]->GetLignes(), 2);
       // printf ("\n G0.4");
        CourbPatronDXF[i]->pts = new Matrix(Xd[i]->GetLignes(), 2);
       // printf ("\n G0.5");
        CourbPatronBack[i]->pts = new Matrix(Xd[i]->GetLignes(), 2);
       // printf ("\n G1");
        for (j = 0; j < Xd[i]->GetLignes(); j++) {
            CourbPatron[i]->pts->SetElement(j, 0, Xd[i]->Element(j, 0));
            CourbPatronDXF[i]->pts->SetElement(j, 0, Xd[i]->Element(j, 0));
            CourbPatronBack[i]->pts->SetElement(Xd[i]->GetLignes() - j - 1, 0, Xd[i]->Element(j, 0));

            CourbPatron[i]->pts->SetElement(j, 1, Yd[i]->Element(j, 0));
            CourbPatronDXF[i]->pts->SetElement(j, 1, Yd[i]->Element(j, 0));
            CourbPatronBack[i]->pts->SetElement(Xd[i]->GetLignes() - j - 1, 1, Yd[i]->Element(j, 0));
        }
       // printf ("\n G2");
        addCourbe(*AxePatronP, CourbPatron[i]);
        distance = Ones(CourbPatron[i]->pts->GetLignes(), 1);
        for (j = 0; j < distance->GetLignes(); j++) distance->MultiplyElement(j, 0, Marge[i] / 100.0f);
      // printf ("\n CC Marge/Marge/MargeBack");
        CourbMarge[i] = new Courbe("Marge");
        CourbMargeDXF[i] = new Courbe("Marge");
        CourbMargeBack[i] = new Courbe("MargeBack");
       // printf ("\n G3");
        if (i == 0) {
            CourbMarge[i]->pts = calcContour(CourbPatron[i]->pts, distance, -1);
            CourbMargeDXF[i]->pts = calcContour(CourbPatron[i]->pts, distance, -1);
            CourbMargeBack[i]->pts = calcContour(CourbPatron[i]->pts, distance, -1);
        } else {
            CourbMarge[i]->pts = calcContour(CourbPatron[i]->pts, distance, +1);
            CourbMargeDXF[i]->pts = calcContour(CourbPatron[i]->pts, distance, +1);
            CourbMargeBack[i]->pts = calcContour(CourbPatron[i]->pts, distance, +1);
        }
       // printf ("\n G4");
        CourbMarge[i]->points = OFF;
        CourbMarge[i]->symX = OFF;
        CourbMargeDXF[i]->points = OFF;
        CourbMargeDXF[i]->symX = OFF;
        CourbMargeBack[i]->points = OFF;
        CourbMargeBack[i]->symX = OFF;
       // printf ("\n G5");
        for (j = 0; j < CourbMarge[i]->pts->GetLignes(); j++) {
            CourbMargeBack[i]->pts->SetElement(CourbMarge[i]->pts->GetLignes() - j - 1, 0, CourbMarge[i]->pts->Element(j, 0));
            CourbMargeBack[i]->pts->SetElement(CourbMarge[i]->pts->GetLignes() - j - 1, 1, CourbMarge[i]->pts->Element(j, 1));
        }
       // printf ("\n G6");
        addCourbe(*AxePatronP, CourbMarge[i]);
        delete(distance);
    }
   // printf ("\n G7");
    addCourbe(*AxePatronDXFP, CourbPatronDXF[0]);
    if ((Xd[0]->Element(0, 0) != Xd[1]->Element(0, 0)) //test points Av cote 1&2 confondus
            || (Yd[0]->Element(0, 0) != Yd[1]->Element(0, 0))) {
       // if (debug) printf ("\n GenerateCourbe Avant");
        CourbAv = new Courbe("Avant");
        CourbAvBack = new Courbe("AvantBack");
        CourbAv->pts = new Matrix(2, 2);
        CourbAvBack->pts = new Matrix(2, 2);

        CourbAv->pts->SetElement(0, 0, Xd[0]->Element(0, 0));
        CourbAvBack->pts->SetElement(1, 0, Xd[0]->Element(0, 0));

        CourbAv->pts->SetElement(1, 0, Xd[1]->Element(0, 0));
        CourbAvBack->pts->SetElement(0, 0, Xd[1]->Element(0, 0));

        CourbAv->pts->SetElement(0, 1, Yd[0]->Element(0, 0));
        CourbAvBack->pts->SetElement(1, 1, Yd[0]->Element(0, 0));

        CourbAv->pts->SetElement(1, 1, Yd[1]->Element(0, 0));
        CourbAvBack->pts->SetElement(0, 1, Yd[1]->Element(0, 0));

        CourbAv->points = OFF;
        CourbAv->symX = OFF;
        CourbAvBack->points = OFF;
        CourbAv->symX = OFF;

        addCourbe(*AxePatronP, CourbAv);
        ajCoAv = 1;

        //marge avant
        distance = Ones(2, 1);
        for (j = 0; j < distance->GetLignes(); j++)
            distance->MultiplyElement(j, 0, margeDeb / 100.0f);

        CourbMargeAv = new Courbe("MargeAV");
        CourbMargeAvBack = new Courbe("MargeAvBack");
        if (debug) printf ("\n CC MargeAv/MargeAvBack");
        CourbMargeAv->pts = calcContour(CourbAv->pts, distance, +1);
        CourbMargeAvBack->pts = calcContour(CourbAv->pts, distance, +1);
        for (j = 0; j < CourbMargeAv->pts->GetLignes(); j++) {
            CourbMargeAvBack->pts->SetElement(CourbMargeAv->pts->GetLignes() - j - 1, 0, CourbMargeAv->pts->Element(j, 0));
            CourbMargeAvBack->pts->SetElement(CourbMargeAv->pts->GetLignes() - j - 1, 1, CourbMargeAv->pts->Element(j, 1));
        }

        CourbMargeAv->points = OFF;
        CourbMargeAv->symX = OFF;
        CourbMargeAvBack->points = OFF;
        CourbMargeAvBack->symX = OFF;

        addCourbe(*AxePatronP, CourbMargeAv);
        ajCoMaAv = 1;

        delete(distance);

        if  ( ( abs(CourbMargeAv->pts->Element(0, 0) -  CourbMargeAv->pts->Element(1, 0)) > 0.000001 ) 
				|| ( abs(CourbMargeAv->pts->Element(0, 1) - CourbMargeAv->pts->Element(1, 1)) > 0.000001 )	) {

			for (i = 0; i < 2; i++) {

				CourbCoin1[i] = new Courbe("Coin1");
				CourbCoin = new Courbe("Coin");
				CourbCoin1Back[i] = new Courbe("Coin1Back");

				CourbCoin1[i]->points = OFF;
				CourbCoin->points = OFF;
				CourbCoin1Back[i]->points = OFF;

				CourbCoin1[i]->pts = Zeros(3, 2);
				CourbCoin->pts = Zeros(3, 2);
				CourbCoin1Back[i]->pts = Zeros(3, 2);

				CourbCoin1[i]->pts->SetElement(0, 0, CourbMarge[i]->pts->Element(0, 0));
				CourbCoin->pts->SetElement(0, 0, CourbMarge[i]->pts->Element(0, 0));

				CourbCoin1[i]->pts->SetElement(0, 1, CourbMarge[i]->pts->Element(0, 1));
				CourbCoin->pts->SetElement(0, 1, CourbMarge[i]->pts->Element(0, 1));

				//			CourbCoin1[i]->pts->SetElement(2,0, CourbMargeAv->pts->Element(0,0));
				//			CourbCoin->pts->SetElement(2,0, CourbMargeAv->pts->Element(0,0));
				CourbCoin1[i]->pts->SetElement(2, 0, CourbMargeAv->pts->Element(i, 0));
				CourbCoin->pts->SetElement(2, 0, CourbMargeAv->pts->Element(i, 0));

				//			CourbCoin1[i]->pts->SetElement(2,1, CourbMargeAv->pts->Element(0,1));
				//			CourbCoin->pts->SetElement(2,1, CourbMargeAv->pts->Element(0,1));
				CourbCoin1[i]->pts->SetElement(2, 1, CourbMargeAv->pts->Element(i, 1));
				CourbCoin->pts->SetElement(2, 1, CourbMargeAv->pts->Element(i, 1));

				double x, y;
				Inter2Vecteurs(
						CourbMarge[i]->pts->Element(1, 0), CourbMarge[i]->pts->Element(1, 1),
						CourbMarge[i]->pts->Element(0, 0), CourbMarge[i]->pts->Element(0, 1),
						CourbMargeAv->pts->Element(1, 0), CourbMargeAv->pts->Element(1, 1),
						CourbMargeAv->pts->Element(0, 0), CourbMargeAv->pts->Element(0, 1),
						&x, &y);

				CourbCoin1[i]->pts->SetElement(1, 0, x);
				CourbCoin->pts->SetElement(1, 0, x);

				CourbCoin1[i]->pts->SetElement(1, 1, y);
				CourbCoin->pts->SetElement(1, 1, y);

				for (int _i = 0; _i < 3; _i++) {
					CourbCoin1Back[i]->pts -> SetElement(2 - _i, 0, CourbCoin1[i]->pts->Element(_i, 0));
					CourbCoin1Back[i]->pts -> SetElement(2 - _i, 1, CourbCoin1[i]->pts->Element(_i, 1));
				}
				addCourbe(*AxePatronP, CourbCoin);
				ajCo1[i] = 1;
			}
		}
	}
    /*else //relie les marges cote 1&2 par un segment
    {
        if (debug) printf ("\n GenerateCourbe AV");
        CourbAv = new Courbe("AV");
        CourbAv->points = OFF;
        CourbAvBack = new Courbe("AVBack");
        CourbAvBack->points = OFF;
        CourbAv->pts = Zeros(2, 2);
        CourbAvBack->pts = Zeros(2, 2);
        CourbAv->pts->SetElement(0, 0, CourbMarge[0]->pts->Element(0, 0));
        CourbAvBack->pts->SetElement(1, 0, CourbMarge[0]->pts->Element(0, 0));

        CourbAv->pts->SetElement(0, 1, CourbMarge[0]->pts->Element(0, 1));
        CourbAvBack->pts->SetElement(1, 1, CourbMarge[0]->pts->Element(0, 1));

        CourbAv->pts->SetElement(1, 0, CourbMarge[1]->pts->Element(0, 0));
        CourbAvBack->pts->SetElement(0, 0, CourbMarge[1]->pts->Element(0, 0));

        CourbAv->pts->SetElement(1, 1, CourbMarge[1]->pts->Element(0, 1));
        CourbAvBack->pts->SetElement(0, 1, CourbMarge[1]->pts->Element(0, 1));

        addCourbe(*AxePatronP, CourbAv);
        // 1
        ajCoAv = 0;
    }*/
    int n0 = Xd[0]->GetLignes() - 1;
    int n1 = Xd[1]->GetLignes() - 1;
    if ((Xd[0]->Element(n0, 0) != Xd[1]->Element(n1, 0)) //test points Ar cote 1&2 confondus
            || (Yd[0]->Element(n0, 0) != Yd[1]->Element(n1, 0))) {
        if (debug) printf ("\n GenerateCourbe AR");
        CourbAr = new Courbe("AR");
        CourbArDXF = new Courbe("AR");

        CourbAr->pts = new Matrix(2, 2);
        CourbArDXF->pts = new Matrix(2, 2);

        CourbAr->pts->SetElement(0, 0, Xd[0]->Element(n0, 0));
        CourbArDXF->pts->SetElement(0, 0, Xd[0]->Element(n0, 0));

        CourbAr->pts->SetElement(1, 0, Xd[1]->Element(n1, 0));
        CourbArDXF->pts->SetElement(1, 0, Xd[1]->Element(n1, 0));

        CourbAr->pts->SetElement(0, 1, Yd[0]->Element(n0, 0));
        CourbArDXF->pts->SetElement(0, 1, Yd[0]->Element(n0, 0));

        CourbAr->pts->SetElement(1, 1, Yd[1]->Element(n1, 0));
        CourbArDXF->pts->SetElement(1, 1, Yd[1]->Element(n1, 0));

        CourbAr->points = OFF;
        CourbAr->symX = OFF;
        CourbArDXF->points = OFF;
        CourbAr->symX = OFF;

        addCourbe(*AxePatronP, CourbAr);
        ajCoAr = 1;

        distance = Ones(2, 1);

        for (j = 0; j < distance->GetLignes(); j++)
            distance->MultiplyElement(j, 0, margeFin / 100.0f);

        CourbMargeAr = new Courbe("MargeAR");
        CourbMargeArDXF = new Courbe("MargeAR");
        if (debug) printf ("\n CC MargeAr/MargeArDXF");
        CourbMargeAr->pts = calcContour(CourbAr->pts, distance, -1);
        CourbMargeArDXF->pts = calcContour(CourbAr->pts, distance, -1);

        CourbMargeAr->points = OFF;
        CourbMargeAr->symX = OFF;
        CourbMargeArDXF->points = OFF;
        CourbMargeArDXF->symX = OFF;

        addCourbe(*AxePatronP, CourbMargeAr);
        ajCoMaAr = 1;

        delete(distance);
        if  ( 	( abs(CourbMargeAr->pts->Element(0, 0) -  CourbMargeAr->pts->Element(1, 0)) > 0.000001 ) 
				|| 
				( abs(CourbMargeAr->pts->Element(0, 1) - CourbMargeAr->pts->Element(1, 1)) > 0.000001 )	) 
			{

			for (i = 0; i < 2; i++) {
				CourbCoin2[i] = new Courbe("Coin2");
				CourbCoin = new Courbe("Coin2");
				CourbCoin2Back[i] = new Courbe("Coin2Back");

				CourbCoin2[i]->points = OFF;
				CourbCoin->points = OFF;
				CourbCoin2Back[i]->points = OFF;

				CourbCoin2[i]->pts = Zeros(3, 2);
				CourbCoin->pts = Zeros(3, 2);
				CourbCoin2Back[i]->pts = Zeros(3, 2);
				if (debug) printf ("\n GenerateCourbe CourbMarge(%d", i);
				n = CourbMarge[i]->pts->GetLignes() - 1;

				CourbCoin2[i]->pts->SetElement(0, 0, CourbMarge[i]->pts->Element(n, 0));
				CourbCoin->pts->SetElement(0, 0, CourbMarge[i]->pts->Element(n, 0));

				CourbCoin2[i]->pts->SetElement(0, 1, CourbMarge[i]->pts->Element(n, 1));
				CourbCoin->pts->SetElement(0, 1, CourbMarge[i]->pts->Element(n, 1));

				//CourbCoin2[i]->pts->SetElement(2,0, CourbMargeAr->pts->Element(0,0));
				//CourbCoin->pts->SetElement(2,0, CourbMargeAr->pts->Element(0,0));
				CourbCoin2[i]->pts->SetElement(2, 0, CourbMargeAr->pts->Element(i, 0));
				CourbCoin->pts->SetElement(2, 0, CourbMargeAr->pts->Element(i, 0));

				//CourbCoin2[i]->pts->SetElement(2,1, CourbMargeAr->pts->Element(0,1));
				//CourbCoin->pts->SetElement(2,1, CourbMargeAr->pts->Element(0,1));
				CourbCoin2[i]->pts->SetElement(2, 1, CourbMargeAr->pts->Element(i, 1));
				CourbCoin->pts->SetElement(2, 1, CourbMargeAr->pts->Element(i, 1));

				double x, y;
				Inter2Vecteurs(
						CourbMarge[i]->pts->Element(n, 0), CourbMarge[i]->pts->Element(n, 1),
						CourbMarge[i]->pts->Element(n - 1, 0), CourbMarge[i]->pts->Element(n - 1, 1),
						CourbMargeAr->pts->Element(1, 0), CourbMargeAr->pts->Element(1, 1),
						CourbMargeAr->pts->Element(0, 0), CourbMargeAr->pts->Element(0, 1),
						&x, &y);

				if ( x < CourbMarge[i]->pts->Element(n, 0)) {
					if (debug) printf ("\n\n x < ... %f\n\n", 1000.0f* (CourbMarge[i]->pts->Element(n, 0)-x));
					CourbMarge[i]->pts->SetElement(n, 0, x);
					CourbMarge[i]->pts->SetElement(n, 1, y);
					CourbMargeDXF[i]->pts->SetElement(n, 0, x);
					CourbMargeDXF[i]->pts->SetElement(n, 1, y);
					CourbMargeBack[i]->pts->SetElement(0, 0, x);
					CourbMargeBack[i]->pts->SetElement(0, 1, y);

					CourbCoin2[i]->pts->SetElement(0, 0, CourbMarge[i]->pts->Element(n, 0));
					CourbCoin2[i]->pts->SetElement(0, 1, CourbMarge[i]->pts->Element(n, 1));
					CourbCoin->pts->SetElement(0, 0, CourbMarge[i]->pts->Element(n, 0));
					CourbCoin->pts->SetElement(0, 1, CourbMarge[i]->pts->Element(n, 1));
				}


				CourbCoin2[i]->pts->SetElement(1, 0, x);
				CourbCoin->pts->SetElement(1, 0, x);

				CourbCoin2[i]->pts->SetElement(1, 1, y);
				CourbCoin->pts->SetElement(1, 1, y);

				for (int _i = 0; _i < 3; _i++) {
					CourbCoin2Back[i]->pts -> SetElement(2 - _i, 0, CourbCoin2[i]->pts->Element(_i, 0));
					CourbCoin2Back[i]->pts -> SetElement(2 - _i, 1, CourbCoin2[i]->pts->Element(_i, 1));
				}

				addCourbe(*AxePatronP, CourbCoin);
				ajCo2[i] = 1;
			}
		}
    } /*else //relie les marges cote 1&2 par un segment
    {
        CourbAr = new Courbe("AR");
        CourbAr->points = OFF;
        CourbArDXF = new Courbe("AR");
        CourbArDXF->points = OFF;

        CourbAr->pts = Zeros(2, 2);
        CourbArDXF->pts = Zeros(2, 2);

        CourbAr->pts->SetElement(0, 0, CourbMarge[0]->pts->Element(n0, 0));
        CourbArDXF->pts->SetElement(0, 0, CourbMarge[0]->pts->Element(n0, 0));

        CourbAr->pts->SetElement(0, 1, CourbMarge[0]->pts->Element(n0, 1));
        CourbArDXF->pts->SetElement(0, 1, CourbMarge[0]->pts->Element(n0, 1));

        CourbAr->pts->SetElement(1, 0, CourbMarge[1]->pts->Element(n1, 0));
        CourbArDXF->pts->SetElement(1, 0, CourbMarge[1]->pts->Element(n1, 0));

        CourbAr->pts->SetElement(1, 1, CourbMarge[1]->pts->Element(n1, 1));
        CourbArDXF->pts->SetElement(1, 1, CourbMarge[1]->pts->Element(n1, 1));

        addCourbe(*AxePatronP, CourbAr);
        // 1
        ajCoAr = 0;
    }*/
    if (ajCoAr) addCourbe(*AxePatronDXFP, CourbArDXF);

    addCourbe(*AxePatronDXFP, CourbPatronBack[1]);
    if (ajCoAv) addCourbe(*AxePatronDXFP, CourbAvBack);


    addCourbe(*AxeMarginDXFP, CourbMargeDXF[0]);
    if (ajCo2[0]) addCourbe(*AxeMarginDXFP, CourbCoin2[0]);

    if (ajCoMaAr) addCourbe(*AxeMarginDXFP, CourbMargeArDXF);
    if (ajCo2[1]) addCourbe(*AxeMarginDXFP, CourbCoin2Back[1]);


    addCourbe(*AxeMarginDXFP, CourbMargeBack[1]);
    if (ajCo1[1]) addCourbe(*AxeMarginDXFP, CourbCoin1[1]);
    if (ajCoMaAv) addCourbe(*AxeMarginDXFP, CourbMargeAvBack);
    if (ajCo1[0]) addCourbe(*AxeMarginDXFP, CourbCoin1Back[0]);

    if (debug) printf ("\n GC go in Reper()");
    int ventisave=0;
    double xk0, xk1, xk2, yk0, yk1, yk2, xkf, ykf;
    if (makeRep)
    for (i = 0; i < 2; i++) {
        if (debug) printf ("\n i=%d\n", i);
            for (j = 0;j < 5; j++) {
                if (debug) printf (" set j=%d ", j);
                xSuspente[i][j] = -100000.0f;
                ySuspente[i][j] = -100000.0f;
            }

			Matrix* m = 0;

			if (gfd->CorrectRepPoints && isPince && (coeff[i] != -1)) {
				if (debug) printf ("\n coeff=%f", coeff[i]);
				if (debug) printf ("\n 1");
				m = getReperPointsPince (gfd, Xd0[i], Yd0[i], coeff[i], 
										Xd[i], Yd[i], P[i], nerv[i], Deb[i], FaceDeb[i], Fin[i], FaceFin[i], (nerv1==nerv2));
				if (debug) printf ("\n ...1");
			} else {
				m = getReperPoints(gfd, Xd[i], Yd[i], P[i], nerv[i], Deb[i], FaceDeb[i], Fin[i], FaceFin[i], (nerv1==nerv2));
				if (coeff[i] == -1) 
					printf ("\nlength XY=%6.2f", 1000*calcCourbeLength(Xd[i], Yd[i]));
			} 
			if (debug) printf ("\n ...--");

            if (debug) for (int _i=0; _i<m->GetLignes(); _i++) {
                printf ("\n grp %d (%f, %f, %f)", _i, m->Element(_i,0), m->Element(_i,1), m->Element(_i,2));
            }
            //add croix au graphe
            for (j = 0;j<m->GetLignes();j++) {
                if (debug)    printf ("\nj3=%d", j);
                    int mark = m->Element(j,2);
                    double xSus = m->Element(j,0);
                    double ySus = m->Element(j,1);

                    if (((mark == REP_TRIANGLE) || (mark == REP_V)) && (m->Element(j, 3) > 0.0f)) {
                        ventisave = (int) m->Element(j, 3) - 1;
                        xSuspente[i][ventisave] = xSus;
                        ySuspente[i][ventisave] = ySus;
                    }
                    if ((i==0) && (nerv1==nerv2)) {
                        if (mark == REP_0) {
                            xk0 = xSus;
                            yk0 = ySus;
                            //printf ("\n REP_0 (%f, %f)", xk0, yk0);
                        }
                        if (mark == REP_VENTDEB) {
                            xk1 = xSus;
                            yk1 = ySus;
                            //printf ("\n REP_VENTDEB (%f, %f)", xk1, yk1);
                        }
                        if (mark == REP_VENTFIN) {
                            xk2 = xSus;
                            yk2 = ySus;
                            //printf ("\n REP_VENTFIN (%f, %f)", xk2, yk2);
                        }
                        if (mark == REP_KLAPAN_FIN) {
                            xkf = xSus;
                            ykf = ySus;
                            //printf ("\n REP_KLAPAN_FIN (%f, %f)", xkf, ykf);
                        }
                    }


                    if ((mark == REP_CROSS) && (gfd->ReperesProfile[i])) {
                            CourbRep = new Courbe("Repers points");
                            CourbRep->points = OFF;
                            CourbRepDXF = new Courbe("Repers points");
                            CourbRepDXF->points = OFF;
                            CourbRep->pts = Zeros(2, 2);
                            CourbRepDXF->pts = Zeros(2, 2);
                            CourbRep->pts->SetElement(0, 0, xSus - gfd->XMashtab * 0.01f);
                            CourbRep->pts->SetElement(0, 1, ySus - gfd->XMashtab * 0.01f);
                            CourbRep->pts->SetElement(1, 0, xSus + gfd->XMashtab * 0.01f);
                            CourbRep->pts->SetElement(1, 1, ySus + gfd->XMashtab * 0.01f);
                            int _i, _j;
                            for (_i = 0; _i < 2; _i++)
                                for (_j = 0; _j < 2; _j++)
                                    CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                            addCourbe(*AxeRepDXFP, CourbRepDXF);
                            CourbRep = new Courbe("Repers points");
                            CourbRep->points = OFF;
                            CourbRepDXF = new Courbe("Repers points");
                            CourbRepDXF->points = OFF;
                            CourbRep->pts = Zeros(2, 2);
                            CourbRepDXF->pts = Zeros(2, 2);
                            CourbRep->pts->SetElement(0, 0, xSus - gfd->XMashtab * 0.01f);
                            CourbRep->pts->SetElement(0, 1, ySus + gfd->XMashtab * 0.01f);
                            CourbRep->pts->SetElement(1, 0, xSus + gfd->XMashtab * 0.01f);
                            CourbRep->pts->SetElement(1, 1, ySus - gfd->XMashtab * 0.01f);
                            for (_i = 0; _i < 2; _i++)
                                for (_j = 0; _j < 2; _j++)
                                    CourbRepDXF->pts->SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                            addCourbe(*AxeRepDXFP, CourbRepDXF);
                        } 
                        if ((mark==REP_LINE) || ((gfd->ReperesProfile[i]) && (mark == REP_MIDDLE_LINE))) {
                                char courbName[255];
                                if (mark==REP_MIDDLE_LINE) 
                                    strcpy(courbName, "Repers points");
                                else
                                    strcpy(courbName, "Repers diag");
                                //perpendiular
                                int direction = -1;
                                if (i == 1) direction = 1;
                                double _x0 = xSus;
                                double _y0 = ySus;
                                int ix, icol;
                                Ind(_x0, Xd[i], &ix, &icol);
                                //Ind(ySuspente[i][j], interpYSuspente, &iy, &icol);
                                double xc, yc;
                                double _l = 3 * gfd->XMashtab * 0.01f;
                                if (mark == REP_MIDDLE_LINE) {
                                    _l = 0.4 * _l;
                                }
                                if ((Xd[i]->Element(ix, 0) > _x0) && (ix != 0)) ix--;
                                if ((Xd[i]->Element(ix,0) == _x0) && (ix != 0)) ix--;
                                if (ix != 0)
                                    calcVecteurNormal(Xd[i]->Element(ix,0), Yd[i]->Element(ix, 0),_x0, _y0, &xc, &yc, _l, direction);
                                else
                                    calcVecteurNormal(_x0, _y0,Xd[i]->Element(1,0), Yd[i]->Element(1, 0), &xc, &yc, _l, direction);
                                double _x1 = _x0;
                                double _y1 = _y0;
                                double _x2 = xc;
                                double _y2 = yc;
                                if (ix == 0) {
                                    _x2 = _x0 + (xc - Xd[i]->Element(1,0));
                                    _y2 = _y0 + (yc - Yd[i]->Element(1,0));
                                }
                                int _i, _j;
                                CourbRep = new Courbe(courbName);
                                CourbRep->points = OFF;
                                CourbRepDXF = new Courbe(courbName);
                                CourbRepDXF->points = OFF;
                                CourbRep->pts = Zeros(2, 2);
                                CourbRepDXF->pts = Zeros(2, 2);
                                CourbRep->pts->SetElement(0, 0, _x1);
                                CourbRep->pts->SetElement(0, 1, _y1);
                                CourbRep->pts->SetElement(1, 0, _x2);
                                CourbRep->pts->SetElement(1, 1, _y2);
                                for (_i = 0; _i < 2; _i++)
                                    for (_j = 0; _j < 2; _j++)
                                        CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                                addCourbe(*AxeRepDXFP, CourbRepDXF);
                            }
                            if ((gfd -> ReperesSuspentes[i]) && ((mark == REP_TRIANGLE) || (mark==REP_V)) ) {
                                int direction = -1;
                                if (i == 1) direction = 1;
                                double _x0 = xSus;
                                double _y0 = ySus;
                                double _l = 2*gfd->XMashtab * 0.01f;
                                double _x1 = _x0;
                                double _y1 = _y0;

                                int ix, icol;
                                Ind(xSus, Xd[i], &ix, &icol);
                                double xc, yc;
                                if ((Xd[i]->Element(ix,0) > _x0) && (ix != 0 )) ix--;
                                if ((Xd[i]->Element(ix,0) == _x0) && (ix != 0) ) ix--;
                                if  (ix != 0)
                                    calcVecteurNormal(Xd[i]->Element(ix,0), Yd[i]->Element(ix, 0), _x0, _y0, &xc, &yc, _l, direction);
                                else
                                    calcVecteurNormal( _x0, _y0, Xd[i]->Element(1,0), Yd[i]->Element(1, 0), &xc, &yc, _l, direction);
                                if (ix == 0) ix=1;

                                double nx = (Xd[i]->Element(ix, 0)-_x0);
                                double ny = (Yd[i]->Element(ix, 0)-_y0);
                                double rasst = sqrt (sqr(nx) + sqr(ny));
                                nx=_l/2.0f*nx/rasst;
                                ny=_l/2.0f*ny/rasst;

                                double _x2 = xc + nx;
                                double _y2 = yc + ny;

                                double _x3 = xc - nx;
                                double _y3 = yc - ny;
                                int _i, _j;
                                CourbRep = new Courbe("Repers suspente");
                                CourbRep->points = OFF;
                                CourbRepDXF = new Courbe("Repers suspente");
                                CourbRepDXF->points = OFF;
                                CourbRep->pts = Zeros(2, 2);
                                CourbRepDXF->pts = Zeros(2, 2);
                                CourbRep->pts->SetElement(0, 0, _x1);
                                CourbRep->pts->SetElement(0, 1, _y1);
                                CourbRep->pts->SetElement(1, 0, _x2);
                                CourbRep->pts->SetElement(1, 1, _y2);
                                for (_i = 0; _i < 2; _i++)
                                    for (_j = 0; _j < 2; _j++)
                                        CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));

                                addCourbe(*AxeRepDXFP, CourbRepDXF);

                                if (mark != REP_V) {
                                    CourbRep = new Courbe("Repers suspente");
                                    CourbRep->points = OFF;
                                    CourbRepDXF = new Courbe("Repers suspente");
                                    CourbRepDXF->points = OFF;
                                    CourbRep->pts = Zeros(2, 2);
                                    CourbRepDXF->pts = Zeros(2, 2);
                                    CourbRep->pts->SetElement(0, 0, _x2);
                                    CourbRep->pts->SetElement(0, 1, _y2);
                                    CourbRep->pts->SetElement(1, 0, _x3);
                                    CourbRep->pts->SetElement(1, 1, _y3);
                                    for (_i = 0; _i < 2; _i++)
                                        for (_j = 0; _j < 2; _j++)
                                            CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                                    addCourbe(*AxeRepDXFP, CourbRepDXF);
                                }
                                CourbRep = new Courbe("Repers suspente");
                                CourbRep->points = OFF;
                                CourbRepDXF = new Courbe("Repers suspente");
                                CourbRepDXF->points = OFF;
                                CourbRep->pts = Zeros(2, 2);
                                CourbRepDXF->pts = Zeros(2, 2);
                                CourbRep->pts->SetElement(0, 0, _x3);
                                CourbRep->pts->SetElement(0, 1, _y3);
                                CourbRep->pts->SetElement(1, 0, _x1);
                                CourbRep->pts->SetElement(1, 1, _y1);
                                for (_i = 0; _i < 2; _i++)
                                    for (_j = 0; _j < 2; _j++)
                                        CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                                addCourbe(*AxeRepDXFP, CourbRepDXF);
                            }
            }

            if ((i==0) && (nerv1==nerv2) && (gfd->LayoutKlapans) && (gfd->VentHoles)) {
                //printf ("\n go make shov!");
                bool match = false;
				int _i, _j;
                if ((nerv1==-
					1) || (nerv1==0)){
                        if (gfd->VentCentralNerv) match = true;
                }
                else {
                        if ((gfd->noNervVH[nerv1] == 1) || (gfd->noNervVH[nerv1 - 1]==1)) match = true;
                }

                if (match) {
                    //printf (" \n MATCH true!");
                    if (gfd->VentHolesDouble) {
                        // xk0, yk0 -> xkfin, ykfin
                        CourbRep = new Courbe("Klapan line");
                        CourbRep->points = OFF;
                        CourbRepDXF = new Courbe("Klapan line");
                        CourbRepDXF->points = OFF;
                        CourbRep->pts = Zeros(2, 2);
                        CourbRepDXF->pts = Zeros(2, 2);
                        CourbRep->pts->SetElement(0, 0, xk0);
                        CourbRep->pts->SetElement(0, 1, yk0);
                        CourbRep->pts->SetElement(1, 0, xkf);
                        CourbRep->pts->SetElement(1, 1, 0.0f);
                        for (_i = 0; _i < 2; _i++)
                            for (_j = 0; _j < 2; _j++)
                                CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                        addCourbe(*AxeRepDXFP, CourbRepDXF);
                    }
                    // xk1, yk1 -> xkfin, ykfin
                    CourbRep = new Courbe("Klapan line");
                    CourbRep->points = OFF;
                    CourbRepDXF = new Courbe("Klapan line");
                    CourbRepDXF->points = OFF;
                    CourbRep->pts = Zeros(2, 2);
                    CourbRepDXF->pts = Zeros(2, 2);
                    CourbRep->pts->SetElement(0, 0, xk1);
                    CourbRep->pts->SetElement(0, 1, yk1);
                    CourbRep->pts->SetElement(1, 0, xkf);
                    CourbRep->pts->SetElement(1, 1, 0.0f);
                    for (_i = 0; _i < 2; _i++)
                        for (_j = 0; _j < 2; _j++)
                            CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                    addCourbe(*AxeRepDXFP, CourbRepDXF);
                    // xk2, yk2 -> xkfin, ykfin
                    CourbRep = new Courbe("Klapan line");
                    CourbRep->points = OFF;
                    CourbRepDXF = new Courbe("Klapan line");
                    CourbRepDXF->points = OFF;
                    CourbRep->pts = Zeros(2, 2);
                    CourbRepDXF->pts = Zeros(2, 2);
                    CourbRep->pts->SetElement(0, 0, xk2);
                    CourbRep->pts->SetElement(0, 1, yk2);
                    CourbRep->pts->SetElement(1, 0, xkf);
                    CourbRep->pts->SetElement(1, 1, 0.0f);
                    for (_i = 0; _i < 2; _i++)
                        for (_j = 0; _j < 2; _j++)
                            CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                    addCourbe(*AxeRepDXFP, CourbRepDXF);
                } else {
                   // printf (" \n NOOOOT MATCH!");
                }

            }

    }

    if (debug) printf ("\n... in Reper()");

    if (Ventilation == 1) {
      if (debug)   printf ("\nif Ventil");
        for (i = 1; i < 4; i++) {
            if ((xSuspente[0][i] != -100000.0f)
                    && (xSuspente[1][i] != -100000.0f)
                    && (xSuspente[0][i - 1] != -100000.0f)
                    && (xSuspente[1][i - 1] != -100000.0f)) {

                CourbCercle = new Courbe("Cercle");
                CourbCercleDXF = new Courbe("Cercle");

                CourbCercle->points = OFF;
                CourbCercleDXF->points = OFF;

                CourbCercle->pts = Cercle(
                        (xSuspente[0][i] + xSuspente[1][i] + xSuspente[0][i - 1] + xSuspente[1][i - 1]) / 4.0,
                        (ySuspente[0][i] + ySuspente[1][i] + ySuspente[0][i - 1] + ySuspente[1][i - 1]) / 4.0,
                        (ySuspente[0][i] + ySuspente[0][i - 1] - ySuspente[1][i] - ySuspente[1][i - 1]) / 8.0,
                        36);
                CourbCercleDXF->pts = Cercle(
                        (xSuspente[0][i] + xSuspente[1][i] + xSuspente[0][i - 1] + xSuspente[1][i - 1]) / 4.0,
                        (ySuspente[0][i] + ySuspente[1][i] + ySuspente[0][i - 1] + ySuspente[1][i - 1]) / 4.0,
                        (ySuspente[0][i] + ySuspente[0][i - 1] - ySuspente[1][i] - ySuspente[1][i - 1]) / 8.0,
                        36);

                addCourbe(*AxePatronP, CourbCercle);
                addCourbe(*AxeCercleDXFP, CourbCercleDXF);
            }

        }
      if (debug)   printf ("\n...if Ventil");
    }

    //if (Numerotation == 1) {
    xText = (Xd[0]->Element(0, 0) + Xd[0]->Element(Xd[0]->GetLignes() - 1, 0)
            + Xd[1]->Element(0, 0) + Xd[1]->Element(Xd[1]->GetLignes() - 1, 0)) / 4.0f;

    double x0len = Xd[0]->Element(Xd[0]->GetLignes() - 1, 0) - Xd[0]->Element(0, 0);
    double x1len = Xd[1]->Element(Xd[1]->GetLignes() - 1, 0) - Xd[1]->Element(0, 0);
    double xlen = (x0len + x1len) / 2.0f;

    //		printf ("\nXd[0]->Element(0,0) %f",Xd[0]->Element(0,0));
    //		printf ("\nXd[0]->Element(Xd[0]->GetLignes()-1,0) %f", Xd[0]->Element(Xd[0]->GetLignes()-1,0));
    //		printf ("\nXd[1]->Element(0,0) %f", Xd[1]->Element(0,0));
    //		printf ("\nXd[1]->Element(Xd[1]->GetLignes()-1,0) %f", Xd[1]->Element(Xd[1]->GetLignes()-1,0));

    xText = gfd->textX * xlen;
    yText = (Yd[0]->Element(0, 0) + Yd[0]->Element(Yd[0]->GetLignes() - 1, 0)
            + Yd[1]->Element(0, 0) + Yd[1]->Element(Yd[1]->GetLignes() - 1, 0)) / 4.0f;

    double y0len = Yd[1]->Element(0, 0) - Yd[0]->Element(0, 0);
    double y1len = Yd[1]->Element(Yd[1]->GetLignes() - 1, 0) - Yd[0]->Element(Yd[0]->GetLignes() - 1, 0);
    double ylen = (y0len + y1len) / 2.0f;

    yText = Yd[1]->Element(0, 0) - gfd->textY*ylen;

    //		printf ("\nYd[0]->Element(0,0) %f",Yd[0]->Element(0,0));
    //		printf ("\nYd[0]->Element(Yd[0]->GetLignes()-1,0) %f", Yd[0]->Element(Yd[0]->GetLignes()-1,0));
    //		printf ("\nYd[1]->Element(0,0) %f", Yd[1]->Element(0,0));
    //		printf ("\nYd[1]->Element(Yd[1]->GetLignes()-1,0) %f", Yd[1]->Element(Yd[1]->GetLignes()-1,0));
    //sprintf(texte, "%d%c%03.1f%c%03.1f_%d%c%03.1f%c%03.1f",
    //	NoNerv[0],texteExtInt[FaceDeb[0]-1],Deb[0],texteExtInt[FaceFin[0]-1],Fin[0],
    //	NoNerv[1],texteExtInt[FaceDeb[1]-1],Deb[1],texteExtInt[FaceFin[1]-1],Fin[1]);
    sprintf(texte, "%s", text);
    addTexte(*AxePatronP, texte, 0.02f, 0.0f, xText, yText);
    sprintf(texte, "%s", text);
    addTexte(*AxePatronTextDXFP, texte, 0.02f, 0.0f, xText, yText);
    if (debug) printf ("\n...GenerateCourbe()");
    //}
    //    printf ("\n ...GenerateCourbe()");
}

void fillXYZP(int* nerv, double* DebT, int* FaceDebT, double* FinT,int* FaceFinT,
			  Matrix* xExtProf, Matrix* XExt,Matrix* YExt,Matrix* ZExt,
			  Matrix* xIntProf, Matrix* XInt,Matrix* YInt,Matrix* ZInt,
              Matrix** X0,Matrix** Y0, Matrix** Z0,Matrix** P0,
		      Matrix** X1,Matrix** Y1, Matrix** Z1,Matrix** P1) {
    Matrix * X[2], *Y[2], *Z[2], *P[2];
	//printf ("\n DebT[%f, %f]", DebT[0], DebT[1]);
	//printf ("\n FinT[%f, %f]", FinT[0], FinT[1]);
    int iDeb, iFin, iCol;
	int i, j;
    for (i = 0; i < 2; i++) {
        if ((FaceDebT[i] == 1) && (FaceFinT[i] == 1)) //debut et fin en extrados
        {
            if (DebT[i] > FinT[i]) //inversion Deb<->Fin
            {
                Ind(FinT[i], xExtProf, &iDeb, &iCol);
                Ind(DebT[i], xExtProf, &iFin, &iCol);
            } else {
                Ind(DebT[i], xExtProf, &iDeb, &iCol);
                Ind(FinT[i], xExtProf, &iFin, &iCol);
            }
            X[i] = new Matrix(iFin - iDeb + 1, 1);
            Y[i] = new Matrix(iFin - iDeb + 1, 1);
            Z[i] = new Matrix(iFin - iDeb + 1, 1);
            P[i] = new Matrix(iFin - iDeb + 1, 1);
            for (j = 0; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Y[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Z[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < P[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xExtProf->Element(iDeb + j, 0));
        } else if (FaceDebT[i] == 1) //debut en extrados et fin en intrados
        {
            Ind(DebT[i], xExtProf, &iDeb, &iCol);
            Ind(FinT[i], xIntProf, &iFin, &iCol);
            X[i] = new Matrix(iFin + iDeb - 1, 1);
            Y[i] = new Matrix(iFin + iDeb - 1, 1);
            Z[i] = new Matrix(iFin + iDeb - 1, 1);
            P[i] = new Matrix(iFin + iDeb - 1, 1);
            for (j = 0; j <= iDeb; j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) P[i]->SetElement(j, 0, xExtProf->Element(iDeb - j, 0));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xIntProf->Element(j - iDeb, 0));
        } else if (FaceFinT[i] == 1) //debut en intrados et fin en extrados
        {
            Ind(DebT[i], xIntProf, &iDeb, &iCol);
            Ind(FinT[i], xExtProf, &iFin, &iCol);
			//printf ("\n iDeb=%d iFin=%d", iDeb, iFin);
            X[i] = new Matrix(iFin + iDeb + 1, 1);
            Y[i] = new Matrix(iFin + iDeb + 1, 1);
            Z[i] = new Matrix(iFin + iDeb + 1, 1);
            P[i] = new Matrix(iFin + iDeb + 1, 1);
            for (j = 0; j <= iDeb; j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) P[i]->SetElement(j, 0, xIntProf->Element(iDeb - j, 0));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xExtProf->Element(j - iDeb, 0));
        } else //debut et fin en intrados
        {
            if (DebT[i] > FinT[i]) //inversion Deb<->Fin
            {
                Ind(FinT[i], xIntProf, &iDeb, &iCol);
                Ind(DebT[i], xIntProf, &iFin, &iCol);
            } else {
                Ind(DebT[i], xIntProf, &iDeb, &iCol);
                Ind(FinT[i], xIntProf, &iFin, &iCol);
            }
            X[i] = new Matrix(iFin - iDeb + 1, 1);
            Y[i] = new Matrix(iFin - iDeb + 1, 1);
            Z[i] = new Matrix(iFin - iDeb + 1, 1);
            P[i] = new Matrix(iFin - iDeb + 1, 1);
            for (j = 0; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Y[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Z[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < P[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xIntProf->Element(iDeb + j, 0));
        }
    }
    (*X0) = X[0]; (*Y0) = Y[0]; (*Z0) = Z[0]; (*P0) = P[0];
    (*X1) = X[1]; (*Y1) = Y[1]; (*Z1) = Z[1]; (*P1) = P[1];
}

/*void makeChangePoints(Matrix* XExt,Matrix*  YExt,Matrix*  ZExt,
                      Matrix*  XInt,Matrix*  YInt,Matrix*  ZInt,
        Matrix**  XExt1,Matrix** YExt1, Matrix** ZExt1,
        Matrix** XInt1, Matrix** YInt1, Matrix** ZInt1, int nerv, float pos)
{
    XExt1=CloneMat(XExt); YExt1=CloneMat(YExt); ZExt1=CloneMat(ZExt);
    XInt1=CloneMat(XInt); YInt1=CloneMat(YInt); ZInt1=CloneMat(ZInt);
    int n = XInt->GetColonnes();
    Matrix* interpXInt = new Matrix(n, 2);
    Matrix* interpYInt = new Matrix(n, 2);
    Matrix* interpZInt = new Matrix(n, 2);
    for (int i = 0; i<n;i++) {
        interpXInt->SetElement(i,0,XInt->Element(nerv, i));
        interpXInt->SetElement(i,1,XInt->Element(nerv, i));
    }
    for (int i = 0; i < XExt->GetColonnes(); i++) {
        XExt1 -> SetElement(nerv, i, XInt->Element);
    }
}*/

void calcPatronPosNerv(WindPatternsProject* gfd, int noNerv1, bool sym1, int FaceDeb1, int FaceFin1, double Deb1, double Fin1,
        int noNerv2, bool sym2, int FaceDeb2, int FaceFin2, double Deb2, double Fin2,
        Matrix **Xdev1, Matrix **Ydev1, Matrix **Xdev2, Matrix **Ydev2,
        Matrix **X0, Matrix **Y0, Matrix **Z0, Matrix **P0,
        Matrix **X1, Matrix **Y1, Matrix **Z1, Matrix **P1)
{
    //printf ("\n calcPatronPosNerv()");
    Matrix *XExt, *YExt, *ZExt;
    Matrix *XInt, *YInt, *ZInt;
    Matrix *XExt1, *YExt1, *ZExt1;
    Matrix *XInt1, *YInt1, *ZInt1;
    Matrix *XExt2, *YExt2, *ZExt2;
    Matrix *XInt2, *YInt2, *ZInt2;

    int i, j;
    bool symetrique[2] = {sym1, sym2};
    int nerv[2];
    Matrix *X[2], *Y[2], *Z[2], *P[2];
	Matrix *Xt[2], *Yt[2], *Zt[2], *Pt[2];
    Matrix *xExtProf, *xIntProf;
    int NoNervT[2] = {noNerv1, noNerv2};
    double DebT[2] = {Deb1, Deb2}, FinT[2] = {Fin1, Fin2};
    int FaceDebT[2], FaceFinT[2];
    FaceDebT[0] = 1;
    FaceDebT[1] = 1;
    FaceFinT[0] = 1;
    FaceFinT[1] = 1;

    //Matrix *ExtProfCent1, *IntProfCent1, *ExtProfCent2, *IntProfCent2;
    //Matrix *ExtProfBout1, *IntProfBout1, *ExtProfBout2, *IntProfBout2;
/*
    printf ("\n before makePosProfile()");
    makePosProfile(gfd->ExtProfCent, gfd->IntProfCent, gfd->PosNerv[0], &ExtProfCent1);
    makePosProfile(gfd->ExtProfBout, gfd->IntProfCent, gfd->PosNerv[0], &ExtProfBout1);
    ExtProfCent1->print(0);
    ExtProfCent1->print(1);
    makePosProfile(gfd->ExtProfCent, gfd->IntProfCent, gfd->PosNerv[1], &ExtProfCent2);
    makePosProfile(gfd->ExtProfBout, gfd->IntProfCent, gfd->PosNerv[1], &ExtProfBout2);
    printf ("\n after makePosProfile()");*/
    //printf ("\n before Form0");
    calcForm3D(gfd->Form,0, 0.0f,
            gfd->ExtProfCent, gfd->IntProfCent, gfd->ExtProfBout, gfd->IntProfBout,
            &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);
    //printf ("\n before Form1");
    calcForm3D(gfd->Form, 1, gfd->PosNerv[0],
            gfd->ExtProfCent, gfd->IntProfCent, gfd->ExtProfBout, gfd->IntProfBout,
            &XExt1, &YExt1, &ZExt1, &XInt1, &YInt1, &ZInt1);
    //printf ("\n before Form1");
    calcForm3D(gfd->Form, 1, gfd->PosNerv[1],
            gfd->ExtProfCent, gfd->IntProfCent, gfd->ExtProfBout, gfd->IntProfBout,
            &XExt2, &YExt2, &ZExt2, &XInt2, &YInt2, &ZInt2);
    //printf ("\n after calcForm3D");
    xExtProf = new Matrix(gfd->ExtProfCent->GetLignes(), 1);
    xIntProf = new Matrix(gfd->IntProfCent->GetLignes(), 1);

    for (i = 0; i < gfd->ExtProfCent->GetLignes(); i++)   xExtProf->SetElement(i, 0, gfd->ExtProfCent->Element(i, 0));
    for (i = 0; i < gfd->IntProfCent->GetLignes(); i++)   xIntProf->SetElement(i, 0, gfd->IntProfCent->Element(i, 0));

    for (i = 0; i < 2; i++) {
        nerv[i] = NoNervT[i];
        if (NoNervT[i] == -1) {
            nerv[i] = 0;
            symetrique[i] = true;
        }
    }

    fillXYZP(nerv, DebT, FaceDebT, FinT, FaceFinT,
                    xExtProf, XExt1, YExt1, ZExt1,
                    xIntProf, XInt1, YInt1, ZInt1,
                   &X[0], &Y[0], &Z[0], &P[0],
		   &Xt[1], &Yt[1], &Zt[1], &Pt[1]) ;
    fillXYZP(nerv, DebT, FaceDebT, FinT, FaceFinT,
                    xExtProf, XExt2, YExt2, ZExt2,
                    xIntProf, XInt2, YInt2, ZInt2,
                   &Xt[0], &Yt[0], &Zt[0], &Pt[0],
		   &X[1], &Y[1], &Z[1], &P[1]);
    //printf ("\n after FillXYZP");
    makeEqualSize(&X[0], &Y[0], &Z[0], &P[0],&X[1], &Y[1], &Z[1], &P[1]);
    //printf ("\n after makeEqualSize");
    for (i = 0; i < 2; i++) {
        if (symetrique[i] == true) {
            for (j = 0; j < X[i]->GetLignes(); j++)
                X[i]->MultiplyElement(j, 0, -1.0f);
        }
    }
    calcDeveloppe(
            X[0], Y[0], Z[0], X[1], Y[1], Z[1],
            Xdev1, Ydev1, Xdev2, Ydev2);
    
    *X0 = CloneMat(X[0]);
    *Y0 = CloneMat(Y[0]);
    *Z0 = CloneMat(Z[0]);
    *P0 = CloneMat(P[0]);
    *X1 = CloneMat(X[1]);
    *Y1 = CloneMat(Y[1]);
    *Z1 = CloneMat(Z[1]);
    *P1 = CloneMat(P[1]);
    
    delete(X[0]);
    delete(Y[0]);
    delete(Z[0]);
    delete(P[0]);
    delete(XExt);
    delete(YExt);
    delete(ZExt);
    delete(XInt);
    delete(YInt);
    delete(ZInt);

    delete(XExt1);
    delete(YExt1);
    delete(ZExt1);
    delete(XInt1);
    delete(YInt1);
    delete(ZInt1);

    delete(XExt2);
    delete(YExt2);
    delete(ZExt2);
    delete(XInt2);
    delete(YInt2);
    delete(ZInt2);


    delete(xExtProf);
    delete(xIntProf);
    delete(X[1]);
    delete(Y[1]);
    delete(Z[1]);
    delete(P[1]);

}
/*        calcPatronKlapan(getWindPatternsProject(), nerv[0], symetrique[0], nerv[1], symetrique[1], PosKlapanInt, PosKlapanFin,
                    &Xd[0], &Yd[0], &Xd[1], &Yd[1],
                    &X[0], &Y[0], &Z[0], &P[0],
                    &X[1], &Y[1], &Z[1], &P[1]);*/
void calcPatronKlapan(WindPatternsProject* gfd, int noNerv1, bool sym1, int noNerv2, bool sym2, double posKlapanInt, double posKlapanFin,
        Matrix **Xdev1, Matrix **Ydev1, Matrix **Xdev2, Matrix **Ydev2,
        Matrix **X0, Matrix **Y0, Matrix **Z0, Matrix **P0,
        Matrix **X1, Matrix **Y1, Matrix **Z1, Matrix **P1)
{

    Matrix *XExt, *YExt, *ZExt,*XExt0, *YExt0, *ZExt0;
    Matrix *XInt, *YInt, *ZInt,*XInt0, *YInt0, *ZInt0;
    int i, j;
    bool symetrique[2] = {sym1, sym2};
    int nerv[2], iDebInt, iFinExt, iCol;
    int NoNervT[2] = {noNerv1, noNerv2};

    Matrix * X[2], *Y[2], *Z[2], *P[2];
    Matrix *xExtProf, *xIntProf;

    calcForm3D(gfd->Form,0, 0.0f,
            gfd->ExtProfCent, gfd->IntProfCent, gfd->ExtProfBout, gfd->IntProfBout,
            &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);

    xExtProf = new Matrix(gfd->ExtProfCent->GetLignes(), 1);
    xIntProf = new Matrix(gfd->IntProfCent->GetLignes(), 1);
    for (i = 0; i < gfd->ExtProfCent->GetLignes(); i++)   xExtProf->SetElement(i, 0, gfd->ExtProfCent->Element(i, 0));
    for (i = 0; i < gfd->IntProfCent->GetLignes(); i++)   xIntProf->SetElement(i, 0, gfd->IntProfCent->Element(i, 0));

    Ind(posKlapanInt, xIntProf, &iDebInt, &iCol);

    Matrix* zeroExtProfCent, *zeroExtProfBout;

    zeroExtProfCent = CloneMat(gfd-> ExtProfCent);
    for ( i = 0; i < gfd->ExtProfCent->GetLignes(); i++) zeroExtProfCent->SetElement(i, 1, 0.0f);
    zeroExtProfBout = CloneMat(gfd-> ExtProfBout);
    for ( i = 0; i < gfd->ExtProfBout->GetLignes(); i++) zeroExtProfBout->SetElement(i, 1, 0.0f);

    calcForm3D(gfd->Form,0, 0.0f,
            zeroExtProfCent, gfd->IntProfCent, zeroExtProfBout, gfd->IntProfBout,
            &XExt0, &YExt0, &ZExt0, &XInt0, &YInt0, &ZInt0);
    Ind(posKlapanFin, xExtProf, &iFinExt, &iCol);

    for (i = 0; i < 2; i++) {
        nerv[i] = NoNervT[i];
        if (NoNervT[i] == -1) {
            nerv[i] = 0;
            symetrique[i] = true;
        }
    }

    X[0] = new Matrix(2, 1);
    X[1] = new Matrix(2, 1);
    Y[0] = new Matrix(2, 1);
    Y[1] = new Matrix(2, 1);
    Z[0] = new Matrix(2, 1);
    Z[1] = new Matrix(2, 1);
    P[0] = new Matrix(2, 1);
    P[1] = new Matrix(2, 1);

    for ( i = 0; i < 2; i++) {
        X[i]->SetElement(0,0,XInt->Element(nerv[i], iDebInt));
        Y[i]->SetElement(0,0,YInt->Element(nerv[i], iDebInt));
        Z[i]->SetElement(0,0,ZInt->Element(nerv[i], iDebInt));
        P[i]->SetElement(0,0,posKlapanInt);

        X[i]->SetElement(1,0,XExt0->Element(nerv[i], iFinExt));
        Y[i]->SetElement(1,0,YExt0->Element(nerv[i], iFinExt));
        Z[i]->SetElement(1,0,ZExt0->Element(nerv[i], iFinExt));
        P[i]->SetElement(1,0,posKlapanFin);
    }

    for (i = 0; i < 2; i++) {
        if (symetrique[i] == true) {
            for (j = 0; j < X[i]->GetLignes(); j++)
                X[i]->MultiplyElement(j, 0, -1.0f);
        }
    }

    calcDeveloppe(
            X[0], Y[0], Z[0], X[1], Y[1], Z[1],
            Xdev1, Ydev1, Xdev2, Ydev2);

    *X0 = CloneMat(X[0]);
    *Y0 = CloneMat(Y[0]);
    *Z0 = CloneMat(Z[0]);
    *P0 = CloneMat(P[0]);
    *X1 = CloneMat(X[1]);
    *Y1 = CloneMat(Y[1]);
    *Z1 = CloneMat(Z[1]);
    *P1 = CloneMat(P[1]);
    delete(X[0]);
    delete(Y[0]);
    delete(Z[0]);
    delete(P[0]);
    delete(XExt);
    delete(YExt);
    delete(ZExt);
    delete(XInt);
    delete(YInt);
    delete(ZInt);
    delete(xExtProf);
    delete(xIntProf);
    delete(X[1]);
    delete(Y[1]);
    delete(Z[1]);
    delete(P[1]);
    
}


void calcPatron(WindPatternsProject* gfd, int noNerv1, bool sym1, int FaceDeb1, int FaceFin1, double Deb1, double Fin1,
        int noNerv2, bool sym2, int FaceDeb2, int FaceFin2, double Deb2, double Fin2,
        Matrix **Xdev1, Matrix **Ydev1, Matrix **Xdev2, Matrix **Ydev2,
        Matrix **X0, Matrix **Y0, Matrix **Z0, Matrix **P0,
        Matrix **X1, Matrix **Y1, Matrix **Z1, Matrix **P1)
{

    //printf ("\n calcPatron()");
	//printf ("\n FaceDeb1=%d FaceFin1=%d FaceDeb2=%d FaceFin2=%d", FaceDeb1, FaceFin1, FaceDeb2, FaceFin2);
	//printf ("\n[%f, %f] [%f, %f]", Deb1, Fin1, Deb2, Fin2);

    Matrix *XExt, *YExt, *ZExt;
    Matrix *XInt, *YInt, *ZInt;

    int i, j;
    bool symetrique[2] = {sym1, sym2};
    int nerv[2];

    Matrix * X[2], *Y[2], *Z[2], *P[2];
    //int iDeb, iFin, iCol;

    Matrix *xExtProf, *xIntProf;
    //Matrix *iAv, *iAp, *Res;

    int NoNervT[2] = {noNerv1, noNerv2};
    double DebT[2] = {Deb1, Deb2}, FinT[2] = {Fin1, Fin2};
    int FaceDebT[2], FaceFinT[2];
    FaceDebT[0] = FaceDeb1;
    FaceDebT[1] = FaceDeb2;
    FaceFinT[0] = FaceFin1;
    FaceFinT[1] = FaceFin2;

    calcForm3D(gfd->Form,0, 0.0f,
            gfd->ExtProfCent, gfd->IntProfCent, gfd->ExtProfBout, gfd->IntProfBout,
            &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);
    /*adde points de suspentage*/
    //	if(ViewPtsSuspentes){
    //		addPtsSuspentage(Axe3d, F, IntProfCent, XExt, YExt, ZExt, XInt, YInt, ZInt,
    //			ViewSymetrique, &PtsSuspentes);}

    //------ BEGIN fill X, Y, Z, P [i] ----------
    xExtProf = new Matrix(gfd->ExtProfCent->GetLignes(), 1);
    xIntProf = new Matrix(gfd->IntProfCent->GetLignes(), 1);
    for (i = 0; i < gfd->ExtProfCent->GetLignes(); i++)   xExtProf->SetElement(i, 0, gfd->ExtProfCent->Element(i, 0));
    for (i = 0; i < gfd->IntProfCent->GetLignes(); i++)   xIntProf->SetElement(i, 0, gfd->IntProfCent->Element(i, 0));


    for (i = 0; i < 2; i++) {
        nerv[i] = NoNervT[i];
        if (NoNervT[i] == -1) {
            nerv[i] = 0;
            symetrique[i] = true;
        }
    }

	//printf ("\n FaceDebT=%d FaceFinT=%d", FaceDebT, FaceFinT);
    fillXYZP(nerv, DebT, FaceDebT, FinT, FaceFinT,
				xExtProf, XExt, YExt, ZExt, xIntProf,
				XInt, YInt, ZInt, 
	           &X[0], &Y[0], &Z[0], &P[0],
			   &X[1], &Y[1], &Z[1], &P[1]) ;


/*
    for (i = 0; i < 2; i++) {
        nerv[i] = NoNervT[i];
        if (NoNervT[i] == -1) {
            nerv[i] = 0;
            symetrique[i] = true;
        }

        if (FaceDebT[i] == 1 && FaceFinT[i] == 1) //debut et fin en extrados
        {
            if (DebT[i] > FinT[i]) //inversion Deb<->Fin
            {
                Ind(FinT[i], xExtProf, &iDeb, &iCol);
                Ind(DebT[i], xExtProf, &iFin, &iCol);
            } else {
                Ind(DebT[i], xExtProf, &iDeb, &iCol);
                Ind(FinT[i], xExtProf, &iFin, &iCol);
            }
            X[i] = new Matrix(iFin - iDeb + 1, 1);
            Y[i] = new Matrix(iFin - iDeb + 1, 1);
            Z[i] = new Matrix(iFin - iDeb + 1, 1);
            P[i] = new Matrix(iFin - iDeb + 1, 1);
            for (j = 0; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Y[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Z[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < P[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xExtProf->Element(iDeb + j, 0));
        } else if (FaceDebT[i] == 1) //debut en extrados et fin en intrados
        {
            Ind(DebT[i], xExtProf, &iDeb, &iCol);
            Ind(FinT[i], xIntProf, &iFin, &iCol);
            X[i] = new Matrix(iFin + iDeb - 1, 1);
            Y[i] = new Matrix(iFin + iDeb - 1, 1);
            Z[i] = new Matrix(iFin + iDeb - 1, 1);
            P[i] = new Matrix(iFin + iDeb - 1, 1);
            for (j = 0; j <= iDeb; j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) P[i]->SetElement(j, 0, xExtProf->Element(iDeb - j, 0));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xIntProf->Element(j - iDeb, 0));
        } else if (FaceFinT[i] == 1) //debut en intrados et fin en extrados
        {
            Ind(DebT[i], xIntProf, &iDeb, &iCol);
            Ind(FinT[i], xExtProf, &iFin, &iCol);
            X[i] = new Matrix(iFin + iDeb + 1, 1);
            Y[i] = new Matrix(iFin + iDeb + 1, 1);
            Z[i] = new Matrix(iFin + iDeb + 1, 1);
            P[i] = new Matrix(iFin + iDeb + 1, 1);
            for (j = 0; j <= iDeb; j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) P[i]->SetElement(j, 0, xIntProf->Element(iDeb - j, 0));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xExtProf->Element(j - iDeb, 0));
        } else //debut et fin en intrados
        {
            if (DebT[i] > FinT[i]) //inversion Deb<->Fin
            {
                Ind(FinT[i], xIntProf, &iDeb, &iCol);
                Ind(DebT[i], xIntProf, &iFin, &iCol);
            } else {
                Ind(DebT[i], xIntProf, &iDeb, &iCol);
                Ind(FinT[i], xIntProf, &iFin, &iCol);
            }
            X[i] = new Matrix(iFin - iDeb + 1, 1);
            Y[i] = new Matrix(iFin - iDeb + 1, 1);
            Z[i] = new Matrix(iFin - iDeb + 1, 1);
            P[i] = new Matrix(iFin - iDeb + 1, 1);
            for (j = 0; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Y[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Z[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < P[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xIntProf->Element(iDeb + j, 0));
        }
    }
 */

    //---- END fill X, Y, Z, P [i] ----------
    //	for (int _i=0;_i < P[0]->GetLignes();_i++) {
    //		printf ("\nP[0]->Element(%d, 0)=%f", _i, P[0]->Element(_i, 0));
    //	}
    //---- BEGIN interpolation equal quantity points X[0] <-> X[1]; Y[0] <-> Y[1]
    makeEqualSize(&X[0], &Y[0], &Z[0], &P[0],&X[1], &Y[1], &Z[1], &P[1]);

    /*if (X[0]->GetLignes() > X[1]->GetLignes()) {
        iAv = LinSpace(0.0f, (double) X[1]->GetLignes() - 1, X[1]->GetLignes());
        iAp = LinSpace(0.0f, (double) X[1]->GetLignes() - 1, X[0]->GetLignes());
        Res = new Matrix(X[0]->GetLignes(), 1);
        InterpLinMat(iAv, X[1], iAp, Res);
        delete(X[1]);
        X[1] = Res;
        Res = new Matrix(X[0]->GetLignes(), 1);
        InterpLinMat(iAv, Y[1], iAp, Res);
        delete(Y[1]);
        Y[1] = Res;
        Res = new Matrix(X[0]->GetLignes(), 1);
        InterpLinMat(iAv, Z[1], iAp, Res);
        delete(Z[1]);
        Z[1] = Res;
        Res = new Matrix(X[0]->GetLignes(), 1);
        InterpLinMat(iAv, P[1], iAp, Res);
        delete(P[1]);
        P[1] = Res;
        delete(iAv);
        delete(iAp);
    } else if (X[1]->GetLignes() > X[0]->GetLignes()) {
        iAv = LinSpace(0.0f, (double) X[0]->GetLignes() - 1, X[0]->GetLignes());
        iAp = LinSpace(0.0f, (double) X[0]->GetLignes() - 1, X[1]->GetLignes());
        Res = new Matrix(X[1]->GetLignes(), 1);
        InterpLinMat(iAv, X[0], iAp, Res);
        delete(X[0]);
        X[0] = Res;
        Res = new Matrix(X[1]->GetLignes(), 1);
        InterpLinMat(iAv, Y[0], iAp, Res);
        delete(Y[0]);
        Y[0] = Res;
        Res = new Matrix(X[1]->GetLignes(), 1);
        InterpLinMat(iAv, Z[0], iAp, Res);
        delete(Z[0]);
        Z[0] = Res;
        Res = new Matrix(X[1]->GetLignes(), 1);
        InterpLinMat(iAv, P[0], iAp, Res);
        delete(P[0]);
        P[0] = Res;
        delete(iAv);
        delete(iAp);
    }*/

    //---- END OF interpolation equal quantity points X[0] <-> X[1]; Y[0] <-> Y[1]
    for (i = 0; i < 2; i++) {
        if (symetrique[i] == true) {
            for (j = 0; j < X[i]->GetLignes(); j++)
                X[i]->MultiplyElement(j, 0, -1.0f);
        }
    }

    //--- BEGIN calcate Developpe in Xd[0], Yd[0]   Xd[1], Yd[1]
    //        printf ("\n calcp |X0|=%d |Y0|=%d |Z0|=%d |P0|=%d ", X[0]->GetLignes(), Y[0]->GetLignes(), Z[0]->GetLignes(), P[0]->GetLignes());
    //        printf ("\n calcp |X1|=%d |Y1|=%d |Z1|=%d |P1|=%d ", X[1]->GetLignes(), Y[1]->GetLignes(), Z[1]->GetLignes(), P[1]->GetLignes());

    calcDeveloppe(
            X[0], Y[0], Z[0], X[1], Y[1], Z[1],
            Xdev1, Ydev1, Xdev2, Ydev2);
    //--- END calcate Developpe in Xd[0], Yd[0]   Xd[1], Yd[1]
    *X0 = CloneMat(X[0]);
    *Y0 = CloneMat(Y[0]);
    *Z0 = CloneMat(Z[0]);
    *P0 = CloneMat(P[0]);
    *X1 = CloneMat(X[1]);
    *Y1 = CloneMat(Y[1]);
    *Z1 = CloneMat(Z[1]);
    *P1 = CloneMat(P[1]);
    delete(X[0]);
    delete(Y[0]);
    delete(Z[0]);
    delete(P[0]);
    delete(XExt);
    delete(YExt);
    delete(ZExt);
    delete(XInt);
    delete(YInt);
    delete(ZInt);
    delete(xExtProf);
    delete(xIntProf);
    delete(X[1]);
    delete(Y[1]);
    delete(Z[1]);
    delete(P[1]);
    //printf ("\n ...calcPatron()");
}
