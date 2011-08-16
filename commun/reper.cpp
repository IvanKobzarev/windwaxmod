#include "reper.h"
#include "logger.h"

Matrix* getReperPoints(WindPatternsProject* gfd, Matrix* Xd, Matrix* Yd, Matrix* P, 
						int nerv, float deb, int faceDeb, float fin, int faceFin, bool klapanShov)
{
	bool debug = false;
	if (debug) {
		printf ("\nget Reper Points()");
		printf ("\n P:");
		printM (P);
		printf ("\n...");
	}
    Matrix *interpXSuspente, *interpYSuspente, *interpXSuspente0, *interpYSuspente0;
    bool cerezNos = false;
    int i0 = 0;
    if ((faceDeb != faceFin) && (deb != 0.0f)) {
        cerezNos = true;
        if (debug) printf ("\n cerezNos");
        
		while (P->Element(i0, 0) != 0.0f) {
            if (debug) printf ("\n i0=%d", i0);
            i0++;
        }

		if (debug) printf ("\n res i0: %d  P->GetLignes(): %d P->Element(i0, 0): %f", i0,  P->GetLignes(), P->Element(i0, 0));

		if (i0 >= P->GetLignes()) {
			printf ("\n\n\nERROR! reper.cpp getReperPoints() faceDeb != faceFin, but no 0.00 in P");
		}
    }
    if (debug) printf ("\n...1 i0=%d", i0);
    interpXSuspente = Zeros(P -> GetLignes() - i0, 2);
    interpYSuspente = Zeros(P -> GetLignes() - i0, 2);
    for (int j = i0; j < P->GetLignes(); j++) {
		if (debug) printf ("\n interpXSuspente %d fill", (j - i0));
        interpXSuspente->SetElement(j - i0, 0, P->Element(j, 0));
        interpYSuspente->SetElement(j - i0, 0, P->Element(j, 0));
        interpXSuspente->SetElement(j - i0, 1, Xd->Element(j, 0));
        interpYSuspente->SetElement(j - i0, 1, Yd->Element(j, 0));
    }
	if (debug) printf ("\n interpXSuspente->GetLignes(): %d", interpXSuspente->GetLignes());


    if (debug) printf ("\n ...2");
    int isave = 0;
    if (cerezNos) {
		if (debug) printf ("\n cerezNos");
        interpXSuspente0 = Zeros(i0 + 1, 2);
        interpYSuspente0 = Zeros(i0 + 1, 2);
        for (int j = i0; j >= 0; j--) {
            interpXSuspente0->SetElement(isave, 0, P->Element(j, 0));
            interpYSuspente0->SetElement(isave, 0, P->Element(j, 0));
            interpXSuspente0->SetElement(isave, 1, Xd->Element(j, 0));
            interpYSuspente0->SetElement(isave, 1, Yd->Element(j, 0));
            isave++;
        }
		if (debug) printf ("\n interpYSuspente0->GetLignes()=%d", interpYSuspente0->GetLignes());
    }
    if (debug) printf ("\n ...3");

    int ndn = 0;
    double posda = 0.0f, posdf = 100.0f;
    for (int _d = 0; _d < gfd->quantDiag; _d++) {
		if (debug) printf ("\n_d=%d", _d);
        int dnerv = gfd->noNervD[_d] - 1;
        int f = 1;
        if ((nerv == dnerv) && (faceDeb == f) && (faceFin == f)) {
            ndn = 2;
            posda = gfd->PosDiagNerv1A;
            posdf = gfd->PosDiagNerv1F;
            break;
        }
        dnerv = gfd->noNervD[_d];
		if (debug) printf ("\n1.dnerv=%d", dnerv);
        f = 2;
        if ((nerv == dnerv) && (faceDeb == f) && (faceFin == f)) {
            ndn = 2;
            posda = gfd->PosDiagNerv2A;
            posdf = gfd->PosDiagNerv2F;
            break;

        }
        dnerv = gfd->noNervD[_d] + 1;
		if (debug) printf ("\n2.dnerv=%d", dnerv);
        f = 1;
        if ((nerv == dnerv) && (faceDeb == f) && (faceFin == f)) {
            ndn = 2;
            posda = gfd->PosDiagNerv1A;
            posdf = gfd->PosDiagNerv1F;
            break;
        }
		if (debug) printf ("\n4.dnerv=%d", dnerv);
    }
    if (debug) printf ("\n...4");
    int nrp;
    if (gfd->ReperPointsFromFile) nrp = gfd->ReperPoints[1] -> GetLignes();
    else nrp = P -> GetLignes();
    int n = 20 * nrp + 12;
    Matrix* res = new Matrix(n, 5);
    isave = 0;
    double posSuspente, x, y;
    int mark, MARK_SUSPENTE=REP_TRIANGLE, MARK_REP=REP_CROSS, MARK_DIAG=REP_LINE;

    if (gfd->ReperPointsPlotterFormat) {
        MARK_SUSPENTE = REP_V;
        MARK_REP = REP_MIDDLE_LINE;
        MARK_DIAG = REP_LINE;
    }
    if (debug) printf ("\n bef");
    for (int j = 0; j < 11; j++) {
		if (debug) printf ("\n j=%d");
        switch (j) {
            case 0: posSuspente = gfd->Form->m_pProfils[nerv]->m_fPosA; break;
            case 1: posSuspente = gfd->Form->m_pProfils[nerv]->m_fPosB; break;
            case 2: posSuspente = gfd->Form->m_pProfils[nerv]->m_fPosC; break;
            case 3: posSuspente = gfd->Form->m_pProfils[nerv]->m_fPosD; break;
            case 4: posSuspente = gfd->Form->m_pProfils[nerv]->m_fPosE; break;
            case 5: posSuspente = posda; break;
            case 6: posSuspente = posdf; break;
            case 7: posSuspente = 0.0f; break;
            case 8: posSuspente = gfd->VentHolesDeb; break;
            case 9: posSuspente = gfd->VentHolesFin; break;
            case 10: posSuspente = gfd->PosKlapanFin; break;
        }
        if (debug) printf ("...5.1");
        mark = MARK_SUSPENTE;
        if ((j == 5) || (j == 6)) {
            if (ndn == 0) continue;
            mark = MARK_DIAG;
        }
        if (debug) printf ("\n-- j=%d", j);
        if (j == 7) {
            if (!klapanShov) break;
            mark = REP_0;
            if (debug) printf ("\nposSuspente=%f isave=%d", posSuspente, isave);
        }
        if (j == 8) {
            if (!klapanShov) continue;
            mark=REP_VENTDEB;
            if (debug) printf ("\nposSuspente=%f isave=%d", posSuspente, isave);
        }
        if (j == 9) {
            if (!klapanShov) continue;
            mark=REP_VENTFIN;
            if (debug) printf ("\nposSuspente=%f isave=%d", posSuspente, isave);
        }
        if (j == 10) {
            if (!klapanShov) continue;
            mark=REP_KLAPAN_FIN;
            if (debug) printf ("\nposSuspente=%f isave=%d", posSuspente, isave);
        }
        
        if (debug) printf ("...5.2");

        if ((!cerezNos && (deb <= posSuspente) && (posSuspente <= fin))
                || (cerezNos && (posSuspente <= fin))) {
            x = InterpLinX(interpXSuspente, posSuspente);
            y = InterpLinX(interpYSuspente, posSuspente);
            res->SetElement(isave, 0, x);
            res->SetElement(isave, 1, y);
            res->SetElement(isave, 2, mark);
            if (!cerezNos && (deb <= posSuspente) && (posSuspente <= fin)) {
                res->SetElement(isave, 3, j+1);
            } else {
                res->SetElement(isave, 3, 0);
            }
			res->SetElement(isave, 4, posSuspente);
            isave++;
        }
        if (debug) printf ("\n...5.3 posSuspente=%f deb=%f", posSuspente, deb);
		if (debug) printf ("\n cerezNos: %d", cerezNos);
        //if (cerezNos) printf ("\ncerezNos=1"); else printf ("\ncerezNos=2");

        if ((cerezNos) && (posSuspente <= deb)) {
            x = InterpLinX(interpXSuspente0, posSuspente);
            y = InterpLinX(interpYSuspente0, posSuspente);
            res->SetElement(isave, 0, x);
            res->SetElement(isave, 1, y);
            res->SetElement(isave, 2, mark);
            res->SetElement(isave, 3, 0);
			res->SetElement(isave, 4, posSuspente);
            isave++;
        }
    }
    
    if (debug) printf ("...6");
    if (debug) if (cerezNos) printf ("\n cerezNos"); else printf ("\n NOT cerezNos");
    if (debug) printf ("\nbefore Reper Points from file");
    if (gfd->ReperPointsFromFile) {
        if (debug) printf ("...6.1");
        if (debug) printf ("\nIN Reper Points from file");
        if ((cerezNos) && (interpXSuspente0->GetLignes()>1)) {
            for (int j = 0; j < gfd->ReperPoints[faceDeb]->GetLignes(); j++) {
                posSuspente = gfd->ReperPoints[faceDeb]->Element(j, 0);
                if (debug) printf ("\ncZ: posSuspente=%f deb=%f", posSuspente, deb);
                if (posSuspente <= deb) {
                    if (debug) printf ("\nbefore InterpLinX");
                    x = InterpLinX(interpXSuspente0, posSuspente);
                    y = InterpLinX(interpYSuspente0, posSuspente);
                    if (debug) printf ("...after InterpLinX");
                    res->SetElement(isave, 0, x);
                    res->SetElement(isave, 1, y);
                    res->SetElement(isave, 2, MARK_REP);
                    res->SetElement(isave, 3, 0);
					res->SetElement(isave, 4, posSuspente);
                    isave++;
                }
            }
        }
        if (debug) printf ("\ngo not cerez nos");
		if (debug) printf ("\n gfd->ReperPoints[%d]->GetLignes(): %f", faceFin, gfd->ReperPoints[faceFin]->GetLignes());
        for (int j = 0; j < gfd->ReperPoints[faceFin]->GetLignes(); j++) {
			if (debug) printf ("\n j=%d", j);
            posSuspente = gfd->ReperPoints[faceFin]->Element(j, 0);
			if (debug) printf ("\n posSuspente: %f", posSuspente);
			if (debug) if (cerezNos) printf ("\n (cerezNos) true"); else printf ("\n (cerezNos) false");
            if ((!cerezNos && (deb <= posSuspente) && (posSuspente <= fin))
                    || (cerezNos && (posSuspente <= fin))) {
				if (debug) printf ("\n in if posSuspente: %f", posSuspente);
				if (debug) printf ("\n in if deb: %f fin: %f", deb, fin);

                x = InterpLinX(interpXSuspente, posSuspente);
                y = InterpLinX(interpYSuspente, posSuspente);
                res->SetElement(isave, 0, x);
                res->SetElement(isave, 1, y);
                res->SetElement(isave, 2, MARK_REP);
                res->SetElement(isave, 3, 0);
				res->SetElement(isave, 4, posSuspente);
                isave++;
            }
        }
    } else {
        if (debug) printf ("...6.2");
		if (debug) printf ("\nReper Points NOT FILE");
        if (cerezNos) {
            for (int j = i0; j > 0; j--) {
                posSuspente = P->Element(j, 0);
                x = InterpLinX(interpXSuspente0, posSuspente);
                y = InterpLinX(interpYSuspente0, posSuspente);
                res->SetElement(isave, 0, x);
                res->SetElement(isave, 1, y);
                res->SetElement(isave, 2, MARK_REP);
                res->SetElement(isave, 3, 0);
				res->SetElement(isave, 4, posSuspente);
                isave++;
            }
        }
        for (int j = i0; j < P->GetLignes(); j++) {
            posSuspente = P->Element(j, 0);
            x = InterpLinX(interpXSuspente, posSuspente);
            y = InterpLinX(interpYSuspente, posSuspente);
            res->SetElement(isave, 0, x);
            res->SetElement(isave, 1, y);
            res->SetElement(isave, 2, MARK_REP);
            res->SetElement(isave, 3, 0);
			res->SetElement(isave, 4, posSuspente);
            isave++;
        }
    }
    if (debug) printf ("...7");
    if (debug) printf ("\n isave=%d", isave);
    Matrix* resexit = new Matrix(isave, 5);
    for (int i = 0; i < isave; i++)
        for (int j = 0; j < 5; j++) resexit->SetElement(i, j, res->Element(i, j));
    if (debug) printf ("\n dels");
    delete (res);
    delete (interpXSuspente);
    delete (interpYSuspente);

    if (cerezNos) {
        delete (interpXSuspente0);
        delete (interpYSuspente0);
    }
   if (debug) printf ("\n...END get Reper Points()");
   return resexit;
}


double calcPolyLength(Matrix* X, Matrix* Y, double xrp, double yrp) {
	int i  = 0, res = 0;
	int n = X->GetLignes();
	double sumLength = 0;
	if ( (xrp==X->Element(0, 0)) && (yrp==Y->Element(0, 0)) ) return 0;

	if (X->Element(0,0) > X->Element(X->GetLignes() - 1,0)) getLayoutLogger()->logprintf ("\n calcPolyLength()... Achtung! X->Element(0,0) > X->Element(%d,0)", X->GetLignes() - 1);

	for (i = 0; i < n-1;i++ ) {
		res = pointAtSegment(xrp, yrp, X->Element(i,0), Y->Element(i,0), X->Element(i+1,0), Y->Element(i+1,0));
		if (1 == res ) break;
		sumLength += dist2d (X->Element(i,0), Y->Element(i,0), X->Element(i+1,0), Y->Element(i+1,0));
	}
	
	if (1 == res) 
		return sumLength + dist2d(xrp, yrp, X->Element(i,0), Y->Element(i,0));

	return -1;
}

Matrix* getReperPointsPince(WindPatternsProject* gfd,  Matrix* Xd0, Matrix* Yd0, double coeff, Matrix* Xd, Matrix* Yd, Matrix* P, 
								int nerv, float deb, int faceDeb, float fin, int faceFin, bool klapanShov)
{
	bool debug = true;
	if (debug) {
		printf ("getReperPointsPince()");
		printf ("\n Xd0->GetLignes()=%d", Xd0->GetLignes());
		printf ("\n Yd0->GetLignes()=%d", Yd0->GetLignes());
		printf ("\n Xd->GetLignes()=%d", Xd->GetLignes());
		printf ("\n Yd->GetLignes()=%d", Yd->GetLignes());
		printf ("\n P->GetLignes()=%d", P->GetLignes());
		printM(P);
		printf ("\n coeff=%f", coeff);
		printf ("\n nerv=%d deb=%f (%d)  fin=%f (%d)", nerv, deb, faceDeb, fin, faceFin);
	}

	// printf ("\n\n\n Xd0 Yd0");
	// printXY (Xd0, Yd0);

	// printf ("\n\n\n Xd Yd");
	// printXY (Xd, Yd);

	if (Xd0->GetLignes() != P->GetLignes()) {
		printf ("\n!!!\n!!!\n!!!ACHTUNG! getRepPointsPince() Xd0->GetL[%d] != P->GetL[%d]\n!!!", Xd0->GetLignes(), P->GetLignes());
	}
	if (debug) printf ("\n=0====call getReperPoints()===========================================");
	Matrix *reperPoints0 = getReperPoints(gfd, Xd0, Yd0, P, nerv, deb, faceDeb, fin, faceFin, klapanShov);
	if (debug) printf ("\n=...0====call getReperPoints()===========================================");
	if (debug) printf ("\n=1======call getReperPoints()==============================");
	Matrix *reperPointsP = getReperPoints(gfd, Xd, Yd, P, nerv, deb, faceDeb, fin, faceFin, klapanShov);
	if (debug) printf ("\n=...1======call getReperPoints()==============================");
	if (debug) printf ("\n=2===============================================");
	
	//printf ("\n\n Xd0 Yd0");
    //printXY(Xd0, Yd0);

	if (debug) printf ("\n(Xd0, Yd0).length=%f coeff=%f", 1000*calcCourbeLength(Xd0, Yd0), coeff);

	for (int i = 0; i < reperPoints0->GetLignes(); i++) {
		if (debug) printf ("\n ");
		double xrp = reperPoints0 -> Element(i, 0);
		double yrp = reperPoints0 -> Element(i, 1);

		double xrpP = reperPointsP -> Element(i, 0);
		double yrpP = reperPointsP -> Element(i, 1);

		double l = calcPolyLength(Xd0, Yd0, xrp, yrp);
		if (debug) printf ("\n l=%f l*coeff=%f", l, l*coeff);
		if (debug) printf ("\n xrp, yrp (%f, %f)", xrp, yrp);
		double lP0 = calcPolyLength(Xd, Yd, xrpP, yrpP);
		double newx = 0, newy = 0;
		if (debug) printf ("\n calcRepPointByLength(...coeff=%f", coeff);
		int res = calcRepPointByLength(Xd, Yd, l*coeff, &newx, &newy);
		if (debug) printf ("\n ...calcRepPointByLength(...coeff=%f", coeff);
		if (debug) printf ("\n res=%d newx, newy (%f, %f)", res, newx, newy);

		if (res) {
			if (debug) printf ("\n Ahtung! calcRepPointByLength() == %d ", res);
			if (debug) (" l=%f   coeff=%f", l, coeff);
			double XY0length = calcCourbeLength(Xd0, Yd0);
			double XYlength =  calcCourbeLength(Xd, Yd);
			double realCoeff = XYlength / XY0length;

			if (res == -1) {

			}

			if (res == 1) {
				if (debug) printf ("\n! XY0length=%f  XYlength=%f coeff=%f", XY0length, XYlength, realCoeff);
				if (debug) printf ("\n delta=%f", abs(XYlength-XY0length));
				if (abs(XYlength - XY0length) < 0.0005)  
					if (debug) printf(" OK:)!");
				else 
					if (debug) printf(" BAD:(!");
			}
		}
		
		//getLayoutLogger()->logprintf ("\n%2d R res=%d pS=%6.2f l=%6.2f (%6.2f, %6.2f)", i, res, reperPoints0->Element(i,4), 1000*l, 1000*xrp, 1000*yrp); 

		/*getLayoutLogger()->logprintf ("\n%d R res=%d pS=%f cf=%f newl=%f  dxyrp=%f  dl=%f",
										  i,     res, reperPoints0->Element(i,4), coeff, 1000*l*coeff, 1000*dist2d(xrpP,yrpP,newx,newy), 1000*(lP0-l*coeff));
		getLayoutLogger()->logprintf ("\n l==%f l*coeff=%f (%f, %f)", 1000*l, 1000*l*coeff, xrp, yrp);*/
		reperPoints0->SetElement(i, 0, newx);
		reperPoints0->SetElement(i, 1, newy);
	}

	return reperPoints0;
}


int calcRepPointByLength(Matrix* Xd, Matrix* Yd, double l, double* x, double* y) {
	bool debug = true;
	if (debug) printf ("\n calcRepPointByLength(... l=%f", l);
	if (l < 0.0f) {
		getLayoutLogger()->logprintf("\nERROR: !!! calcRepPointByLength l=%f return -1", l);
		return -1;
	}
	if (l == 0) {
		*x = Xd -> Element(0, 0);
		*y = Yd -> Element(0, 0);
		return 0;
	}
	int i = 0;
	double sumLength = 0;

	while ( (i < (Xd->GetLignes() - 1)) && (sumLength < l) ) {
		sumLength += dist2d (Xd->Element(i, 0), Yd->Element(i, 0), Xd->Element(i + 1, 0), Yd->Element(i + 1, 0));
		i++;
	}

	if ((i == (Xd->GetLignes() - 1)) && (sumLength < l)) {
		//getLayoutLogger()->logprintf("\n calcRepPointByLength l=%f sumLength=%f return 1 (put point to the end:))", l, sumLength);
		*x = Xd -> Element(i, 0);
		*y = Yd -> Element(i, 0);
		return 1;
	}

	double x1 = Xd->Element(i - 1, 0);
	double y1 = Yd->Element(i - 1, 0);
	double x2 = Xd->Element(i, 0);
	double y2 = Yd->Element(i, 0);

	double l2 = sumLength;
	double l1 = sumLength - dist2d (x1, y1, x2, y2);

	*x = x1 + (x2 - x1)*(l - l1) / (l2 - l1);
	*y = y1 + (y2 - y1)*(l - l1) / (l2 - l1);
	return 0;
}





