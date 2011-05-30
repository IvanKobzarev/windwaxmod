#include "reper.h"
#include "logger.h"

Matrix* getReperPoints(WindPatternsProject* gfd, Matrix* Xd, Matrix* Yd, Matrix* P, 
						int nerv, float deb, int faceDeb, float fin, int faceFin, bool klapanShov)
{
    //printf ("\nget Reper Points()");
	//printf ("\n P:");
	//printM (P);
	
	//printf ("\nP");
	//printM(P);
	//printf ("\n...");
    Matrix *interpXSuspente, *interpYSuspente, *interpXSuspente0, *interpYSuspente0;
    bool cerezNos = false;
    int i0 = 0;
    if ((faceDeb != faceFin) && (deb != 0.0f)) {
        cerezNos = true;
      //  printf ("\n cerezNos");
        while (P->Element(i0, 0) != 0.0f) {
        //    printf ("\n i0=%d", i0);
            i0++;
         }
     }
    //printf ("...1 i0=%d", i0);
    interpXSuspente = Zeros(P -> GetLignes() - i0, 2);
    interpYSuspente = Zeros(P -> GetLignes() - i0, 2);
    for (int j = i0; j < P->GetLignes(); j++) {
        interpXSuspente->SetElement(j - i0, 0, P->Element(j, 0));
        interpYSuspente->SetElement(j - i0, 0, P->Element(j, 0));
        interpXSuspente->SetElement(j - i0, 1, Xd->Element(j, 0));
        interpYSuspente->SetElement(j - i0, 1, Yd->Element(j, 0));
    }
    //printf ("...2");
    int isave = 0;
    if (cerezNos) {
		//printf ("\n cerezNos");
        interpXSuspente0 = Zeros(i0 + 1, 2);
        interpYSuspente0 = Zeros(i0 + 1, 2);
        for (int j = i0; j >= 0; j--) {
            interpXSuspente0->SetElement(isave, 0, P->Element(j, 0));
            interpYSuspente0->SetElement(isave, 0, P->Element(j, 0));
            interpXSuspente0->SetElement(isave, 1, Xd->Element(j, 0));
            interpYSuspente0->SetElement(isave, 1, Yd->Element(j, 0));
            isave++;
        }
    }
    //printf ("...3");
    //printf ("\n interpYSuspente0->GetLignes()=%d", interpYSuspente0->GetLignes());
    int ndn = 0;
    double posda = 0.0f, posdf = 100.0f;
    for (int _d = 0; _d < gfd->quantDiag; _d++) {
        int dnerv = gfd->noNervD[_d] - 1;
        int f = 1;
        if ((nerv == dnerv) && (faceDeb == f) && (faceFin == f)) {
            ndn = 2;
            posda = gfd->PosDiagNerv1A;
            posdf = gfd->PosDiagNerv1F;
            break;
        }
        dnerv = gfd->noNervD[_d];
        f = 2;
        if ((nerv == dnerv) && (faceDeb == f) && (faceFin == f)) {
            ndn = 2;
            posda = gfd->PosDiagNerv2A;
            posdf = gfd->PosDiagNerv2F;
            break;

        }
        dnerv = gfd->noNervD[_d] + 1;
        f = 1;
        if ((nerv == dnerv) && (faceDeb == f) && (faceFin == f)) {
            ndn = 2;
            posda = gfd->PosDiagNerv1A;
            posdf = gfd->PosDiagNerv1F;
            break;

        }
    }
    //printf ("...4");
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
    
    for (int j = 0; j < 11; j++) {
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
        //printf ("...5.1");
        mark = MARK_SUSPENTE;
        if ((j == 5) || (j == 6)) {
            if (ndn == 0) continue;
            mark = MARK_DIAG;
        }
        //printf ("\n-- j=%d", j);
        if (j == 7) {
            if (!klapanShov) break;
            mark=REP_0;
            //printf ("\nposSuspente=%f isave=%d", posSuspente, isave);
        }
        if (j == 8) {
            if (!klapanShov) continue;
            mark=REP_VENTDEB;
            //printf ("\nposSuspente=%f isave=%d", posSuspente, isave);
        }
        if (j == 9) {
            if (!klapanShov) continue;
            mark=REP_VENTFIN;
            //printf ("\nposSuspente=%f isave=%d", posSuspente, isave);
        }
        if (j == 10) {
            if (!klapanShov) continue;
            mark=REP_KLAPAN_FIN;
            //printf ("posSuspente=%f isave=%d", posSuspente, isave);
        }
        
        //printf ("...5.2");

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
        //printf ("...5.3 posSuspente=%f deb=%f", posSuspente, deb);
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
    
    //printf ("...6");
    //if (cerezNos) printf ("\n cerezNos"); else printf ("\n NOT cerezNos");
    //printf ("\nbefore Reper Points from file");
    if (gfd->ReperPointsFromFile) {
        //printf ("...6.1");
        //printf ("\nIN Reper Points from file");
        if ((cerezNos) && (interpXSuspente0->GetLignes()>1)) {
            for (int j = 0; j < gfd->ReperPoints[faceDeb]->GetLignes(); j++) {
                posSuspente = gfd->ReperPoints[faceDeb]->Element(j, 0);
                //printf ("\ncZ: posSuspente=%f deb=%f", posSuspente, deb);
                if (posSuspente <= deb) {
                    //printf ("\nbefore InterpLinX");
                    x = InterpLinX(interpXSuspente0, posSuspente);
                    y = InterpLinX(interpYSuspente0, posSuspente);
                    //printf ("...after InterpLinX");
                    res->SetElement(isave, 0, x);
                    res->SetElement(isave, 1, y);
                    res->SetElement(isave, 2, MARK_REP);
                    res->SetElement(isave, 3, 0);
					res->SetElement(isave, 4, posSuspente);
                    isave++;
                }
            }
        }
        //printf ("\nnot cerez nos");
        for (int j = 0; j < gfd->ReperPoints[faceFin]->GetLignes(); j++) {
            posSuspente = gfd->ReperPoints[faceFin]->Element(j, 0);
            if ((!cerezNos && (deb <= posSuspente) && (posSuspente <= fin))
                    || (cerezNos && (posSuspente <= fin))) {
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
        //printf ("...6.2");
      //  printf ("\nReper Points NOT FILE");
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
    //printf ("...7");
    //printf ("\n isave=%d", isave);
    Matrix* resexit = new Matrix(isave, 5);
    for (int i = 0; i < isave; i++)
        for (int j = 0; j < 5; j++) resexit->SetElement(i, j, res->Element(i, j));
    //printf ("\n dels");
    delete (res);
    delete (interpXSuspente);
    delete (interpYSuspente);

    if (cerezNos) {
        delete (interpXSuspente0);
        delete (interpYSuspente0);
    }
   //printf ("\n...END get Reper Points()");
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
	// printf ("\n\n\n Xd0 Yd0");
	// printXY (Xd0, Yd0);

	// printf ("\n\n\n Xd Yd");
	// printXY (Xd, Yd);

	if (Xd0->GetLignes() != P->GetLignes()) {
		printf ("\n!!!\n!!!\n!!!ACHTUNG! getRepPointsPince() Xd0->GetL[%d] != P->GetL[%d]\n!!!", Xd0->GetLignes(), P->GetLignes());
	}

	Matrix *reperPoints0 = getReperPoints(gfd, Xd0, Yd0, P, nerv, deb, faceDeb, fin, faceFin, klapanShov);

	//printf ("\n=1===============================================");
	Matrix *reperPointsP = getReperPoints(gfd, Xd, Yd, P, nerv, deb, faceDeb, fin, faceFin, klapanShov);
	//printf ("\n=2===============================================");
	
	//printf ("\n\n Xd0 Yd0");
    //printXY(Xd0, Yd0);

	printf ("\n(Xd0, Yd0).length=%f coeff=%f", 1000*calcCourbeLength(Xd0, Yd0), coeff);

	for (int i = 0; i < reperPoints0->GetLignes(); i++) {
		// printf ("\n ");
		double xrp = reperPoints0 -> Element(i, 0);
		double yrp = reperPoints0 -> Element(i, 1);

		double xrpP = reperPointsP -> Element(i, 0);
		double yrpP = reperPointsP -> Element(i, 1);

		double l = calcPolyLength(Xd0, Yd0, xrp, yrp);
		// printf ("\n l=%f l*coeff=%f", l, l*coeff);
		// printf ("\n xrp, yrp (%f, %f)", xrp, yrp);
		double lP0 = calcPolyLength(Xd, Yd, xrpP, yrpP);

		double newx = 0, newy = 0;
		
		int res = calcRepPointByLength(Xd, Yd, l*coeff, &newx, &newy);

		// printf ("\n res=%d newx, newy (%f, %f)", res, newx, newy);

		if (res) {
			//getLayoutLogger()->logprintf ("\n Ahtung! calcRepPointByLength() == %d ", res);
			//getLayoutLogger()->logprintf (" l=%f   coeff=%f", l, coeff);
			double XY0length =  calcCourbeLength(Xd0, Yd0);
			double XYlength =  calcCourbeLength(Xd, Yd);
			double realCoeff = XYlength / XY0length;

			if (res == -1) {

			}

			if (res == 1) {
				//getLayoutLogger()->logprintf ("\n! XY0length=%f  XYlength=%f coeff=%f", XY0length, XYlength, realCoeff);
				//getLayoutLogger()->logprintf ("\n delta=%f", abs(XYlength-XY0length));
				/*if (abs(XYlength - XY0length) < 0.0005)  
					getLayoutLogger()->logprintf(" OK:)!");
				else 
					getLayoutLogger()->logprintf(" BAD:(!");*/
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

	if (l < 0) {
		getLayoutLogger()->logprintf("\n!!! calcRepPointByLength l=%f return -1", l);
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





