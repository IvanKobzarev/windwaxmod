#include "profil.h"
#include "form.h"

Profil::Profil()
{
	m_fLength = 0.0f;
	m_fWidth = 0.0f;
	m_fNezX = 0.0f;
	m_fNezY = 0.0f;
	m_fNezZ = 0.0f;
	m_fInclin = 0.0f;
	m_fWash = 0.0f;
	m_fMorph = 0.0f;
	m_fPosA = 0.0f;
	m_fPosB = 0.0f;
	m_fPosC = 0.0f;
	m_fPosD = 0.0f;
	m_fPosE = 0.0f;
}

Profil::Profil(
			double lg,		//longueur
			double width,	//epaisseur relative
			double xnez,		//xnez
			double ynez,		//ynez
			double znez,		//znez
			double incl,		//inclinaison / horizontale
			double wash,		//vrillage / horizontale
			double morph,	//morphing / nervure centrale
			double posA,		//pos relative ligne A
			double posB,		//pos relative ligne B
			double posC,		//pos relative ligne C
			double posD,		//pos relative ligne D
			double posE		//pos relative ligne E
		)
{
	m_fLength = lg;
	m_fWidth = width;
	m_fNezX = xnez;
	m_fNezY = ynez;
	m_fNezZ = znez;
	m_fInclin = incl;
	m_fWash = wash;
	m_fMorph = morph;
	m_fPosA = posA;
	m_fPosB = posB;
	m_fPosC = posC;
	m_fPosD = posD;
	m_fPosE = posE;
}

Profil::~Profil()
{}

void ProfilGeom::print() {
	printf ("\nProfilGeom->print()");
	printf ("\n ExtProf");
	for (int i = 0; i < ExtProf->GetLignes(); i++) {
		printf ("\n%d (%f, %f)", i, ExtProf->Element(i,0), ExtProf->Element(i,1));
	}
	printf ("\n ");
	printf ("\n IntProf");
	for (int i = 0; i < IntProf->GetLignes(); i++) {
		printf ("\n%d (%f, %f)", i, IntProf->Element(i,0), IntProf->Element(i,1));
	}
}

/*ProfilGeom* getHvostByPince(ProfilGeom* pg0, double posF, double pincePower, double Amp, int* ie, int* ii) {
	ProfilGeom* pg = new ProfilGeom();
	int ne = pg0->ExtProf->GetLignes();
	int ni = pg0->IntProf->GetLignes();
	pg -> ExtProf = new Matrix (ne, 2);
	pg -> IntProf = new Matrix (ni, 2);

	

	for (int i = 0; i < ne; i++) {
		double xi0 = pg0 -> ExtProf->Element(i, 0);
		double yi0 = pg0 -> ExtProf->Element(i, 1);
		if (xi0 >= posF) {
			double xi = x0_ + xi0 * l_/l;
			double yi = y0_ + yi0 * wN/w0;
			pg-> ExtProf->SetElement(i, 0, xi);
			pg-> ExtProf->SetElement(i, 1, yi);
		} else {
			pg-> ExtProf->SetElement(i, 0, xi0);
			pg-> ExtProf->SetElement(i, 1, yi0);
		}
	}
	for (int i = 0; i < ni; i++) {
		double xi0 = pg0 -> IntProf->Element(i, 0);
		double yi0 = pg0 -> IntProf->Element(i, 1);
		if (xi0 == posF) xi0F = xi0;
		if (xi0 >= posF) {
			
			double xi = x0_ + xi0 * l_/l;
			double yi = y0_ + yi0 * wN/w0;

			// here need to use pince function
			double amp = Amp * functionPincePower((IntProf->Element(ni - 1, 0) - xi)/(xi0F-xi), (1.0f/pincePower));



			pg-> IntProf->SetElement(i, 0, xi);
			pg-> IntProf->SetElement(i, 1, yi);
		} else {
			pg-> IntProf->SetElement(i, 0, xi0);
			pg-> IntProf->SetElement(i, 1, yi0);
		}
	}
	

	

	return pg;
}
*/

ProfilGeom* getBalloneProfilGeom(ProfilGeom* pg0, double kChord, double kMf, double w0, double wN, double dyw) {
	ProfilGeom* pg = new ProfilGeom();

	int ne = pg0->ExtProf->GetLignes();
	int ni = pg0->IntProf->GetLignes();
	pg -> ExtProf = new Matrix (ne, 2);
	pg -> IntProf = new Matrix (ni, 2);
	
	double l = abs (pg0->ExtProf->Element(ne-1, 0) - pg0->ExtProf->Element(0, 0));

	double l_ = l * kChord;
	double dl = l * (kChord - 1.0);


	double x0_ = pg0 -> ExtProf -> Element(0, 0) - dl * kMf; 
	double wabs0 = EpaisseurRelative(pg0 -> ExtProf, pg0 -> IntProf);
	double y0_ = wabs0 * dyw/w0;

	for (int i = 0; i < ne; i++) {
		double xi0 = pg0 -> ExtProf->Element(i, 0);
		double yi0 = pg0 -> ExtProf->Element(i, 1);
		
		double xi = x0_ + xi0 * l_/l;
		double yi = y0_ + yi0 * wN/w0;

		pg-> ExtProf->SetElement(i, 0, xi);
		pg-> ExtProf->SetElement(i, 1, yi);
	}
	for (int i = 0; i < ni; i++) {
		double xi0 = pg0 -> IntProf->Element(i, 0);
		double yi0 = pg0 -> IntProf->Element(i, 1);

		double xi = x0_ + xi0 * l_/l;
		double yi = y0_ + yi0 * wN/w0;

		pg-> IntProf->SetElement(i, 0, xi);
		pg-> IntProf->SetElement(i, 1, yi);
	}

	return pg;
}

ProfilGeom* getProfilGeomTailDown(ProfilGeom* pg1, ProfilGeom* pg0, double xv, double power) {
	int ne1 = pg1->ExtProf->GetLignes();
	int ni1 = pg1->IntProf->GetLignes();
	int ne0 = pg0->ExtProf->GetLignes();
	int ni0 = pg0->IntProf->GetLignes();
	//printf ("\n ne1=%d ni1=%d", ne1, ni1);
	//printf ("\n ne0=%d ni0=%d", ne0, ni0);

	int iExtTail1 = -1, iIntTail1 = -1, iExtTail0 = -1, iIntTail0 = -1;
 
	// coping pg1
	for (int i = 0; i < ne1; i++) {
		double xi1 = pg1 -> ExtProf->Element(i, 0);
		double yi1 = pg1 -> ExtProf->Element(i, 1);
		if (xi1 <= xv ) iExtTail1 = i;
	}
	for (int i = 0; i < ni1; i++) {
		double xi1 = pg1 -> IntProf->Element(i, 0);
		double yi1 = pg1 -> IntProf->Element(i, 1);
		if (xi1 <= xv ) iIntTail1 = i;
	}

	for (int i = 0; i < ne0; i++) {
		double xi = pg0 -> ExtProf->Element(i, 0);
		double yi = pg0 -> ExtProf->Element(i, 1);
		if (xi <= xv ) iExtTail0 = i;
	}
	for (int i = 0; i < ni0; i++) {
		double xi = pg0 -> IntProf->Element(i, 0);
		double yi = pg0 -> IntProf->Element(i, 1);
		if (xi <= xv ) iIntTail0 = i;
	}

	//printf ("\n iExtTail1=%d, iIntTail1=%d, iExtTail0=%d, iIntTail0=%d",  iExtTail1 , iIntTail1 , iExtTail0 , iIntTail0 );

	// now pg1h equal pg1

	// need to make tail, good tail

	// 1) find index where tail starts [ iExtTail, iIntTail ]
	// make equal quantity tail points
	ProfilGeom* pgtmp = new ProfilGeom();
	pgtmp -> ExtProf = new Matrix (iExtTail1 + ne0 - iExtTail0, 2);
	pgtmp -> IntProf = new Matrix (iIntTail1 + ni0 - iIntTail0, 2);

	for (int i = 0; i <= iExtTail1; i++ ) {
		pgtmp->ExtProf->SetElement(i, 0, pg1->ExtProf->Element(i, 0));
		pgtmp->ExtProf->SetElement(i, 1, pg1->ExtProf->Element(i, 1));
	}

	for (int i = iExtTail0 + 1; i < ne0; i++ ) {
		pgtmp->ExtProf->SetElement(iExtTail1 + i - iExtTail0, 0, pg0->ExtProf->Element(i, 0));
		pgtmp->ExtProf->SetElement(iExtTail1 + i - iExtTail0, 1, pg0->ExtProf->Element(i, 1));
	}

	for (int i = 0; i <= iIntTail1; i++ ) {
		pgtmp->IntProf->SetElement(i, 0, pg1->IntProf->Element(i, 0));
		pgtmp->IntProf->SetElement(i, 1, pg1->IntProf->Element(i, 1));
	}

	for (int i = iIntTail0 + 1; i < ni0; i++ ) {
		pgtmp->IntProf->SetElement(iIntTail1 + i - iIntTail0, 0, pg0->IntProf->Element(i, 0));
		pgtmp->IntProf->SetElement(iIntTail1 + i - iIntTail0, 1, pg0->IntProf->Element(i, 1));
	}

	// constructing Ext part

	// amp ?

    double xn, yn, x , y, ampExt = -1000.0f;
	double x0 = pg1->ExtProf->Element(iExtTail1, 0);
    double y0 = pg1->ExtProf->Element(iExtTail1, 1);
	calcVecteurBissec(pg1->ExtProf->Element(iExtTail1-1, 0), pg1->ExtProf->Element(iExtTail1-1, 1),
					x0, y0,
					pg1->ExtProf->Element(iExtTail1+1, 0), pg1->ExtProf->Element(iExtTail1+1, 1),
					&xn, &yn, 10.0f, +1);
	for (int j = 0; j < ne0; j++) {
		double x1 = pg0->ExtProf->Element(j, 0);
        double y1 = pg0->ExtProf->Element(j, 1);
        double x2 = pg0->ExtProf->Element(j + 1,0);
        double y2 = pg0->ExtProf->Element(j + 1,1);
        Inter2Vecteurs(x0, y0, xn, yn, x1, y1, x2, y2, &x, &y);
        if ((x >= x1) && (x <= x2)) {
            ampExt= dist2d (x,y,x0,y0);
            break;
        }
    }

	//amp est'!
	//printf ("\n ampExt=%f", ampExt);

	double xbegin = pg0->ExtProf->Element(iExtTail0, 0);
	double xend = pg0->ExtProf->Element(ne0-1, 0);
	int ne1tail = iExtTail1 + ne0 - iExtTail0;
	for (int i = iExtTail1 + 1; i < ne1tail-1; i++ ) {
		//pg0->ExtProf->Element(iExtTail0, 0) -> ampExt
		//pg0->ExtProf->Element(i, 0) -> amp? ampExt * pow (xend-xc/xend-xi)
		int i0 = iExtTail0+i-iExtTail1;
		//printf ("\n i=%d i0=%d", i, i0);
		double xc = pg0->ExtProf->Element(i0, 0);
		//pg0->ExtProf->Element(ne0-1, 0) -> 0
		double x = 0, y =0; 
		double amp = ampExt * pow ( (xend-xc)/(xend-xbegin), 1/power);
		calcVecteurBissec(pg0->ExtProf->Element(i0-1, 0), pg0->ExtProf->Element(i0-1, 1),
						pg0->ExtProf->Element(i0, 0), pg0->ExtProf->Element(i0, 1),
						pg0->ExtProf->Element(i0+1, 0), pg0->ExtProf->Element(i0+1, 1),
						&x, &y, amp, +1);
		pgtmp->ExtProf->SetElement(i, 0, x);
		pgtmp->ExtProf->SetElement(i, 1, y);
	}

	pgtmp->ExtProf->SetElement(ne1tail-1, 0, pg0->ExtProf->Element(ne0-1, 0));
	pgtmp->ExtProf->SetElement(ne1tail-1, 1, pg0->ExtProf->Element(ne0-1, 1));
	
	// constructing Int part
    double ampInt = -1000.0f;
	x0 = pg1->IntProf->Element(iIntTail1, 0);
    y0 = pg1->IntProf->Element(iIntTail1, 1);
	calcVecteurBissec(pg1->IntProf->Element(iIntTail1-1, 0), pg1->IntProf->Element(iIntTail1-1, 1),
					x0, y0,
					pg1->IntProf->Element(iIntTail1+1, 0), pg1->IntProf->Element(iIntTail1+1, 1),
					&xn, &yn, 10.0f, +1);
	for (int j = 0; j < ni0; j++) {
		double x1 = pg0->IntProf->Element(j, 0);
        double y1 = pg0->IntProf->Element(j, 1);
        double x2 = pg0->IntProf->Element(j + 1,0);
        double y2 = pg0->IntProf->Element(j + 1,1);
        Inter2Vecteurs(x0, y0, xn, yn, x1, y1, x2, y2, &x, &y);
        if ((x >= x1) && (x <= x2)) {
            ampInt= dist2d (x,y,x0,y0);
            break;
        }
    }
	//amp est'!
	//printf ("\n ampInt=%f", ampInt);

	xbegin = pg0->IntProf->Element(iIntTail0, 0);
	xend = pg0->IntProf->Element(ni0-1, 0);
	int ni1tail = iIntTail1 + ni0 - iIntTail0;
	for (int i = iIntTail1 + 1; i < ni1tail-1; i++ ) {
		//pg0->ExtProf->Element(iExtTail0, 0) -> ampExt
		//pg0->ExtProf->Element(i, 0) -> amp? ampExt * pow (xend-xc/xend-xi)
		int i0 = iIntTail0+i-iIntTail1;
		//printf ("\n i=%d i0=%d", i, i0);
		double xc = pg0->IntProf->Element(i0, 0);
		//pg0->ExtProf->Element(ne0-1, 0) -> 0
		double x = 0, y =0; 
		double amp = ampInt * pow ( (xend-xc)/(xend-xbegin), 1/power);
		calcVecteurBissec(pg0->IntProf->Element(i0-1, 0), pg0->IntProf->Element(i0-1, 1),
						pg0->IntProf->Element(i0, 0), pg0->IntProf->Element(i0, 1),
						pg0->IntProf->Element(i0+1, 0), pg0->IntProf->Element(i0+1, 1),
						&x, &y, amp, -1);
		pgtmp->IntProf->SetElement(i, 0, x);
		pgtmp->IntProf->SetElement(i, 1, y);
	}

	pgtmp->IntProf->SetElement(ni1tail-1, 0, pg0->IntProf->Element(ni0-1, 0));
	pgtmp->IntProf->SetElement(ni1tail-1, 1, pg0->IntProf->Element(ni0-1, 1));



/*
	double dw = pg1 -> ExtProf -> Element (ne-1, 1);
	printf ("\n dw=%f", dw);
	bool found = false;
	double xf = 0.0, yf = 0.0;

	for (int i = 0; i < ne; i++) {
		double xi0 = pg1 -> ExtProf->Element(i, 0);
		double yi0 = pg1 -> ExtProf->Element(i, 1);
		
		if (found) {
			double xi = xi0;
			double yi = (yi0-dw)*yf/(yf-dw);
			pg1h -> ExtProf->SetElement(i, 0, xi);
			pg1h -> ExtProf->SetElement(i, 1, yi);
		} else {
			pg1h -> ExtProf->SetElement(i, 0, xi0);
			pg1h -> ExtProf->SetElement(i, 1, yi0);
		}

		if (xi0 >= xv) {
			found = true;
			xf = xi0;
			yf = yi0;
		}
	}
	found = false;
	for (int i = 0; i < ni; i++) {
		double xi0 = pg1 -> IntProf->Element(i, 0);
		double yi0 = pg1 -> IntProf->Element(i, 1);

		
		if (found) {
			double xi = xi0;
			double yi = yf*(dw-yi0)/(dw + abs(yf));
			pg1h -> IntProf->SetElement(i, 0, xi);
			pg1h -> IntProf->SetElement(i, 1, yi);
		} else {
			pg1h -> IntProf->SetElement(i, 0, xi0);
			pg1h -> IntProf->SetElement(i, 1, yi0);
		}


		if (xi0 >= xv) {
			found = true;
			xf = xi0;
			yf = yi0;
		}

	}
*/
	return pgtmp;
}

void GetProfileXY(WindPatternsProject* gfd, Form* F, int nerv, int face, Matrix** XProf, Matrix** YProf) 
{
	//printf ("\n getProfileXY()");
    Matrix *XExt, *YExt, *ZExt, *XInt, *YInt, *ZInt;
	//printf ("\n 1");
    calcForm3D(F, 0, 0.0f,
		gfd->ExtProfCent, gfd->IntProfCent, gfd->ExtProfBout, gfd->IntProfBout, &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);
	//printf ("\n 1.1");
    double LongNerv = 100.0f;
	//LongNerv = F->m_pProfils[nerv]->m_fLength;
	LongNerv = 100.0f;
	double coeffx = LongNerv/100.0f;
	//coeffx = 1.0;
	//printf ("\n 2");
    double EpaiRel = 0.0f;
    double m = -1.0f;
    if (nerv != -1) {
        EpaiRel = F->m_pProfils[nerv]->m_fWidth;
        m = F->m_pProfils[nerv]->m_fMorph;
    } else {
        EpaiRel = F->m_pProfils[0]->m_fWidth;
        m = 1.0f;
    }
	//printf ("\n 3");
    double EpaiRelProfCent = EpaisseurRelative(gfd->ExtProfCent, gfd->IntProfCent);
    double EpaiRelProfBout = EpaisseurRelative(gfd->ExtProfBout, gfd->IntProfBout);
	//printf ("\n inGetProfileXY: EpaiRelProfCent=%f", EpaiRelProfCent);
	//printf ("\n inGetProfileXY: EpaiRelProfBout=%f", EpaiRelProfBout);
    double coeffyCent = LongNerv * EpaiRel / (EpaiRelProfCent * 100.0f);
    double coeffyBout = LongNerv * EpaiRel / (EpaiRelProfBout * 100.0f);
	//printf ("\n coeffx=%f, coeffyCent=%f, coeffyBout=%f", coeffx, coeffyCent, coeffyBout);
    double xp, yp;
    int n = 0, j = 0;
    if (face == 1) n = gfd->ExtProfCent->GetLignes(); else n = gfd->IntProfCent->GetLignes();
    *XProf = new Matrix(n, 1);
    *YProf = new Matrix(n, 1);
	//printf ("\n 4");
    if (face == 1) {
        for (j = 0; j < gfd->ExtProfCent->GetLignes(); j++) {
            xp = gfd->ExtProfCent->Element(j, 0) * coeffx;
            yp = gfd->ExtProfCent->Element(j, 1) * coeffyCent * m + gfd->ExtProfBout->Element(j, 1) * coeffyBout * (1.0f - m);
            (*XProf)->SetElement(j, 0, xp);
            (*YProf)->SetElement(j, 0, yp);
        }
    } else {
        for (j = 0; j < gfd->IntProfCent->GetLignes(); j++) {
            xp = gfd->IntProfCent->Element(j, 0) * coeffx;
            yp = gfd->IntProfCent->Element(j, 1) * coeffyCent * m + gfd->IntProfBout->Element(j, 1) * coeffyBout * (1.0f - m);
            (*XProf)->SetElement(j, 0, xp);
            (*YProf)->SetElement(j, 0, yp);
        }
    }
	//printf ("\n ..getProfileXY()");
}

ProfilGeom* getProfile(WindPatternsProject* gfd, Form* F, int nerv) {
	//printf ("\n getProfileGeom()");
	ProfilGeom* pg = new ProfilGeom();

	Matrix* XExt, *YExt, *XInt, *YInt;
	GetProfileXY(gfd, F, nerv, 1, &XExt, &YExt);
	int n = XExt->GetLignes();
	pg->ExtProf = new Matrix (n,2);
	for (int i = 0; i < n; i++ ) {
		pg->ExtProf->SetElement(i, 0, XExt->Element(i, 0));
		pg->ExtProf->SetElement(i, 1, YExt->Element(i, 0));
	}

	GetProfileXY(gfd, F, nerv, 2, &XInt, &YInt);
	n = XInt->GetLignes();
	pg->IntProf = new Matrix (n,2);
	for (int i = 0; i < n; i++ ) {
		pg->IntProf->SetElement(i, 0, XInt->Element(i, 0));
		pg->IntProf->SetElement(i, 1, YInt->Element(i, 0));
	}
	//printf ("\n ...getProfileGeom()");
	return pg;
}

void GetMiddleProfile(WindPatternsProject* gfd, Form* F, int nerv1, int nerv2, int face, int realMashtab, Matrix** XProf, Matrix** YProf) {
    /*Matrix *XExt, *YExt, *ZExt, *XInt, *YInt, *ZInt;
    calcForm3D(F, 0, 0.0f,
		gfd->ExtProfCent, gfd->IntProfCent, gfd->ExtProfBout, gfd->IntProfBout, &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);*/
    int i1 = nerv1;
    int i2 = nerv2;
	if (i1 == -1) i1 = 0;
	if (i2 == -1) i2 = 0;

	double LongNerv = 100.0f;
	if (realMashtab == 1) {
		LongNerv = 0.5f*(F->m_pProfils[i1]->m_fLength + F->m_pProfils[i2]->m_fLength);
	} 
    double EpaiRel = 0.0f;
    double m = -1.0f;
    if (i1 != -1) {
        EpaiRel = 0.5f * (F->m_pProfils[i1]->m_fWidth + F->m_pProfils[i2]->m_fWidth);
        m = 0.5f * (F->m_pProfils[i1]->m_fMorph + F->m_pProfils[i2]->m_fMorph);
    } else {
        EpaiRel = F->m_pProfils[0]->m_fWidth;
        m = 1.0f;
    }
    double EpaiRelProfCent = EpaisseurRelative(gfd->ExtProfCent, gfd->IntProfCent);
    double EpaiRelProfBout = EpaisseurRelative(gfd->ExtProfBout, gfd->IntProfBout);
    double coeffyCent = LongNerv * EpaiRel / (EpaiRelProfCent * 100.0f);
    double coeffyBout = LongNerv * EpaiRel / (EpaiRelProfBout * 100.0f);
	double coeffx = LongNerv/100.0f;
	//printf ("\n in getMiddleProfile() LongNerv=%f", LongNerv);
    double xp, yp;
    int n = 0, j = 0;
    if (face == 1) n = gfd->ExtProfCent->GetLignes(); else n = gfd->IntProfCent->GetLignes();
    *XProf = new Matrix(n, 1);
    *YProf = new Matrix(n, 1);
    if (face == 1) {
        for (j = 0; j < gfd->ExtProfCent->GetLignes(); j++) {
            xp = gfd->ExtProfCent->Element(j, 0) * coeffx;
            yp = gfd->ExtProfCent->Element(j, 1) * coeffyCent * m + gfd->ExtProfBout->Element(j, 1) * coeffyBout * (1.0f - m);
            (*XProf)->SetElement(j, 0, xp);
            (*YProf)->SetElement(j, 0, yp);
        }
    } else {
        for (j = 0; j < gfd->IntProfCent->GetLignes(); j++) {
            xp = gfd->IntProfCent->Element(j, 0) * coeffx;
            yp = gfd->IntProfCent->Element(j, 1) * coeffyCent * m + gfd->IntProfBout->Element(j, 1) * coeffyBout * (1.0f - m);
            (*XProf)->SetElement(j, 0, xp);
            (*YProf)->SetElement(j, 0, yp);
        }
    }

}

void GetMiddleProfileBal(WindPatternsProject* gfd, Form* F, int nerv1, int nerv2, int face,  int realMashtab, Matrix** XProf, Matrix** YProf) {
	double EpaiRel, xp, yp;
	double EpaiRelProfCent, EpaiRelProfBout;
	double coeffx, coeffy, coeffyCent, coeffyBout;
	int i1, i2, j;
	i1 = nerv1;
	if (i1 == -1) i1=0;
	i2 = nerv2;
	if (i2 == -1) i2=0;
    //bool isCenterPanel = (1 & forme->NbCaiss);
	double LongNerv = 100.0f;
	if (realMashtab == 1) {
		LongNerv = (F->m_pProfils[i1]->m_fLength + F->m_pProfils[i2]->m_fLength) * 0.5f;
	} 
	EpaiRel = (F->m_pProfils[i1]->m_fWidth + F->m_pProfils[i2]->m_fWidth) * 0.5f;
    if (i1 != -1) {
        EpaiRel = 0.5f * (F->m_pProfils[i1]->m_fWidth + F->m_pProfils[i2]->m_fWidth);
    } else {
        EpaiRel = F->m_pProfils[0]->m_fWidth;
    }
	ProfilGeom* pgCur = getProfile(gfd, F, i1);
	ProfilGeom* pgCurBal = getBalloneProfilGeom(pgCur, gfd->ballonement->kChord->Element(i1, 0), gfd->ballonement->kMf->Element(i1, 0), EpaiRel, gfd->ballonement->wN->Element(i1, 0), gfd->ballonement->dyw->Element(i1, 0));
	double l = abs (pgCur->ExtProf->Element(pgCur->ExtProf->GetLignes() - 1, 0) - pgCur->ExtProf->Element(0, 0));
	double xv = l * (100.0f - (gfd->PosPinceBF[0])) * 0.01f;
	ProfilGeom* pg = getProfilGeomTailDown(pgCurBal, pgCur, xv, gfd->ballonement->powerTail->Element(i1, 0));
	coeffx = LongNerv/100.0f;
	double EpaiRelCur = EpaisseurRelative(pgCur->ExtProf, pgCur->IntProf);
	coeffy = LongNerv*EpaiRel/(EpaiRelCur*100.0f);

    int n = 0;
    if (face == 1) n = pg->ExtProf->GetLignes(); else n = pg->IntProf->GetLignes();
    *XProf = new Matrix(n, 1);
    *YProf = new Matrix(n, 1);
	if (face == 1) {
		for (j=0; j < pg->ExtProf->GetLignes(); j++)
		{
			xp = pg->ExtProf->Element(j, 0) * coeffx;
			yp = pg->ExtProf->Element(j, 1) * coeffy;
			(*XProf)->SetElement(j, 0, xp);
			(*YProf)->SetElement(j, 0, yp);
		}
	} else {
		for (j=0; j < pg->IntProf->GetLignes(); j++)
		{
			xp = pg->IntProf->Element(j, 0) * coeffx;
			yp = pg->IntProf->Element(j, 1) * coeffy;
			(*XProf)->SetElement(j, 0, xp);
			(*YProf)->SetElement(j, 0, yp);
		}
	}
}
void makePosProfile(Matrix* Ext, Matrix* Int, double percent, Matrix** ExtRes) {
    //printf ("\n makePosProfile()");
    //IntRes = CloneMat(Int);
    (*ExtRes) = CloneMat(Ext);
    double perc = percent*0.01f;
    double yInt, yExt;
    for (int i = 0; i < Ext->GetLignes(); i++) {
        yExt = Ext->Element(i, 1);
        yInt = InterpLinX(Int, Ext->Element(i,0));
        (*ExtRes)->SetElement(i, 1, yInt + perc*(yExt-yInt));
        //printf ("\n%d -> (%f, %f) => (%f,%f)", i, Ext->Element(i,0), Ext->Element(i,1), (*ExtRes)->Element(i,0), (*ExtRes)->Element(i,1));
    }
    //printf ("\n ...makePosProfile()");
}
