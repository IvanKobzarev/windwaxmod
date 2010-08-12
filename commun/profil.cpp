#include "profil.h"

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
	pg -> ExtProf = new Matrice (ne, 2);
	pg -> IntProf = new Matrice (ni, 2);

	

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
	pg -> ExtProf = new Matrice (ne, 2);
	pg -> IntProf = new Matrice (ni, 2);
	
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
	pgtmp -> ExtProf = new Matrice (iExtTail1 + ne0 - iExtTail0, 2);
	pgtmp -> IntProf = new Matrice (iIntTail1 + ni0 - iIntTail0, 2);

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
	CalculVecteurBissec(pg1->ExtProf->Element(iExtTail1-1, 0), pg1->ExtProf->Element(iExtTail1-1, 1),
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
		CalculVecteurBissec(pg0->ExtProf->Element(i0-1, 0), pg0->ExtProf->Element(i0-1, 1),
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
	CalculVecteurBissec(pg1->IntProf->Element(iIntTail1-1, 0), pg1->IntProf->Element(iIntTail1-1, 1),
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
		CalculVecteurBissec(pg0->IntProf->Element(i0-1, 0), pg0->IntProf->Element(i0-1, 1),
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
