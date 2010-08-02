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

ProfilGeom* getProfilGeomTailDown(ProfilGeom* pg1, double xv) {
	ProfilGeom* pg1h = new ProfilGeom();

	int ne = pg1->ExtProf->GetLignes();
	int ni = pg1->IntProf->GetLignes();
	
	pg1h -> ExtProf = new Matrice (ne, 2);
	pg1h -> IntProf = new Matrice (ni, 2);

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

	return pg1h;
}
