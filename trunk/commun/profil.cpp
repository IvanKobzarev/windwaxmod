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

ProfilGeom* getBalloneProfilGeom(ProfilGeom* pg0, double kChord, double kMf, double w0, double wN, double dyw) {
	ProfilGeom* pg = new ProfilGeom();

	int ne = pg0->ExtProf->GetLignes();
	int ni = pg0->IntProf->GetLignes();
	pg -> ExtProf = new Matrice (ne, 2);
	pg -> IntProf = new Matrice (ni, 2);
	
	double l = abs (pg0->ExtProf->Element(ne-1, 0) - pg0->ExtProf->Element(0, 0));

	double l_ = l * kChord;
	double dl = l * (kChord-1.0);


	double x0_ = pg0->ExtProf->Element(0, 0) - dl * kMf; 
	double wabs0 = EpaisseurRelative(pg0->ExtProf, pg0->IntProf);
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
