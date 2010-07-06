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

