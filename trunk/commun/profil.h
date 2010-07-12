#ifndef __PROFIL_H__
#define __PROFIL_H__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "matrice.h"
#include "geom.h"

class ProfilGeom;

class Profil
{
public:
	Profil();
	Profil(
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
		);
	virtual ~Profil();

	double m_fLength;	//longueur
	double m_fWidth;		//epaisseur relative
	double m_fNezX;		//xnez
	double m_fNezY;		//ynez
	double m_fNezZ;		//znez
	double m_fInclin;	//inclinaison / horizontale
	double m_fWash;		//vrillage / horizontale
	double m_fMorph;		//morphing / nervure centrale
	double m_fPosA;		//pos relative ligne A
	double m_fPosB;		//pos relative ligne B
	double m_fPosC;		//pos relative ligne C
	double m_fPosD;		//pos relative ligne D
	double m_fPosE;		//pos relative ligne D
};

class ProfilGeom
{
public:
	Matrice* ExtProf;
	Matrice* IntProf;
	double LongProf;
	void print();
};

ProfilGeom* getBalloneProfilGeom(ProfilGeom* pg0, double kChord, double kMf, double w0, double wN, double dyw);
#endif