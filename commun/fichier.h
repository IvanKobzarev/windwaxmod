#pragma once

#define NO_VERSION 0.1f
#define DEG2RAD	(3.1415926f/180.0f)
#define RAD2DEG	(180.0f/3.1415926f)

#include <string>

#include "plot.h"
#include "pince.h"
#include "matrice.h"
#include "rasklad.h"
#include "patternsproject.h"

/* using namespace std; */
class Matrice;
class Profil;
class WindPatternsProject;
class Ballonement;

class Forme
{

public:
	Forme();
	virtual ~Forme();
	

	//parametres principaux
	double CoeffProgGeom;
	double CoeffExp;
	double EpaiRelCent;	//epaisseur relative du profil central 
	int NbCaiss;

	//points de controle = courbes de 4 pts(x,y)
	Matrice *mCtrlNez, *mCtrlFui, *mCtrlA, *mCtrlB, *mCtrlC, *mCtrlD, *mCtrlE;
	Matrice *mCourbNez, *mCourbFui, *mCourbA, *mCourbB, *mCourbC, *mCourbD, *mCourbE;
	Matrice *mCtrlDiedre, *mCtrlMorphing, *mCtrlVrillage, *mCtrlEpaiRel;
	Matrice *mCourbDiedre, *mCourbMorphing, *mCourbVrillage, *mCourbEpaiRel;


	// used only for Serialize
    Matrice *ExtProfCent, *IntProfCent, *ExtProfBout, *IntProfBout;

	bool courbInput;

	//nom des profils central/extremite

	string m_strNomProfilCent;
	string m_strNomProfilBout;

	//tableau de forme

	Profil** m_pProfils;
	int m_nbProfils;

	Ballonement *ballon;

/*	void SerializeToWpa(QDataStream &ar);
	void SerializeFoil(QDataStream &ar, Matrice* ExtProf, Matrice* IntProf);
	void SerializeWing (QDataStream &ar);
	void WritePolars (QDataStream &ar);*/
	void DeleteProfils();
	void AllocateProfils( int i );
	void Validate();
};

class Forme3D
{
public:
	Forme3D();
	virtual ~Forme3D();
	
	Matrice *XExt, *YExt, *ZExt;
	Matrice *XInt, *YInt, *ZInt;
};

class FormeProjection
{
public:
	FormeProjection();
	virtual ~FormeProjection();

	Matrice *X, *Y;
};

class Ballonement
{
	//kChord, kMf, wN, dyw
	public:
		Ballonement();
		virtual ~Ballonement();
		Matrice* kChord;
		Matrice* kMf;
		Matrice* wN;
		Matrice* dyw;
		Matrice* powerTail;
		void loadFromFile(const char* fileName);
};

typedef struct InfoForme TInfoForme;

struct InfoForme
{
	double surface, envergure;
	double surfaceProj, envergureProj;
	double allongement, allongementProj;
	double cordeMin, cordeMax;

	double largMin, largMax;

};

void LectureFichierProfil(const char* NomProf, Matrice** extrados, Matrice** intrados);

Ballonement* readBallonementFromFile(char* NomFic);

Forme* LectureFichierForme(char* NomFic);

Ballonement* readBallonementFromFile(char* NomFic);

Forme* LectureFichierForme2(char* NomFic);

void EcritureFichierForme(char *NomFichier, Forme *f);

void EcritureFichierForme2(char *NomFichier, Forme *f);

bool TrouveMotDansFichierTexte(FILE* fid, char* Mot);

void EcritureFichierDXF(char *NomFichier, TAxe *axe);

bool EcritureFichierWpa(char *NomFichier, Forme *forme);

void EcritureFichierPolyDXF(char *NomFichier, TAxe *axe, TAxe *axe2, int rep, TAxe *axeR, int vent, TAxe *axeC, int num, TAxe *axeT);

void EcritureFichierPolyDXFDelta(FILE *fid, TAxe *axe, TAxe *axe2, int rep, TAxe *axeR, int vent, TAxe *axeC, int num, TAxe *axeT, double dx, double dy, double dz, int n );

void EcritureFichierFGen(char *NomFichier, Forme *foil);

void CalculInfoForme( Forme* F, TInfoForme* info );

void AfficheInfoForme( TInfoForme info );

void GenerateCourbe(WindPatternsProject* gfd, Matrice *Xd1,
					Matrice *Yd1, Matrice *P1,
					int nerv1, double deb1, int faceDeb1, double fin1, int faceFin1,
					Matrice *Xd2, Matrice *Yd2, Matrice *P2,
					int nerv2, double deb2, int faceDeb2, double fin2, int faceFin2, char *text,
			        TAxe **AxePatronP, TAxe **AxePatronDXFP, TAxe **AxePatronTextDXFP, TAxe **AxeMarginDXFP, TAxe **AxeCercleDXFP, TAxe **AxeRepDXFP, int Ventilation,
                    double marge1, double marge2, double margeDeb, double margeFin,bool makeRep, bool debug,
					bool isPince=false, Matrice *Xd01=0, Matrice *Yd01=0,double coeff1=0, Matrice *Xd02=0, Matrice *Yd02=0,double coeff2=0);

WindPatternsProject* LectureWindPatternsProject(char* NomFic);

int* LectureFichierVentHoles(char* NomFic, int* quant, int* central);

int* LectureFichierDiagNervs(char* NomFic, int* quant);

void EcritureManyFichierPolyDXF(char *NomFichier, int np, int n, TAxe **axe, TAxe **axe2, int rep, TAxe **axeR, int vent, TAxe **axeC, int num, TAxe **axeT, double* W, double* H);

void EcritureManyFichierPolyDXF2(char *NomFichier, int np, int n, TAxe **axe, TAxe **axe2, int rep, TAxe **axeR, int vent, TAxe **axeC, int num, TAxe **axeT, double* W, double* H, int* numncon);

void EcritureWindPatternsProject(char *NomFichier, WindPatternsProject *wpp);
