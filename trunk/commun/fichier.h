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


class Matrix;
class Profil;
class WindPatternsProject;
class Ballonement;
class KiteDesign;
class KiteDesignElement;
class Line;

class Form
{

public:
	Form();
	virtual ~Form();

	//parametres principaux
	double CoeffProgGeom;
	double CoeffExp;
	double EpaiRelCent;	//epaisseur relative du profil central 
	int NbCaiss;

	//points de controle = courbes de 4 pts(x,y)
	Matrix *mCtrlNez, *mCtrlFui, *mCtrlA, *mCtrlB, *mCtrlC, *mCtrlD, *mCtrlE;
	Matrix *mCourbNez, *mCourbFui, *mCourbA, *mCourbB, *mCourbC, *mCourbD, *mCourbE;
	Matrix *mCtrlDiedre, *mCtrlMorphing, *mCtrlVrillage, *mCtrlEpaiRel;
	Matrix *mCourbDiedre, *mCourbMorphing, *mCourbVrillage, *mCourbEpaiRel;

	// used only for Serialize
    Matrix *ExtProfCent, *IntProfCent, *ExtProfBout, *IntProfBout;
	double EpaiRelProfCent;
	double EpaiRelProfBout;

	bool courbInput;
	std::string m_strNomProfilCent;
	std::string m_strNomProfilBout;

	Matrix* getExtProf(int nerv, bool realSize);
	Matrix* getIntProf(int nerv, bool realSize);

	//tableau de forme
	Profil** m_pProfils;
	int m_nbProfils;

	Ballonement *ballon;

/*	void SerializeToWpa(QDataStream &ar);
	void SerializeFoil(QDataStream &ar, Matrix* ExtProf, Matrix* IntProf);
	void SerializeWing (QDataStream &ar);
	void WritePolars (QDataStream &ar);*/
	void DeleteProfils();
	void AllocateProfils( int i );
	void Validate();
};

class Form3D
{
public:
	Form3D();
	virtual ~Form3D();
	
	Matrix *XExt, *YExt, *ZExt;
	Matrix *XInt, *YInt, *ZInt;
	Form* forme;
	//Matrix *ExtProfCent, *IntProfCent, *ExtProfBout, *IntProfBout;
};

class FormProjection
{
public:
	FormProjection();
	virtual ~FormProjection();

	Matrix *X, *Y;
};

class Ballonement
{
	//kChord, kMf, wN, dyw
	public:
		Ballonement();
		virtual ~Ballonement();
		Matrix* kChord;
		Matrix* kMf;
		Matrix* wN;
		Matrix* dyw;
		Matrix* powerTail;
		void loadFromFile(const char* fileName);
};

typedef struct InfoForm TInfoForm;

struct InfoForm
{
	double surface, envergure;
	double surfaceProj, envergureProj;
	double allongement, allongementProj;
	double cordeMin, cordeMax;

	double largMin, largMax;

};

void LectureFichierProfil(const char* NomProf, Matrix** extrados, Matrix** intrados);

Ballonement* readBallonementFromFile(char* NomFic);

Form* LectureFichierForm(char* NomFic);

Ballonement* readBallonementFromFile(char* NomFic);

Form* LectureFichierForm2(char* NomFic);

void EcritureFichierForm(char *fileName, Form *f);

void EcritureFichierForm2(char *fileName, Form *f);

bool TrouveMotDansFichierTexte(FILE* fid, char* Mot);

void EcritureFichierDXF(char *fileName, TAxe *axe);

bool EcritureFichierWpa(char *fileName, Form *forme);

void EcritureFichierPolyDXF(char *fileName, TAxe *axe, TAxe *axe2, int rep, TAxe *axeR, int vent, TAxe *axeC, int num, TAxe *axeT);

void EcritureFichierPolyDXFDelta(FILE *fid, TAxe *axe, TAxe *axe2, int rep, TAxe *axeR, int vent, TAxe *axeC, int num, TAxe *axeT, double dx, double dy, double dz, int n );

void EcritureFichierFGen(char *fileName, Form *foil);

void calcInfoForm( Form* F, TInfoForm* info );

void AfficheInfoForm( TInfoForm info );

void GenerateCourbe(WindPatternsProject* gfd, Matrix *Xd1,
					Matrix *Yd1, Matrix *P1,
					int nerv1, double deb1, int faceDeb1, double fin1, int faceFin1,
					Matrix *Xd2, Matrix *Yd2, Matrix *P2,
					int nerv2, double deb2, int faceDeb2, double fin2, int faceFin2, char *text,
			        TAxe **AxePatronP, TAxe **AxePatronDXFP, TAxe **AxePatronTextDXFP, TAxe **AxeMarginDXFP, TAxe **AxeCercleDXFP, TAxe **AxeRepDXFP, int Ventilation,
                    double marge1, double marge2, double margeDeb, double margeFin,bool makeRep, bool debug,
					bool isPince=false, Matrix *Xd01=0, Matrix *Yd01=0,double coeff1=0, Matrix *Xd02=0, Matrix *Yd02=0,double coeff2=0);

WindPatternsProject* LectureWindPatternsProject(char* NomFic);

int* LectureFichierVentHoles(char* NomFic, int* quant, int* central);

int* LectureFichierDiagNervs(char* NomFic, int* quant);

void EcritureManyFichierPolyDXF(char *fileName, int np, int n, TAxe **axe, TAxe **axe2, int rep, TAxe **axeR, int vent, TAxe **axeC, int num, TAxe **axeT, double* W, double* H);

void EcritureManyFichierPolyDXF2(char *fileName, int np, int n, TAxe **axe, TAxe **axe2, int rep, TAxe **axeR, int vent, TAxe **axeC, int num, TAxe **axeT, double* W, double* H, int* numncon);

void EcritureWindPatternsProject(char *fileName, WindPatternsProject *wpp);

KiteDesign* readKiteDesignFromFile(const char* FilePath);