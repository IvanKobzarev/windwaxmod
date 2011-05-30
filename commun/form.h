#pragma once

#include <string>

#include "profil.h"
#include "plot.h"
#include "pince.h"
#include "matrice.h"
#include "layout.h"
#include "patternsproject.h"


class Matrix;
class Profil;
class WindPatternsProject;
class Ballonement;
class KiteDesign;
class KiteDesignElement;
class Line;
class Layout;
class Form3D;

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

