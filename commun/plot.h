#pragma once

#define DIR_NOSE_UP 1
#define DIR_NOSE_DOWN 0

#define ON 1
#define OFF 0

#define NON_DEF		0
#define	CTRL_POINT	1

#define PROJ_ORTHOGONALE 0 /*orthogonale*/
#define PROJ_PERSPECTIVE 1 /*perspective*/

#define DEG2RAD	(3.1415926f/180.0f)
#define RAD2DEG	(180.0f/3.1415926f)

#define EXT_SIDE 1
#define INT_SIDE 2


#include <string>
#include "profil.h"

//using namespace std;

class Matrix;
class Form;
class Form3D;
class FormProjection;
class Courbe;
class KiteDesign;
class KiteDesignElement;
class Line;

class Courbe
{
public:
	Courbe(const char* name);
	~Courbe();

public:
	int type; /*type de courbe: 0: non defini, 1: point de controle*/
	int points, segments;	/*1 pour visible, 0 sinon*/
	int symX; /*trace courbe sym�trique /axe X: 0:OFF, 1:ON*/
	int symY; /*trace courbe sym�trique /axe Y: 0:OFF, 1:ON*/
	int symZ; /*trace courbe sym�trique /axe Y: 0:OFF, 1:ON si axe3d*/
	int alt; /*alterne segment et vide au trac�*/
	double colorSegments[4];
	double colorPoints[4];
	Matrix *pts;		/*matrice des points de la courbe*/
	Courbe *CourbSuiv;	/*pointeur vers la courbe suivante de l'axe courant*/
	bool visible;
	std::string name;
};

/**************************************/
/*structure de mesh (maillage) chain�e*/
/**************************************/

typedef struct Mesh TMesh;
struct Mesh
{
	int points, segments, faces;	/*1 pour visible, 0 sinon*/
	int symX; /*trace courbe sym�trique /axe X: 0:OFF, 1:ON*/
	int symY; /*trace courbe sym�trique /axe Y: 0:OFF, 1:ON*/
	int symZ; /*trace courbe sym�trique /axe Y: 0:OFF, 1:ON si axe3d*/
	double colorFaces[4];
	double colorSegments[4];
	double colorPoints[4];
	int InvNormales;	/*inversion normale ON ou OFF*/
	int side;
	Matrix *x, *y, *z;
	TMesh *MeshSuiv;	/*pointeur vers la courbe suivante de l'axe courant*/
};

/**************************/
/*structure d'axe 2d ou 3d*/
/**************************/

struct Axe
{
	int fenetre;	/*numero fenetre contenant l'axe*/
	int axe3d;		/*0: axe 2d, 1:axe 3d*/
	double xmin,xmax,ymin,ymax,zmin,zmax;	/*valeur min et max des axes*/
	double colorFond[4];
	double colorGrid[4];
	double colorRep[4];
	int XAuto, YAuto, ZAuto;	/*auto-ajustement des axes*/
	int Norme;					/*axes norm�s*/
	int XGrid, YGrid, ZGrid;	/*0: pas de graduation, 1: graduation*/
	double xcam, ycam, zcam, twist, incidence, azimuth; /*position camera dans la vue 3d*/
	int eclairage;	/*parametres d'�clairage OpenGL: 0:d�sactiv�, 1:actif*/
	int proj;		/*type de projection*/

	Courbe *Courb;	/*liste des courbes de l'axe*/
	TMesh *Mesh;	/*liste des mesh de l'axe*/
};

typedef struct Axe  TAxe;

/************/
/*procedures*/
/************/
TAxe* createAxe(int fen);
TMesh* createMesh(void);

void addCourbe(TAxe *axe, Courbe *courbe);
void addCCourbe(Courbe *courbeBegin, Courbe *courbe);
void showCourbe(Courbe *courbe);
void addMesh(TAxe *axe, TMesh *mesh);
void ViewAxe(TAxe *axe);
void clearAxe(TAxe *axe);
void clearMesh(TMesh *mesh);
void clearCourbesAxe(TAxe *axe);
void clearMeshsAxe(TAxe *axe);
void addTexte(TAxe *axe, char *texte, double taille, double orientation, double posx, double posy);

void addForm3D( TAxe *Axe3d, 
			 Matrix *XExt, Matrix *YExt, Matrix *ZExt,
			 Matrix *XInt, Matrix *YInt, Matrix *ZInt,
			 int mesh, int symetrie);

void addPtsSuspentage(
						TAxe *Axe3d, Form *forme, Matrix *IntProfCent,
						Matrix *XExt, Matrix *YExt, Matrix *ZExt,
						Matrix *XInt, Matrix *YInt, Matrix *ZInt,
						int symetrie, Matrix **PosPtsSuspente);

