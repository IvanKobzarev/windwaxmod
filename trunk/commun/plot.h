#pragma once
/* header de plot.c */

/*pour activer ou desactiver une fonction*/

#define ON 1

#define OFF 0

/*definition des diff�rents type de courbe*/

#define NON_DEF		0

#define	CTRL_POINT	1



/*def des type de projection*/

#define PROJ_ORTHOGONALE 0 /*orthogonale*/

#define PROJ_PERSPECTIVE 1 /*perspective*/



#define DEG2RAD	(3.1415926f/180.0f)

#define RAD2DEG	(180.0f/3.1415926f)


#include <string>
using namespace std;

class Matrice;
class Forme;
class Courbe;

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

	double CouleurSegments[3];

	double CouleurPoints[3];

	Matrice *pts;		/*matrice des points de la courbe*/

	Courbe *CourbSuiv;	/*pointeur vers la courbe suivante de l'axe courant*/

	bool visible;
	string name;
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

	double CouleurFaces[3];

	double CouleurSegments[3];

	double CouleurPoints[3];

	int InvNormales;	/*inversion normale ON ou OFF*/

	Matrice *x, *y, *z;

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

	double CouleurFond[3];

	double CouleurGrid[3];

	double CouleurRep[3];

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

TAxe* CreerAxe(int fen);

TMesh* CreerMesh(void);



void AjoutCourbe(TAxe *axe, Courbe *courbe);
void AjoutCCourbe(Courbe *courbeBegin, Courbe *courbe);

void showCourbe(Courbe *courbe);

void AjoutMesh(TAxe *axe, TMesh *mesh);

void VisuAxe(TAxe *axe);



void LibererAxe(TAxe *axe);
void LibererMesh(TMesh *mesh);



void LibererCourbesAxe(TAxe *axe);

void LibererMeshsAxe(TAxe *axe);



void AjoutTexte(TAxe *axe, char *texte, double taille, double orientation, double posx, double posy);



void AjoutForme3D( TAxe *Axe3d, 

			 Matrice *XExt, Matrice *YExt, Matrice *ZExt,

			 Matrice *XInt, Matrice *YInt, Matrice *ZInt,

			 int mesh, int symetrie);



void AjoutPtsSuspentage(

						TAxe *Axe3d, Forme *forme, Matrice *IntProfCent,

						Matrice *XExt, Matrice *YExt, Matrice *ZExt,

						Matrice *XInt, Matrice *YInt, Matrice *ZInt,

						int symetrie, Matrice **PosPtsSuspente);



