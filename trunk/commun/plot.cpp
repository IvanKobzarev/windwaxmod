#pragma warning(disable:4514)
#pragma warning(disable:4505)

#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>

#include "matrice.h"
#include "geom.h"
#include "plot.h"
#include "fichier.h"
#include "profil.h"
#include "design.h"
#include "GL/glut.h"

#ifndef DEBUG
#define DEBUG false
#endif

#define EXT_SIDE 1
#define INT_SIDE 2
/********************************/
/* constantes pour tracer texte */
/********************************/
	const char TabXYch[38][18]={
		{'0',	4,	15,	25,	34,	32,	21,	11,	2,	4,	32,	0,	0,	0,	0,	0,	0,	0},
		{'1',	12,	3,	33,	32,	34,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'2',	11,	6,	2,	4,	10,	15,	31,	35,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'3',	6,	2,	4,	10,	15,	19,	18,	19,	25,	30,	34,	32,	26,	0,	0,	0,	0},
		{'4',	3,	21,	25,	24,	14,	34,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'5',	5,	1,	16,	19,	25,	30,	34,	31,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'6',	10,	4,	2,	6,	26,	32,	34,	30,	25,	19,	16,	0,	0,	0,	0,	0,	0},
		{'7',	1,	5,	31,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'8',	2,	4,	10,	15,	19,	17,	21,	26,	32,	34,	30,	25,	19,	17,	11,	6,	2},
		{'9',	20, 17,	11,	6,	2,	4,	10,	30,	34,	32,	26,	0,	0,	0,	0,	0,	0},
		{'A',	31, 16,	7,	3,	9,	20,	16,	20,	35,	0,	0,	0,	0,	0,	0,	0,	0},
		{'B',	1,	4,	10,	15,	19,	16,	19,	25,	30,	34,	31,	1,	0,	0,	0,	0,	0},
		{'C',	10, 4,	2,	6,	26,	32,	34,	30,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'D',	1,	4,	10,	30,	34,	31,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'E',	5,	1,	16,	19,	16,	31,	35,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'F',	5,	1,	16,	19,	16,	31,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'G',	10, 4,	2,	6,	26,	32,	34,	30,	20,	18,	0,	0,	0,	0,	0,	0,	0},
		{'H',	1,	31,	16,	20,	5,	35,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'I',	2,	4,	3,	33,	32,	34,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'J',	3,	5,	4,	29,	33,	32,	26,	21,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'K',	5,	16,	1,	31,	16,	35,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'L',	1,	31,	35,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'M',	31, 1,	18,	5,	35,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'N',	31, 1,	35,	5,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'O',	11, 2,	4,	15,	25,	34,	32,	21,	11,	0,	0,	0,	0,	0,	0,	0,	0},
		{'P',	31, 1,	4,	10,	20,	24,	21,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'Q',	35, 23,	29,	33,	32,	26,	6,	2,	4,	10,	25,	29,	0,	0,	0,	0,	0},
		{'R',	31, 1,	4,	10,	20,	24,	21,	24,	35,	0,	0,	0,	0,	0,	0,	0,	0},
		{'S',	10, 4,	2,	6,	11,	17,	19,	25,	30,	34,	32,	26,	0,	0,	0,	0,	0},
		{'T',	1,	5,	3,	33,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'U',	1,	26,	32,	34,	30,	5,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'V',	1,	33,	5,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'W',	1,	32,	13,	34,	5,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'X',	1,	35,	18,	31,	5,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'Y',	1,	18,	33,	18,	5,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'Z',	1,	5,	31,	35,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'.',	22, 24,	34,	32,	22,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
		{'-',	16, 20,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0},
	};

	const double TabXYpoint[40][2]={
		{0.0f, 7.0f},{ 1.0f, 7.0f},{ 2.0f, 7.0f},{ 3.0f, 7.0f},{ 4.0f, 7.0f},
		{0.0f, 6.0f},{ 1.0f, 6.0f},{ 2.0f, 6.0f},{ 3.0f, 6.0f},{ 4.0f, 6.0f},
		{0.0f, 5.0f},{ 1.0f, 5.0f},{ 2.0f, 5.0f},{ 3.0f, 5.0f},{ 4.0f, 5.0f},
		{0.0f, 4.0f},{ 1.0f, 4.0f},{ 2.0f, 4.0f},{ 3.0f, 4.0f},{ 4.0f, 4.0f},
		{0.0f, 3.0f},{ 1.0f, 3.0f},{ 2.0f, 3.0f},{ 3.0f, 3.0f},{ 4.0f, 3.0f},
		{0.0f, 2.0f},{ 1.0f, 2.0f},{ 2.0f, 2.0f},{ 3.0f, 2.0f},{ 4.0f, 2.0f},
		{0.0f, 1.0f},{ 1.0f, 1.0f},{ 2.0f, 1.0f},{ 3.0f, 1.0f},{ 4.0f, 1.0f},
		{0.0f, 0.0f},{ 1.0f, 0.0f},{ 2.0f, 0.0f},{ 3.0f, 0.0f},{ 4.0f, 0.0f}
	};

TAxe* CreerAxe(int fen)
{
	TAxe *axe;
	/*allocation memoire*/
	axe = (TAxe*)malloc(sizeof(TAxe));
	/*init parametres par defaut*/
	axe->axe3d = OFF;
	axe->XAuto = ON; axe->YAuto = ON; axe->ZAuto = ON;
	axe->Norme = ON;
	axe->Courb = NULL;
	axe->Mesh = NULL;
	axe->fenetre = fen;
	axe->XGrid = ON; axe->YGrid = ON; axe->ZGrid = ON;
	axe->xmin = -1.0f; axe->xmax = +1.0f;
	axe->ymin = -1.0f; axe->ymax = +1.0f;
	axe->zmin = -1.0f; axe->zmax = +1.0f;
	/*couleur de fond: gris*/
	axe->CouleurFond[0] = 0.5f;	axe->CouleurFond[1] = 0.5f;	axe->CouleurFond[2] = 0.5f;
	/*couleur grille: blanc*/
	axe->CouleurGrid[0] = 1.0f;	axe->CouleurGrid[1] = 1.0f; axe->CouleurGrid[2] = 1.0f;
	/*couleur repere: rouge*/
	axe->CouleurRep[0] = 1.0f; axe->CouleurRep[1] = 0.0f; axe->CouleurRep[2] = 0.0f;
	/*position camera dans vue 3d*/
	axe->xcam = 0.0; axe->ycam = 0.0;	axe->zcam = 10.0;
	axe->twist = 0.0; axe->incidence = 0.0;	axe->azimuth = 0.0;
	/*eclairage*/
	axe->eclairage = OFF;
	/*type de projection 3d*/
	axe->proj = PROJ_ORTHOGONALE;
	return axe;
}

Courbe::Courbe(const char* myname)
{
	name = myname;

	/*init parametres par defaut*/
	type = NON_DEF;
	CourbSuiv = NULL;
	pts = NULL;
	points = ON;
	segments = ON;
	symX = OFF;
	symY = OFF;
	symZ = OFF;
	alt = OFF;

	/*jaune*/
	CouleurPoints[0] = 1.0f;
	CouleurPoints[1] = 1.0f;
	CouleurPoints[2] = 0.0f;

	/*blanc*/
	CouleurSegments[0] = 1.0f;
	CouleurSegments[1] = 1.0f;
	CouleurSegments[2] = 1.0f;
	visible = true;
}

Courbe::~Courbe()
{
	if ( pts != NULL )
		delete pts;
	pts = NULL;
}

TMesh* CreerMesh(void)
{
	TMesh *mesh;
	/*allocation memoire*/
	mesh = (TMesh*)malloc(sizeof(TMesh));
	/*init parametres par defaut*/
	mesh->MeshSuiv = NULL;
	mesh->x = NULL;
	mesh->y = NULL;
	mesh->z = NULL;
	mesh->points = OFF;
	mesh->segments = ON;
	mesh->faces = OFF;
	mesh->symX = OFF;
	mesh->symY = OFF;
	mesh->symZ = OFF;
	mesh->InvNormales = OFF;
	/*jaune*/
	mesh->CouleurPoints[0] = 1.0f;
	mesh->CouleurPoints[1] = 1.0f;
	mesh->CouleurPoints[2] = 0.0f;
	/*blanc*/
	mesh->CouleurSegments[0] = 1.0f;
	mesh->CouleurSegments[1] = 1.0f;
	mesh->CouleurSegments[2] = 1.0f;
	/*jaune*/
	mesh->CouleurFaces[0] = 1.0f;
	mesh->CouleurFaces[1] = 1.0f;
	mesh->CouleurFaces[2] = 0.0f;
	return mesh;
}

void VisuAxe(TAxe *axe)
{
	int i,l,c, w, h, nbAff;
	double xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,ratio;
	Courbe* CourbCour;
	TMesh* MeshCour;
	double mx,my,mz; /*coeff -1 ou 1 pour affichage sym�trique*/
	/*double cap, roulis, tangage;*/
	Matrice *rotTwist, *rotIncidence, *rotAzimuth, *rotFinal, *rotInt; /*matrice de rotation*/
	Matrice *uni, *uniRot;	/*vecteur unitaire*/
	double sina, cosa; /*pour pre-calcul sinus/cosinus*/
	double xaxe=0.0, yaxe=0.0, zaxe=0.0;	/*position graduation dans vue 3d*/
	Matrice *x,*y,*z;	/*pour simplifier ecriture*/
	double xv1,yv1,zv1,xv2,yv2,zv2;
//	double light_position[]={0.0, 0.0, 1.0, 0.0};
	/*selection de la fenetre de cet axe*/
	glutSetWindow(axe->fenetre);
	/***************/
	/*axe auto en X*/
	/***************/
	if((axe->XAuto == ON)&&((axe->Courb != NULL)||(axe->Mesh != NULL)))
	{
		xmin=10000.0f; xmax=-10000.0f;

		/*recherche des min/max des courbes*/
		CourbCour=axe->Courb;

		while (CourbCour != NULL)
		{	
			if (CourbCour->pts != NULL && CourbCour->visible)
			{
				for(l=0; l<CourbCour->pts->GetLignes(); l++)
				{
					/*mise a jour xmin,...*/
					if (xmin>CourbCour->pts->Element(l,0)) xmin=CourbCour->pts->Element(l,0);
					if (xmax<CourbCour->pts->Element(l,0)) xmax=CourbCour->pts->Element(l,0);

					/*test symetrie X*/
					if ((CourbCour->symX==ON)&&(xmin>-xmax)) xmin=-xmax;
					if ((CourbCour->symX==ON)&&(xmin<-xmax)) xmax=-xmin;
				}
			}
			/*passage � la courbe suivante*/
			CourbCour=CourbCour->CourbSuiv;	
		}

		/*recherche des min/max des mesh*/

		MeshCour=axe->Mesh;

		while (MeshCour != NULL)
		{	
			if (MeshCour->x != NULL)
			{
				for(l=0; l<MeshCour->x->GetLignes(); l++)
				{
					for(c=0; c<MeshCour->x->GetColonnes(); c++)
					{
						/*mise a jour xmin,...*/
						if (xmin>MeshCour->x->Element(l,c)) xmin=MeshCour->x->Element(l,c);
						if (xmax<MeshCour->x->Element(l,c)) xmax=MeshCour->x->Element(l,c);

						/*test symetrie X*/
						if ((MeshCour->symX==ON)&&(xmin>-xmax)) xmin=-xmax;
						if ((MeshCour->symX==ON)&&(xmin<-xmax)) xmax=-xmin;
					}
				}
			}
			/*passage au Mesh suivant*/
			MeshCour=MeshCour->MeshSuiv;
		}

		/*bordure supplementaire de 5%*/
		xmin=xmin-(xmax-xmin)*0.05f; xmax=xmax+(xmax-xmin)*0.05f;

		/*mise a jour axe*/
		axe->xmin = xmin; axe->xmax = xmax;
	}

	/***************/
	/*axe auto en Y*/
	/***************/

	if ((axe->YAuto == ON)&&((axe->Courb != NULL)||(axe->Mesh != NULL)))
	{
		ymin=10000.0f; ymax=-10000.0f;
		/*recherche des min/max des courbes*/
		CourbCour=axe->Courb;
		while (CourbCour != NULL)
		{	
			if (CourbCour->pts != NULL && CourbCour->visible)
			{
				for(l=0; l<CourbCour->pts->GetLignes(); l++)
				{
					/*mise a jour ymin,...*/
					if (ymin>CourbCour->pts->Element(l,1)) ymin=CourbCour->pts->Element(l,1);
					if (ymax<CourbCour->pts->Element(l,1)) ymax=CourbCour->pts->Element(l,1);
					/*test symetrie Y*/
					if ((CourbCour->symY==ON)&&(ymin>-ymax)) ymin=-ymax;
					if ((CourbCour->symY==ON)&&(ymin<-ymax)) ymax=-ymin;
				}
			}
			/*passage � la courbe suivante*/
			CourbCour=CourbCour->CourbSuiv;
		}

		/*recherche des min/max des mesh*/
		MeshCour=axe->Mesh;
		while (MeshCour != NULL)
		{	
			if (MeshCour->y != NULL)
			{
				for(l=0; l<MeshCour->y->GetLignes(); l++)
				{
					for(c=0; c<MeshCour->y->GetColonnes(); c++)
					{
						/*mise a jour ymin,...*/
						if (ymin>MeshCour->y->Element(l,c)) ymin=MeshCour->y->Element(l,c);
						if (ymax<MeshCour->y->Element(l,c)) ymax=MeshCour->y->Element(l,c);
						/*test symetrie Y*/
						if ((MeshCour->symY==ON)&&(ymin>-ymax)) ymin=-ymax;
						if ((MeshCour->symY==ON)&&(ymin<-ymax)) ymax=-ymin;
					}
				}
			}
			/*passage au Mesh suivant*/
			MeshCour=MeshCour->MeshSuiv;
		}

		/*bordure supplementaire de 5%*/
		ymin=ymin-(ymax-ymin)*0.05f; ymax=ymax+(ymax-ymin)*0.05f;

		/*mise a jour axe*/
		axe->ymin = ymin; axe->ymax = ymax;
	}


	/***************/
	/*axe auto en Z*/
	/***************/
	if ((axe->axe3d == ON)&&((axe->Courb != NULL)||(axe->Mesh != NULL)))
	{
		zmin=10000.0f; zmax=-10000.0f;
		/*recherche des min/max des courbes*/
		CourbCour=axe->Courb;
		while (CourbCour != NULL)
		{	
			if (CourbCour->pts != NULL && CourbCour->visible)
				for(l=0; l<CourbCour->pts->GetLignes(); l++)
				{
					/*mise a jour zmin,...*/
					if (zmin>CourbCour->pts->Element(l,2)) zmin=CourbCour->pts->Element(l,2);
					if (zmax<CourbCour->pts->Element(l,2)) zmax=CourbCour->pts->Element(l,2);
					/*test symetrie Z*/
					if ((CourbCour->symZ==ON)&&(zmin>-zmax)) zmin=-zmax;
					if ((CourbCour->symZ==ON)&&(zmin<-zmax)) zmax=-zmin;
				}
				/*passage � la courbe suivante*/
				CourbCour=CourbCour->CourbSuiv;
		}
		/*recherche des min/max des mesh*/
		MeshCour=axe->Mesh;
		while (MeshCour != NULL)
		{	
			if (MeshCour->z != NULL)
				for(l=0; l<MeshCour->z->GetLignes(); l++)
				for(c=0; c<MeshCour->z->GetColonnes(); c++)
				{
					/*mise a jour zmin,...*/
					if (zmin>MeshCour->z->Element(l,c)) zmin=MeshCour->z->Element(l,c);
					if (zmax<MeshCour->z->Element(l,c)) zmax=MeshCour->z->Element(l,c);
					/*test symetrie Z*/
					if ((MeshCour->symZ==ON)&&(zmin>-zmax)) zmin=-zmax;
					if ((MeshCour->symZ==ON)&&(zmin<-zmax)) zmax=-zmin;
				}
				/*passage au Mesh suivant*/
				MeshCour=MeshCour->MeshSuiv;
		}
		/*bordure supplementaire de 5%*/
		zmin=zmin-(zmax-zmin)*0.05f; zmax=zmax+(zmax-zmin)*0.05f;
		/*mise a jour axe*/
		axe->zmin = zmin; axe->zmax = zmax;
	}

	/*****************/
	/*axe x/y Norm�s */
	/*****************/
	if (axe->Norme == ON)
	{
		w = glutGet(GLUT_WINDOW_WIDTH);
		h = glutGet(GLUT_WINDOW_HEIGHT);
		/*test fenetre en icone ...*/
		if(h!=0)
		{
			dx = axe->xmax - axe->xmin;
			dy = axe->ymax - axe->ymin;
			ratio = (dy/dx)/((double)h/(double)w);
			if (ratio>1.0f)
			{
				axe->xmin -= dx*(ratio-1.0f)/2.0f;
				axe->xmax += dx*(ratio-1.0f)/2.0f;
			}
			else
			{
				axe->ymin -= dy*(1.0f/ratio-1.0f)/2.0f;
				axe->ymax += dy*(1.0f/ratio-1.0f)/2.0f;
			}
			/*
			printf("\nxmin=%1.3f, xmax=%1.3f, ymin=%1.3f, ymax=%1.3f, r=%1.3f, h=%d, w=%d    ",
				axe->xmin,axe->xmax,axe->ymin,axe->ymax,ratio,h,w); 
			*/
		}
	}

	/**********************************************/
	/*couleur d'effacement et parametre eclairage */
	/**********************************************/
	glClearColor(axe->CouleurFond[0], axe->CouleurFond[1], axe->CouleurFond[2], 0.0);
	if(axe->axe3d == ON)
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		if(axe->eclairage == OFF)
			glDisable(GL_LIGHTING);
		else
			glEnable(GL_LIGHTING);
	}
	else
	{
		glClear(GL_COLOR_BUFFER_BIT);
	}
	
	/**********************/
	/*modele de projection*/
	/**********************/
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (axe->axe3d == OFF)
	{
		glOrtho(axe->xmin, axe->xmax, axe->ymin, axe->ymax, -5.0, 5.0);
	}
	else
	{
		if(axe->proj == PROJ_ORTHOGONALE)
			glOrtho(axe->xmin, axe->xmax, axe->ymin, axe->ymax, 0.1, 100.0);
		else /*PROJ_PERSPECTIVE*/
			gluPerspective(50.0, 1.5, 0.1, 100.0);
	}

	/*************************/
	/*modele de visualisation*/
	/*************************/
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	if (axe->axe3d == ON)
	{
		/*cap = 30; roulis = 30; tangage =30;*/
		glTranslatef(axe->xcam, axe->ycam, -axe->zcam);
		glRotated(axe->twist,0,0,1);
		glRotated(axe->incidence,1,0,0);
		glRotated(axe->azimuth,0,1,0);
	}
	
	/*******************************************/
	/* si axe 3d trac� d'une boite anglobante  */
	/* et calcul de la position des graduations*/
	/*******************************************/
	if((axe->axe3d == ON)&&(axe->ZGrid))
	{
		/*trac� boite englobante*/
		glLineStipple (1, 0xFFFF);
		glColor3dv(axe->CouleurRep);
		glBegin(GL_LINE_STRIP);
		glVertex3f(axe->xmin,axe->ymin,axe->zmin);
		glVertex3f(axe->xmax,axe->ymin,axe->zmin);
		glVertex3f(axe->xmax,axe->ymax,axe->zmin);
		glVertex3f(axe->xmin,axe->ymax,axe->zmin);
		glVertex3f(axe->xmin,axe->ymin,axe->zmin);
		glVertex3f(axe->xmin,axe->ymin,axe->zmax);
		glVertex3f(axe->xmax,axe->ymin,axe->zmax);
		glVertex3f(axe->xmax,axe->ymin,axe->zmin);
		glVertex3f(axe->xmax,axe->ymin,axe->zmax);
		glVertex3f(axe->xmax,axe->ymax,axe->zmax);
		glVertex3f(axe->xmax,axe->ymax,axe->zmin);
		glVertex3f(axe->xmax,axe->ymax,axe->zmax);
		glVertex3f(axe->xmin,axe->ymax,axe->zmax);
		glVertex3f(axe->xmin,axe->ymax,axe->zmin);
		glVertex3f(axe->xmin,axe->ymax,axe->zmax);
		glVertex3f(axe->xmin,axe->ymin,axe->zmax);
		glEnd();

		/*calcul matrice des 3 rotations*/
		/*pour d�terminer quelle faces se trouvent en arri�re plan*/
		/*afin d'y tracer les graduations ... */
		rotTwist=Zeros(3,3); rotIncidence=Zeros(3,3); rotAzimuth=Zeros(3,3);
		sina=sin(axe->twist*DEG2RAD); cosa=cos(axe->twist*DEG2RAD);
		rotTwist->SetElement(0,0,cosa);
		rotTwist->SetElement(0,1,-sina);

		rotTwist->SetElement(1,0,sina); 
		rotTwist->SetElement(1,1,cosa);

		rotTwist->SetElement(2,2,1.0f);

		sina=sin(axe->incidence*DEG2RAD); 
		cosa=cos(axe->incidence*DEG2RAD);

		rotIncidence->SetElement(0,0,1.0f);

		rotIncidence->SetElement(1,1,cosa); 
		rotIncidence->SetElement(1,2,-sina);

		rotIncidence->SetElement(2,1,sina); 
		rotIncidence->SetElement(2,2,cosa);

		sina=sin(axe->azimuth*DEG2RAD);
		cosa=cos(axe->azimuth*DEG2RAD);

		rotAzimuth->SetElement(0,0,cosa); 
		rotAzimuth->SetElement(0,2,sina);

		rotAzimuth->SetElement(1,1,1.0f);

		rotAzimuth->SetElement(2,0,-sina); 
		rotAzimuth->SetElement(2,2,cosa);

		/*calcul matrice resultante*/
		rotInt=MultMat(rotTwist,rotIncidence);
		rotFinal=MultMat(rotInt,rotAzimuth);
		/*visu pour debug*/
		/*VoirMat(rotTwist);
		VoirMat(rotIncidence);
		VoirMat(rotAzimuth);
		VoirMat(rotInt);
		VoirMat(rotFinal);*/
		delete(rotInt);
		/*calcul position x des graduations*/
		uni=Zeros(3,1); 
		uni->SetElement(0,0,1);	
		uniRot=MultMat(rotFinal, uni);

		xaxe = (uniRot->Element(2,0)>0) ? axe->xmin : axe->xmax ;

		delete(uniRot);

		/*calcul position y des graduations*/

		uni=Zeros(3,1); 
		uni->SetElement(1,0,1);	
		uniRot=MultMat(rotFinal, uni);

		yaxe = (uniRot->Element(2,0)>0) ? axe->ymin : axe->ymax;

		delete(uniRot);

		/*calcul position z des graduations*/

		uni=Zeros(3,1); 
		uni->SetElement(2,0,1);	
		uniRot=MultMat(rotFinal, uni);

		zaxe = (uniRot->Element(2,0)>0) ? axe->zmin : axe->zmax;

		delete(uniRot);

	}


	/*******************/
	/*trac� grille en X*/
	/*******************/
	if(axe->XGrid == ON)
	{
		/*dx = axe->xmax - axe->xmin;*/
		for(i=ceil(axe->xmin); i<ceil(axe->xmax); i++)
		{
			if(i!=0)
			{
				glLineStipple (1, 0x0F0F);
				glColor3dv(axe->CouleurGrid);
			}
			else
			{
				glLineStipple (1, 0xFFFF);
				glColor3dv(axe->CouleurRep);
			}
			if(axe->axe3d == OFF)
			{
				glBegin(GL_LINES);
				glVertex2f((double)i,axe->ymin);
				glVertex2f((double)i,axe->ymax);
				glEnd();
			}
			else
			{
				glBegin(GL_LINES);
				glVertex3f((double)i,axe->ymin, zaxe);
				glVertex3f((double)i,axe->ymax, zaxe);
				glEnd();
				glBegin(GL_LINES);
				glVertex3f((double)i, yaxe, axe->zmin);
				glVertex3f((double)i, yaxe, axe->zmax);
				glEnd();
			}
		}
		/*retour en ligne continu par defaut*/
		glLineStipple (1, 0xFFFF);
	}

	/*******************/
	/*trac� grille en Y*/
	/*******************/
	if(axe->YGrid == ON)
	{
		/*dy = axe->ymax - axe->ymin;*/
		for(i=ceil(axe->ymin); i<ceil(axe->ymax); i++)
		{
			if(i!=0)
			{
				glLineStipple (1, 0x0F0F);
				glColor3dv(axe->CouleurGrid);
			}
			else
			{
				glLineStipple (1, 0xFFFF);
				glColor3dv(axe->CouleurRep);
			}
			if(axe->axe3d == OFF)
			{
				glBegin(GL_LINES);
				glVertex2f(axe->xmin,(double)i);
				glVertex2f(axe->xmax,(double)i);
				glEnd();
			}
			else
			{
				glBegin(GL_LINES);
				glVertex3f(axe->xmin,(double)i, zaxe);
				glVertex3f(axe->xmax,(double)i, zaxe);
				glEnd();
				glBegin(GL_LINES);
				glVertex3f(xaxe, (double)i, axe->zmin);
				glVertex3f(xaxe, (double)i, axe->zmax);
				glEnd();
			}
		}
		/*retour en ligne continu par defaut*/
		glLineStipple (1, 0xFFFF);
	}

	/*******************/
	/*trac� grille en Z*/
	/*******************/
	if((axe->ZGrid == ON)&&(axe->axe3d == ON))
	{

		/*dz = axe->zmax - axe->zmin;*/
		for(i=ceil(axe->zmin); i<ceil(axe->zmax); i++)
		{
			if(i!=0)
			{
				glLineStipple (1, 0x0F0F);
				glColor3dv(axe->CouleurGrid);
			}
			else
			{
				glLineStipple (1, 0xFFFF);
				glColor3dv(axe->CouleurRep);
			}
			glBegin(GL_LINES);
			glVertex3f(axe->xmin,yaxe,(double)i);
			glVertex3f(axe->xmax,yaxe,(double)i);
			glEnd();
			glBegin(GL_LINES);
			glVertex3f(xaxe,axe->ymin,(double)i);
			glVertex3f(xaxe,axe->ymax,(double)i);
			glEnd();
		}
		/*retour en ligne continu par defaut*/
		glLineStipple (1, 0xFFFF);
	}

	/*******************/
	/*trac� des courbes*/
	/*******************/
	CourbCour=axe->Courb;
	while (CourbCour!=NULL)
	{
		if (CourbCour->pts != NULL && CourbCour->visible)
		{
			/*determine si affichage sym�trique*/
			nbAff=1;
			if ((CourbCour->symX==ON)||(CourbCour->symX==ON)||(CourbCour->symX==ON)) nbAff=2;
			/*boucle affichage normal et sym�trique*/
			for(i=0; i<nbAff; i++)
			{
				/*mise a jour coeff mx, my, mz a 1.0 par defaut*/
				mx =1.0f; my=1.0f, mz=1.0f;
				if(i==1)
				{
					if (CourbCour->symX==ON) mx=-1.0f;
					if (CourbCour->symY==ON) my=-1.0f;
					if (CourbCour->symZ==ON) mz=-1.0f;
				}
				
				/*affichage des points*/
				if (CourbCour->points == ON)
				{
					glPointSize(5.0);
					glColor3dv(CourbCour->CouleurPoints);
					glBegin(GL_POINTS);
					for(l=0; l<CourbCour->pts->GetLignes(); l++)
						if(axe->axe3d == OFF)
							glVertex2f(mx*CourbCour->pts->Element(l,0), my*CourbCour->pts->Element(l,1));
						else
							glVertex3f(mx*CourbCour->pts->Element(l,0), my*CourbCour->pts->Element(l,1),
							mz*CourbCour->pts->Element(l,2));
						glEnd();
				}
				
				/*affichage des segments*/
				if (CourbCour->segments == ON)
				{
					glColor3dv(CourbCour->CouleurSegments);
					if(CourbCour->alt == OFF)
					{
						glBegin(GL_LINE_STRIP);
						for(l=0; l<CourbCour->pts->GetLignes(); l++)
							if(axe->axe3d == OFF)
								glVertex2f(mx*CourbCour->pts->Element(l,0), my*CourbCour->pts->Element(l,1));
							else
								glVertex3f(mx*CourbCour->pts->Element(l,0), my*CourbCour->pts->Element(l,1),
								mz*CourbCour->pts->Element(l,2));
					}
					else
					{						
						glBegin(GL_LINES);
						for(l=0; l<CourbCour->pts->GetLignes(); l+=2)
							if(axe->axe3d == OFF)
							{
								glVertex2f(mx*CourbCour->pts->Element(l,0), my*CourbCour->pts->Element(l,1));
								glVertex2f(mx*CourbCour->pts->Element(l+1,0), my*CourbCour->pts->Element(l+1,1));
							}
							else
							{
								glVertex3f(mx*CourbCour->pts->Element(l,0), my*CourbCour->pts->Element(l,1),
									mz*CourbCour->pts->Element(l,2));
								glVertex3f(mx*CourbCour->pts->Element(l+1,0), my*CourbCour->pts->Element(l+1,1),
									mz*CourbCour->pts->Element(l+1,2));
							}
					}
					glEnd();
				}/*test affichage segments*/
			}/*boucle affichage normal/symetrique*/
		}/*test courbe non vide*/
		/*passage � la courbe suivante*/
		CourbCour=CourbCour->CourbSuiv;
	}/*boucle sur les courbes*/

	/*******************/
	/* trac� des meshs */
	/*******************/
	MeshCour=axe->Mesh;
	while (MeshCour!=NULL)
	{
		/*pour simplifier ecriture*/
		x=MeshCour->x; y=MeshCour->y; z=MeshCour->z;
		/*test matrices non nulles*/
		if ( (MeshCour->x != NULL)&&(MeshCour->y != NULL)
			&&((axe->axe3d == OFF)||(MeshCour->z != NULL)) )
		{
			/*determine si affichage sym�trique*/
			nbAff=1;
			if ((MeshCour->symX==ON)||(MeshCour->symX==ON)||(MeshCour->symX==ON)) nbAff=2;
			/*boucle affichage normal et sym�trique*/
			for(i=0; i<nbAff; i++)
			{
				/*mise a jour coeff mx, my, mz a 1.0 par defaut*/
				mx =1.0f; my=1.0f, mz=1.0f;
				if(i==1)
				{
					if (MeshCour->symX==ON) mx=-1.0f;
					if (MeshCour->symY==ON) my=-1.0f;
					if (MeshCour->symZ==ON) mz=-1.0f;
				}
				
				/*affichage des points*/
				if (MeshCour->points == ON)
				{
					glPointSize(5.0);
					glColor3dv(MeshCour->CouleurPoints);
					glBegin(GL_POINTS);
					for(l=0; l<MeshCour->x->GetLignes(); l++)
						for(c=0; c<MeshCour->x->GetColonnes(); c++)
							if(axe->axe3d == OFF)
								glVertex2f(mx*MeshCour->x->Element(l,c), my*MeshCour->y->Element(l,c));
							else
								glVertex3f(mx*MeshCour->x->Element(l,c), my*MeshCour->y->Element(l,c),
								mz*MeshCour->z->Element(l,c));
					glEnd();
				}
				
				/*affichage des segments ou faces*/
				if ((MeshCour->segments == ON)||(MeshCour->faces == ON))
				{
					/*couleur et mode d'affichage*/
					if (MeshCour->segments == ON)
					{
						glColor3dv(MeshCour->CouleurSegments);
						glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
					}
					else
					{
						glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
						glColor3dv(MeshCour->CouleurFaces);
					}
					/*boucle de trac�*/
					for(l=0; l<MeshCour->x->GetLignes()-1; l++)
					{
						glBegin(GL_QUAD_STRIP);
						for(c=0; c<MeshCour->x->GetColonnes(); c++)
						{
							if(axe->axe3d == OFF)
							{
								glVertex2f(mx*MeshCour->x->Element(l,c), my*MeshCour->y->Element(l,c));
								glVertex2f(mx*MeshCour->x->Element(l+1,c), my*MeshCour->y->Element(l+1,c));
							}
							else
							{
								/*calcul normale*/
								if(c < x->GetColonnes() - 1)
								{
									xv1 = x->Element(l+1,c+1) - x->Element(l,c+1);
									yv1 = y->Element(l+1,c+1) - y->Element(l,c+1);
									zv1 = z->Element(l+1,c+1) - z->Element(l,c+1);
									xv2 = x->Element(l+1,c+1) - x->Element(l+1,c);
									yv2 = y->Element(l+1,c+1) - y->Element(l+1,c);
									zv2 = z->Element(l+1,c+1) - z->Element(l+1,c);
									/*inversion normale si demand�*/
									if(MeshCour->InvNormales == OFF)
										glNormal3f(-(yv1*zv2-zv1*yv2)*mx,
											(xv1*zv2-zv1*xv2)*my,
											(yv1*xv2-xv1*yv2)*mz);
									else
										glNormal3f((yv1*zv2-zv1*yv2)*mx,
											(zv1*xv2-xv1*zv2)*my,
											(xv1*yv2-yv1*xv2)*mz);
								}
								/*envoi vecteurs*/
								glVertex3f(mx*MeshCour->x->Element(l,c), my*MeshCour->y->Element(l,c),
									mz*MeshCour->z->Element(l,c));
								glVertex3f(mx*MeshCour->x->Element(l+1,c), my*MeshCour->y->Element(l+1,c),
									mz*MeshCour->z->Element(l+1,c));
							}
						}
						glEnd();
					}
				}/*test affichage segments*/
			}/*boucle affichage normal/symetrique*/
		}/*test courbe non vide*/
		/*passage � la courbe suivante*/
		MeshCour=MeshCour->MeshSuiv;
	}/*boucle sur les Mesh*/

}



/***************/

/* AjoutCourbe */

/***************/

void AjoutCourbe(TAxe *axe, Courbe *courbe)

{

	Courbe* CourbCour;

	/*test s'il y a deja une courbe sur l'axe*/

	if (axe->Courb==NULL)

	{

		axe->Courb=courbe;

	}

	else

	{

		/*recherche derni�re courbe de la liste*/

		CourbCour=axe->Courb;

		int kk=0;
		while(CourbCour->CourbSuiv!=NULL) {

			CourbCour=CourbCour->CourbSuiv;
			kk++;
//			printf(" -> %d ", kk);
		}

		/*ajout courbe*/

		CourbCour->CourbSuiv=courbe;

	}

}

void showCourbe(Courbe *courbe)

{
	printf ("\nshowCourbe()");
	Courbe* CourbCour;

	/*test s'il y a deja une courbe sur l'axe*/

	if (courbe==NULL)

	{

		printf (" courbe==NULL");

	}

	else

	{

		/*recherche derni�re courbe de la liste*/

		CourbCour=courbe;

		int kk=0;
		while(CourbCour->CourbSuiv!=NULL) {

			CourbCour=CourbCour->CourbSuiv;
			kk++;
			printf ("showCourbe(): - %d ", kk);
		}

		/*ajout courbe*/

		CourbCour->CourbSuiv=courbe;

	}

}


void AjoutCCourbe(Courbe *head, Courbe *courbe)

{

	Courbe* CourbCour;

	/*test s'il y a deja une courbe sur l'axe*/
	printf ("\nAjoutCC");
	if (head==NULL)

	{
		printf ("\nAjoutCC courbeBegin == NULL");
		head = courbe;
		if (courbe==NULL)

		{
			printf (" /courbe == NULL");

		}
		if (head==NULL)

		{
			printf (" /courbeBegin == NULL");

		} else {
			printf (" /courbeBegin != NULL");

		} 


	}

	else

	{
		printf ("\nAjoutCC courbeBegin != NULL");
		/*recherche derni�re courbe de la liste*/

		CourbCour = head;

		int kk=0;
		while(CourbCour->CourbSuiv!=NULL) {

			CourbCour=CourbCour->CourbSuiv;
			kk++;
			printf(" -cc- -> %d ", kk);
		}

		/*ajout courbe*/

		CourbCour->CourbSuiv=courbe;

	}

}


/***************/

/* AjoutMesh   */

/***************/

void AjoutMesh(TAxe *axe, TMesh *mesh)

{

	TMesh* MeshCour;

	/*test s'il y a deja une courbe sur l'axe*/

	if (axe->Mesh==NULL)

	{

		axe->Mesh=mesh;

	}

	else

	{

		/*recherche derni�re courbe de la liste*/

		MeshCour=axe->Mesh;

		while(MeshCour->MeshSuiv!=NULL)

			MeshCour=MeshCour->MeshSuiv;

		/*ajout courbe*/

		MeshCour->MeshSuiv=mesh;

	}

}



/***************/

/* LibererMesh */

/***************/

void LibererMesh(TMesh *mesh)

{

	if (mesh != NULL)

	{

		delete(mesh->x);

		delete(mesh->y);

		delete(mesh->z);

		free(mesh);

	}

}



/*********************/

/* LibererCourbesAxe */

/*********************/

void LibererCourbesAxe(TAxe *axe)

{

	Courbe *CourbCour, *CourbSuiv;

	if (axe != NULL)

	{

		CourbCour = axe->Courb;

		while(CourbCour!=NULL)

		{

			CourbSuiv = CourbCour->CourbSuiv;

			delete CourbCour;

			CourbCour=CourbSuiv;

		}

		axe->Courb = NULL;

	}

}



/*******************/

/* LibererMeshsAxe */

/*******************/

void LibererMeshsAxe(TAxe *axe)

{

	TMesh *MeshCour, *MeshSuiv;

	if(axe != NULL)

	{

		MeshCour = axe->Mesh;

		while(MeshCour!=NULL)

		{

			MeshSuiv = MeshCour->MeshSuiv;

			LibererMesh(MeshCour);

			MeshCour=MeshSuiv;

		}

		axe->Mesh = NULL;

	}

}



/**************/

/* LibererAxe */

/**************/

void LibererAxe(TAxe *axe)

{

	LibererCourbesAxe(axe);

	LibererMeshsAxe(axe);

	free(axe);

}



/**************/

/* AjoutTexte */

/**************/

void AjoutTexte(TAxe *axe, char *texte, double taille, double orientation, double posx, double posy)

{

	char *TEXTE;

	unsigned int i,j, ichar, izero=0;

	Matrice *X,*Y,*T,*R;

	double xdec, ydec;

	Courbe *CourbChar;

	double TabXYpt[40][2];



	//mise en majuscule du texte

	TEXTE = _strupr( _strdup( texte ) );

	//mise a la taille

	for(i=0; i<40; i++){for(j=0; j<2; j++){TabXYpt[i][j] = TabXYpoint[i][j] * taille/7.0f;}}

	//rotation

	X=new Matrice(40,2); Y=new Matrice(40,2);

	for(i=0; i<40; i++) 
	{
		X->SetElement(i,0,TabXYpt[i][0]); 
		Y->SetElement(i,0,TabXYpt[i][1]);
	}

	Cart2Pol(X,Y,&T,&R); delete(X); delete(Y); 

	for(i=0; i<40; i++) 
	{
		T->AddElement(i,0, orientation*3.1415926f/180.0f);
	}

	Pol2Cart(T,R,&X,&Y);

	//Calcul decalage entre chiffre

	xdec = (double)cos(orientation*3.1415926f/180.0f)*taille*6.0f/7.0f;

	ydec = (double)sin(orientation*3.1415926f/180.0f)*taille*6.0f/7.0f;

	

	//boucle sur les chiffres

	for(i=0; i<strlen(TEXTE); i++)

	{

		//recherche indice du caractere courant

		ichar=100;

		for(j=0; j<38; j++)

		{

			if(TEXTE[i]==TabXYch[j][0])

			{

				ichar=j;

				break;

			}

		}

		//si trouve ...

		if(ichar!=100)

		{

			//recherche colonne du premier zero

			for(j=1; j<18; j++)

			{

				izero=j;

				if(TabXYch[ichar][j]==0) break;

			}

			//trac� du caract�re

			CourbChar=new Courbe("Text"); CourbChar->points = OFF;

			CourbChar->pts=Zeros(izero-1, 2);

			for(j=0; j<izero-1; j++)

			{

				CourbChar->pts->SetElement(j,0, X->Element(TabXYch[ichar][j+1]-1,0) + posx + xdec*i);

				CourbChar->pts->SetElement(j,1, Y->Element(TabXYch[ichar][j+1]-1,0) + posy + ydec*i);

			}

			AjoutCourbe(axe, CourbChar);

		}
	}
	delete(X); delete(Y); delete(T); delete(R);
}

/**************************/
/* AjoutForme3DKiteDesign */
/**************************/
void AjoutForme3DKiteDesign( TAxe *Axe3d, Forme3D* f3d, KiteDesign* kdExt, KiteDesign* kdInt, int mesh, int symetric) 
{
	for (int i = 0; i < kdExt->n_elements; i++) {
		KiteDesignElement* kde = kdExt->kiteDesignElements[i];
		kde -> ajoutCourbesToAxe3d(Axe3d, f3d, kdExt, EXT_SIDE, symetric);
	}

	for (int i = 0; i < kdInt->n_elements; i++) {
		KiteDesignElement* kde = kdInt->kiteDesignElements[i];
		kde -> ajoutCourbesToAxe3d(Axe3d, f3d, kdInt, INT_SIDE, symetric);
	}


}

void AjoutForme3D( TAxe *Axe3d, 
			 Matrice *XExt, Matrice *YExt, Matrice *ZExt,
			 Matrice *XInt, Matrice *YInt, Matrice *ZInt,
			 int mesh, int symetrie)
{
        
	Courbe *CourbCour;
	TMesh *MeshExt, *MeshInt;
	int NbNerv, i, j;
	int nbPtsInt, nbPtsExt;

	/*init variables pour simplifier ecriture*/
	nbPtsInt = XInt->GetColonnes(); nbPtsExt = XExt->GetColonnes();	NbNerv = XInt->GetLignes();
    /*creation courbes profils*/
	if(mesh == 0)
	{
        //printf ("\nbefore Axe3d->eclairage = OFF;");
		Axe3d->eclairage = OFF;
        //printf ("\nafter Axe3d->eclairage = OFF;");
        //printf ("\nNbNerv=%d", NbNerv);
		for (i=0; i<NbNerv; i++)
		{
            //      printf ("\ni=%d", i);
			/*extrados*/
			CourbCour = new Courbe("Extrados");
			CourbCour->points = OFF;
			CourbCour->pts = Zeros(XExt->GetColonnes(),3);
            //printf ("\n after Zeros()");
			for(j=0; j<XExt->GetColonnes(); j++)
			{
                // printf ("\nj=%d",j);
				CourbCour->pts->SetElement(j,0,XExt->Element(i,j));
				CourbCour->pts->SetElement(j,1,YExt->Element(i,j));
				CourbCour->pts->SetElement(j,2,ZExt->Element(i,j));
			}
            //  printf ("\nbefore 1");
			AjoutCourbe(Axe3d, CourbCour);
            //printf ("\nafter 1");
			if(symetrie==1) CourbCour->symX = ON;
			/*intrados*/
			CourbCour=new Courbe("Intrados");
			CourbCour->points = OFF;
			CourbCour->pts = Zeros(XInt->GetColonnes(),3);
			for(j=0; j<XInt->GetColonnes(); j++)
			{
				CourbCour->pts->SetElement(j,0,XInt->Element(i,j));
				CourbCour->pts->SetElement(j,1,YInt->Element(i,j));
				CourbCour->pts->SetElement(j,2,ZInt->Element(i,j));
			}
			AjoutCourbe(Axe3d, CourbCour);
			if(symetrie==1) CourbCour->symX = ON;
		}

		/* front line 'NOSE' */
		CourbCour=new Courbe("BA");
		CourbCour->points = OFF;
		CourbCour->pts = Zeros(NbNerv,3);
		for (i=0; i<NbNerv; i++)
		{
			CourbCour->pts->SetElement(i,0,XExt->Element(i,0));
			CourbCour->pts->SetElement(i,1,YExt->Element(i,0));
			CourbCour->pts->SetElement(i,2,ZExt->Element(i,0));
		}

		AjoutCourbe(Axe3d, CourbCour);

		if(symetrie==1) CourbCour->symX = ON;
		/* back line 'TALE'*/
		CourbCour=new Courbe("BF");
		CourbCour->points = OFF;
		CourbCour->pts = Zeros(NbNerv,3);
		for (i=0; i<NbNerv; i++)
		{
			int nXExt = XExt->GetColonnes();
			CourbCour->pts->SetElement(i,0,XExt->Element(i, nXExt - 1));
			CourbCour->pts->SetElement(i,1,YExt->Element(i, nXExt - 1));
			CourbCour->pts->SetElement(i,2,ZExt->Element(i, nXExt - 1));
		}

		AjoutCourbe(Axe3d, CourbCour);
		if(symetrie==1) CourbCour->symX = ON;
        }

	/*creation mesh*/

	else
	{
		Axe3d->eclairage = ON;
		MeshExt=CreerMesh(); MeshInt=CreerMesh();
		MeshInt->InvNormales=ON;
		MeshExt->segments=OFF; MeshExt->faces=ON;
		MeshInt->segments=OFF; MeshInt->faces=ON;
		MeshExt->x=new Matrice(nbPtsExt, NbNerv);
		MeshExt->y=new Matrice(nbPtsExt, NbNerv);
		MeshExt->z=new Matrice(nbPtsExt, NbNerv);
		MeshInt->x=new Matrice(nbPtsInt, NbNerv);
		MeshInt->y=new Matrice(nbPtsInt, NbNerv);
		MeshInt->z=new Matrice(nbPtsInt, NbNerv);
		if(symetrie==1)
		{
			MeshExt->symX = ON;
			MeshInt->symX = ON;
		}

                /*boucle sur les points du profil en extrados*/
		for (i=0; i<nbPtsExt; i++)
		{
			/*boucle � partir de la 1ere nervure du centre vers l'Intr�mit�*/
			for (j=0; j<NbNerv; j++)
			{
				MeshExt->x->SetElement(i,j, XExt->Element(j,i));
				MeshExt->y->SetElement(i,j, YExt->Element(j,i));
				MeshExt->z->SetElement(i,j, ZExt->Element(j,i));
			}
		}

                /*boucle sur les points du profil en intrados*/
		for (i=0; i<nbPtsInt; i++)
		{
			/*boucle � partir de la 1ere nervure du centre vers l'Intr�mit�*/
			for (j=0; j<NbNerv; j++)
			{
				MeshInt->x->SetElement(i,j, XInt->Element(j,i));
				MeshInt->y->SetElement(i,j, YInt->Element(j,i));
				MeshInt->z->SetElement(i,j, ZInt->Element(j,i));
			}
		}

		/*ajout mesh a l'axe3d*/
		AjoutMesh(Axe3d, MeshExt); AjoutMesh(Axe3d, MeshInt);
	}
}

void getCourbeFromProfilGeom(ProfilGeom* pg, Courbe** courbeExt, Courbe** courbeInt){
	*courbeExt = new Courbe("ProfileExt");
    (*courbeExt)->points = OFF;
    (*courbeExt)->symX = OFF;

	int nExt = pg->ExtProf->GetLignes();
	int nInt = pg->IntProf->GetLignes();
	(*courbeExt)->pts = new Matrice(nExt,2);
	for (int i = 0; i < nExt; i++) {
		(*courbeExt)->pts->SetElement(i, 0, pg->ExtProf->Element(i, 0));
		(*courbeExt)->pts->SetElement(i, 1, pg->ExtProf->Element(i, 1));
	}
	(*courbeInt) = new Courbe("ProfileInt");
    (*courbeInt)->points = OFF;
    (*courbeInt)->symX = OFF;
	(*courbeInt)->pts = new Matrice(nInt,2);
	for (int i = 0; i < nInt; i++) {
		(*courbeInt)->pts->SetElement( i, 0, pg->IntProf->Element(i, 0));
		(*courbeInt)->pts->SetElement( i, 1, pg->IntProf->Element(i, 1));
	}
}

void ajoutFormeProjectionCourbesToAxe(TAxe* axe, FormeProjection* fp, KiteDesign* kd, int symetric, double dy, int dir) {
	printf ("\n ajoutFormeProjectionCourbesToAxe()");
	double ymult = 1;
	if (dir == DIR_NOSE_DOWN) ymult=-1;
	int n = fp->X->GetLignes();
	int m = fp->X->GetColonnes();
	//printf ("\n n=%d, m=%d", n, m);
	
	Courbe* courbe0 = new Courbe("CourbeUp");
	courbe0->points = OFF;
	if (symetric) courbe0->symX = ON;
	courbe0->pts = new Matrice(n, 2);

	Courbe* courbe1 = new Courbe("CourbeDown");
	courbe1->points = OFF;
	if (symetric) courbe1->symX = ON;
	courbe1->pts = new Matrice(n, 2);
	

	double ymin = 1000000;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			double y = fp->Y->Element(i, j);
			if (y*ymult < ymin) ymin = y*ymult;
		}
	}

	for (int i = 0; i < n; i++) {
		//printf ("\n i=%d", i);
		Courbe* courbe = new Courbe("Courbe");
		courbe->points = OFF;
		courbe->symX = OFF;
		if (symetric) courbe->symX = ON;
		courbe0->pts->SetElement(i, 0, fp->X->Element(i, 0));
		courbe0->pts->SetElement(i, 1, dy + ymult*fp->Y->Element(i, 0)- ymin);
		courbe->pts = new Matrice(m, 2);
		for (int j = 0; j < m; j++) {
			double x = fp->X->Element(i, j);
			courbe->pts->SetElement(j, 0, x);
			double y = fp->Y->Element(i, j);
			courbe->pts->SetElement(j, 1, dy + (ymult*y - ymin));
			//printf ("\n (%d, %d)    %f, %f", i, j, fp->X->Element(i, j), fp->Y->Element(i, j));
		}
		courbe1->pts->SetElement(i, 0, fp->X->Element(i, m-1) );
		courbe1->pts->SetElement(i, 1, dy + ymult*fp->Y->Element(i, m-1)-ymin);

		//printf ("\n AjoutCourbe()");
		AjoutCourbe(axe, courbe);
	}
	
	//calculate kite design courbes

	for (int i = 0; i < kd->n_elements; i++) {
		KiteDesignElement* kde = kd->kiteDesignElements[i];
		kde -> ajoutCourbesToAxe(axe, fp, symetric, dy, ymult, ymin);
	}

	AjoutCourbe(axe, courbe0);
	AjoutCourbe(axe, courbe1);
}

/**********************/
/* AjoutPtsSuspentage */
/**********************/

void AjoutPtsSuspentage(
						TAxe *Axe3d, Forme *forme, Matrice *IntProfCent,
						Matrice* /*XExt*/, Matrice* /*YExt*/, Matrice* /*ZExt*/,
						Matrice *XInt, Matrice *YInt, Matrice *ZInt,
						int symetrie, Matrice **PosSuspentes)

{
	Courbe *CourbCour = NULL;
	Matrice *interpSuspente = NULL;
	int NbNerv, i, j;



	/*init variables pour simplifier ecriture*/

	NbNerv = XInt->GetLignes();

	//creation courbe des pts de suspentage
	CourbCour=new Courbe("Points de suspentage");
	CourbCour->points = ON;
	CourbCour->segments =OFF;
	CourbCour->pts = Zeros(5*NbNerv,3); //par defaut 4 pts de suspentage
	if(symetrie==1) CourbCour->symX = ON;

	//boucle nervures
	for (i=0; i<NbNerv; i++)
	{
		//creation matrice intermediaire pour interpolation points de suspentage
		interpSuspente=Zeros(IntProfCent->GetLignes(),2);
		for (j=0; j<IntProfCent->GetLignes(); j++) 
			interpSuspente->SetElement(j,0,IntProfCent->Element(j,0));
		//interpolation X
		for (j=0; j<IntProfCent->GetLignes(); j++) 
			interpSuspente->SetElement(j,1,XInt->Element(i,j));

		CourbCour->pts->SetElement(i*5,  0,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosA));
		CourbCour->pts->SetElement(i*5+1,0,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosB));
		CourbCour->pts->SetElement(i*5+2,0,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosC));
		CourbCour->pts->SetElement(i*5+3,0,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosD));
                CourbCour->pts->SetElement(i*5+4,0,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosE));
		//interpolation Y

		for (j=0; j<IntProfCent->GetLignes(); j++) 
			interpSuspente->SetElement(j,1,YInt->Element(i,j));

		CourbCour->pts->SetElement(i*5,  1,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosA));
		CourbCour->pts->SetElement(i*5+1,1,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosB));
		CourbCour->pts->SetElement(i*5+2,1,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosC));
		CourbCour->pts->SetElement(i*5+3,1,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosD));
		CourbCour->pts->SetElement(i*5+4,1,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosE));
		//interpolation Z

		for (j=0; j<IntProfCent->GetLignes(); j++) 
			interpSuspente->SetElement(j,1,ZInt->Element(i,j));

		CourbCour->pts->SetElement(i*5,  2,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosA));
		CourbCour->pts->SetElement(i*5+1,2,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosB));
		CourbCour->pts->SetElement(i*5+2,2,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosC));
		CourbCour->pts->SetElement(i*5+3,2,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosD));
		CourbCour->pts->SetElement(i*5+4,2,InterpLinX(interpSuspente, forme->m_pProfils[i]->m_fPosE));
	}

	AjoutCourbe(Axe3d, CourbCour);

	*PosSuspentes = CourbCour->pts;

	delete(interpSuspente);
}