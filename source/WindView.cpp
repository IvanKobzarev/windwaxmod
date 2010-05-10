/*  WindView.c			
*  visualise forme 3d  ...

* issue de WindShape ou généré à partir d'une autre application

* mais respectant le format WindShape !
*/

#pragma warning(disable:4786)
#pragma warning(disable:4514)


/***********/

/* include */

/***********/


#include <stdlib.h>

#include <stdio.h>

#include <string.h>

#include <math.h>

#include <GL/glut.h>


#include <afx.h>		//class CString

#include <afxdlgs.h>	//class CFileDialog



#include "../commun/matrice.h"

#include "../commun/geom.h"

#include "../commun/plot.h"

#include "../commun/fichier.h"

#include "GL/glui.h"



/**********/

/* define */

/**********/



#define AUTEUR_DATE "(Thierry Pebayle - 18/08/2002)\n[Export FGen + Refactoring Olivier REMAUD 2007]"



/**********************/

/* variables globales */

/**********************/



/*live-variable GLUI*/

int VisuFace=0;		/*choix visu nervures(0) ou surfaces(1)*/

int ProjOrthoPers=0;	/*choix projection orthogonale ou perspective*/

int	VisuSymetrique=1;	/*visualise les 2 ou 1 moitié de l'aile en 3D*/

int VisuPtsSuspentes=1; /*visualise pts de suspentage*/



/*global pour control rotation, ...*/

int xSouris, ySouris;		/*coordonnées fenetre de la souris*/

int Rotation3d = OFF;		/*flag rotation 3d -> bouton gauche souris*/

int Translation3d = OFF;	/*flag translation 3d -> bouton droit souris*/

int ZoomTwist3d = OFF;		/*flag twist 3d -> bouton gauche et droit souris*/



/*fenetres OpenGL*/

int Fenetre3d;		/*visu 3D*/



/* axes graphiques*/

TAxe *Axe3d;		

TAxe *AxeSel; 



/*forme*/

Forme *F;



/*divers tableaux*/

Matrice *IntProfCent, *ExtProfCent, *IntProfBout ,*ExtProfBout;



char NomFichierForme[255];

GLUI_StaticText *FicForme;



//pour display des infos de la forme

TInfoForme info;

GLUI_StaticText *textSurface, *textEnvergure, *textAllongement, *textCorde, *textLarg;



/***********************************************/

/*liste des procedures definies dans ce fichier*/

/***********************************************/



void display(void);

void reshape(int w, int h);

void motion(int x, int y);

void BoutonSouris(int button, int state, int x, int y);

void PositionneFenetres(void);

void keyboard(unsigned char key, int x, int y);

void InitFenetre(void);

void InitLumiere(void);



void CalculVue3d(void);



void ModifProjection3d( int control );

void ModifVisu3d( int control );

void ModifVisuSymetrique( int control );



void Quitter( int control );



void ChargerFichierForme( int control );

void SauverFichierForme( int control );



void InterpoleProfilBout(Matrice** XYBout, Matrice* XYCent);



/***********/

/* MajInfo */

/***********/

void MajInfo(void)

{

	char surf[100], env[100], all[100], cord[100], larg[100];

	//maj chaine de caractère

	sprintf(surf,"\n\tSurface: plat=%2.2fm2, proj=%2.2fm2, r=%2.1f%%",

		info.surface, info.surfaceProj, info.surfaceProj/info.surface*100.0f);

	sprintf(env,"\n\tEnvergure: plat=%2.2fm, proj=%2.2fm, r=%2.1f%%",

		info.envergure, info.envergureProj, info.envergureProj/info.envergure*100.0f);

	sprintf(all,"\n\tAllongement: plat=%2.2f, proj=%2.2f, r=%2.1f%%",

		info.allongement, info.allongementProj, info.allongementProj/info.allongement*100.0f);

	sprintf(cord,"\n\tCorde mini = %2.3fm, maxi = %2.3fm", info.cordeMin, info.cordeMax);

	sprintf(larg,"\n\tLarg caisson mini = %2.3fm, maxi = %2.3fm", info.largMin, info.largMax);

	//maj interface

	textSurface->set_text(surf);

	textEnvergure->set_text(env);

	textAllongement->set_text(all);

	textCorde->set_text(cord);

	textLarg->set_text(larg);

}



/********************/

/* SauverFichierDXF */

/********************/

void SauverFichierDXF( int /*control*/ )

{

	CString NomFichier;

	LPTSTR PtrNomFichier;

	

	//ouverture boite de dialogue

	CFileDialog DlgOpen(FALSE, NULL, "*.dxf", OFN_OVERWRITEPROMPT, NULL, NULL);

	if (DlgOpen.DoModal()==IDOK)

	{

		//recupere nom de fichier

		NomFichier=DlgOpen.GetPathName();

		PtrNomFichier = NomFichier.GetBuffer(1);



		//ecriture fichier DXF

		EcritureFichierDXF(PtrNomFichier, Axe3d);

	}

}



/***********************/

/* ChargerFichierForme */

/***********************/

void ChargerFichierForme( int /*control*/ )

{

	CString NomFichier;

	char* PtrNomFichier;

	Forme* oldF;

	

	//ouverture boite de dialogue

	CFileDialog DlgOpen(TRUE, NULL, "*.txt", OFN_OVERWRITEPROMPT, NULL, NULL);

	if (DlgOpen.DoModal()==IDOK)

	{

		//recupere nom de fichier

		NomFichier=DlgOpen.GetPathName();

		PtrNomFichier=NomFichier.GetBuffer(1);

		strcpy(NomFichierForme, PtrNomFichier);

		

		//lecture fichier de forme

		//LibererForme(F);

		oldF = F;

		F = LectureFichierForme(NomFichierForme);

		

		//libere memoire ancienne forme

		delete(oldF);



		//maj interface GLUI

		FicForme->set_text(NomFichierForme);

		

		//chargement profils

		LectureFichierProfil(F->m_strNomProfilCent.c_str(), &ExtProfCent, &IntProfCent);

		LectureFichierProfil(F->m_strNomProfilBout.c_str(), &ExtProfBout, &IntProfBout);

		InterpoleProfilBout(&ExtProfBout, ExtProfCent);

		InterpoleProfilBout(&IntProfBout, IntProfCent);

	}

	/*mise a jour affichage*/

	CalculVue3d();

	display();



	/*maj info*/

	CalculInfoForme(F, &info);

	MajInfo();

}





/*********************/

/* ModifProjection3d */

/*********************/

void ModifProjection3d( int /*control*/ )

{

	/*test type de proj.*/

	if(ProjOrthoPers == 0)

		Axe3d->proj = PROJ_ORTHOGONALE;

	else

		Axe3d->proj = PROJ_PERSPECTIVE;

	/*retrace uniquement la vue 3d*/

	VisuAxe(Axe3d); glutSwapBuffers();

}



/***************/

/* ModifVisu3d */

/***************/

void ModifVisu3d( int /*control*/ )

{

	/*recalcul et retrace la vue 3d*/

	CalculVue3d(); VisuAxe(Axe3d); glutSwapBuffers();

}



/***********************/

/* ModifVisuSymetrique */

/***********************/

void ModifVisuSymetrique( int /*control*/ )

{

	/*recalcul et retrace la vue 3d*/

	CalculVue3d(); VisuAxe(Axe3d); glutSwapBuffers();

}



/***********/

/* Quitter */

/***********/

void Quitter( int /*control*/ )

{
	if ( MessageBox( NULL, "Voulez-vous vraiment quitter l'application ?", "Confirmation", MB_YESNO ) == IDYES )
	{
		/*liberation espace memoire*/

		LibererAxe(Axe3d);

		delete(F);

		delete(IntProfCent);

		delete(ExtProfCent);

		delete(IntProfBout);

		delete(ExtProfBout);

		/*quitte interface*/

		exit(0);
	}
}



/***************/

/* CalculVue3d */

/***************/

void CalculVue3d(void)

{

	Matrice *XExt,*YExt,*ZExt; //coordonnées 3D extrados

	Matrice *XInt,*YInt,*ZInt; //coordonnées 3D intrados

	Matrice *PtsSuspentes; //pour recuperer la position 3D des pts de suspentage

	

	/*par defaut destruction des courbes ou mesh de l'axe 3D*/

	LibererCourbesAxe(Axe3d); LibererMeshsAxe(Axe3d);

	

	/*calcul forme*/

	CalculForme3D(F,

		ExtProfCent, IntProfCent, ExtProfBout, IntProfBout,

		&XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);

	

	/*ajoute visu 3D de la forme*/

	AjoutForme3D(Axe3d, XExt, YExt, ZExt, XInt, YInt, ZInt, VisuFace, VisuSymetrique);



	/*ajoute points de suspentage*/

	if(VisuPtsSuspentes){

		AjoutPtsSuspentage(Axe3d, F, IntProfCent, XExt, YExt, ZExt, XInt, YInt, ZInt,

			VisuSymetrique, &PtsSuspentes);}



	/*liberation matrices intermediaires !!!*/

	delete(XExt); delete(YExt); delete(ZExt);

	delete(XInt); delete(YInt); delete(ZInt);

}



/****************/

/* display      */

/****************/

void display(void)
{
	//CalculVue3d();

	VisuAxe(Axe3d);

	glutSwapBuffers();

}


/*****************/

/* reshape       */

/*****************/

void reshape(int w, int h)
{
	glViewport(0, 0, w, h);

	display();

}



/*****************/

/* motion        */

/*****************/

void motion(int x, int y)

{



	/*message pour test*/

	//printf("\n<%d,%d>",x,y);

	

	/*axe 3d -> rotation*/

	if (Rotation3d)

	{

		Axe3d->azimuth = Axe3d->azimuth + (x - xSouris);

		Axe3d->incidence = Axe3d->incidence + (y - ySouris);

		xSouris = x; ySouris = y;

	}

	

	/*axe 3d -> translation*/

	else if (Translation3d)

	{

		Axe3d->xcam = Axe3d->xcam + (float)(x - xSouris)/100.0;

		Axe3d->ycam = Axe3d->ycam - (float)(y - ySouris)/100.0;

		xSouris = x; ySouris = y;

	}

	

	/*axe 3d -> zoom*/

	else if (ZoomTwist3d)

	{

		Axe3d->twist = Axe3d->twist - (float)(x - xSouris);

		Axe3d->zcam = Axe3d->zcam + (float)(y - ySouris)/100.0;

		xSouris = x; ySouris = y;

	}

	

	/*affichage*/

	//printf("\nx=%d, y=%d, xs=%f, ys=%f, dist=%f       ",x,y,xs,ys,distmin);

	

	glutPostRedisplay();

	//glutSwapBuffers();

}



/****************/

/* BoutonSouris */

/****************/

void BoutonSouris(int button, int state, int x, int y)

{

	int Fenetre;

	

	/*memorise position de la souris dans la fenetre*/

	xSouris =x; ySouris = y;

	

	/*par defaut pas de point ni d'axe selectionné 

	et pas de rotation/translation/zoom 3d

	et ne pas redessinner tous les axes*/

	AxeSel = NULL;

	

	/*Rotation3d = OFF;

	Translation3d = OFF;

	ZoomTwist3d = OFF;*/

	

	/*recupere axe correspondant à la fenetre courante*/

	Fenetre = glutGetWindow();

	if (Fenetre == Fenetre3d) {AxeSel=Axe3d; /*printf("Axe3d ->");*/}

	

	/*test debut/fin, rotation/translation/zoom pour Axe3d*/

	if (AxeSel == Axe3d)

	{

		if (state == GLUT_DOWN)

		{

			switch(button)

			{

			case GLUT_LEFT_BUTTON :

				if (Translation3d == ON) {Translation3d = OFF; ZoomTwist3d = ON;}

				else Rotation3d = ON;

				break;

			case GLUT_RIGHT_BUTTON :

				if (Rotation3d == ON) {Rotation3d = OFF; ZoomTwist3d = ON;}

				else Translation3d = ON;

				break;

			case GLUT_MIDDLE_BUTTON :

				ZoomTwist3d = ON;

				break;

			default: printf("\n type de bouton souris inconnu !!!");

			}

		}

		else /*state = GLUT_UP*/

		{

			Rotation3d = OFF;

			Translation3d = OFF;

			ZoomTwist3d = OFF;

		}

	}	

}



/*****************/

/* keyboard      */

/*****************/

void keyboard(unsigned char key, int /*x*/, int /*y*/)

{

	switch (key)

	{

	case 27: Quitter(0); break;

	default:;

	}

}



/*****************/

/* InitFenetre   */

/*****************/

void InitFenetre(void)

{

	glEnable(GL_LINE_STIPPLE);

	glutDisplayFunc(&display);

	glutReshapeFunc(&reshape);

	glutKeyboardFunc(&keyboard);

	glutMotionFunc(&motion);

	glutMouseFunc(&BoutonSouris);

}



/*****************/

/* InitLumiere   */

/*****************/

void InitLumiere(void)

{

	float mat_specular[]={1.0, 1.0, 1.0, 1.0};

	float mat_shininess[]={50.0};

	float light_position[]={0.0, 0.0, 1.0, 0.0};

	float white_light[]={1.0, 1.0, 0.0, 1.0};

	

	glShadeModel(GL_SMOOTH);

	/*glClearColor(0.0, 0.0, 0.0, 0.0);*/

	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);

	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

	

	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);

	glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);

	

	glEnable(GL_LIGHTING);

	glEnable(GL_LIGHT0);

	glEnable(GL_DEPTH_TEST);

	glEnable(GL_NORMALIZE);

}



/*****************/

/* M A I N       */

/*****************/

int main(int argc, char** argv)
{
	int ws, hs;

	float EpaiRelProfCent, EpaiRelProfBout;



	/*message*/

	printf("\nWindView ");

	printf(AUTEUR_DATE);

	printf("\n");



	/*init glut*/
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize (500, 500);
	glutInitWindowPosition (100, 100);

	

	/* creation fenetre/axe/courbes Axe3d*/

	Fenetre3d=glutCreateWindow ("Visualisation 3D");

	InitFenetre(); InitLumiere();

	Axe3d=CreerAxe(Fenetre3d); Axe3d->axe3d=ON;

	

	/*recupere taille de l'ecran en pixel*/

	ws = glutGet( GLUT_SCREEN_WIDTH );

	hs = glutGet( GLUT_SCREEN_HEIGHT );

	

	/*redefini taille et position fenetre a partir de coordonnées normalisées*/

	glutSetWindow( Fenetre3d );

	glutPositionWindow((int)(0.1*(float)ws), (int)(0.1*(float)hs));

	glutReshapeWindow((int)(0.8*(float)ws), (int)(0.8*(float)hs));

	

	//lecture fichier de forme

	strcpy(NomFichierForme,"FormeDefaut.txt");

	F = LectureFichierForme(NomFichierForme);

	
	if ( F == NULL )
		return 1;

	//chargement profils

	LectureFichierProfil(F->m_strNomProfilCent.c_str(), &ExtProfCent, &IntProfCent);

	LectureFichierProfil(F->m_strNomProfilBout.c_str(), &ExtProfBout, &IntProfBout);



	//calcul epaisseur relative ?

	EpaiRelProfCent=EpaisseurRelative(ExtProfCent, IntProfCent);

	EpaiRelProfBout=EpaisseurRelative(ExtProfBout, IntProfBout);



	//applique distribution des x du profil du bout comme au centre

	InterpoleProfilBout(&ExtProfBout, ExtProfCent);

	InterpoleProfilBout(&IntProfBout, IntProfCent);

	

	/****************************************/

	/*         Here's the GLUI code         */

	/****************************************/

	

	/*affichage No de version et creation fenetre*/

	printf( "Ce programme utilise GLUI version: %3.2f\n", GLUI_Master.get_version() );

	GLUI *glui = GLUI_Master.create_glui( "Boite de dialogue",0,0,0);

	

	/****** rollout fichiers *****/

	GLUI_Rollout *RolloutFichiers =	glui->add_rollout("Fichiers",false);

	

	

	/*panel forme*/

	GLUI_Panel *panel_forme = glui->add_panel_to_panel(RolloutFichiers,"Forme");

	/*nom fichier forme*/

	FicForme = glui->add_statictext_to_panel( panel_forme, "???" );

	FicForme->set_text(NomFichierForme);

	glui->add_column_to_panel(panel_forme, false);

	GLUI_Button *bouton = glui->add_button_to_panel( panel_forme, "Load", 0, &ChargerFichierForme );

	bouton->set_w(10);

	/*fichier DXF*/

	bouton = glui->add_button_to_panel( RolloutFichiers, "3D->DXF", 0, &SauverFichierDXF );

	bouton->set_w(10);

	

	/****** rollout visu ******/

	//glui->add_column(false);

	GLUI_Panel *RolloutVisu = glui->add_rollout( "Visualisation",false);

	/*visualisation nervures/surface*/

	GLUI_Panel *panel_visu = glui->add_panel_to_panel( RolloutVisu, "Type" );

	GLUI_RadioGroup *radio_visu =

		glui->add_radiogroup_to_panel(panel_visu,&VisuFace,0,&ModifVisu3d);

	glui->add_radiobutton_to_group( radio_visu, "Nervures" );

	glui->add_radiobutton_to_group( radio_visu, "Surfaces" );

	/*choix projection orthogonale ou perspective*/

	GLUI_Panel *panel_proj = glui->add_panel_to_panel( RolloutVisu, "Projection" );

	GLUI_RadioGroup *radio_proj =

		glui->add_radiogroup_to_panel(panel_proj,&ProjOrthoPers,0,&ModifProjection3d);

	glui->add_radiobutton_to_group( radio_proj, "Orthogonale" );

	glui->add_radiobutton_to_group( radio_proj, "Perspective" );

	/*choix visu symétriqe ou non*/

	glui->add_checkbox_to_panel(

		RolloutVisu,"symetrique",&VisuSymetrique,0,&ModifVisuSymetrique);

	/*choix visu des points de suspentage*/

	glui->add_checkbox_to_panel(

		RolloutVisu,"pts suspentes",&VisuPtsSuspentes,0,&ModifVisuSymetrique);

	

	/**** pour affichage info ****/

	GLUI_Panel *RolloutInfo = glui->add_rollout( "Info forme",false);

	textSurface = glui->add_statictext_to_panel(RolloutInfo, "surface..."); 

	textEnvergure = glui->add_statictext_to_panel(RolloutInfo, "envergure..."); 

	textAllongement = glui->add_statictext_to_panel(RolloutInfo, "allongement..."); 

	textCorde = glui->add_statictext_to_panel(RolloutInfo, "corde..."); 

	textLarg = glui->add_statictext_to_panel(RolloutInfo, "largeur caisson...");



	/**** bouton quitter ****/

	glui->add_button( "Quitter", 0, Quitter );

	

	/**** Link windows to GLUI, and register idle callback ******/

	glui->set_main_gfx_window( Fenetre3d );

	

	/* We register the idle callback with GLUI, not with GLUT */

	GLUI_Master.set_glutIdleFunc( NULL );

	

	/*premiers calculs et trace des axes*/

	AxeSel = Axe3d;

	CalculVue3d();

	display();



	/*maj info*/

	CalculInfoForme(F, &info);

	MajInfo();



	/*boucle glut*/

	glutMainLoop();

	return 0;
}
