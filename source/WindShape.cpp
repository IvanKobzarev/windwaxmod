/*  WindShape.c			
*  modeleur de forme pour ailes de traction ...
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
#include <afx.h>		//class CString
#include <afxdlgs.h>	//class CFileDialog


#include "../commun/profil.h"
#include "../commun/matrice.h"
#include "../commun/geom.h"
#include "../commun/plot.h"
#include "../commun/rasklad.h"
#include "../commun/pince.h"
#include "../commun/fichier.h"
#include "../commun/logger.h"
#include "../commun/patternsproject.h"

#include "GL/glui.h"
#include "GL/glut.h"
/**********/
/* define */
/**********/
#define BTN_SLING 1
#define BTN_EPAIREL	2
#define BTN_MORPHING 3
#define BTN_VRILLAGE 4
#define BTN_DIEDRE 5

#define AUTEUR_DATE "(Thierry Pebayle - 18/08/2002)\n[Export FGen + Refactoring Olivier REMAUD 2007]"
#define sqr(f1) ((f1)*(f1))
#define NBP_BEZ		20	/*nombre de points pour les courbes de bezier*/

/*pour identification callback GLUI*/
#define TAILLE_PLUS_1	1
#define TAILLE_PLUS_5	2
#define TAILLE_MOINS_1	3
#define TAILLE_MOINS_5	4
#define PROF_CENT	1
#define PROF_BOUT	2

/**********************/
/* variables globales */
/**********************/

/** These are the live variables passed into GLUI ***/
float CoeffProgGeom=1.0f;  /*coeff progression g�om�trique*/
float CoeffExp=0.0f;
float EpaiRelCent=17.0;	/*epaisseur relative du profil central*/ 
int NbCaiss=20;			/*nombre de caissons au total*/
int VisuFace=0;		/*choix visu nervures(0) ou surfaces(1)*/
int ProjOrthoPers=0;	/*choix projection orthogonale ou perspective*/
int pinceType=0;
float VrillageCent; /*vrillage en degres au bout*/
float VrillageBout; /*vrillage en degres au bout*/
int	VisuSymetrique=1;	/*visualise les 2 ou 1 moiti� de l'aile en 3D*/
int	RaskladSymetrique=0;
int VisuPtsSuspentes=1; /*visualise pts de suspentage*/
int AxesAuto=1; /*Virer le dimensionnement auto*/

/*divers*/
double *PointSel;	/*ptr sur point selectionn�*/
TAxe *AxeSel;		/*ptr sur axe selectionn�*/
int xSouris, ySouris;		/*coordonn�es fenetre de la souris*/
int VisuTousLesAxes = 0;	/*flag pour demande de dessin de tous les axes*/
int Rotation3d = OFF;		/*flag rotation 3d -> bouton gauche souris*/
int Translation3d = OFF;	/*flag translation 3d -> bouton droit souris*/
int ZoomTwist3d = OFF;		/*flag twist 3d -> bouton gauche et droit souris*/
int NbNerv=11;		/*nombre de nervures de la demi-aile*/	
float EpaiRelProfCent, EpaiRelProfBout; /*epaisseur relative des profils*/

/*fenetres OpenGL*/
int FenetreForme;	/*forme a plat*/
int FenetreDiedre;	/*diedre*/
int FenetreMorphing;/*profil evolutif*/
int FenetreVrillage;/*vrillage de l'aile*/
int FenetreEpaiRel;	/*Epaisseur relative*/
int Fenetre3d;		/*visu 3D*/

/*axes correspondant aux fenetres*/
TAxe *AxeForme;
TAxe *AxeDiedre;	
TAxe *AxeMorphing;
TAxe *AxeVrillage;
TAxe *Axe3d;		
TAxe *AxeEpaiRel;

/*points de controle = courbes de 4 pts(x,y)*/
Courbe *CtrlNez, *CtrlFui, *CtrlA, *CtrlB, *CtrlC, *CtrlD, *CtrlE;
Courbe *CtrlDiedre, *CtrlMorphing, *CtrlVrillage, *CtrlEpaiRel;

/*Forme!*/
Forme *F=NULL;

/*courbes resultantes*/
Courbe *CourbNez, *CourbFui, *CourbA, *CourbB, *CourbC, *CourbD, *CourbE;
Courbe *CourbNerv, *CourbDiedre, *CourbMorphing;
Courbe  *CourbVrillage, *CourbEpaiRel;

/*divers tableaux*/
Matrice *IntProfCent, *ExtProfCent, *IntProfBout ,*ExtProfBout;
Matrice *XNerv = NULL;

//char NomProfilCent[255], NomProfilBout[255];
char NomFichierForme[255];
char NomFichierProject[255];
char NomRepPoints[255];
//interface GLUI
GLUI *glui;
GLUI_StaticText *ProfCent, *ProfBout, *FicForme,*FicProject, *FicRepPoints;
//pour display des infos de la forme
TInfoForme info;
GLUI_StaticText *textSurface, *textTmp, *textEnvergure, *textAllongement, *textCorde, *textLarg;
bool wasLoad = false;
bool loadConfig = false;

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
void CalculCourbesCtrl(void);
void CalculVue3d(void);
void ModifNbrAlveoles( int control );
void ModifTaille( int control );
void ModifProjection3d( int control );
void ModifVisu3d( int control );
void ModifVisuSymetrique( int control );
void ModifVisuPtsSuspentes( int control );
void ModifAxesAuto( int control );
void ModifProgGeom( int control );
void ModifVrillage( int control );

void Quitter( int control );
void ExportTemplate ( int control );
void LoadCourbe ( int control );
void ChargerFichierForme( int control );
void SauverFichierForme( int control );
void ChargerFichierProfil( int control );
void ChargerFichierConfig(char *NomFicConfig);
void SauverFichierConfig(char *NomFicConfig);
void CalculTableauForme(void);

/***********/
/* MajInfo */
/***********/
void MajInfo(void)
{
	char surf[100], env[100], all[100], cord[100], larg[100];
	//maj chaine de caract�re
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

/**********************/
/* CalculTableauForme */
/**********************/

void CalculTableauForme(void)
{
	double L;
	int i;
	/*reset tableau*/
	F->AllocateProfils(NbNerv);
	/*boucle � partir de la 1ere nervure du centre vers l'extr�mit�*/
	for (i=0; i<NbNerv; i++)
	{	
		//longueur nervure
		L = CourbNez->pts->Element(i,1)-CourbFui->pts->Element(i,1);
		F->m_pProfils[i]->m_fLength = L; 
		//epaisseur relative
		F->m_pProfils[i]->m_fWidth = EpaiRelCent*CourbEpaiRel->pts->Element(i,1);
		//xnez	//ynez	//znez
		F->m_pProfils[i]->m_fNezX = CourbDiedre->pts->Element(i,0);	
		F->m_pProfils[i]->m_fNezY = CourbDiedre->pts->Element(i,1);
		F->m_pProfils[i]->m_fNezZ = CourbNez->pts->Element(i,1);	
		
		//inclinaison de la nervure par rapport a l'horizontale en radian
		if(i==0)
			F->m_pProfils[i]->m_fInclin = 1.5708f; // pi/2 au centre par defaut
		else if(i==NbNerv-1)
			F->m_pProfils[i]->m_fInclin = 1.5708f + 
			atan2(CourbDiedre->pts->Element(i,1)-CourbDiedre->pts->Element(i-1,1),
			CourbDiedre->pts->Element(i,0)-CourbDiedre->pts->Element(i-1,0));
		else
			F->m_pProfils[i]->m_fInclin = 1.5708f + 
			atan2(CourbDiedre->pts->Element(i+1,1)-CourbDiedre->pts->Element(i-1,1),
			CourbDiedre->pts->Element(i+1,0)-CourbDiedre->pts->Element(i-1,0));
			//( atan2(CourbDiedre->pts->Element(i+1,1)-CourbDiedre->pts->Element(i,1),
			//CourbDiedre->pts->Element(i+1,0)-CourbDiedre->pts->Element(i,0))
			//+ atan2(CourbDiedre->pts->Element(i,1)-CourbDiedre->pts->t[i-1][1],
			//CourbDiedre->pts->Element(i,0)-CourbDiedre->pts->t[i-1][0]) )/2.0f;
		
		//vrillage / horizontale en radian
		F->m_pProfils[i]->m_fWash = CourbVrillage->pts->Element(i,1)*DEG2RAD;

		//morphing / nervure centrale
		F->m_pProfils[i]->m_fMorph = CourbMorphing->pts->Element(i,1);

		//pos relative ligne A, B ,C ,D
		F->m_pProfils[i]->m_fPosA = (CourbNez->pts->Element(i,1)-CourbA->pts->Element(i,1))*100.0f/L;
		F->m_pProfils[i]->m_fPosB = (CourbNez->pts->Element(i,1)-CourbB->pts->Element(i,1))*100.0f/L;
		F->m_pProfils[i]->m_fPosC = (CourbNez->pts->Element(i,1)-CourbC->pts->Element(i,1))*100.0f/L;
		F->m_pProfils[i]->m_fPosD = (CourbNez->pts->Element(i,1)-CourbD->pts->Element(i,1))*100.0f/L;
		F->m_pProfils[i]->m_fPosE = (CourbNez->pts->Element(i,1)-CourbE->pts->Element(i,1))*100.0f/L;
	}
}

/************************/
/* ChargerFichierConfig */
/***********************$*/
void ChargerFichierConfig(char *NomFicConfig)
{
	int ws,hs;
	float xf,yf,wf,hf;
	int i;
	FILE *fid;
	char texte[50];
	int f[6];
	/*recupere taille de l'ecran*/
	ws = glutGet( GLUT_SCREEN_WIDTH );
	hs = glutGet( GLUT_SCREEN_HEIGHT );
	/*init tableau des no de fenetre*/
	f[0]=FenetreForme; f[1]=FenetreDiedre; f[2]=FenetreMorphing;
	f[3]=FenetreVrillage; f[4]=FenetreEpaiRel; f[5]=Fenetre3d;
	/*ouverture fichier en lecture*/
	printf("\nLecture fichier de Config: '%s' ", NomFicConfig);
	if( (fid = fopen( NomFicConfig, "rt" )) == NULL )
	{
      printf( "\nErreur ouverture fichier '%s'", NomFicConfig);
	  exit(0);
	}
	/*lecture nom du fichier de forme*/
	fscanf(fid,"%s %s",texte, NomFichierForme);
	/*boucle de lecture des positions/tailles normalis�es des fenetres*/
	for (i=0; i<6; i++)
	{
		glutSetWindow( f[i] );
		fscanf(fid,"%s %f %f %f %f\n",texte, &xf, &yf, &wf, &hf);
        glutPositionWindow((int)(xf*(double)ws), (int)(yf*(double)hs));
		glutReshapeWindow((int)(wf*(double)ws), (int)(hf*(double)hs));
	}
	/*fermeture fichier*/
	if(fclose(fid))
	{
		printf("\nProbleme � la fermeture du fichier");
		exit(0);
	}
	//lecture fichier de forme 
    F = LectureFichierForme(NomFichierForme);
	//met a jour matrice des points de controle
	//met a jour matrice des points de controle
	CtrlNez->pts = F->mCtrlNez;
	CtrlFui->pts = F->mCtrlFui;
	CtrlA->pts = F->mCtrlA;
	CtrlB->pts = F->mCtrlB;
	CtrlC->pts = F->mCtrlC;
	CtrlD->pts = F->mCtrlD;
	CtrlE->pts = F->mCtrlE;
	CtrlDiedre->pts = F->mCtrlDiedre;
	CtrlMorphing->pts = F->mCtrlMorphing;
	CtrlVrillage->pts = F->mCtrlVrillage;
	CtrlEpaiRel->pts = F->mCtrlEpaiRel;
	loadConfig = true;
	//maj interface GLUI
	CoeffProgGeom = F->CoeffProgGeom;
        CoeffExp = F->CoeffExp;
	EpaiRelCent = F->EpaiRelCent;
	NbCaiss = F->NbCaiss;
	if(NbCaiss%2==0)
		NbNerv=NbCaiss/2+1;
	else
		NbNerv=(NbCaiss+1)/2;
	VrillageCent = CtrlVrillage->pts->Element(0,1);
	VrillageBout = CtrlVrillage->pts->Element(3,1);
	//chargement profils
	LectureFichierProfil(F->m_strNomProfilCent.c_str(), &ExtProfCent, &IntProfCent);
	LectureFichierProfil(F->m_strNomProfilBout.c_str(), &ExtProfBout, &IntProfBout);
	EpaiRelProfCent=EpaisseurRelative(ExtProfCent, IntProfCent);
	EpaiRelProfBout=EpaisseurRelative(ExtProfBout, IntProfBout);
	InterpoleProfilBout(&ExtProfBout, ExtProfCent);
	InterpoleProfilBout(&IntProfBout, IntProfCent);
}

/***********************/
/* SauverFichierConfig */
/***********************/

void SauverFichierConfig(char *NomFicConfig)
{
	int x,y,w,h,ws,hs,i,f[6];
	FILE *fid;
	char* texte[6]=
	{"FENETRE_FORME","FENETRE_DIEDRE","FENETRE_MORPHING",
	"FENETRE_VRILLAGE","FENETRE_EPAISSEUR_RELATIVE","FENETRE_3D"};
	/*init tableau des no de fenetre*/
	f[0]=FenetreForme; f[1]=FenetreDiedre; f[2]=FenetreMorphing;
	f[3]=FenetreVrillage; f[4]=FenetreEpaiRel; f[5]=Fenetre3d;

	/*recupere taille de l'ecran*/
	ws = glutGet( GLUT_SCREEN_WIDTH );
	hs = glutGet( GLUT_SCREEN_HEIGHT );

	/*ouverture fichier en ecriture*/
	if( (fid = fopen( NomFicConfig, "wt" )) == NULL )
	{
      printf( "\nErreur ouverture fichier '%s'", NomFicConfig );
	  exit(0);
	}

        /*ecriture fichier de forme*/
	fprintf(fid,"FORME %s",	NomFichierForme);
	/*boucle ecriture des positions/tailles des fenetres normalis�s*/
	for (i=0; i<6; i++)
	{
		glutSetWindow( f[i] );
		x = glutGet( GLUT_WINDOW_X );
		y = glutGet( GLUT_WINDOW_Y );
		w = glutGet( GLUT_WINDOW_WIDTH );
		h = glutGet( GLUT_WINDOW_HEIGHT );
		fprintf(fid,"\n%s %1.4f %1.4f %1.4f %1.4f",	texte[i],
			(double)x/(double)ws, y/(double)hs, (double)w/(double)ws, (double)h/(double)hs);
	}
	/*fermeture fichier*/
	if(fclose(fid))
	{
		printf("\nProbleme � la fermeture du fichier");
		exit(0);
	}
}

/************************/
/* ChargerFichierProfil */
/************************/

void ChargerFichierProfil( int control )
{
	CString NomFichier;
	char* PtrNomFichier;
	//ouverture boite de dialogue
	CFileDialog DlgOpen(TRUE, NULL, "*.*", OFN_OVERWRITEPROMPT, NULL, NULL);
	if (DlgOpen.DoModal()==IDOK)
	{
        	NomFichier=DlgOpen.GetPathName();
		PtrNomFichier = NomFichier.GetBuffer(1);
		//test profil centre ou bout ?
		if(control==PROF_CENT)
		{
			LectureFichierProfil(PtrNomFichier, &ExtProfCent, &IntProfCent);
			EpaiRelProfCent=EpaisseurRelative(ExtProfCent, IntProfCent);
			ProfCent->set_text(PtrNomFichier);
			F->m_strNomProfilCent = PtrNomFichier;
		}
		else
		{
			LectureFichierProfil(PtrNomFichier, &ExtProfBout, &IntProfBout);
			EpaiRelProfBout=EpaisseurRelative(ExtProfBout, IntProfBout);
			ProfBout->set_text(PtrNomFichier);
			F->m_strNomProfilBout = PtrNomFichier;
		}
		InterpoleProfilBout(&ExtProfBout, ExtProfCent);
		InterpoleProfilBout(&IntProfBout, IntProfCent);
	}
	/*mise a jour axes*/
	VisuTousLesAxes=1; display();
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
		oldF = F;
		F = LectureFichierForme(NomFichierForme);
		//met a jour matrice des points de controle
		CtrlNez->pts = F->mCtrlNez;
		CtrlFui->pts = F->mCtrlFui;
		CtrlA->pts = F->mCtrlA;
		CtrlB->pts = F->mCtrlB;
		CtrlC->pts = F->mCtrlC;
		CtrlD->pts = F->mCtrlD;
		CtrlE->pts = F->mCtrlE;
		CtrlDiedre->pts = F->mCtrlDiedre;
		CtrlMorphing->pts = F->mCtrlMorphing;
		CtrlVrillage->pts = F->mCtrlVrillage;
		CtrlEpaiRel->pts = F->mCtrlEpaiRel;
		//libere memoire ancienne forme
		delete(oldF);
		oldF = NULL;
		//maj interface GLUI
		FicForme->set_text(NomFichierForme);
		ProfCent->set_text(F->m_strNomProfilCent.c_str());
		ProfBout->set_text(F->m_strNomProfilBout.c_str());
		//maj parametres (live-variables)
		CoeffProgGeom = F->CoeffProgGeom;
                CoeffExp = F->CoeffExp;
		EpaiRelCent = F->EpaiRelCent;
		NbCaiss = F->NbCaiss;
		VrillageCent = CtrlVrillage->pts->Element(0,1);
		VrillageBout = CtrlVrillage->pts->Element(3,1);
		//chargement profils !!!
		LectureFichierProfil(F->m_strNomProfilCent.c_str(), &ExtProfCent, &IntProfCent);
		LectureFichierProfil(F->m_strNomProfilBout.c_str(), &ExtProfBout, &IntProfBout);
		EpaiRelProfCent=EpaisseurRelative(ExtProfCent, IntProfCent);
		EpaiRelProfBout=EpaisseurRelative(ExtProfBout, IntProfBout);
		InterpoleProfilBout(&ExtProfBout, ExtProfCent);
		InterpoleProfilBout(&IntProfBout, IntProfCent);
	}
	/*synchronisation variables GLUI*/
	glui->sync_live();
	/*mise a jour axes*/
	VisuTousLesAxes=1; display();
}


/***********************/
/* ChargerFichierForme2 */
/***********************/

void ChargerFichierForme2( int /*control*/ )
{
	printf ("\nChargerFichierForme2\n");
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
		oldF = F;
		F = LectureFichierForme2(NomFichierForme);
		wasLoad = true;
		//met a jour matrice des points de controle
		CtrlNez->pts = F->mCtrlNez;
		CtrlFui->pts = F->mCtrlFui;
		CtrlA->pts = F->mCtrlA;
		CtrlB->pts = F->mCtrlB;
		CtrlC->pts = F->mCtrlC;
		CtrlD->pts = F->mCtrlD;
		CtrlE->pts = F->mCtrlE;
		CtrlDiedre->pts = F->mCtrlDiedre;
		CtrlMorphing->pts = F->mCtrlMorphing;
		CtrlVrillage->pts = F->mCtrlVrillage;
		CtrlEpaiRel->pts = F->mCtrlEpaiRel;
		CourbNez->pts = F->mCourbNez;
		CourbFui->pts = F->mCourbFui;
		CourbA->pts = F->mCourbA;
		CourbB->pts = F->mCourbB;
		CourbC->pts = F->mCourbC;
		CourbD->pts = F->mCourbD;
		CourbE->pts = F->mCourbE;
		CourbDiedre->pts = F->mCourbDiedre;
		CourbMorphing->pts = F->mCourbMorphing;
		CourbVrillage->pts = F->mCourbVrillage;
		CourbEpaiRel->pts = F->mCourbEpaiRel;
		//libere memoire ancienne forme
		delete(oldF);
		oldF = NULL;
		//maj interface GLUI
		FicForme->set_text(NomFichierForme);
		ProfCent->set_text(F->m_strNomProfilCent.c_str());
		ProfBout->set_text(F->m_strNomProfilBout.c_str());
		//maj parametres (live-variables)
		CoeffProgGeom = F->CoeffProgGeom;
                CoeffExp = F->CoeffExp;
		EpaiRelCent = F->EpaiRelCent;
		NbCaiss = F->NbCaiss;
		VrillageCent = CtrlVrillage->pts->Element(0,1);
		VrillageBout = CtrlVrillage->pts->Element(3,1);
		//chargement profils !!!
		LectureFichierProfil(F->m_strNomProfilCent.c_str(), &ExtProfCent, &IntProfCent);
		LectureFichierProfil(F->m_strNomProfilBout.c_str(), &ExtProfBout, &IntProfBout);
		EpaiRelProfCent=EpaisseurRelative(ExtProfCent, IntProfCent);
		EpaiRelProfBout=EpaisseurRelative(ExtProfBout, IntProfBout);
		InterpoleProfilBout(&ExtProfBout, ExtProfCent);
		InterpoleProfilBout(&IntProfBout, IntProfCent);
	}
	/*synchronisation variables GLUI*/
	glui->sync_live();
	/*mise a jour axes*/
	VisuTousLesAxes=1; 
	display();
}

void ExportTemplateSling( int control )
{
	Courbe *CtrlCour=NULL, *CourbCour=NULL; 
	CString NomFichier;
	LPTSTR PtrNomFichier;
	printf("\nExport template %d", control);
	CFileDialog DlgOpen(FALSE, NULL, "*.txt", OFN_OVERWRITEPROMPT, NULL, NULL);
	if (DlgOpen.DoModal()==IDOK)
	{
		NomFichier=DlgOpen.GetPathName();
		PtrNomFichier = NomFichier.GetBuffer(1);
        	FILE *fid;
		if( (fid = fopen( NomFichier, "wt" )) == NULL )
		{
			printf( "\nError open '%s'", NomFichier);
			exit(0);
		}
		for (int _i = 0; _i < 7; _i++) {
			switch( _i )
			{
        			case 0 : CourbCour = CourbNez;
					break;
				case 1 : CourbCour = CourbA;
					break;
				case 2 : CourbCour = CourbB;
					break;
				case 3 : CourbCour = CourbC;
        				break;
        			case 4 : CourbCour = CourbD;
					break;
				case 5 : CourbCour = CourbE;
					break;
				case 6 : CourbCour = CourbFui;
					break;
			}
			fprintf(fid, "%s %d\n", CourbCour ->name.c_str(), NbNerv);
			printf("\n%s %d", CourbCour ->name.c_str(), NbNerv);
			for (int i = 0; i < NbNerv; i++) {
				fprintf(fid,"%f %f\n", CourbCour->pts->Element(i,0), CourbCour->pts->Element(i,1));
				printf("\n%f %f", CourbCour->pts->Element(i,0), CourbCour->pts->Element(i,1));
			}
		}
		if(fclose(fid))
		{
			printf("\nError close %s", NomFichier);
			exit(0);
		}
	}
}

void LoadCourbeSling( int control )
{
	Courbe *CourbCour=NULL; 
	CString NomProf;
	//char* PtrNomFichier;
	CFileDialog DlgOpen(TRUE, NULL, "*.*", OFN_OVERWRITEPROMPT, NULL, NULL);
	if (DlgOpen.DoModal()==IDOK)
	{
		NomProf=DlgOpen.GetPathName();
		//PtrNomFichier = NomFichier.GetBuffer(1);
		FILE *fid;
		if( (fid = fopen( NomProf, "rt" )) == NULL )
		{
			printf( "\nErreur ouverture fichier: '%s'", NomProf );
			//fflush(stdin); getch();
		}
		else
		{
			for (int _i = 0; _i < 7; _i++)
		
			{

				switch( _i )

				{

					case 0 : CourbCour = CourbNez;
						break;

					case 1 : CourbCour = CourbA;
						break;

					case 2 : CourbCour = CourbB;
						break;

					case 3 : CourbCour = CourbC;
						break;

					case 4 : CourbCour = CourbD;
						break;

					case 5 : CourbCour = CourbE;
						break;

					case 6 : CourbCour = CourbFui;
						break;

				}

				char temp[1024];
				int dNbNerv;
				fscanf(fid,"%s %d", temp, &dNbNerv); 

				printf ("\nLoading %s %d(%d)", temp, dNbNerv, NbNerv);

				for(int i = 0; i < NbNerv; i++)
				{
					double f1, f2;
					fscanf( fid, "%f %f", &f1, &f2 );
					printf ("\nLoading %f %f", f1, f2);
					CourbCour -> pts -> SetElement(i,0,f1);
					CourbCour -> pts -> SetElement(i,1,f2);

				}

			}

			if(fclose(fid))

			{

				printf("\nProblem close '%s'", NomProf);

				fflush(stdin); //getch();

			}


		}

	}

	wasLoad = true;

	display();

}

void ExportTemplate( int control )
{

	Courbe *CtrlCour=NULL, *CourbCour=NULL; 

	CString NomFichier;

	LPTSTR PtrNomFichier;


	switch( control )
	{
		case BTN_DIEDRE : CourbCour = CourbDiedre;
			break;

		case BTN_MORPHING : CourbCour = CourbMorphing;
			break;

		case BTN_VRILLAGE : CourbCour = CourbVrillage;
			break;

		case BTN_EPAIREL : CourbCour = CourbEpaiRel;
			break;
	}

	printf("\nExport template %d", control);

	CFileDialog DlgOpen(FALSE, NULL, "*.txt", OFN_OVERWRITEPROMPT, NULL, NULL);

	if (DlgOpen.DoModal()==IDOK)

	{

		NomFichier=DlgOpen.GetPathName();

		PtrNomFichier = NomFichier.GetBuffer(1);


		FILE *fid;

		if( (fid = fopen( NomFichier, "wt" )) == NULL )

		{

			printf( "\nError open '%s'", NomFichier);

			exit(0);

		}

		fprintf(fid, "%s %d", CourbCour ->name.c_str(), NbNerv);
	
		printf("\n%s %d", CourbCour ->name.c_str(), NbNerv);

		for (int i = 0; i < NbNerv; i++) {

			fprintf(fid,"\n%f %f", CourbCour->pts->Element(i,0), CourbCour->pts->Element(i,1));

			printf("\n%f %f", CourbCour->pts->Element(i,0), CourbCour->pts->Element(i,1));
		}
	
		if(fclose(fid))
		{

			printf("\nError close %s", NomFichier);

			exit(0);

		}



	}

}

void LoadCourbe( int control )
{

	Courbe *CtrlCour=NULL, *CourbCour=NULL; 

	switch( control )
	{
		case BTN_DIEDRE : CourbCour = CourbDiedre;
			break;

		case BTN_MORPHING : CourbCour = CourbMorphing;
			break;

		case BTN_VRILLAGE : CourbCour = CourbVrillage;
			break;

		case BTN_EPAIREL : CourbCour = CourbEpaiRel;
			break;
	}


	CString NomProf;

	//char* PtrNomFichier;

	CFileDialog DlgOpen(TRUE, NULL, "*.*", OFN_OVERWRITEPROMPT, NULL, NULL);

	if (DlgOpen.DoModal()==IDOK)

	{

		NomProf=DlgOpen.GetPathName();

		//PtrNomFichier = NomFichier.GetBuffer(1);

		FILE *fid;

		if( (fid = fopen( NomProf, "rt" )) == NULL )

		{

			printf( "\nErreur ouverture fichier: '%s'", NomProf );

			//fflush(stdin); getch();

		}

		else

		{
				
				char temp[1024];

				int dNbNerv;

				fscanf(fid,"%s %d", temp, &dNbNerv); 

				printf ("\nLoad %s %d", temp, NbNerv);

				for(int i=0; i<NbNerv; i++)
				{
					double f1, f2;
					fscanf( fid, "%f %f", &f1, &f2 );
					printf( "\n%f %f", f1, f2 );
					CourbCour -> pts -> SetElement(i,0,f1);
					CourbCour -> pts -> SetElement(i,1,f2);

				}
	
				if(fclose(fid))

				{

					printf("\nProbleme � la fermeture du fichier '%s'", NomProf);

					fflush(stdin); //getch();

				}


		}

	}

	wasLoad = true;

	display();
}

/**********************/

/* SauverFichierForme */

/**********************/

void SauverFichierForme( int /*control*/ )
{
	CString NomFichier;
	LPTSTR PtrNomFichier;
	//ouverture boite de dialogue
	CFileDialog DlgOpen(FALSE, NULL, "*.txt", OFN_OVERWRITEPROMPT, NULL, NULL);
	if (DlgOpen.DoModal()==IDOK)
	{
		//recupere nom de fichier
		NomFichier=DlgOpen.GetPathName();
		PtrNomFichier = NomFichier.GetBuffer(1);
		//met a jour parametre de forme
		F->CoeffProgGeom = CoeffProgGeom;
        F->CoeffExp = CoeffExp;
		F->EpaiRelCent = EpaiRelCent;
		F->NbCaiss = NbCaiss;
		//ecriture fichier
		EcritureFichierForme(PtrNomFichier, F);
		//maj nom affich� dans l'interface
		strcpy(NomFichierForme, PtrNomFichier);
		FicForme->set_text(NomFichierForme);
	}
}

/**********************/
/* SauverFichierForme2 */
/**********************/

void SauverFichierForme2( int /*control*/ )
{
	printf ("\nSauverFichierForme2\n");
	CString NomFichier;
	LPTSTR PtrNomFichier;
	//ouverture boite de dialogue
	CFileDialog DlgOpen(FALSE, NULL, "*.txt", OFN_OVERWRITEPROMPT, NULL, NULL);
	if (DlgOpen.DoModal()==IDOK)
	{
		//recupere nom de fichier
		NomFichier=DlgOpen.GetPathName();
		PtrNomFichier = NomFichier.GetBuffer(1);
			F->mCourbA = CourbA->pts;
			F->mCourbB = CourbB->pts;
			F->mCourbC = CourbC->pts;
			F->mCourbD = CourbD->pts;
			F->mCourbE = CourbE->pts;
			F->mCourbDiedre = CourbDiedre->pts;
			F->mCourbEpaiRel = CourbEpaiRel->pts;
			F->mCourbMorphing = CourbMorphing->pts;
			F->mCourbFui = CourbFui->pts;
			F->mCourbNez = CourbNez->pts;
			F->mCourbVrillage = CourbVrillage->pts;
		//met a jour parametre de forme
		F->CoeffProgGeom = CoeffProgGeom;
                F->CoeffExp = CoeffExp;
		F->EpaiRelCent = EpaiRelCent;
		F->NbCaiss = NbCaiss;
		//ecriture fichier
		EcritureFichierForme2(PtrNomFichier, F);
		//maj nom affich� dans l'interface
		strcpy(NomFichierForme, PtrNomFichier);
		FicForme->set_text(NomFichierForme);
	}
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

/********************/
/* SauverFichierFGen */
/********************/

void SauverFichierFGen( int /*control*/ )
{
	CString NomFichier;
	LPTSTR PtrNomFichier;
	//ouverture boite de dialogue
	CFileDialog DlgOpen(FALSE, NULL, "*.fgen", OFN_OVERWRITEPROMPT, NULL, NULL);
	if (DlgOpen.DoModal()==IDOK)
	{
		//recupere nom de fichier
		NomFichier=DlgOpen.GetPathName();
		PtrNomFichier = NomFichier.GetBuffer(1);
		//ecriture fichier FGen
		EcritureFichierFGen(PtrNomFichier, F);
	}
}

/*****************/
/* ModifProgGeom */
/*****************/
void ModifProgGeom( int /*control*/ )
{
	/*mise a jour axes*/
	VisuTousLesAxes=1; display();
}

/*****************/
/* ModifVrillage */
/*****************/

void ModifVrillage( int /*control*/ )
{
	/*modif vrillage*/
	CtrlVrillage->pts->SetElement(0,1,VrillageCent);
	CtrlVrillage->pts->SetElement(3,1,VrillageBout);
	/*mise a jour axes*/
	VisuTousLesAxes=1; display();
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
/* ModifRaskladSymetrique */
/***********************/
void ModifRaskladSymetrique( int /*control*/ )
{
    //
    //
    //
}


/***********************/
/* ModifVisuSymetrique */
/***********************/
void ModifVisuSymetrique( int /*control*/ )
{
	/*recalcul et retrace la vue 3d*/
	CalculVue3d(); VisuAxe(Axe3d); glutSwapBuffers();
}


/***********************/
/* ModifVisuPtsSuspentes */
/***********************/
void ModifVisuPtsSuspentes( int /*control*/ )
{
	bool value = (VisuPtsSuspentes == 1);

	CtrlA->visible = value;
	CtrlB->visible = value;
	CtrlC->visible = value;
	CtrlD->visible = value;
	CtrlE->visible = value;

	CourbA->visible = value;
	CourbB->visible = value;
	CourbC->visible = value;
	CourbD->visible = value;
	CourbE->visible = value;

	/*recalcul et retrace la vue 3d*/
	CalculVue3d(); VisuAxe(Axe3d); glutSwapBuffers();
}

/***********************/
/* ModifVisuPtsSuspentes */
/***********************/
void ModifAxesAuto( int /*control*/ )
{
	AxeForme->XAuto = AxesAuto;
	AxeForme->YAuto = AxesAuto;
	AxeForme->ZAuto = AxesAuto;

	AxeDiedre->XAuto = AxesAuto;
	AxeDiedre->YAuto = AxesAuto;
	AxeDiedre->ZAuto = AxesAuto;

	AxeMorphing->XAuto = AxesAuto;
	AxeMorphing->YAuto = AxesAuto;
	AxeMorphing->ZAuto = AxesAuto;

	AxeVrillage->XAuto = AxesAuto;
	AxeVrillage->YAuto = AxesAuto;
	AxeVrillage->ZAuto = AxesAuto;

	AxeEpaiRel->XAuto = AxesAuto;
	AxeEpaiRel->YAuto = AxesAuto;
	AxeEpaiRel->ZAuto = AxesAuto;

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
		/*sauvegarde derni�re config*/

		//SauverFichierConfig("config.txt");

		/*libere m�moire axe et courbes associ�es*/

		LibererAxe(AxeForme);

		LibererAxe(AxeDiedre);	

		LibererAxe(AxeMorphing);

		LibererAxe(AxeVrillage);

		LibererAxe(Axe3d);		

		LibererAxe(AxeEpaiRel);

		/*libere m�moire forme*/

		/*attention les pts de controle ont deja �t� lib�r� !!!*/

		/* -> on ne peut donc pas utiliser LibererForme...*/

		delete(F->m_pProfils);
		F->m_pProfils = NULL;
		/*libere m�moire diverses matrices*/
		delete(IntProfCent);
		delete(ExtProfCent);
		delete(IntProfBout);
		delete(ExtProfBout);
		delete(XNerv);
		/*malgres tout �a le debugger dit qu'il y a encore
		plein de truc qui traine en m�moire !? */
		/*sortie*/
		exit(0);
	}
}



/********************/
/* ModifNbrAlveoles */
/********************/

void ModifNbrAlveoles( int /*control*/ )

{
	int AncienNbNerv;
	AncienNbNerv = NbNerv;
	/*calcul du nouveau nombre de nervures*/
	if(NbCaiss%2==0) 
		NbNerv=NbCaiss/2+1;
	else
		NbNerv=(NbCaiss+1)/2;
	/*modif points de controle des axes morphing, vrillage, epai-rel*/
	CtrlMorphing->pts->SetElement(3,0, NbNerv-1); 
	CtrlMorphing->pts->AddElement(2,0, NbNerv - AncienNbNerv); 
	CtrlVrillage->pts->SetElement(3,0, NbNerv-1); 
	CtrlVrillage->pts->AddElement(2,0, NbNerv - AncienNbNerv); 
	CtrlEpaiRel->pts->SetElement(3,0, NbNerv-1); 
	CtrlEpaiRel->pts->AddElement(2,0, NbNerv - AncienNbNerv); 
	/*mise a jour tous les axes*/
	VisuTousLesAxes=1; display();
}



/***************/

/* ModifTaille */

/***************/

void ModifTaille( int control )

{

	double coef;

	int l,c;

	/*test decremente ou incremente*/

	switch( control)
	{
		case TAILLE_MOINS_5 : coef=0.95f;
			break;

		case TAILLE_MOINS_1 : coef=0.99f;
			break;

		case TAILLE_PLUS_5 : coef=1.05f;
			break;

		case TAILLE_PLUS_1 : coef=1.01f;
			break;

		default : coef = 0.0f;
	}

	/*modif des points de controle*/

	for(l=0; l<CtrlNez->pts->GetLignes(); l++)

	{

		for(c=0; c<CtrlNez->pts->GetColonnes(); c++)

		{

			CtrlNez->pts->MultiplyElement(l,c,  coef);

			CtrlFui->pts->MultiplyElement(l,c,  coef);

			CtrlA->pts->MultiplyElement(l,c,  coef);

			CtrlB->pts->MultiplyElement(l,c,  coef);

			CtrlC->pts->MultiplyElement(l,c,  coef);

			CtrlD->pts->MultiplyElement(l,c,  coef);

			CtrlE->pts->MultiplyElement(l,c,  coef);

		}

	}

	/*mise a jour axes*/

	VisuTousLesAxes=1; display();

}





/***************/

/* CalculVue3d */

/***************/

void CalculVue3d(void)
{
	Matrice *XExt,*YExt,*ZExt; //coordonn�es 3D extrados
	Matrice *XInt,*YInt,*ZInt; //coordonn�es 3D intrados
	Matrice *PtsSuspentes; //pour recuperer la position 3D des pts de suspentage
	/*par defaut destruction des courbes ou mesh de l'axe 3D*/
	LibererCourbesAxe(Axe3d); LibererMeshsAxe(Axe3d);
	/*calcul forme*/
	CalculTableauForme();
	CalculForme3D(F, 0, 0.0f,
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
	/*maj info de la forme*/
	CalculInfoForme(F, &info);
	MajInfo();
}

/***************/
/* CalculXNerv */
/***************/
void CalculXNerv(void)
{
	double CumulProg,MultG,LargCent;
	int i;

	/*calcul du nombre de nervures en fonctions du nombre de caissons*/
	/* et mise a jour position relative 1ere nervure*/
	if(NbCaiss%2==0) 
	{
		NbNerv=NbCaiss/2+1;
		CumulProg=0.0f;
	}
	else
	{
		NbNerv=(NbCaiss+1)/2;
		CumulProg=0.5f;
	}
	/*creation matrice XNerv*/
	delete(XNerv);
	XNerv = new Matrice(NbNerv, 1);
	MultG=1.0f;
	double MultE=0.0f;
	XNerv->SetElement(0,0, CumulProg);
	for(i=1; i<NbNerv; i++)
	{
		//
		//CumulProg += Mult;
		MultG *= CoeffProgGeom;
		MultE = pow(2.718281828f,-CoeffExp*(i-1));
		CumulProg += MultG*MultE; 
		XNerv->SetElement(i,0, CumulProg);
	}
	/*calcul largeur caisson central*/
	LargCent=CtrlNez->pts->Element(3,0)/CumulProg;
	/*calcul position des nervures*/
	for(i=0; i<NbNerv; i++)	
		XNerv->MultiplyElement(i,0, LargCent);
}

/********************/
/* CalculCourbesCtrl */
/********************/

void CalculCourbesCtrl()
{
	Courbe *CtrlCour=NULL, *CourbCour=NULL; 
	Matrice *Bez, *XBez, *YBez;	/*courbe de bezier*/
	Matrice *DevBez;
	Matrice *XDiedre, *YDiedre;
	Matrice *Xi, *Yi;
	double Ratio;
	int i,c;
	/************************/
	/*courbes de l'axe forme*/
	/************************/
	for(c=0; c<7; c++)
	{
		/*affectation ptr Courb et Ctrl courant*/
		switch(c)
		{
		case 0: CtrlCour = CtrlNez; CourbCour = CourbNez; break;
		case 1: CtrlCour = CtrlFui; CourbCour = CourbFui; break;
		case 2: CtrlCour = CtrlA; CourbCour = CourbA; break;
		case 3: CtrlCour = CtrlB; CourbCour = CourbB; break;
		case 4: CtrlCour = CtrlC; CourbCour = CourbC; break;
		case 5: CtrlCour = CtrlD; CourbCour = CourbD; break;
		case 6: CtrlCour = CtrlE; CourbCour = CourbE; break;
		default:;
		}
		/*calcul courde de bezier*/
		Bez = MonBezier(CtrlCour->pts, NBP_BEZ);
		delete(CourbCour->pts);
		CourbCour->pts = new Matrice(NbNerv,2);
		/*interpolation lin�aire*/
		for (i=0; i<NbNerv; i++)
		{
			CourbCour->pts->SetElement(i,0, XNerv->Element(i,0));
			CourbCour->pts->SetElement(i,1, InterpLinX(Bez, XNerv->Element(i,0)));
		}
		delete(Bez);
	}
	/*creation courbe pour visualiser les nervures*/
	delete(CourbNerv->pts);
	CourbNerv->pts =new Matrice(NbNerv*2,2);
	for(i=0; i<NbNerv; i++)
	{
		CourbNerv->pts->SetElement(2*i,0, CourbNez->pts->Element(i,0));
		CourbNerv->pts->SetElement(2*i,1, CourbNez->pts->Element(i,1));
		CourbNerv->pts->SetElement(2*i+1,0, CourbFui->pts->Element(i,0));
		CourbNerv->pts->SetElement(2*i+1,1, CourbFui->pts->Element(i,1));
	}
	/************************/
	/*courbe de l'axe diedre*/
	/************************/
	/*calcul de la courbe de bezier*/
	Bez = MonBezier(CtrlDiedre->pts, NBP_BEZ);
	/*calcul du developp� du diedre*/
	DevBez = Zeros(NBP_BEZ,1);
	for(i=0; i<NBP_BEZ-1; i++)
		DevBez->SetElement(i+1,0, DevBez->Element(i,0) + 
			sqrt(sqr(Bez->Element(i+1,0)-Bez->Element(i,0))+sqr(Bez->Element(i+1,1)-Bez->Element(i,1))));
	/*mise a l'echelle de la courbe de bezier du diedre 
	 et des points de controle */
	Ratio = CtrlNez->pts->Element(3,0) / DevBez->Element(NBP_BEZ-1,0);
	for (i=0; i<NBP_BEZ; i++)
        {
		Bez->MultiplyElement(i,0,Ratio); 
		Bez->MultiplyElement(i,1,Ratio); 
		DevBez->MultiplyElement(i,0,Ratio);
	}
	for (i=0; i<4; i++)
	{
		CtrlDiedre->pts->MultiplyElement(i,0,Ratio); 
		CtrlDiedre->pts->MultiplyElement(i,1,Ratio);
	} 
	/*interpolation des positions de nervure sur le diedre*/
	XDiedre=Zeros(NbNerv,1); YDiedre=Zeros(NbNerv,1);
	XBez=Zeros(NBP_BEZ,1); YBez=Zeros(NBP_BEZ,1);
	for(i=0; i<NBP_BEZ; i++)
        {
		XBez->SetElement(i,0,Bez->Element(i,0)); 
		YBez->SetElement(i,0,Bez->Element(i,1));
	}
	InterpLinMat(DevBez, XBez, XNerv, XDiedre); 
	InterpLinMat(DevBez, YBez, XNerv, YDiedre);
	/*affectation courbe diedre*/
	delete(CourbDiedre->pts);
	CourbDiedre->pts = new Matrice(NbNerv,2);
	for(i=0; i<NbNerv; i++)
	{
		CourbDiedre->pts->SetElement(i,0,XDiedre->Element(i,0)); 
		CourbDiedre->pts->SetElement(i,1,YDiedre->Element(i,0));
	}
	/*destruction variable intermedaire*/
	delete(XDiedre); delete(YDiedre);
	delete(Bez); delete(XBez); delete(YBez);
	/*************************************/
	/* courbe morphing/vrillage/epai-rel */
	/*************************************/
	for(c=0; c<3; c++)
	{
		/*affectation ptr Courb et Ctrl courant*/
		switch(c)
		{
		case 0: CtrlCour = CtrlMorphing; CourbCour = CourbMorphing; break;
		case 1: CtrlCour = CtrlVrillage; CourbCour = CourbVrillage; break;
		case 2: CtrlCour = CtrlEpaiRel; CourbCour = CourbEpaiRel; break;
		default:;
		}
		/*calcul de la courbe de bezier*/
		Bez = MonBezier(CtrlCour->pts, NBP_BEZ);
		/*interpolation des Y correspondants aux No de nervure*/
		Xi=Zeros(NbNerv,1); Yi=Zeros(NbNerv,1);
		XBez=Zeros(NBP_BEZ,1); YBez=Zeros(NBP_BEZ,1);
		for(i=0; i<NBP_BEZ; i++)
		{
			XBez->SetElement(i,0,Bez->Element(i,0)); 
			YBez->SetElement(i,0,Bez->Element(i,1));
		}
		for(i=0; i<NbNerv; i++) 
			Xi->SetElement(i,0,i);
		InterpLinMat(XBez, YBez, Xi, Yi); 
		/*affectation courbe morphing*/
		delete(CourbCour->pts);
		//CourbMorphing->pts = Bez;
		CourbCour->pts = new Matrice(NbNerv,2);
		for(i=0; i<NbNerv-1; i++)
		{
			CourbCour->pts->SetElement(i,0,Xi->Element(i,0));
			CourbCour->pts->SetElement(i,1,Yi->Element(i,0));
		}
		/*correction en extremite*/
		CourbCour->pts->SetElement(NbNerv-1,0,CtrlCour->pts->Element(3,0));
		CourbCour->pts->SetElement(NbNerv-1,1,CtrlCour->pts->Element(3,1));
		/*destruction variable intermedaire*/
		delete(Xi); delete(Yi);	delete(Bez);
		delete(XBez); delete(YBez);
	}
	if (loadConfig) {
		loadConfig = false;
	}
}
/****************/
/* display      */
/****************/

void display(void)
{
		CalculXNerv();
		//if (!wasLoad) 
		CalculCourbesCtrl();
		//if (AxeSel == Axe3d) CalculVue3d();
		/*visualise tous les axes*/
		if(VisuTousLesAxes == 1)
		{
			CalculVue3d();
			VisuAxe(AxeForme); glutSwapBuffers();
			VisuAxe(AxeDiedre);	glutSwapBuffers();
			VisuAxe(AxeMorphing); glutSwapBuffers();
			VisuAxe(AxeVrillage); glutSwapBuffers();
			VisuAxe(AxeEpaiRel); glutSwapBuffers();
			VisuAxe(Axe3d);	glutSwapBuffers();
			VisuTousLesAxes = 0;
		}
		else /*visualise uniquement l'axe s�lectionn�*/
		{
			if (AxeSel != NULL)	
				VisuAxe(AxeSel); glutSwapBuffers();
		}
//	}
	//printf("\nerreur: %d",glGetError());
}


/*****************/
/* reshape       */
/*****************/
void reshape(int w, int h)
{
	glViewport(0, 0, w, h);
	VisuTousLesAxes=1;
	display();
	//glutPostRedisplay();
}

/*****************/
/* motion        */
/*****************/
void motion(int x, int y)
{
	double xs,ys; /*coordonn�es souris dans l'axe*/
	double xmin,xmax,ymin,ymax; /*min max de l'axe selectionn�*/
	int h,w; /*taille de la fenetre*/
	int i,j; /*indice pour centre ou extremite*/
	double RatioY, RatioX, DeltaX, DeltaY, OldY, DeltaYNezFui, YNezOuFui;
	int n; /*indice de boucle*/
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
		Axe3d->xcam = Axe3d->xcam + (double)(x - xSouris)/100.0;
		Axe3d->ycam = Axe3d->ycam - (double)(y - ySouris)/100.0;
		xSouris = x; ySouris = y;
	}
	/*axe 3d -> zoom*/
	else if (ZoomTwist3d)
	{
		Axe3d->twist = Axe3d->twist - (double)(x - xSouris);
		Axe3d->zcam = Axe3d->zcam + (double)(y - ySouris)/100.0;
		xSouris = x; ySouris = y;
        }
	/*test d'un point selectionn�*/
	else if (PointSel!=NULL)
	{
		/*conversion en coordonn�es repere*/
		w = glutGet(GLUT_WINDOW_WIDTH);
		h = glutGet(GLUT_WINDOW_HEIGHT);
		xmin=AxeSel->xmin; xmax=AxeSel->xmax;
		ymin=AxeSel->ymin; ymax=AxeSel->ymax;
		xs=xmin+(double)x*(xmax-xmin)/(double)w;
		ys=ymin+((double)h-(double)y)*(ymax-ymin)/h;
		/*calcul deplacement en X et Y*/
		DeltaY = ys - PointSel[1]; DeltaX = xs - PointSel[0];
		/*message pour test*/
		//printf(",<%2.3f,%2.3f>",xs,ys);
		/*affectation position souris au point selectionn�*/
		/*et actions sur les autres points suivant le type de point...*/
		/* nez ou fuite, au centre ou en extremite*/
		if((PointSel == CtrlFui->pts->GetLigne(0))
        		||(PointSel == CtrlNez->pts->GetLigne(0))
			||(PointSel == CtrlFui->pts->GetLigne(3))
			||(PointSel == CtrlNez->pts->GetLigne(3)))
		{
			/*modification des Y*/
			/*maj des indices i et j en fonction de centre ou extremite*/
			/*utilisation de ces 2 indices pour �viter de r�p�ter 2x le code suivant*/
			if ((PointSel == CtrlFui->pts->GetLigne(0))
			||(PointSel == CtrlNez->pts->GetLigne(0)))
			{i = 0; j = 1;} /*centre*/
			else
			{i = 3; j = 2;} /*extremite*/
			/*calcul larg au centre*/
			DeltaYNezFui=(CtrlNez->pts->Element(i,1)-CtrlFui->pts->Element(i,1));
			if(DeltaYNezFui!=0.0f)
			{
				/*maj Nez Ou Fuite*/
				if(PointSel == CtrlFui->pts->GetLigne(i))
				{
        				YNezOuFui = CtrlNez->pts->Element(i,1);
					CtrlFui->pts->SetElement(j,1, CtrlFui->pts->Element(j,1) + (ys - CtrlFui->pts->Element(i,1)));
					CtrlFui->pts->SetElement(i,1, ys);
					RatioY = (YNezOuFui-ys)/DeltaYNezFui;
				}
				else
				{
					YNezOuFui = CtrlFui->pts->Element(i,1);
					CtrlNez->pts->SetElement(j,1, CtrlNez->pts->Element(j,1) + (ys - CtrlNez->pts->Element(i,1)));
					CtrlNez->pts->SetElement(i,1, ys);
					RatioY = -(YNezOuFui-ys)/DeltaYNezFui;
				}
				/*applique le meme ratio aux lignes de suspente*/
				/*ligne A*/
				OldY=CtrlA->pts->Element(i,1);
				CtrlA->pts->SetElement(i,1, YNezOuFui + (CtrlA->pts->Element(i,1)-YNezOuFui)*RatioY);
				CtrlA->pts->SetElement(j,1, CtrlA->pts->Element(j,1) + (CtrlA->pts->Element(i,1)-OldY));
				/*ligne B*/
				OldY=CtrlB->pts->Element(i,1);
				CtrlB->pts->SetElement(i,1, YNezOuFui + (CtrlB->pts->Element(i,1)-YNezOuFui)*RatioY);
				CtrlB->pts->SetElement(j,1, CtrlB->pts->Element(j,1) + (CtrlB->pts->Element(i,1)-OldY));
				/*ligne C*/
				OldY=CtrlC->pts->Element(i,1);
				CtrlC->pts->SetElement(i,1, YNezOuFui + (CtrlC->pts->Element(i,1)-YNezOuFui)*RatioY);
				CtrlC->pts->SetElement(j,1, CtrlC->pts->Element(j,1) + (CtrlC->pts->Element(i,1)-OldY));
				/*ligne D*/
				OldY=CtrlD->pts->Element(i,1);
				CtrlD->pts->SetElement(i,1, YNezOuFui + (CtrlD->pts->Element(i,1)-YNezOuFui)*RatioY);
				CtrlD->pts->SetElement(j,1, CtrlD->pts->Element(j,1) + (CtrlD->pts->Element(i,1)-OldY));
				/*ligne E*/
				OldY=CtrlE->pts->Element(i,1);
				CtrlE->pts->SetElement(i,1, YNezOuFui + (CtrlE->pts->Element(i,1)-YNezOuFui)*RatioY);
				CtrlE->pts->SetElement(j,1, CtrlE->pts->Element(j,1) + (CtrlE->pts->Element(i,1)-OldY));
			}
			/*modification des X*/
			if( ((PointSel == CtrlFui->pts->GetLigne(3))
			||(PointSel == CtrlNez->pts->GetLigne(3)))
			&&(xs != 0) ) /*pour eviter une mechante division par 0*/
			{
				/*RatioX = xs / CtrlNez->pts->Element(3,0);*/
				RatioX = xs / PointSel[0];
				CtrlNez->pts->SetElement(3,0, xs);
				CtrlFui->pts->SetElement(3,0, xs);
				CtrlA->pts->SetElement(3,0, xs);
				CtrlB->pts->SetElement(3,0, xs);
				CtrlC->pts->SetElement(3,0, xs);
				CtrlD->pts->SetElement(3,0, xs);
				CtrlE->pts->SetElement(3,0, xs);
				for(n=1; n<3; n++)
				{
					CtrlNez->pts->SetElement(n,0, CtrlNez->pts->Element(n,0)*RatioX);
					CtrlFui->pts->SetElement(n,0, CtrlFui->pts->Element(n,0)*RatioX);
					CtrlA->pts->SetElement(n,0, CtrlA->pts->Element(n,0)*RatioX);
					CtrlB->pts->SetElement(n,0, CtrlB->pts->Element(n,0)*RatioX);
					CtrlC->pts->SetElement(n,0, CtrlC->pts->Element(n,0)*RatioX);
					CtrlD->pts->SetElement(n,0, CtrlD->pts->Element(n,0)*RatioX);
					CtrlE->pts->SetElement(n,0, CtrlE->pts->Element(n,0)*RatioX);
				}
			}
		}
		//lignes A->D en central
		else if(PointSel == CtrlA->pts->GetLigne(0)){
			CtrlA->pts->AddElement(1,1, DeltaY); PointSel[1] = ys; PointSel[0]=0.0;}
		else if(PointSel == CtrlB->pts->GetLigne(0)){
			CtrlB->pts->AddElement(1,1, DeltaY); PointSel[1] = ys; PointSel[0]=0.0;}
		else if(PointSel == CtrlC->pts->GetLigne(0)){
			CtrlC->pts->AddElement(1,1, DeltaY); PointSel[1] = ys; PointSel[0]=0.0;}
		else if(PointSel == CtrlD->pts->GetLigne(0)){
			CtrlD->pts->AddElement(1,1, DeltaY); PointSel[1] = ys; PointSel[0]=0.0;}
		else if(PointSel == CtrlE->pts->GetLigne(0)){
			CtrlE->pts->AddElement(1,1, DeltaY); PointSel[1] = ys; PointSel[0]=0.0;}
		//vrillage/morphing/EpaiRel/diedre en central
		//else if(PointSel == &CtrlVrillage->pts->Element(0,0)){
		//	CtrlVrillage->pts->Element(1,1) += DeltaY; PointSel[1] = ys; PointSel[0]=0.0;}
		else if(PointSel == CtrlDiedre->pts->GetLigne(0)){
			CtrlDiedre->pts->AddElement(1,1, DeltaY); PointSel[1] = ys; PointSel[0]=0.0;}
		//lignes A->D en extremite
		else if(PointSel == CtrlA->pts->GetLigne(3)){
			CtrlA->pts->AddElement(2,1, DeltaY); PointSel[1] = ys; PointSel[0]=CtrlNez->pts->Element(3,0);}
		else if(PointSel == CtrlB->pts->GetLigne(3)){
			CtrlB->pts->AddElement(2,1, DeltaY); PointSel[1] = ys; PointSel[0]=CtrlNez->pts->Element(3,0);}
		else if(PointSel == CtrlC->pts->GetLigne(3)){
			CtrlC->pts->AddElement(2,1, DeltaY); PointSel[1] = ys; PointSel[0]=CtrlNez->pts->Element(3,0);}
		else if(PointSel == CtrlD->pts->GetLigne(3)){
			CtrlD->pts->AddElement(2,1, DeltaY); PointSel[1] = ys; PointSel[0]=CtrlNez->pts->Element(3,0);}
		else if(PointSel == CtrlE->pts->GetLigne(3)){
			CtrlE->pts->AddElement(2,1, DeltaY); PointSel[1] = ys; PointSel[0]=CtrlNez->pts->Element(3,0);}
		//vrillage/morphing/EpaiRel/diedre en extremite
		//else if(PointSel == &CtrlVrillage->pts->Element(3,0)){
		//	CtrlVrillage->pts->Element(2,1) += DeltaY; PointSel[1] = ys; PointSel[0]=NbNerv-1;}
		else if(PointSel == CtrlDiedre->pts->GetLigne(3)){
			CtrlDiedre->pts->AddElement(2,0, DeltaX); CtrlDiedre->pts->AddElement(2,1, DeltaY);
			PointSel[0] = xs; PointSel[1] = ys;}

		//autres points, sauf points centre/extr�mit� de morphing/epairel/vrillage

		else if( (PointSel != CtrlMorphing->pts->GetLigne(0))

			&&(PointSel != CtrlMorphing->pts->GetLigne(3))

			&&(PointSel != CtrlEpaiRel->pts->GetLigne(0))

			&&(PointSel != CtrlEpaiRel->pts->GetLigne(3))

			&&(PointSel != CtrlVrillage->pts->GetLigne(0))

			&&(PointSel != CtrlVrillage->pts->GetLigne(3)) )

		{

			DeltaY = ys - PointSel[1];

			PointSel[0]=xs;	PointSel[1]=ys;

		}

		

		//impose certaines limites au tangentes

		if( (PointSel == CtrlVrillage->pts->GetLigne(2))

			||(PointSel == CtrlMorphing->pts->GetLigne(2))

			||(PointSel == CtrlEpaiRel->pts->GetLigne(2)) )

			if (PointSel[0]>NbNerv-1) PointSel[0]=NbNerv-1;	

			



		//impose position au vrillage central/extremite ...

		CtrlVrillage->pts->SetElement(0,1,VrillageCent);

		CtrlVrillage->pts->SetElement(3,1,VrillageBout);

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
	int i,imin;
	double dist,distmin;
	double xs,ys,xmin,xmax,ymin,ymax;
	int h,w;
	Courbe* courbe;
	int Fenetre;
	/*memorise position de la souris dans la fenetre*/
	xSouris =x; ySouris = y;
	/*par defaut pas de point ni d'axe selectionn� 
	et pas de rotation/translation/zoom 3d
	et ne pas redessinner tous les axes*/
	PointSel = NULL;
	AxeSel = NULL;
	VisuTousLesAxes = 0;
	/*Rotation3d = OFF;
	Translation3d = OFF;
	ZoomTwist3d = OFF;*/
	/*recupere axe correspondant � la fenetre courante*/
	Fenetre = glutGetWindow();
	if (Fenetre == FenetreForme) {AxeSel=AxeForme; /*printf("AxeForme ->");*/}
	if (Fenetre == FenetreDiedre) {AxeSel=AxeDiedre; /*printf("AxeDiedre ->");*/}
	if (Fenetre == FenetreMorphing) {AxeSel=AxeMorphing; /*printf("AxeMorphing ->");*/}
	if (Fenetre == FenetreVrillage) {AxeSel=AxeVrillage; /*printf("AxeVrillage ->");*/}
	if (Fenetre == FenetreEpaiRel) {AxeSel=AxeEpaiRel; /*printf("AxeEpaiRel ->");*/}
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
	/*autres axes = selection d'un point*/
	else if ((AxeSel != NULL)&&(button == GLUT_LEFT_BUTTON)&&(state == GLUT_DOWN))
	{
		//printf("\nbouton ->");
		/*recupere taille de la fenetre courante*/
		w = glutGet(GLUT_WINDOW_WIDTH);
		h = glutGet(GLUT_WINDOW_HEIGHT);
		/*conversion en coordonn�es de l'axe des coordonn�es souris*/
		xmin=AxeSel->xmin; xmax=AxeSel->xmax;
		ymin=AxeSel->ymin; ymax=AxeSel->ymax;
		xs=xmin+(double)x*(xmax-xmin)/(double)w;
		ys=ymin+((double)h-(double)y)*(ymax-ymin)/h;
		/*recherche du point le plus pret*/
		distmin=100000.0F; imin=0;
		courbe=AxeSel->Courb;
		while(courbe != NULL)
		{
			if ((courbe->type == CTRL_POINT)&&(courbe->pts != NULL)&&(courbe->points == 1))
			{
				for(i=0; i<courbe->pts->GetLignes(); i++)
				{
					dist = sqrt(sqr(xs-courbe->pts->Element(i,0)) + sqr(ys-courbe->pts->Element(i,1)));
					if (dist < distmin)
					{
						distmin=dist;
						PointSel=courbe->pts->GetLigne(i);					
					}
				}
			}
			courbe=courbe->CourbSuiv;
		}
		//printf("dist=%f",distmin);
		if (distmin>0.5) PointSel=NULL;
	}
	/*fin de selection/deplacement d'un point*/
	else if ((AxeSel != NULL)&&(button == GLUT_LEFT_BUTTON)&&(state == GLUT_UP))
	{
		/*redessinne tous les axes*/
		VisuTousLesAxes = 1;
		display();
	}
}

/***************************/
/* positionne les fenetres */
/***************************/
void PositionneFenetres(void)
{
    printf ("\nPositionneFenetres(void)");
	int ws,hs;
	double xf,yf,wf,hf;
	int i;
	FILE *fid;
	char texte[50];
	int f[6];
	/*recupere taille de l'ecran*/
	ws = glutGet( GLUT_SCREEN_WIDTH );
	hs = glutGet( GLUT_SCREEN_HEIGHT );
	/*init tableau des no de fenetre*/
	f[0]=FenetreForme; f[1]=FenetreDiedre; f[2]=FenetreMorphing;
	f[3]=FenetreVrillage; f[4]=FenetreEpaiRel; f[5]=Fenetre3d;
	/*ouverture fichier en lecture*/
	if( (fid = fopen( "PosFenetres.txt", "rt" )) == NULL )
	{
      printf( "\nErreur ouverture fichier 'PosFenetres.txt'" );
	  exit(0);
	}
	/*boucle de lecture des positions/tailles normalis�es des fenetres*/
	for (i=0; i<6; i++)
	{
		glutSetWindow( f[i] );
		fscanf(fid,"%s %f %f %f %f\n",texte, &xf, &yf, &wf, &hf);
		glutPositionWindow((int)(xf*(double)ws), (int)(yf*(double)hs));
		glutReshapeWindow((int)(wf*(double)ws), (int)(hf*(double)hs));
	}
	/*fermeture fichier*/
	if(fclose(fid))
	{
		printf("\nProbleme � la fermeture du fichier");
		exit(0);
	}
}

/*****************/
/* keyboard      */
/*****************/
void keyboard(unsigned char key, int x, int y)
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
	/*message*/
	printf("\nWindShape ");
	printf(AUTEUR_DATE);
	printf("\n");

	/*init glut*/
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize (500, 500);
	glutInitWindowPosition (100, 100);

	/* creation  fenetre/axe/courbes Forme*/
	FenetreForme=glutCreateWindow ("Forme a plat");
	InitFenetre();
	AxeForme=CreerAxe(FenetreForme);
	CourbNez=new Courbe("Nez");
	CtrlNez=new Courbe("Control Nez");
	CourbFui=new Courbe("Fuite"); CtrlFui=new Courbe("Control Fuite");
	CourbA=new Courbe("A"); CtrlA=new Courbe("Control A"); 
	CourbB=new Courbe("B"); CtrlB=new Courbe("Control B"); 
	CourbC=new Courbe("C"); CtrlC=new Courbe("Control C");
	CourbD=new Courbe("D"); CtrlD=new Courbe("Control D");
	CourbE=new Courbe("E"); CtrlE=new Courbe("Control E"); 
	CourbNerv=new Courbe("Nervure");
	/* creation fenetre/axe/courbes Diedre*/
	FenetreDiedre=glutCreateWindow ("Diedre");
	InitFenetre();
	AxeDiedre=CreerAxe(FenetreDiedre);
	CtrlDiedre=new Courbe("Di�dre control"); CourbDiedre=new Courbe("Di�dre");
	/* creation fenetre/axe/courbes Morphing*/
	FenetreMorphing=glutCreateWindow ("Morphing des profils");
	InitFenetre();
	AxeMorphing=CreerAxe(FenetreMorphing);
	AxeMorphing->Norme = OFF;
	//AxeMorphing->XGrid = OFF;
	CtrlMorphing=new Courbe("Morphing Control"); CourbMorphing=new Courbe("Morphing");
	/* creation fenetre/axe/courbes Vrillage*/
	FenetreVrillage=glutCreateWindow ("Vrillage de l'aile");
	InitFenetre();
	AxeVrillage=CreerAxe(FenetreVrillage);
	AxeVrillage->Norme = OFF;
	//AxeVrillage->XGrid = OFF;
	CtrlVrillage=new Courbe("Vrillage Control"); CourbVrillage=new Courbe("Vrillage");
	/* creation fenetre/axe/courbes EpaiRel*/
	FenetreEpaiRel=glutCreateWindow ("Epaisseurs Relatives");
	InitFenetre();
	AxeEpaiRel=CreerAxe(FenetreEpaiRel);
	AxeEpaiRel->Norme = OFF;
	//AxeEpaiRel->XGrid = OFF;
	CtrlEpaiRel=new Courbe("Epaisseur Control"); CourbEpaiRel=new Courbe("Epaisseur");
	/* creation fenetre/axe/courbes Axe3d*/
	Fenetre3d=glutCreateWindow ("Visualisation 3D");
	InitFenetre(); InitLumiere();
	Axe3d=CreerAxe(Fenetre3d); Axe3d->axe3d=ON;
	/*charge fichier de config par default*/
	ChargerFichierConfig("ConfigDefaut.txt");
	/*proprietes*/
	CtrlNez->alt=ON; CtrlFui->alt=ON;
	CtrlA->alt=ON; CtrlB->alt=ON; CtrlC->alt=ON; CtrlD->alt=ON; CtrlE->alt=ON;
	CtrlVrillage->alt=ON; CtrlDiedre->alt=ON;
	CtrlMorphing->alt=ON; CtrlEpaiRel->alt=ON;
	CtrlNez->type=CTRL_POINT; CtrlFui->type=CTRL_POINT;
	CtrlA->type=CTRL_POINT; CtrlB->type=CTRL_POINT;
	CtrlC->type=CTRL_POINT; CtrlD->type=CTRL_POINT; CtrlE->type=CTRL_POINT;
	CtrlVrillage->type=CTRL_POINT; CtrlDiedre->type=CTRL_POINT;
	CtrlMorphing->type=CTRL_POINT;
	CtrlEpaiRel->type=CTRL_POINT;
	/*proprietes*/
	CourbNez->points=OFF; CourbNez->symX=ON;
	CourbFui->points=OFF; CourbFui->symX=ON;
	CourbA->points=OFF; CourbA->symX=ON;
	CourbB->points=OFF; CourbB->symX=ON;
        CourbC->points=OFF; CourbC->symX=ON;
	CourbD->points=OFF; CourbD->symX=ON;
	CourbE->points=OFF; CourbE->symX=ON;
	CourbNerv->points=OFF; CourbNerv->symX=ON; CourbNerv->alt=ON;
	CourbDiedre->points=ON; CourbDiedre->symX=ON;
	/*chainage Axe -> Courbe*/
	AjoutCourbe(AxeForme, CtrlNez); AjoutCourbe(AxeForme, CourbNez);
	AjoutCourbe(AxeForme, CtrlFui); AjoutCourbe(AxeForme, CourbFui);
	AjoutCourbe(AxeForme, CtrlA); AjoutCourbe(AxeForme, CourbA);
	AjoutCourbe(AxeForme, CtrlB); AjoutCourbe(AxeForme, CourbB);
	AjoutCourbe(AxeForme, CtrlC); AjoutCourbe(AxeForme, CourbC);
	AjoutCourbe(AxeForme, CtrlD); AjoutCourbe(AxeForme, CourbD);
	AjoutCourbe(AxeForme, CtrlE); AjoutCourbe(AxeForme, CourbE);
	AjoutCourbe(AxeForme, CourbNerv);
        AjoutCourbe(AxeDiedre, CtrlDiedre); AjoutCourbe(AxeDiedre, CourbDiedre);
	AjoutCourbe(AxeMorphing, CtrlMorphing); AjoutCourbe(AxeMorphing, CourbMorphing);
	AjoutCourbe(AxeVrillage, CtrlVrillage); AjoutCourbe(AxeVrillage, CourbVrillage);
	AjoutCourbe(AxeEpaiRel, CtrlEpaiRel); AjoutCourbe(AxeEpaiRel, CourbEpaiRel);

	/****************************************/
	/*         Here's the GLUI code         */
	/****************************************/
	/*affichage No de version et creation fenetre*/
	printf( "\nCe programme utilise GLUI version: %3.2f",
		GLUI_Master.get_version() );
	glui = GLUI_Master.create_glui( "Boite de dialogue",0,0,0);
		//glutGet( GLUT_SCREEN_HEIGHT )-200);

	/****** rollout fichiers *****/
	GLUI_Rollout *RolloutFichiers =	glui->add_rollout("Fichiers",false);

	/*panel profil centre*/
	GLUI_Panel *panel_profil_centre =
		glui->add_panel_to_panel(RolloutFichiers,"Profil centre" );

	/*nom profil du centre*/
	ProfCent = glui->add_statictext_to_panel( panel_profil_centre, "???" );
	ProfCent->set_text(F->m_strNomProfilCent.c_str());
	/*bouton de chargement*/
	glui->add_column_to_panel(panel_profil_centre, false);
	GLUI_Button *bouton = glui->add_button_to_panel(
		panel_profil_centre, "Load", PROF_CENT, ChargerFichierProfil );
	bouton->set_w(10);
	/*panel profil bout*/
	GLUI_Panel *panel_profil_bout =
		glui->add_panel_to_panel(RolloutFichiers,"Profil bout" );
	/*nom profil du bout*/
	ProfBout = glui->add_statictext_to_panel( panel_profil_bout, "???" );
	ProfBout->set_text(F->m_strNomProfilBout.c_str());
	/*boutons de chargement*/
	glui->add_column_to_panel(panel_profil_bout, false);
	bouton = glui->add_button_to_panel(
		panel_profil_bout, "Load", PROF_BOUT, &ChargerFichierProfil );
	bouton->set_w(10);
	/*panel forme*/
	GLUI_Panel *panel_forme = glui->add_panel_to_panel(RolloutFichiers,"Forme");
	/*nom fichier forme*/
	FicForme = glui->add_statictext_to_panel( panel_forme, "???" );
	FicForme->set_text(NomFichierForme);
	glui->add_column_to_panel(panel_forme, false);
	bouton = glui->add_button_to_panel( panel_forme, "Load", 0, &ChargerFichierForme );
	bouton->set_w(10);
	bouton = glui->add_button_to_panel( panel_forme, "Save", 0, &SauverFichierForme );
	bouton->set_w(10);
	glui->add_column_to_panel(panel_forme, false);
	bouton = glui->add_button_to_panel( panel_forme, "Load2", 0, &ChargerFichierForme2 );
	bouton->set_w(10);
	bouton = glui->add_button_to_panel( panel_forme, "Save2", 0, &SauverFichierForme2 );
	bouton->set_w(10);
	/*fichier DXF*/
        bouton = glui->add_button_to_panel( RolloutFichiers, "3D->DXF", 0, &SauverFichierDXF );
	bouton->set_w(10);
	/*fichier FGen*/
	bouton = glui->add_button_to_panel( RolloutFichiers, "3D->FGen", 0, &SauverFichierFGen );
	bouton->set_w(10);

	/**** rollout parametres ****/
	//glui->add_column(false);
	GLUI_Panel *RolloutParametres = glui->add_rollout( "Parametres", false);
	/*Nombre de caissons*/
	GLUI_Spinner *SpinNbCaiss  =
		glui->add_spinner_to_panel( RolloutParametres,
		"Nbr caissons", GLUI_SPINNER_INT, &(NbCaiss), 0, &ModifNbrAlveoles);
	SpinNbCaiss -> set_int_limits( 4, 1000 );
	/*epaisseur relative au centre*/
	GLUI_Spinner *SpinEpaiRel =
		glui->add_spinner_to_panel( RolloutParametres,
		"Epai. Rel. cent.:", GLUI_SPINNER_FLOAT, &(EpaiRelCent));
	SpinEpaiRel -> set_float_limits( 5.0, 30.0 );



	GLUI_Panel *panel_dwidth = glui->add_panel_to_panel( RolloutParametres, "D.Width" );
	/*progression g�om�trique*/
	GLUI_Spinner *SpinProgGeom =
		glui->add_spinner_to_panel( panel_dwidth,
		"Prog. geom.:", GLUI_SPINNER_FLOAT, &(CoeffProgGeom), 0, &ModifProgGeom );
	SpinProgGeom -> set_float_limits( 0.5, 1.0 );

	GLUI_Spinner *SpinCoeffExp =
		glui->add_spinner_to_panel( panel_dwidth,
		"Exp k:", GLUI_SPINNER_FLOAT, &(CoeffExp), 0, &ModifProgGeom );
	SpinCoeffExp -> set_float_limits( 0.0, 20.0 );

	/*panel vrillage*/
	GLUI_Panel *panel_vrillage = glui->add_panel_to_panel( RolloutParametres, "Vrillage" );
	/*vrillage au centre*/
	GLUI_Spinner *SpinVrillCent =
		glui->add_spinner_to_panel( panel_vrillage,
		"Centre:", GLUI_SPINNER_FLOAT, &VrillageCent, 0, &ModifVrillage );
	SpinVrillCent -> set_float_limits( -10.0, 30.0 );
        /*vrillage au bout*/
	GLUI_Spinner *SpinVrillBout =
    		glui->add_spinner_to_panel( panel_vrillage,
		"Bout:", GLUI_SPINNER_FLOAT, &VrillageBout, 0, &ModifVrillage );
	SpinVrillBout -> set_float_limits( -10.0, 30.0 );
	/*panel taille*/
	GLUI_Panel *panel_taille = glui->add_panel_to_panel( RolloutParametres, "Taille" );
	/*taille -1%*/
	bouton = glui->add_button_to_panel(panel_taille, "-1%", TAILLE_MOINS_1, &ModifTaille);
	bouton->set_w(10);
	/*taille -5%*/
	glui->add_column_to_panel(panel_taille,false);
	bouton = glui->add_button_to_panel(panel_taille, "-5%", TAILLE_MOINS_5, &ModifTaille);
	bouton->set_w(10);
	/*taille +1%*/
	glui->add_column_to_panel(panel_taille,false);
	bouton = glui->add_button_to_panel(panel_taille, "+1%", TAILLE_PLUS_1, &ModifTaille);
	bouton->set_w(10);

	/*taille +5%*/

	glui->add_column_to_panel(panel_taille,false);
	bouton = glui->add_button_to_panel(panel_taille, "+5%", TAILLE_PLUS_5, &ModifTaille);
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
	/*choix visu sym�triqe ou non*/
	glui->add_checkbox_to_panel(
		RolloutVisu,"symetrique",&VisuSymetrique,0,&ModifVisuSymetrique);
	/*choix visu des points de suspentage*/
	glui->add_checkbox_to_panel(
		RolloutVisu,"pts suspentes",&VisuPtsSuspentes,0,&ModifVisuPtsSuspentes);
	glui->add_checkbox_to_panel(
		RolloutVisu,"Axes auto",&AxesAuto,0,&ModifAxesAuto);
	//pour affichage info
	GLUI_Panel *RolloutInfo = glui->add_rollout( "Info forme",false);
	textSurface = glui->add_statictext_to_panel(RolloutInfo, "surface..."); 
	textEnvergure = glui->add_statictext_to_panel(RolloutInfo, "envergure..."); 
	textAllongement = glui->add_statictext_to_panel(RolloutInfo, "allongement..."); 
	textCorde = glui->add_statictext_to_panel(RolloutInfo, "corde..."); 
	textLarg = glui->add_statictext_to_panel(RolloutInfo, "largeur caisson...");
	// Vanya's hack
        GLUI_Panel *RolloutDirectLoad = glui->add_rollout( "Direct load",false);
	// Sling rows points
	GLUI_Panel *panel_load_sling_points = glui->add_panel_to_panel(RolloutDirectLoad,"Nez Fui Sling rows");
	glui->add_column_to_panel(panel_load_sling_points, false);
	bouton = glui->add_button_to_panel( panel_load_sling_points, "Export template", BTN_SLING, &ExportTemplateSling );
	bouton->set_w(10);
	glui->add_column_to_panel(panel_load_sling_points, false);
	bouton = glui->add_button_to_panel( panel_load_sling_points, "Load", BTN_SLING, &LoadCourbeSling );
	bouton->set_w(10);
	// Diedre
	GLUI_Panel *panel_load_diedre = glui->add_panel_to_panel(RolloutDirectLoad,"Diedre points");
	glui->add_column_to_panel(panel_load_diedre, false);
	bouton = glui->add_button_to_panel( panel_load_diedre, "Export template", BTN_DIEDRE, &ExportTemplate );
	bouton->set_w(10);
	glui->add_column_to_panel(panel_load_diedre, false);
	bouton = glui->add_button_to_panel( panel_load_diedre, "Load", BTN_DIEDRE, &LoadCourbe );
	bouton->set_w(10);
	// EpaiRel
	GLUI_Panel *panel_load_epairel = glui->add_panel_to_panel(RolloutDirectLoad,"EpaiRel");
	glui->add_column_to_panel(panel_load_epairel, false);
	bouton = glui->add_button_to_panel( panel_load_epairel, "Export template", BTN_EPAIREL, &ExportTemplate );
	bouton->set_w(10);
	glui->add_column_to_panel(panel_load_epairel, false);
	bouton = glui->add_button_to_panel( panel_load_epairel, "Load", BTN_EPAIREL, &LoadCourbe );
	bouton->set_w(10);
	// morphing
	GLUI_Panel *panel_load_morphing = glui->add_panel_to_panel(RolloutDirectLoad,"Morphing");
	glui->add_column_to_panel(panel_load_morphing, false);
	bouton = glui->add_button_to_panel( panel_load_morphing, "Export template", BTN_MORPHING, &ExportTemplate );
	bouton->set_w(10);
	glui->add_column_to_panel(panel_load_morphing, false);
	bouton = glui->add_button_to_panel( panel_load_morphing, "Load", BTN_MORPHING, &LoadCourbe );
	bouton->set_w(10);
	// vrillage
	GLUI_Panel *panel_load_vrillage = glui->add_panel_to_panel(RolloutDirectLoad,"Vrillage");
        glui->add_column_to_panel(panel_load_vrillage, false);
	bouton = glui->add_button_to_panel( panel_load_vrillage, "Export template", BTN_VRILLAGE, &ExportTemplate );
	bouton->set_w(10);
	glui->add_column_to_panel(panel_load_vrillage, false);
	bouton = glui->add_button_to_panel( panel_load_vrillage, "Load", BTN_VRILLAGE, &LoadCourbe );
	bouton->set_w(10);
	/**** bouton sauver / info quitter ****/
	glui->add_button( "Quitter", 0, Quitter );
	/**** Link windows to GLUI, and register idle callback ******/
	glui->set_main_gfx_window( FenetreForme );
	/* We register the idle callback with GLUI, not with GLUT */
	GLUI_Master.set_glutIdleFunc( NULL );
	/*premiers calculs et trace des axes*/
	AxeSel = Axe3d;
	VisuTousLesAxes = 1;
	display();
	//maj info
	CalculInfoForme(F, &info);
	MajInfo();
	/*boucle glut*/
	glutMainLoop();
	return 0;
}
