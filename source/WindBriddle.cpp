/*  WindBriddle.c			
* calcul suspentage de la forme 3d  ...
* issue de WindShape ou gÃ¯Â¿Â½nÃ¯Â¿Â½rÃ¯Â¿Â½ Ã¯Â¿Â½ partir d'une autre application
* mais respectant le format WindShape !
*/

#pragma warning(disable:4786)
#pragma warning(disable:4514)

/***********/
/* include */
/***********/
#include <afx.h>	//class CString
#include <afxdlgs.h>	//class CFileDialog

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "GL/glut.h"

#include "../commun/plot.h"
#include "../commun/matrice.h"
#include "../commun/geom.h"
#include "../commun/fichier.h"
#include "../commun/profil.h"
#include "../commun/pince.h"
#include "../commun/rasklad.h"
#include "../commun/logger.h"

#include "GL/glui.h"

/**********/
/* define */
/**********/

#define AUTEUR_DATE "(Thierry Pebayle - 24/11/2002)\n[Export FGen + Refactoring Olivier REMAUD 2007]"
#define sqr(f1) ((f1)*(f1))

/**********************/
/* variables globales */
/**********************/

/*live-variable GLUI*/

GLUI *glui = NULL;

int VisuFace=0;		/*choix visu nervures(0) ou surfaces(1)*/
int ProjOrthoPers=0;	/*choix projection orthogonale ou perspective*/
int VisuSymetrique=0;	/*visualise les 2 ou 1 moitiÃ¯Â¿Â½ de l'aile en 3D*/
int VisuPtsSuspentes=1; /*visualise pts de suspentage*/
int DebRangee=1, FinRangee=3, nbNervGroupe=2, DebNerv=1, SauteNerv=0;

float LongSec=70.0f, Incidence=4.0f, PosXYZ[3]={0.0f,0.0f,0.0f};
float profondeur=65.0f;

/*global pour control rotation, ...*/

int xSouris, ySouris;		/*coordonnÃ¯Â¿Â½es fenetre de la souris*/
int Rotation3d = OFF;		/*flag rotation 3d -> bouton gauche souris*/
int Translation3d = OFF;	/*flag translation 3d -> bouton droit souris*/
int ZoomTwist3d = OFF;		/*flag twist 3d -> bouton gauche et droit souris*/

/*fenetres OpenGL*/

int Fenetre3d;		/*visu 3D*/

/* axes graphiques*/

TAxe *Axe3d = NULL;		
TAxe *AxeSel = NULL; 

/*forme*/

Forme *F = NULL;

/*divers tableaux*/

Matrice *IntProfCent = NULL;
Matrice *ExtProfCent = NULL;
Matrice *IntProfBout = NULL;
Matrice *ExtProfBout = NULL;

char NomFichierForme[255];
GLUI_StaticText *FicForme = NULL;

//pour visu pression sur le profil central
Matrice *PtsSuspentes = NULL; //pour recuperer la position 3D des pts de suspentage
int FenetrePression;
TAxe *AxePression = NULL;
Courbe *CourbPression = NULL, *CourbProfil = NULL;
Matrice *P = NULL, *TabBari = NULL, *TabRenc = NULL;
int recalculPressions = FALSE, recalculCentrage = FALSE;
float finesse;
Courbe *CourbBari = NULL, *CourbPrim = NULL, *CourbSec = NULL; 


//pour display des infos de la forme

TInfoForme info;

GLUI_StaticText *textSurface = NULL, *textEnvergure = NULL, *textAllongement = NULL, *textCorde = NULL, *textLarg = NULL;

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



void CalculSuspentage( int control );

void CalculPressions( void );

void AjoutSuspentage( void );



void ModifIncidence( int control );



/***********/
/* MajInfo */
/***********/

void MajInfo(void)
{
	char surf[100], env[100], all[100], cord[100], larg[100];
	//maj chaine de caractÃ¯Â¿Â½re
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
/* CalculCentrageAuto */
/**********************/

void CalculCentrageAuto( void )
{
	int i,j,k;
	float somme[3]={0.0f,0.0f,0.0f}, sommePoids=0.0f, bari[3];
	float assiette, longueur;
	//calcul baricentre de tous les points d'ancrage
	for(i=0; i<P->GetLignes(); i++)
	{
		for(j=0; j<P->GetColonnes(); j++)
		{
			sommePoids += P->Element(i,j);
			for(k=0; k<3; k++) somme[k] += P->Element(i,j) * PtsSuspentes->Element(i*5+j,k);
		}
	}

	for(k=0; k<3; k++) bari[k] = somme[k]/sommePoids;

	//calcul assiette et longueur du cone de suspentage
	assiette = (float)atan(1/finesse) - Incidence*DEG2RAD;
	longueur = profondeur * F->m_pProfils[F->m_nbProfils-1]->m_fNezX / 50.0f;

	//maj "centroid"
	PosXYZ[0] = bari[0]; //x inchangÃ¯Â¿Â½
	PosXYZ[1] = bari[1] - longueur * (float)cos(assiette);
	PosXYZ[2] = bari[2] + longueur * (float)sin(assiette);
	glui->sync_live();
}

/***************/
/* InitFinesse */
/***************/
void InitFinesse( void )
{
	//calcul info pour rÃ¯Â¿Â½cupÃ¯Â¿Â½rer allongement projetÃ¯Â¿Â½
	CalculInfoForme(F, &info);
        //formule "trÃ¯Â¿Â½s" empirique Ã¯Â¿Â½ partir de l'allongement projetÃ¯Â¿Â½
	finesse = info.allongementProj + 1.0f;
	glui->sync_live();
}

/******************/
/* ModifIncidence */
/******************/

void ModifIncidence( int /*control*/ )
{
	recalculPressions = TRUE;
}

/******************/
/* AppliquerAuto */
/******************/

void AppliquerAuto( int /*control*/ )
{
	recalculCentrage = TRUE;
}

/*******************/
/* CalculPressions */
/*******************/

void CalculPressions( void )
{
	//pour appel LVFoil
	Matrice *XB = NULL, *YB = NULL, *X = NULL, *Y = NULL, *Cp = NULL;
	int i, j, nbNerv;
	float v, m, coeffx, coeffyCent, coeffyBout;
	float LongNerv, EpaiRel, EpaiRelProfCent, EpaiRelProfBout, IncidenceCourante;
	float PPan, ratio, xsus[4];
	Matrice *ptsP = NULL, *ptsM = NULL; //pour visu pressions
	Matrice *T = NULL, *R = NULL, *XBpiv = NULL, *YBpiv = NULL; //pour pivotement profil
        float largMin; //pour pondÃ¯Â¿Â½ration par largeur du caisson
	float pas, pond; //pour ponderation effet d'extremitÃ¯Â¿Â½
	//init tableaux & variables
	nbNerv = F->m_nbProfils;
	if ( P != NULL )
		delete(P); 
	P=Zeros(nbNerv,4);
	EpaiRelProfCent = EpaisseurRelative(ExtProfCent, IntProfCent);
	EpaiRelProfBout = EpaisseurRelative(ExtProfBout, IntProfBout);

	//boucle Ã¯Â¿Â½ partir de la 1ere nervure du centre vers l'extrÃ¯Â¿Â½mitÃ¯Â¿Â½
	printf("\n\n");

	for (i = 0; i < nbNerv; i++)
	{
		///////////////////////////
		// calcul profil courant
		///////////////////////////
		//longueur nervure courante
		LongNerv = F->m_pProfils[i]->m_fLength;
		//epaisseur relative
		EpaiRel = F->m_pProfils[i]->m_fWidth;
		if (EpaiRel < 2.0f) EpaiRel = 2.0f; //valeur mini pour le calcul ...
		//angle vrillage
		v = F->m_pProfils[i]->m_fWash;
		//coeff morphing
		m = F->m_pProfils[i]->m_fMorph;
		//calcul coeffx et coeffy des points du profil en
		//fonction de l'Ã¯Â¿Â½paisseur relative et de la longueur de nervure
		coeffx = LongNerv/100.0f;
		coeffyCent = LongNerv*EpaiRel/(EpaiRelProfCent*100.0f);
		coeffyBout = LongNerv*EpaiRel/(EpaiRelProfBout*100.0f);

		///////////////////////////////////////////
		//calcul pression le long du profil courant
		///////////////////////////////////////////

		delete(XB); delete(YB);
		XB=Zeros(ExtProfCent->GetLignes() + IntProfCent->GetLignes() -1,1);
		YB=Zeros(ExtProfCent->GetLignes() + IntProfCent->GetLignes() -1,1);

		//boucle intrados
		for(j=0; j<IntProfCent->GetLignes(); j++)
		{
			XB->SetElement(j,0, IntProfCent->Element(IntProfCent->GetLignes()-j-1,0)*coeffx);
			YB->SetElement(j,0, IntProfCent->Element(IntProfCent->GetLignes()-j-1,1)*coeffyCent*m
				+ IntProfBout->Element(IntProfBout->GetLignes()-j-1,1)*coeffyBout*(1.0f-m));
		}

		//boucle extrados
		for(j=1; j<ExtProfCent->GetLignes(); j++)
		{
			XB->SetElement(j-1+IntProfCent->GetLignes(),0, ExtProfCent->Element(j,0)*coeffx);
			YB->SetElement(j-1+IntProfCent->GetLignes(),0, ExtProfCent->Element(j,1)*coeffyCent*m
        			+ ExtProfBout->Element(j,1)*coeffyBout*(1.0f-m));
		}

		//affichage
		IncidenceCourante = Incidence + (v - F->m_pProfils[0]->m_fWash) * RAD2DEG;
		printf("\rCalcul distribution de pression de la nervure No %02d, incidence=%1.2f",
			i+1, IncidenceCourante);

		delete(X); delete(Y); delete(Cp);
		LVFoil(XB, YB, IncidenceCourante, &X, &Y, &Cp );

		//VoirMat(XB); VoirMat(YB); VoirMat(Cp);
		//////////////////////////////////////////////////
		//calcul pression associÃ¯Â¿Â½e a chaque pts d'attache
		//////////////////////////////////////////////////
		//position en x des points de suspentage

		xsus[0]=coeffx*F->m_pProfils[i]->m_fPosA;
		xsus[1]=coeffx*F->m_pProfils[i]->m_fPosB;
		xsus[2]=coeffx*F->m_pProfils[i]->m_fPosC;
		xsus[3]=coeffx*F->m_pProfils[i]->m_fPosD;

		//boucle sur les panneaux du profil
		for(j=0; j<X->GetLignes(); j++)
		{
			//calcul "portance verticale" du panneau courant
			PPan = Cp->Element(j,0)*(float)fabs(XB->Element(j,0) - XB->Element(j+1,0));
			if (j>=IntProfCent->GetLignes()) PPan = -PPan; // -> positif = porteur 
			//affectation de la pression du panneau courant aux points d'attache
			if		(X->Element(j,0) < xsus[0])
				P->AddElement(i,0, PPan);
			else if (X->Element(j,0) < xsus[1])
			{
				ratio = (X->Element(j,0) - xsus[0])/(xsus[1] - xsus[0]);
				P->AddElement(i,0, (1.0f-ratio)*PPan);
				P->AddElement(i,1, ratio*PPan);
			}
			else if (X->Element(j,0) < xsus[2])
			{
				ratio = (X->Element(j,0) - xsus[1])/(xsus[2] - xsus[1]);
				P->AddElement(i,1, (1.0f-ratio)*PPan);
				P->AddElement(i,2, ratio*PPan);
			}
			else if (X->Element(j,0) < xsus[3])
			{
				ratio = (X->Element(j,0) - xsus[2])/(xsus[3] - xsus[2]);
				P->AddElement(i,2, (1.0f-ratio)*PPan);
				P->AddElement(i,3, ratio*PPan);
			}
			else
				P->AddElement(i,3, PPan);
		}

        	///////////////////////////////////////////////
		//visu des pressions le long du profil central
		///////////////////////////////////////////////	
		if(i==0)
		{	
			//pivotement profil
			Cart2Pol(XB,YB,&T,&R);
			for(j=0; j<T->GetLignes(); j++) 
        			T->RemoveElement(j,0,Incidence*DEG2RAD);

			Pol2Cart(T,R,&XBpiv,&YBpiv);

			//init pts de la courbe profil
			delete(CourbProfil->pts);

			CourbProfil->pts=Zeros(XBpiv->GetLignes(),2);
			for(j=0; j<CourbProfil->pts->GetLignes(); j++)
			{
				CourbProfil->pts->SetElement(j,0,XBpiv->Element(j,0));
				CourbProfil->pts->SetElement(j,1,YBpiv->Element(j,0));
			}
			//recalcul milieu des panels du profil pivotÃ¯Â¿Â½

			for(j=0; j<X->GetLignes(); j++)

			{
				X->SetElement(j,0,(XBpiv->Element(j,0)+XBpiv->Element(j+1,0))/2.0f);
				Y->SetElement(j,0,(YBpiv->Element(j,0)+YBpiv->Element(j+1,0))/2.0f);
			}
			//calcul pts de la courbe pression

			ptsM=Zeros(X->GetLignes(),2);
			for(j=0; j<ptsM->GetLignes(); j++) 
			{ 
				ptsM->SetElement(j,0,X->Element(j,0)); 
				ptsM->SetElement(j,1,Y->Element(j,0)); 
			}
			ptsP = CalculContour(ptsM, Cp, -1);
			delete(CourbPression->pts);
			CourbPression->pts=Zeros(X->GetLignes()*2,2);
			for(j=0; j<X->GetLignes(); j++)
			{
				CourbPression->pts->SetElement(j*2,0,X->Element(j,0));
				CourbPression->pts->SetElement(j*2,1,Y->Element(j,0));
				CourbPression->pts->SetElement(j*2+1,0,ptsP->Element(j,0));
				CourbPression->pts->SetElement(j*2+1,1,ptsP->Element(j,1));
			}

			//visu
			//VisuAxe(AxePression);
			//glutSwapBuffers();
			//libere matrices intermedaires

			if ( T != NULL )
				delete(T);
			if ( R != NULL )
				delete(R);
			if ( XBpiv != NULL )
				delete(XBpiv);
			if ( YBpiv != NULL )
				delete(YBpiv);
			if ( ptsM != NULL )
				delete(ptsM);
			if ( ptsP != NULL )
				delete(ptsP);
			T = R = XBpiv = YBpiv = ptsM = ptsP = NULL;
		}	

	}


        for (i = 0; i < nbNerv; i++) {
            printf ("\nP%2d: %7.3f %7.3f %7.3f %7.3f", i, P->Element(i, 0), P->Element(i, 1), P->Element(i, 2), P->Element(i, 3));
        }
	/////////////////////////////////////////////////////
	//pondÃ¯Â¿Â½ration en fonction de la largeur des caissons
	/////////////////////////////////////////////////////
	//determine largeur mini

	for (largMin=10000.0f, i=0; i<F->m_nbProfils-1; i++)
	{
		//larg = sqrt(sqr(TabForme->Element(i,2)-TabForme->Element(i+1,2))
		//+ sqr(TabForme->Element(i,3)-TabForme->Element(i+1,3)));
		float larg = F->m_pProfils[i]->m_fNezX - F->m_pProfils[i+1]->m_fNezX; //largeur projetÃ¯Â¿Â½e !!!
		if (largMin > larg) largMin = larg;
	}

	//pondere nervure

	float larg = 0.0f;
	for (i=0; i<nbNerv-1; i++)
	{
		//larg = sqrt(sqr(TabForme->Element(i,2)-TabForme->Element(i+1,2))
		//+ sqr(TabForme->Element(i,3)-TabForme->Element(i+1,3)));
		larg = F->m_pProfils[i]->m_fNezX - F->m_pProfils[i+1]->m_fNezX; //largeur projetÃ¯Â¿Â½e !!!
		if (F->m_pProfils[i]->m_fNezX==0.0f) larg*=0.5f;
		for(j=0; j<4; j++) 
        		P->MultiplyElement(i,j,larg/largMin);
	}

	//cas de la derniere nervure ...

	for(j=0; j<4; j++) 
		P->MultiplyElement(i,j,(larg/largMin)/2.0f);



	////////////////////////////////////////////////////////////

	//pondÃ¯Â¿Â½ration pour prendre en compte les effets d'extrÃ¯Â¿Â½mitÃ¯Â¿Â½

	////////////////////////////////////////////////////////////
	pas=1.0f/(nbNerv-1);
	for (i=0; i<nbNerv; i++)
	{
		if(i==nbNerv-1)
			pond=0.4f;
		else
			pond=0.4f+0.6f*sqrt(1-sqr(i*pas));
		for(j=0; j<4; j++) 
			P->MultiplyElement(i,j,pond);
	}
	//VoirMat(P);
	if ( XB != NULL )
		delete(XB);
	if ( YB != NULL )
		delete(YB);

	if ( X != NULL )
		delete(X);
	if ( Y != NULL )
		delete(Y);
	if ( Cp != NULL )
		delete(Cp);

}

/********************/
/* CalculSuspentage */
/********************/

void CalculSuspentage( int /*control*/ )
{
	int nbNerv, nbGroupes, i, j, g, k, nerv;
	float dist, num[3], den;

	//calcul pressions si besoin (modif incidence)
	if (recalculPressions == TRUE) 
	{
		recalculPressions = FALSE;
		CalculPressions();
	}

	//calcul "centrage" si besoin
	if (recalculCentrage == TRUE)
	{
		recalculCentrage = FALSE;
		CalculCentrageAuto();
	}

	//init variables et tableaux de coordonnÃ¯Â¿Â½es 
	nbNerv = F->m_nbProfils;
    // printf ("\n nbNerv==%d", nbNerv);
	//nbGroupes = (nbNerv-DebNerv+1)/nbNervGroupe;
	nbGroupes = (nbNerv-DebNerv+1+SauteNerv)/(nbNervGroupe*(SauteNerv+1));

    // printf ("\n nbGroupes==%d", nbGroupes);
	if ( TabBari != NULL )
		delete(TabBari); 
	TabBari=Zeros(nbGroupes, 3); //coordonnÃ¯Â¿Â½es XYZ des baricentres

	if ( TabRenc != NULL )
		delete(TabRenc); 
	TabRenc=Zeros(nbGroupes, 3); //coordonnÃ¯Â¿Â½es XYZ pts rencontre prim.

	//boucle groupes
	for(g=0; g<nbGroupes; g++)
	{
		//calcul baricentre du groupe de suspentes considÃ¯Â¿Â½rÃ¯Â¿Â½
		num[0]=0.0f; num[1]=0.0f; num[2]=0.0f; den=0.0f;
		//boucle nervures
		//nerv=DebNerv-1 + g*nbNervGroupe;
		nerv=DebNerv-1 + g*(nbNervGroupe*(SauteNerv+1));
		//for(i=nerv; i<nerv+nbNervGroupe; i++)
		for(i=nerv; i<nerv+nbNervGroupe*(SauteNerv+1); i+=1) // i+=SauteNerv+1)
		{
			//boucle rangÃ¯Â¿Â½e
			for(j=DebRangee-1; j<FinRangee; j++)
			{
				//boucle coordonnÃ¯Â¿Â½es x,y,z
				for(k=0; k<3; k++) num[k] += P->Element(i,j)*PtsSuspentes->Element(i*5+j,k);
				den += P->Element(i,j);
			}
		}
		//boucle coordonnÃ¯Â¿Â½es x,y,z
		for(k=0; k<3; k++) 
			TabBari->SetElement(g,k,num[k]/den);
		//calcul point de rencontre des suspentes primaires
		dist=0.0f;
		//pour secondaire de longueur constante
        	//for(k=0; k<3; k++) dist+=sqr(PosXYZ[k]-TabBari->Element(g,k)); dist=(float)sqrt(dist);
		//for(k=0; k<3; k++) TabRenc->Element(g,k)=PosXYZ[k]+(TabBari->Element(g,k)-PosXYZ[k])*LongSec/dist;
		//pour secondaire en %
		//boucle coordonnÃ¯Â¿Â½es x,y,z
		for(k=0; k<3; k++) 
			TabRenc->SetElement(g,k,PosXYZ[k]+(TabBari->Element(g,k)-PosXYZ[k])*LongSec/100.0f);

	}

}



/*******************/
/* AjoutSuspentage */
/*******************/

void AjoutSuspentage( void )
{
	int i, j, g, k;
	int nbNerv, nbGroupes, nbPtsGroupe, nbRangGroupe;
	int ind, indPts;
	//init variables et tableaux de coordonnÃ¯Â¿Â½es 
	nbNerv = F->m_nbProfils; //nombre de nervure de la demi-aile
	nbGroupes = (nbNerv-DebNerv+1+SauteNerv)/(nbNervGroupe*(SauteNerv+1));
	nbRangGroupe = FinRangee-DebRangee+1; //nombre de rangees par groupe
	nbPtsGroupe = nbNervGroupe*nbRangGroupe; //nombre de pts d'attache par groupe

	//init courbes
	CourbBari = new Courbe("Bari"); CourbPrim = new Courbe("Prim"); CourbSec = new Courbe("Sec");
	CourbBari->segments = OFF; CourbBari->pts=Zeros(nbGroupes,3);
	CourbBari->CouleurPoints[0] = 1.0f; 
	CourbBari->CouleurPoints[1] = 0.0f; 
	CourbBari->CouleurPoints[2] = 0.0f; 
	CourbPrim->alt = ON; CourbPrim->pts=Zeros(nbGroupes*nbPtsGroupe*2,3); 
	CourbSec->alt = ON; CourbSec->pts=Zeros(nbGroupes*2,3); 



	//test symÃ¯Â¿Â½trique

	if(VisuSymetrique == 1)
	{
		CourbBari->symX = ON;
		CourbPrim->symX = ON;
		CourbSec->symX = ON;
	}

	//boucle groupes
	for(g=0; g<nbGroupes; g++)
	{
		//boucle nervures
		for(i=0; i<nbNervGroupe; i++)
		{
			//boucle rangÃ¯Â¿Â½e
			for(j=0; j<nbRangGroupe; j++)
			{
				ind = (g*nbPtsGroupe + i*nbRangGroupe + j)*2;
				indPts = (DebNerv-1+(g*nbNervGroupe+i)*(SauteNerv+1))*5+(j+DebRangee-1);
				//boucle coordonnÃ¯Â¿Â½es x,y,z
				for(k=0; k<3; k++)
				{
					//suspentes primaires
					CourbPrim->pts->SetElement(ind,k, PtsSuspentes->Element(indPts,k));
					CourbPrim->pts->SetElement(ind+1,k, TabRenc->Element(g,k));
				}
			}
		}
		//suspentes secondaires et baricentres
		for(k=0; k<3; k++)
		{
			CourbSec->pts->SetElement(g*2,k, TabRenc->Element(g,k));
			CourbSec->pts->SetElement(g*2+1,k, PosXYZ[k]);
			CourbBari->pts->SetElement(g,k, TabBari->Element(g,k));
		}
	}

	//ajout courbes a l'axe 3D
	AjoutCourbe(Axe3d, CourbBari);
	AjoutCourbe(Axe3d, CourbPrim);
	AjoutCourbe(Axe3d, CourbSec);
}



/*******************/

/* SauveSuspentage */

/*******************/

void SauveSuspentage( int /*control*/ )

{
	int i, j, g, k, ind;
	int nbNerv, nbGroupes, nbPtsGroupe, nbRangGroupe;
	float dist;

	CString NomFichier;
	LPTSTR PtrNomFichier;
	FILE *fid;

	//ouverture boite de dialogue
	CFileDialog DlgOpen(FALSE, NULL, "*.txt", OFN_OVERWRITEPROMPT, NULL, NULL);
	if (DlgOpen.DoModal()==IDOK)
	{
		//recupere nom de fichier
		NomFichier=DlgOpen.GetPathName();
		PtrNomFichier = NomFichier.GetBuffer(1);

		//message
		printf("\nEcriture fichier de suspentage: '%s'",NomFichier);
		//ouverture fichier en ecriture

		if( (fid = fopen( NomFichier, "wt" )) == NULL )
		{
			printf( "\nErreur ouverture fichier '%s'", NomFichier);
			exit(0);
		}



		//init variables et tableaux de coordonnÃ¯Â¿Â½es 
		nbNerv = F->m_nbProfils; //nombre de nervure de la demi-aile
		nbGroupes = (nbNerv-DebNerv+1+SauteNerv)/(nbNervGroupe*(SauteNerv+1));
		nbRangGroupe = FinRangee-DebRangee+1; //nombre de rangees par groupe
		nbPtsGroupe = nbNervGroupe*nbRangGroupe; //nombre de pts d'attache par groupe

		//boucle groupes
		for(g=0; g<nbGroupes; g++)
		{
			fprintf(fid,"\n\nGroupe de suspentes No %02d:",g+1);
			fprintf(fid,"\n\tLongueurs des suspentes primaires (nervure/rangee):");
			//boucle nervures
			for(i=0; i<nbNervGroupe; i++)
			{
				//boucle rangÃ¯Â¿Â½e
				for(j=0; j<nbRangGroupe; j++)
				{
					//calcul longueur suspente primaire courante
					ind = (g*nbPtsGroupe + i*nbRangGroupe + j)*2;
					dist = 0.0f;
					for(k=0; k<3; k++)
						dist += sqr(CourbPrim->pts->Element(ind,k) - CourbPrim->pts->Element(ind+1,k));					
					dist = (float)sqrt(dist);
					//fprintf(fid,"\n\t\t(%02d,%02d) %2.3fm", i+1, j+1, dist);
					fprintf(fid,"\n\t\t(%02d,%02d) %2.3fm",
						DebNerv+(g*nbNervGroupe+i)*(SauteNerv+1),
						j+DebRangee,
						dist);
				}
			}

			//suspentes secondaires et baricentres
			dist = 0.0f;
			for(k=0; k<3; k++)
				dist += sqr(CourbSec->pts->Element(g*2,k) - CourbSec->pts->Element(g*2+1,k));					
			dist = (float)sqrt(dist);
			fprintf(fid,"\n\tSuspente secondaire: %2.3fm", dist);
		}

		//fermeture fichier
		if(fclose(fid))
		{
			printf("\nProbleme Ã¯Â¿Â½ la fermeture du fichier");
			exit(0);
		}

	}//fichier selectionnÃ¯Â¿Â½

}

/***********************/
/* ChargerFichierForme */
/***********************/

void ChargerFichierForme( int /*control*/ )
{
	CString NomFichier;

	char* PtrNomFichier = NULL;

	Forme* oldF = NULL;

	

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

		if ( oldF != NULL )
			delete(oldF);
		oldF = NULL;


		//maj interface GLUI

		FicForme->set_text(NomFichierForme);

		

		//chargement profils

		LectureFichierProfil(F->m_strNomProfilCent.c_str(), &ExtProfCent, &IntProfCent);

		LectureFichierProfil(F->m_strNomProfilBout.c_str(), &ExtProfBout, &IntProfBout);

		InterpoleProfilBout(&ExtProfBout, ExtProfCent);

		InterpoleProfilBout(&IntProfBout, IntProfCent);

	}

	/*mise a jour affichage*/

	AxeSel = Axe3d;

	PtsSuspentes = NULL;

	recalculPressions = TRUE;

	recalculCentrage = TRUE;

	InitFinesse();

	display();

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

		if ( F != NULL )
			delete(F);

		if ( IntProfCent != NULL )
			delete(IntProfCent);

		if ( ExtProfCent != NULL )
			delete(ExtProfCent);

		if ( IntProfBout != NULL )
			delete(IntProfBout);

		if ( ExtProfBout != NULL )
			delete(ExtProfBout);

		exit(0);
	}
}



/***************/

/* CalculVue3d */

/***************/

void CalculVue3d(void)
{
	Matrice *XExt=NULL,*YExt=NULL,*ZExt=NULL; //coordonnÃ¯Â¿Â½es 3D extrados
	Matrice *XInt=NULL,*YInt=NULL,*ZInt=NULL; //coordonnÃ¯Â¿Â½es 3D intrados

	/*par defaut destruction des courbes ou mesh de l'axe 3D*/
	LibererCourbesAxe(Axe3d); LibererMeshsAxe(Axe3d);

	/*calcul forme*/
	CalculForme3D( F, 0, 0.0f,
		ExtProfCent, IntProfCent, ExtProfBout, IntProfBout,
		&XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);

	/*ajoute visu 3D de la forme*/
	AjoutForme3D(Axe3d, XExt, YExt, ZExt, XInt, YInt, ZInt, VisuFace, VisuSymetrique);

	/*ajoute points de suspentage*/
	//delete(PtsSuspentes); inutile car dÃ¯Â¿Â½ja fait avec LibererCourbesAxes ...

	if(VisuPtsSuspentes){
		AjoutPtsSuspentage(Axe3d, F, IntProfCent, XExt, YExt, ZExt, XInt, YInt, ZInt,
			VisuSymetrique, &PtsSuspentes);}

	/*ajoute le suspentage*/
	CalculSuspentage(0);
	AjoutSuspentage();

	/*liberation matrices intermediaires !!!*/
	delete(XExt); delete(YExt); delete(ZExt);
	delete(XInt); delete(YInt); delete(ZInt);
}

/****************/
/* display      */
/****************/

void display(void)
{
	CalculVue3d();
	VisuAxe(Axe3d);
	glutSwapBuffers();
	VisuAxe(AxePression);
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

	/*par defaut pas de point ni d'axe selectionnÃ¯Â¿Â½ 
	et pas de rotation/translation/zoom 3d
	et ne pas redessinner tous les axes*/

	AxeSel = NULL;

	/*Rotation3d = OFF;
	Translation3d = OFF;
	ZoomTwist3d = OFF;*/

	/*recupere axe correspondant Ã¯Â¿Â½ la fenetre courante*/
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
	printf("\nWindBriddle ");
	printf(AUTEUR_DATE);
	printf("\n");
	/*init glut*/
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize (500, 500);
	glutInitWindowPosition (100, 100);

	/* pour affichage pressions ... */
	// creation fenetre/axe/courbes Pression

        FenetrePression=glutCreateWindow ("Distributions des Pressions en central");
	InitFenetre();
	AxePression=CreerAxe(FenetrePression);
	CourbPression=new Courbe("Pressions");
	CourbProfil=new Courbe("Profil");
	AjoutCourbe(AxePression, CourbPression);
	AjoutCourbe(AxePression, CourbProfil);
	CourbPression->alt = ON; CourbPression->points = OFF;

	/* creation fenetre/axe/courbes Axe3d*/
	Fenetre3d=glutCreateWindow ("Visualisation 3D");
	InitFenetre(); InitLumiere();
	Axe3d=CreerAxe(Fenetre3d); Axe3d->axe3d=ON;

	/*recupere taille de l'ecran en pixel*/
	ws = glutGet( GLUT_SCREEN_WIDTH );
	hs = glutGet( GLUT_SCREEN_HEIGHT );
	
	/*redefini taille et position fenetre a partir de coordonnÃ¯Â¿Â½es normalisÃ¯Â¿Â½es*/
	glutSetWindow( Fenetre3d );
	glutPositionWindow((int)(0.1*(float)ws), (int)(0.1*(float)hs));
	glutReshapeWindow((int)(0.8*(float)ws), (int)(0.8*(float)hs));

	//lecture fichier de forme
	strcpy(NomFichierForme,"10test.txt");
	F = LectureFichierForme(NomFichierForme);

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
	glui = GLUI_Master.create_glui( "Boite de dialogue",0,0,0);

	/****** rollout fichiers *****/

        GLUI_Rollout *RolloutForme =	glui->add_rollout("Forme",false);

	
	/*nom fichier forme*/

	FicForme = glui->add_statictext_to_panel( RolloutForme, "???" );
	FicForme->set_text(NomFichierForme);
	glui->add_column_to_panel(RolloutForme, false);
	GLUI_Button *bouton = glui->add_button_to_panel( RolloutForme, "Load", 0, &ChargerFichierForme );
	bouton->set_w(10);

	/****** rollout visu ******/
	//glui->add_column(false);

	GLUI_Panel *RolloutVisu = glui->add_rollout( "Visualisation",false);
	/*visualisation nervures/surface*/

	GLUI_Panel *panel_visu = glui->add_panel_to_panel( RolloutVisu, "Type" );
	GLUI_RadioGroup *radio_visu =
		glui->add_radiogroup_to_panel(panel_visu,&VisuFace,0,ModifVisu3d);

	glui->add_radiobutton_to_group( radio_visu, "Nervures" );
	glui->add_radiobutton_to_group( radio_visu, "Surfaces" );

	/*choix projection orthogonale ou perspective*/

	GLUI_Panel *panel_proj = glui->add_panel_to_panel( RolloutVisu, "Projection" );
	GLUI_RadioGroup *radio_proj =
		glui->add_radiogroup_to_panel(panel_proj,&ProjOrthoPers,0,&ModifProjection3d);
	glui->add_radiobutton_to_group( radio_proj, "Orthogonale" );
	glui->add_radiobutton_to_group( radio_proj, "Perspective" );
	/*choix visu symÃ¯Â¿Â½triqe ou non*/
	glui->add_checkbox_to_panel(
		RolloutVisu,"symetrique",&VisuSymetrique,0,&ModifVisuSymetrique);

	/*choix visu des points de suspentage*/
	//fait planter si actif !!!
	//glui->add_checkbox_to_panel(
	//	RolloutVisu,"pts suspentes",&VisuPtsSuspentes,0,ModifVisuSymetrique);

	/****** rollout suspentage ******/

	GLUI_Panel *RolloutSuspentage = glui->add_rollout( "Suspentage",false);

	/*groupe de rangÃ¯Â¿Â½e: debut et fin*/

	GLUI_Panel *panel_rangee = glui->add_panel_to_panel( RolloutSuspentage, "Grouper les rangees" );
	GLUI_Spinner *SpinDebRangee =
		glui->add_spinner_to_panel( panel_rangee, "de ", GLUI_SPINNER_INT, &(DebRangee));
	GLUI_Spinner *SpinFinRangee =
		glui->add_spinner_to_panel( panel_rangee, " a ", GLUI_SPINNER_INT, &(FinRangee));
	SpinDebRangee -> set_int_limits( 1, 4 );
	SpinFinRangee -> set_int_limits( 1, 4 );

	/*groupe de nervures*/

	GLUI_Spinner *SpinnbNervGroupe =
		glui->add_spinner_to_panel( RolloutSuspentage,
		"Grouper les nervures par ", GLUI_SPINNER_INT, &(nbNervGroupe));

	SpinnbNervGroupe -> set_int_limits( 1, 10 );
	/*commencement du suspentage*/

	GLUI_Spinner *SpinDebSuspentage =
		glui->add_spinner_to_panel( RolloutSuspentage,
		"Commencer a partir de la nervure ", GLUI_SPINNER_INT, &(DebNerv));
	SpinDebSuspentage -> set_int_limits( 1, F->NbCaiss );

	/*sauter nervures*/

	GLUI_Spinner *SpinSauteNervure =
		glui->add_spinner_to_panel( RolloutSuspentage,
		"Sauter nervures ", GLUI_SPINNER_INT, &(SauteNerv));
	SpinSauteNervure -> set_int_limits( 0, 3 );

	/*longueur suspente secondaire*/

	GLUI_Spinner *SpinLongSec =	glui->add_spinner_to_panel( RolloutSuspentage,
		"Longueur secondaires (%)= ", GLUI_SPINNER_FLOAT, &(LongSec));

	SpinLongSec -> set_float_limits( 0.0, 90.0 );
	/*incidence centrale*/

	GLUI_EditText *EditIncidence =	glui->add_edittext_to_panel( RolloutSuspentage,
		"Incidence nervure centrale ", GLUI_EDITTEXT_FLOAT, &(Incidence), 0, &ModifIncidence);

	EditIncidence -> set_float_limits( -5.0, 10.0 );
	/*position point de rencontre des suspentes en X,Y,Z*/

	GLUI_Panel *panel_pts = glui->add_panel_to_panel( RolloutSuspentage, "Point de rencontre" );
	GLUI_Spinner *SpinPosX =
		glui->add_spinner_to_panel( panel_pts, "X = ", GLUI_SPINNER_FLOAT, &(PosXYZ[0]));

	GLUI_Spinner *SpinPosY =
		glui->add_spinner_to_panel( panel_pts, "Y = ", GLUI_SPINNER_FLOAT, &(PosXYZ[1]));

	GLUI_Spinner *SpinPosZ =
		glui->add_spinner_to_panel( panel_pts, "Z = ", GLUI_SPINNER_FLOAT, &(PosXYZ[2]));

	SpinPosX -> set_float_limits( -10.0, 10.0 );
	SpinPosY -> set_float_limits( -10.0, 10.0 );
	SpinPosZ -> set_float_limits( -10.0, 10.0 );

	/*positionnement automatique*/

	GLUI_Panel *panel_auto = glui->add_panel_to_panel( panel_pts, "Automatique" );
	GLUI_Spinner *SpinFinesse =
		glui->add_spinner_to_panel( panel_auto, "Finesse = ", GLUI_SPINNER_FLOAT, &finesse);
	SpinFinesse -> set_float_limits( 1.0, 10.0 );
	GLUI_Spinner *SpinProfondeur =
		glui->add_spinner_to_panel( panel_auto, "Profondeur (% envergure) = ",
		GLUI_SPINNER_FLOAT, &profondeur);
	SpinProfondeur -> set_float_limits( 10.0, 200.0 );
	glui->add_button_to_panel( panel_auto ,"Appliquer", 0, &AppliquerAuto );
	/*fichier suspentage*/
	glui->add_button_to_panel( RolloutSuspentage, "SAUVER SUSPENTAGE", 0, &SauveSuspentage );

	/**** affichage info forme ****/
	GLUI_Panel *RolloutInfo = glui->add_rollout( "Info forme",false);
	textSurface = glui->add_statictext_to_panel(RolloutInfo, "surface..."); 
	textEnvergure = glui->add_statictext_to_panel(RolloutInfo, "envergure..."); 
	textAllongement = glui->add_statictext_to_panel(RolloutInfo, "allongement..."); 
	textCorde = glui->add_statictext_to_panel(RolloutInfo, "corde..."); 
	textLarg = glui->add_statictext_to_panel(RolloutInfo, "largeur caisson...");
	CalculInfoForme(F, &info);
	MajInfo();

	/**** bouton quitter ****/
	glui->add_button( "Quitter", 0, &Quitter );

	/**** Link windows to GLUI, and register idle callback ******/
	glui->set_main_gfx_window( Fenetre3d );

	/* We register the idle callback with GLUI, not with GLUT */
	GLUI_Master.set_glutIdleFunc( NULL );

	/*premiers calculs et trace des axes*/
	AxeSel = Axe3d;
	PtsSuspentes = NULL;
	recalculPressions = TRUE;
	recalculCentrage = TRUE;
	InitFinesse();
	display();

	

	/*boucle glut*/

	glutMainLoop();

	return 0;
}
