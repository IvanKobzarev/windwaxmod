/**************************************

* gestion chargement divers fichiers

* forme ,profil, patrons ...

***************************************/



#include <stdlib.h>

#include <stdio.h>

#include <string.h>

#include <conio.h>

#include <math.h>

#include <afx.h>		//class CString

#include <afxdlgs.h>	//class CFileDialog


#include <fstream>
#include <iomanip>
using namespace std;

#define sqr(f1) ((f1)*(f1))


#include "plot.h"
#include "matrice.h"
#include "fichier.h"
#include "profil.h"

/////////////////////////////////////////////////////////////

// NOM TrouveMotDansFichierTexte

// Role: Recherche un mot dans un fichier texte

// Entree: Ptr sur le fichier Texte et Mot a rechercher

// Sortie: Ptr sur le fichier place apres le mot a rechercher

//		(pour pouvoir lire le texte suivant ce	mot)

//		ou code d'erreur si pas trouve:

//			TRUE si trouve, FALSE sinon

//////////////////////////////////////////////////////////////
Forme::Forme() 
	: m_pProfils( NULL ), m_nbProfils(0)
{
	//valeurs par default
	CoeffProgGeom = 0.95f;
	EpaiRelCent = 17.0f;

	NbCaiss = 21;

	//init ptr tableaux a NULL
	mCtrlNez = NULL;
	mCtrlFui = NULL;
	mCtrlA = NULL;
	mCtrlB = NULL;
	mCtrlC = NULL;
	mCtrlD = NULL;
	mCtrlE = NULL;
	mCtrlDiedre = NULL;
	mCtrlEpaiRel = NULL;
	mCtrlMorphing = NULL;
	mCtrlVrillage = NULL;
}

void Forme::DeleteProfils()
{
	if ( m_pProfils != NULL )
	{
		for( int i = 0; i<m_nbProfils; i++ )
		{
			if ( m_pProfils[i] != NULL )
				delete m_pProfils[i];
			m_pProfils[i] = NULL;
		}
		delete m_pProfils; 
	}
	m_pProfils = NULL;
	m_nbProfils = 0;
}

void Forme::AllocateProfils( int nb )
{
	DeleteProfils();
	m_nbProfils = nb;
	m_pProfils = new Profil*[nb];
	for( int i = 0; i < m_nbProfils; i++ )
	{
		m_pProfils[i] = new Profil();
	}
}

Forme::~Forme() 
{ 
	DeleteProfils();

	if ( mCtrlNez != NULL ) 
		delete(mCtrlNez);
	mCtrlNez = NULL;

	if ( mCtrlFui != NULL ) 
		delete(mCtrlFui);
	mCtrlFui = NULL;

	if ( mCtrlA != NULL ) 
		delete(mCtrlA);
	mCtrlA = NULL;

	if ( mCtrlB != NULL ) 
		delete(mCtrlB);
	mCtrlB = NULL;

	if ( mCtrlC != NULL ) 
		delete(mCtrlC);
	mCtrlC = NULL;

	if ( mCtrlD != NULL ) 
		delete(mCtrlD);
	mCtrlD = NULL;

	if ( mCtrlE != NULL ) 
		delete(mCtrlE);
	mCtrlE = NULL;

	if ( mCtrlDiedre != NULL ) 
		delete(mCtrlDiedre);
	mCtrlDiedre = NULL;

	if ( mCtrlEpaiRel != NULL ) 
		delete(mCtrlEpaiRel);
	mCtrlEpaiRel = NULL;

	if ( mCtrlMorphing != NULL ) 
		delete(mCtrlMorphing);
	mCtrlMorphing = NULL;

	if ( mCtrlVrillage != NULL ) 
		delete(mCtrlVrillage);
	mCtrlVrillage = NULL;
}

bool TrouveMotDansFichierTexte(FILE* fid, char* Mot)

{

	//var locale

	char Chaine[100];

	//retour au debut du fichier

	rewind(fid);

	//recherche du mot dans le fichier texte

	while ( ((fscanf(fid,"%s", Chaine))!= EOF) &&

		(strcmp(Chaine,Mot)!=0) ) {}

	//a-t-on trouve le mot ?

	return strcmp(Chaine,Mot) == 0 ;
}



/************************/

/* LectureFichierProfil */

/************************/

void LectureFichierProfil(const char* NomProf, Matrice** extrados, Matrice** intrados)

{

	FILE *fid;

	char line[255];

	int NbInt, NbExt, i;

	float x,y,coeff;

	//rajoute pour format Selig ...

	int nbData, iNez;

	Matrice *data = NULL;



	/* ouverture fichier profil */

	printf("\nLecture fichier profil: '%s'", NomProf);

	if( (fid = fopen( NomProf, "rt" )) == NULL )

	{

		printf( "\nErreur ouverture fichier: '%s'", NomProf );

		//fflush(stdin); getch();

	}

	else

	{

		/* lecture entete du fichier */

		if( fgets( line, 255, fid ) == NULL)

		{

			printf( "\nfgets error" );

			rewind(fid);

		}

		else 

			printf( "\n\tdescription: %s", line);



		//test format SELIG ou NACA ...

		fscanf(fid,"%f",&x);

		if(x==1.0f)

		{

			///////////////

			//FORMAT SELIG 

			///////////////

			

			//1ere lecture pour determiner le nombre de data

			nbData = 0;

			while( fscanf(fid,"%f %f", &x, &y) != EOF ) nbData++;

			//lecture des datas

			rewind(fid);

			fgets( line, 255, fid ); //premiere ligne

			data = new Matrice(nbData, 2);

			for(i=0; i<nbData; i++)
			{
				float f1, f2;
				fscanf( fid, "%f %f", &f1, &f2 );
				data->SetElement(i,0,f1);
				data->SetElement(i,1,f2);
			}
			//determine position du nez !

			iNez = 0; x=10000.0f;

			for(i=0; i<nbData; i++)

			{

				if(x > data->Element(i,0))

				{

					x=data->Element(i,0);

					iNez=i;

				}

			}

			//recopie data dans le tableau d'extrados

			delete(*extrados); NbExt = iNez+1;

			*extrados=new Matrice(NbExt,2);

			for(i=0; i<NbExt; i++)

			{

				(*extrados)->SetElement(NbExt-i-1,0,data->Element(i,0));

				(*extrados)->SetElement(NbExt-i-1,1,data->Element(i,1));

			}

			//recopie data dans le tableau d'intrados

			delete(*intrados); NbInt = nbData-iNez;

			*intrados=new Matrice(NbInt,2);

			for(i=0; i<NbInt; i++)

			{

				(*intrados)->SetElement(i,0,data->Element(i+NbExt-1,0));

				(*intrados)->SetElement(i,1,data->Element(i+NbExt-1,1));

			}

			//liberation memoire

			if ( data != NULL )
				delete(data);
			data = NULL;

		}

		else

		{

			//////////////

			//FORMAT NACA

			//////////////

			

			// lecture nombre de points en intrados et extrados

			rewind(fid);

			fgets( line, 255, fid ); //premiere ligne

			fscanf(fid,"%f %f", &x, &y);

			NbExt=(int)floor(x);	NbInt=(int)floor(y);		

			// lecture pts Extrados

			delete(*extrados);

			*extrados=new Matrice(NbExt,2);

			for(i=0; i<NbExt; i++)
			{
				float f1, f2;
				fscanf(fid,"%f %f", &f1, &f2 );
				(*extrados)->SetElement(i,0,f1);
				(*extrados)->SetElement(i,1,f2);
			}

			// lecture pts Intrados

			delete(*intrados);

			*intrados=new Matrice(NbInt,2);

			for(i=0; i<NbInt; i++)
			{
				float f1, f2;
				fscanf(fid,"%f %f", &f1, &f2 );
				(*intrados)->SetElement(i,0,f1);
				(*intrados)->SetElement(i,1,f2);
			}

		}

		

		//fermeture fichier

		if(fclose(fid))

		{

			printf("\nProbleme � la fermeture du fichier '%s'", NomProf);

			fflush(stdin); getch();

		}



		//met bord d'attaque en x=0.0f

		if((*extrados)->Element(0,0) != 0.0f){

			for(i=0; i<NbExt; i++) (*extrados)->RemoveElement(i,0,(*extrados)->Element(0,0));}

		if((*intrados)->Element(0,0) != 0.0f){

			for(i=0; i<NbInt; i++) (*intrados)->RemoveElement(i,0, (*intrados)->Element(0,0));}

		

		//met bord de fuite en x=100.0f

		if((*extrados)->Element(NbExt-1,0) != 100.0f)

		{

			coeff = 100.0f/(*extrados)->Element(NbExt-1,0);

			for(i=0; i<NbExt; i++)

			{

				(*extrados)->MultiplyElement(i,0, coeff);

				(*extrados)->MultiplyElement(i,1, coeff);

			}

		}

		if((*intrados)->Element(NbInt-1,0) != 100.0f)

		{

			coeff = 100.0f/(*intrados)->Element(NbInt-1,0);

			for(i=0; i<NbInt; i++)

			{

				(*intrados)->MultiplyElement(i,0, coeff);

				(*intrados)->MultiplyElement(i,1, coeff);

			}

		}

	}

}



/***********************/

/* LectureFichierForme */

/***********************/

Forme* LectureFichierForme(char* NomFic)

{

	int i,j,n;

	FILE *fid;

	Matrice *m=NULL;

	float noVer;

	Forme *f=NULL;

	//liste des mots clefs du fichier

	char* motsClef[18] = {"VERSION",

		"NB_ALVEOLES", "COEFF_PROG_GEOM", "EPAI_REL_CENTRE",

		"BORD_ATTAQUE_FORME","BORD_DE_FUITE_FORME","SUSPENTAGE_LIGNE_A",

		"SUSPENTAGE_LIGNE_B","SUSPENTAGE_LIGNE_C","SUSPENTAGE_LIGNE_D", "SUSPENTAGE_LIGNE_E",

		"DIEDRE","MORPHING","VRILLAGE","EPAISSEUR_RELATIVE",

		"PROFIL_CENTRE", "PROFIL_BOUT", "TABLEAU_FORME" };

	

	/**** message ****/

	printf("\nLecture fichier de Forme: '%s'",NomFic);

	

	/**** ouverture fichier en lecture ****/

	if( (fid = fopen( NomFic, "rt" )) == NULL )

	{

		printf( "\nErreur ouverture fichier '%s'", NomFic);

	}

	else

	{

		/**** allocation memoire ****/

		f = new Forme();

		/**** Lecture des points de controle ****/

		/*Boucle de lecture*/

		for (i=0; i<18; i++)

		{

			//recherche mot clef

			if(!TrouveMotDansFichierTexte(fid, motsClef[i]))

				printf("\n\tWarning: mot clef '%s' non trouv� !", motsClef[i]);

			else

			{

				//action en fonction du mot clef

				switch(i)

				{

				case 0: fscanf(fid,"%f", &noVer); 

					//if(noVer==NO_VERSION)

					break;

				case 1: fscanf(fid,"%d", &(f->NbCaiss));break;

				case 2: fscanf(fid,"%f", &(f->CoeffProgGeom)); break;

				case 3: fscanf(fid,"%f", &(f->EpaiRelCent)); break;

				case 4:

				case 5:

				case 6:

				case 7:

				case 8:

				case 9:

				case 10:

				case 11:

				case 12:

				case 13:

				case 14: {

					m = new Matrice(4,2);

					float f1, f2, f3, f4, f5, f6, f7, f8;

					fscanf(fid,"%f %f %f %f %f %f %f %f", &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8 );

					m->SetElement(0,0,f1); m->SetElement(0,1,f2);
					m->SetElement(1,0,f3); m->SetElement(1,1,f4);
					m->SetElement(2,0,f5); m->SetElement(2,1,f6);
					m->SetElement(3,0,f7); m->SetElement(3,1,f8);

					switch(i)

					{

					case 4: f->mCtrlNez = m; break;

					case 5: f->mCtrlFui = m; break;

					case 6: f->mCtrlA = m; break;

					case 7: f->mCtrlB = m; break;

					case 8: f->mCtrlC = m; break;

					case 9: f->mCtrlD = m; break;

					case 10: f->mCtrlE = m; break;

					case 11: f->mCtrlDiedre = m; break;

					case 12: f->mCtrlMorphing = m; break;

					case 13: f->mCtrlVrillage = m; break;

					case 14: f->mCtrlEpaiRel = m; break;

					default:;

					}

					break;}

				case 15: {
					char temp[1024];
							fscanf(fid,"%s", temp); 
							f->m_strNomProfilCent = temp;
							break;
						 }

				case 16: {
							char temp[1024];
							fscanf(fid,"%s", temp); 
							f->m_strNomProfilBout = temp;
							break;
						 }

				case 17:

					f->DeleteProfils();

					fscanf(fid,"\n%d",&n);

						//fscanf(fid,"%d",&n);

						f->AllocateProfils(n);

						//boucle � partir de la 1ere nervure du centre vers l'extr�mit�

						for (j=0; j<f->m_nbProfils; j++)
						{
							float f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13;

							fscanf(fid,"\n%d %f %f %f %f %f %f %f %f %f %f %f %f %f",
										  &n, //no Nervure
										  &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11, &f12, &f13 );

								f->m_pProfils[j]->m_fLength = f1;	//longueur
								f->m_pProfils[j]->m_fWidth  = f2;	//epaisseur relative
								f->m_pProfils[j]->m_fNezX	= f3;	//xnez
								f->m_pProfils[j]->m_fNezY	= f4;	//ynez
								f->m_pProfils[j]->m_fNezZ	= f5;	//znez
								f->m_pProfils[j]->m_fInclin	= f6 * DEG2RAD;	//inclinaison / horizontale
								f->m_pProfils[j]->m_fWash	= f7 * DEG2RAD;	//vrillage / horizontale
								f->m_pProfils[j]->m_fMorph	= f8;	//morphing / nervure centrale
								f->m_pProfils[j]->m_fPosA	= f9;	//pos relative ligne A
								f->m_pProfils[j]->m_fPosB	= f10;	//pos relative ligne B
								f->m_pProfils[j]->m_fPosC	= f11;	//pos relative ligne C
								f->m_pProfils[j]->m_fPosD	= f12;	//pos relative ligne D
								f->m_pProfils[j]->m_fPosE	= f13;	//pos relative ligne E
						}//boucle lecture tableau de forme

						break;

					default:;

				}//switch mot clef

			}//test trouve mot clef

		}//boucle

		

		/**** fermeture fichier ****/

		if(fclose(fid))

		{

			printf("\nProbleme � la fermeture du fichier");

		}

	}//test ouverture fichier ok

	return f;

}

/***********************/

/* LectureFichierForme2 */

/***********************/

Forme* LectureFichierForme2(char* NomFic)

{

	int i,j,n;

	FILE *fid;

	Matrice *m=NULL;

	float noVer;

	Forme *f=NULL;

	//liste des mots clefs du fichier

/*
	char* texte[11]={"BORD_ATTAQUE_FORME_COURBE","BORD_DE_FUITE_FORME_COURBE","SUSPENTAGE_LIGNE_A_COURBE",

		"SUSPENTAGE_LIGNE_B_COURBE","SUSPENTAGE_LIGNE_C_COURBE","SUSPENTAGE_LIGNE_D_COURBE", "SUSPENTAGE_LIGNE_E_COURBE",

		"DIEDRE_COURBE","MORPHING_COURBE","VRILLAGE_COURBE","EPAISSEUR_RELATIVE_COURBE"};

	char* texte2[11]={"BORD_ATTAQUE_FORME_CTRL","BORD_DE_FUITE_FORME_CTRL","SUSPENTAGE_LIGNE_A_CTRL",

		"SUSPENTAGE_LIGNE_B_CTRL","SUSPENTAGE_LIGNE_C_CTRL","SUSPENTAGE_LIGNE_D_CTRL", "SUSPENTAGE_LIGNE_E_CTRL",

		"DIEDRE_CTRL","MORPHING_CTRL","VRILLAGE_CTRL","EPAISSEUR_RELATIVE_CTRL"};

*/

	char* motsClef[29] = {"VERSION",

		"NB_ALVEOLES", "COEFF_PROG_GEOM", "EPAI_REL_CENTRE",

		"BORD_ATTAQUE_FORME_COURBE","BORD_DE_FUITE_FORME_COURBE","SUSPENTAGE_LIGNE_A_COURBE",

		"SUSPENTAGE_LIGNE_B_COURBE","SUSPENTAGE_LIGNE_C_COURBE","SUSPENTAGE_LIGNE_D_COURBE", "SUSPENTAGE_LIGNE_E_COURBE",

		"DIEDRE_COURBE","MORPHING_COURBE","VRILLAGE_COURBE","EPAISSEUR_RELATIVE_COURBE",

		"BORD_ATTAQUE_FORME_CTRL","BORD_DE_FUITE_FORME_CTRL","SUSPENTAGE_LIGNE_A_CTRL",

		"SUSPENTAGE_LIGNE_B_CTRL","SUSPENTAGE_LIGNE_C_CTRL","SUSPENTAGE_LIGNE_D_CTRL", "SUSPENTAGE_LIGNE_E_CTRL",

		"DIEDRE_CTRL","MORPHING_CTRL","VRILLAGE_CTRL","EPAISSEUR_RELATIVE_CTRL",


		"PROFIL_CENTRE", "PROFIL_BOUT", "TABLEAU_FORME" };

	

	/**** message ****/

	printf("\nLecture fichier de Forme: '%s'",NomFic);

	

	/**** ouverture fichier en lecture ****/

	if( (fid = fopen( NomFic, "rt" )) == NULL )

	{

		printf( "\nErreur ouverture fichier '%s'", NomFic);

	}

	else

	{

		/**** allocation memoire ****/

		f = new Forme();

		f -> courbInput = true;

		int NbNerv = 0;

		/**** Lecture des points de controle ****/

		/*Boucle de lecture*/

		for (i = 0; i < 29; i++)

		{

			//recherche mot clef

			if(!TrouveMotDansFichierTexte(fid, motsClef[i]))

				printf("\n\tWarning: mot clef '%s' non trouv� !", motsClef[i]);

			else

			{

				//action en fonction du mot clef

				switch(i)

				{

				case 0: fscanf(fid,"%f", &noVer); 

					//if(noVer==NO_VERSION)

					break;

				case 1: 
					{ 
						fscanf(fid,"%d", &(f->NbCaiss));

						if(f->NbCaiss%2==0) 

								NbNerv=f->NbCaiss/2+1;

							else

								NbNerv=(f->NbCaiss+1)/2;

						break;

					}

				case 2: fscanf(fid,"%f", &(f->CoeffProgGeom)); break;

				case 3: fscanf(fid,"%f", &(f->EpaiRelCent)); break;

				case 4:

				case 5:

				case 6:

				case 7:

				case 8:

				case 9:

				case 10:

				case 11:

				case 12:

				case 13:

				case 14:
					
				case 15:

				case 16:

				case 17:

				case 18:

				case 19:

				case 20:

				case 21:

				case 22:

				case 23:

				case 24:

				case 25: {

					int NN = 0;

					if (i < 15) NN = NbNerv; else NN=4;

					m = new Matrice(NN,2);

					float f1, f2;

					for (int _i = 0; _i < NN; _i++) {
						fscanf(fid,"%f %f", &f1, &f2);
						m->SetElement(_i, 0, f1); m->SetElement(_i, 1, f2);
					}
					
					switch(i)

					{

					case 4: f->mCourbNez = m; break;

					case 5: f->mCourbFui = m; break;

					case 6: f->mCourbA = m; break;

					case 7: f->mCourbB = m; break;

					case 8: f->mCourbC = m; break;

					case 9: f->mCourbD = m; break;

					case 10: f->mCourbE = m; break;

					case 11: f->mCourbDiedre = m; break;

					case 12: f->mCourbMorphing = m; break;

					case 13: f->mCourbVrillage = m; break;

					case 14: f->mCourbEpaiRel = m; break;

					case 15: f->mCtrlNez = m; break;

					case 16: f->mCtrlFui = m; break;

					case 17: f->mCtrlA = m; break;

					case 18: f->mCtrlB = m; break;

					case 19: f->mCtrlC = m; break;

					case 20: f->mCtrlD = m; break;

					case 21: f->mCtrlE = m; break;

					case 22: f->mCtrlDiedre = m; break;

					case 23: f->mCtrlMorphing = m; break;

					case 24: f->mCtrlVrillage = m; break;

					case 25: f->mCtrlEpaiRel = m; break;

					default:;

					}

					break;}

				case 26: {
					char temp[1024];
							fscanf(fid,"%s", temp); 
							f->m_strNomProfilCent = temp;
							break;
						 }

				case 27: {
							char temp[1024];
							fscanf(fid,"%s", temp); 
							f->m_strNomProfilBout = temp;
							break;
						 }

				case 28:

					f->DeleteProfils();

					fscanf(fid,"\n%d",&n);

						//fscanf(fid,"%d",&n);

						f->AllocateProfils(n);

						//boucle � partir de la 1ere nervure du centre vers l'extr�mit�

						for (j=0; j<f->m_nbProfils; j++)
						{
							float f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13;

							fscanf(fid,"\n%d %f %f %f %f %f %f %f %f %f %f %f %f %f",
										  &n, //no Nervure
										  &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11, &f12, &f13 );

								f->m_pProfils[j]->m_fLength = f1;	//longueur
								f->m_pProfils[j]->m_fWidth  = f2;	//epaisseur relative
								f->m_pProfils[j]->m_fNezX	= f3;	//xnez
								f->m_pProfils[j]->m_fNezY	= f4;	//ynez
								f->m_pProfils[j]->m_fNezZ	= f5;	//znez
								f->m_pProfils[j]->m_fInclin	= f6 * DEG2RAD;	//inclinaison / horizontale
								f->m_pProfils[j]->m_fWash	= f7 * DEG2RAD;	//vrillage / horizontale
								f->m_pProfils[j]->m_fMorph	= f8;	//morphing / nervure centrale
								f->m_pProfils[j]->m_fPosA	= f9;	//pos relative ligne A
								f->m_pProfils[j]->m_fPosB	= f10;	//pos relative ligne B
								f->m_pProfils[j]->m_fPosC	= f11;	//pos relative ligne C
								f->m_pProfils[j]->m_fPosD	= f12;	//pos relative ligne D
								f->m_pProfils[j]->m_fPosE	= f13;	//pos relative ligne E
						}//boucle lecture tableau de forme

						break;

					default:;

				}//switch mot clef

			}//test trouve mot clef

		}//boucle

		

		/**** fermeture fichier ****/

		if(fclose(fid))

		{

			printf("\nProbleme � la fermeture du fichier");

		}

	}//test ouverture fichier ok

	return f;

}


/************************/

/* EcritureFichierForme */

/************************/

void EcritureFichierForme(char *NomFichier, Forme *f)

{

	int i;

	FILE *fid;

	char* texte[11]={"BORD_ATTAQUE_FORME","BORD_DE_FUITE_FORME","SUSPENTAGE_LIGNE_A",

		"SUSPENTAGE_LIGNE_B","SUSPENTAGE_LIGNE_C","SUSPENTAGE_LIGNE_D", "SUSPENTAGE_LIGNE_E",

		"DIEDRE","MORPHING","VRILLAGE","EPAISSEUR_RELATIVE"};

	Matrice* m[10];



	/**** message ****/

	printf("\nEcriture fichier de Forme: '%s'",NomFichier);



	/**** ouverture fichier en ecriture ****/

	if( (fid = fopen( NomFichier, "wt" )) == NULL )

	{

		printf( "\nErreur ouverture fichier '%s'", NomFichier);

		exit(0);

	}

	

	/**** Ecriture No de version ****/

	

	fprintf(fid,"\nVERSION %1.2f", NO_VERSION);

	

	/**** Ecriture des parametres ****/

	

	fprintf(fid,"\n\nNB_ALVEOLES %d", f->NbCaiss);

	fprintf(fid,"\nCOEFF_PROG_GEOM %1.5f", f->CoeffProgGeom);

	fprintf(fid,"\nEPAI_REL_CENTRE %2.3f", f->EpaiRelCent);

	

	/**** Ecriture des points de controle ****/

	

	/*init pointeur des matrices de coordonn�es des points de controle*/

	m [0]=f->mCtrlNez; m [1]=f->mCtrlFui;

	m [2]=f->mCtrlA; m [3]=f->mCtrlB; m [4]=f->mCtrlC; m [5]=f->mCtrlD;

	m [6]=f->mCtrlE;

	m [7]=f->mCtrlDiedre; m [8]=f->mCtrlMorphing;

	m [9]=f->mCtrlVrillage; m [10]=f->mCtrlEpaiRel;

//	m [6]=f->mCtrlDiedre; m [7]=f->mCtrlMorphing;

//	m [8]=f->mCtrlVrillage; m [9]=f->mCtrlEpaiRel;

	

	/*boucle ecriture des points de controle*/

	for (i=0; i<11; i++)

		fprintf(fid,"\n%s %1.3f %1.3f  %1.3f %1.3f  %1.3f %1.3f  %1.3f %1.3f",

		texte[i],

		m [i]->Element(0,0), m [i]->Element(0,1),

		m [i]->Element(1,0), m [i]->Element(1,1),

		m [i]->Element(2,0), m [i]->Element(2,1),

		m [i]->Element(3,0), m [i]->Element(3,1));

	

	/**** Ecriture noms des profils ****/



	fprintf(fid,"\n\nPROFIL_CENTRE %s", f->m_strNomProfilCent.c_str());

	fprintf(fid,"\nPROFIL_BOUT %s", f->m_strNomProfilBout.c_str());



	/**** Ecriture du tableau de forme ****/

	

	fprintf(fid,"\n\nTABLEAU_FORME %d",f->m_nbProfils);

	/*boucle � partir de la 1ere nervure du centre vers l'extr�mit�*/

	for (i=0; i<f->m_nbProfils; i++)

	{

		fprintf(fid,

			"\n%d %1.3f %2.2f %2.3f %2.3f %2.3f %1.2f %1.2f %1.2f %3.1f %3.1f %3.1f %3.1f %3.1f",

			i+1,							//no Nervure
			f->m_pProfils[i]->m_fLength,
			f->m_pProfils[i]->m_fWidth,
			f->m_pProfils[i]->m_fNezX,	
			f->m_pProfils[i]->m_fNezY,	
			f->m_pProfils[i]->m_fNezZ,	
			f->m_pProfils[i]->m_fInclin*RAD2DEG,
			f->m_pProfils[i]->m_fWash*RAD2DEG,
			f->m_pProfils[i]->m_fMorph,
			f->m_pProfils[i]->m_fPosA,
			f->m_pProfils[i]->m_fPosB,	
			f->m_pProfils[i]->m_fPosC,	
			f->m_pProfils[i]->m_fPosD,
			f->m_pProfils[i]->m_fPosE
			);

	}

	

	/**** fermeture fichier ****/

	if(fclose(fid))

	{

		printf("\nProbleme � la fermeture du fichier");

		exit(0);

	}

}

/************************/

/* EcritureFichierForme2 */

/************************/

void EcritureFichierForme2(char *NomFichier, Forme *f)

{
	printf ("EcritureFichierForme2\n");

	int i;

	FILE *fid;

	char* texte[11]={"BORD_ATTAQUE_FORME_COURBE","BORD_DE_FUITE_FORME_COURBE","SUSPENTAGE_LIGNE_A_COURBE",

		"SUSPENTAGE_LIGNE_B_COURBE","SUSPENTAGE_LIGNE_C_COURBE","SUSPENTAGE_LIGNE_D_COURBE", "SUSPENTAGE_LIGNE_E_COURBE",

		"DIEDRE_COURBE","MORPHING_COURBE","VRILLAGE_COURBE","EPAISSEUR_RELATIVE_COURBE"};

	char* texte2[11]={"BORD_ATTAQUE_FORME_CTRL","BORD_DE_FUITE_FORME_CTRL","SUSPENTAGE_LIGNE_A_CTRL",

		"SUSPENTAGE_LIGNE_B_CTRL","SUSPENTAGE_LIGNE_C_CTRL","SUSPENTAGE_LIGNE_D_CTRL", "SUSPENTAGE_LIGNE_E_CTRL",

		"DIEDRE_CTRL","MORPHING_CTRL","VRILLAGE_CTRL","EPAISSEUR_RELATIVE_CTRL"};


	Matrice* m[10];

	Matrice* m2[10];



	/**** message ****/

	printf("\nEcriture fichier de Forme: '%s'",NomFichier);



	/**** ouverture fichier en ecriture ****/

	if( (fid = fopen( NomFichier, "wt" )) == NULL )

	{

		printf( "\nErreur ouverture fichier '%s'", NomFichier);

		exit(0);

	}

	

	/**** Ecriture No de version ****/

	

	fprintf(fid,"\nVERSION %1.2f", NO_VERSION);

	

	/**** Ecriture des parametres ****/

	

	fprintf(fid,"\n\nNB_ALVEOLES %d", f->NbCaiss);

	fprintf(fid,"\nCOEFF_PROG_GEOM %1.5f", f->CoeffProgGeom);

	fprintf(fid,"\nEPAI_REL_CENTRE %2.3f", f->EpaiRelCent);

	

	/**** Ecriture des points de controle ****/

	

	/*init pointeur des matrices de coordonn�es des points de controle*/

	m [0]=f->mCourbNez; m [1]=f->mCourbFui;

	m [2]=f->mCourbA; m [3]=f->mCourbB; m [4]=f->mCourbC; m [5]=f->mCourbD;	m [6]=f->mCourbE;

	m [7]=f->mCourbDiedre; m [8]=f->mCourbMorphing;

	m [9]=f->mCourbVrillage; m [10]=f->mCourbEpaiRel;

//	m [6]=f->mCtrlDiedre; m [7]=f->mCtrlMorphing;

//	m [8]=f->mCtrlVrillage; m [9]=f->mCtrlEpaiRel;

	m2 [0]=f->mCtrlNez; m2 [1]=f->mCtrlFui;

	m2 [2]=f->mCtrlA; m2 [3]=f->mCtrlB; m2 [4]=f->mCtrlC; m2 [5]=f->mCtrlD;	m2 [6]=f->mCtrlE;

	m2 [7]=f->mCtrlDiedre; m2 [8]=f->mCtrlMorphing;

	m2 [9]=f->mCtrlVrillage; m2 [10]=f->mCtrlEpaiRel;
	

	/*boucle ecriture des points de controle*/
	int NbNerv=0;

	if(f->NbCaiss%2==0) 

		NbNerv=f->NbCaiss/2+1;

	else

		NbNerv=(f->NbCaiss+1)/2;

	for (i=0; i<11; i++)
	{

		fprintf(fid,"\n%s", texte[i]);
		for (int j = 0; j < NbNerv; j++) {
			fprintf(fid, " %1.3f %1.3f", m [i]->Element(j,0), m [i]->Element(j,1));
		}
	}

	for (i=0; i<11; i++)
	{

		fprintf(fid,"\n%s", texte2[i]);
		for (int j = 0; j < 4; j++) {
			fprintf(fid, " %1.3f %1.3f", m2 [i]->Element(j,0), m2 [i]->Element(j,1));
		}
	}
	

	/**** Ecriture noms des profils ****/



	fprintf(fid,"\n\nPROFIL_CENTRE %s", f->m_strNomProfilCent.c_str());

	fprintf(fid,"\nPROFIL_BOUT %s", f->m_strNomProfilBout.c_str());



	/**** Ecriture du tableau de forme ****/

	

	fprintf(fid,"\n\nTABLEAU_FORME %d",f->m_nbProfils);

	/*boucle � partir de la 1ere nervure du centre vers l'extr�mit�*/

	for (i=0; i<f->m_nbProfils; i++)

	{

		fprintf(fid,

			"\n%d %1.3f %2.2f %2.3f %2.3f %2.3f %1.2f %1.2f %1.2f %3.1f %3.1f %3.1f %3.1f %3.1f",

			i+1,							//no Nervure
			f->m_pProfils[i]->m_fLength,
			f->m_pProfils[i]->m_fWidth,
			f->m_pProfils[i]->m_fNezX,	
			f->m_pProfils[i]->m_fNezY,	
			f->m_pProfils[i]->m_fNezZ,	
			f->m_pProfils[i]->m_fInclin*RAD2DEG,
			f->m_pProfils[i]->m_fWash*RAD2DEG,
			f->m_pProfils[i]->m_fMorph,
			f->m_pProfils[i]->m_fPosA,
			f->m_pProfils[i]->m_fPosB,	
			f->m_pProfils[i]->m_fPosC,	
			f->m_pProfils[i]->m_fPosD,
			f->m_pProfils[i]->m_fPosE
			);

	}

	

	/**** fermeture fichier ****/

	if(fclose(fid))

	{

		printf("\nProbleme � la fermeture du fichier");

		exit(0);

	}

}


/********************/

/* CalculInfoForme */

/********************/

void CalculInfoForme( Forme* F, TInfoForme* info )

{

	float xNerv, xNervPrec, lNerv, lNervPrec, larg;

	int i, nNerv;

	Matrice *tabXNerv = NULL;



	/*calcul position nervures a plat*/

	nNerv = F->m_nbProfils;

	tabXNerv = new Matrice(nNerv,1);

	tabXNerv->SetElement(0,0, F->m_pProfils[0]->m_fNezX);

	for(i=1; i<nNerv; i++)
	{
		tabXNerv->SetElement(i,0,
			tabXNerv->Element(i-1,0)
			+ (float)sqrt(
				sqr(F->m_pProfils[i]->m_fNezX - F->m_pProfils[i-1]->m_fNezX)
				+ sqr(F->m_pProfils[i]->m_fNezY - F->m_pProfils[i-1]->m_fNezY)
			)
		);
	}

	/*envergure*/

	info->envergure =  2.0f * tabXNerv->Element(nNerv-1,0);

	/*surface*/

	info->surface = 0.0f; xNervPrec = 0.0f;

	lNervPrec = F->m_pProfils[0]->m_fLength;
	for(i=0; i<nNerv; i++)
	{
		lNerv = F->m_pProfils[i]->m_fLength;
		xNerv = tabXNerv->Element(i,0);
		info->surface += (lNerv + lNervPrec) * (xNerv - xNervPrec);
		lNervPrec = lNerv; xNervPrec = xNerv;
	}

	/*allongement*/

	info->allongement=sqr(info->envergure)/info->surface;



	/*envergure projet�e*/

	info->envergureProj =  2.0f * F->m_pProfils[nNerv-1]->m_fNezX;

	/*surface projet�e*/

	info->surfaceProj = 0.0f; xNervPrec = 0.0f;

	lNervPrec = F->m_pProfils[0]->m_fLength;
	for(i=0; i<nNerv; i++)
	{
		lNerv = F->m_pProfils[i]->m_fLength;
		xNerv = F->m_pProfils[i]->m_fNezX;
		info->surfaceProj += (lNerv + lNervPrec) * (xNerv - xNervPrec);
		lNervPrec = lNerv; xNervPrec = xNerv;
	}

	/*allongement projet�*/

	info->allongementProj=sqr(info->envergureProj)/info->surfaceProj;



	/*corde max et min*/

	info->cordeMin = 100000.0f; info->cordeMax = 0.0f;

	for(i=0; i<nNerv; i++)
	{
		lNerv = F->m_pProfils[i]->m_fLength;
		if(info->cordeMin > lNerv) info->cordeMin = lNerv;
		if(info->cordeMax < lNerv) info->cordeMax = lNerv;
	}


	/*largeur caisson min et max*/

	info->largMin = 10000.0f;

	info->largMax = tabXNerv->Element(0,0) * 2.0f; //par defaut = largeur en central

	for(i=1; i<nNerv; i++)

	{

		larg = tabXNerv->Element(i,0) - tabXNerv->Element(i-1,0);

		if(info->largMin > larg) info->largMin = larg;

		if(info->largMax < larg) info->largMax = larg;

	}



	/*liberation memoire*/

	if ( tabXNerv != NULL )
		delete tabXNerv;
}



/********************/

/* AfficheInfoForme */

/********************/

void AfficheInfoForme( TInfoForme info )

{

	printf("\n\nINFO:");

	printf("\n\tSurface: plat=%2.2fm2, proj=%2.2fm2, r=%2.1f%%",

		info.surface, info.surfaceProj, info.surfaceProj/info.surface*100.0f);

	printf("\n\tEnvergure: plat=%2.2fm, proj=%2.2fm, r=%2.1f%%",

		info.envergure, info.envergureProj, info.envergureProj/info.envergure*100.0f);

	printf("\n\tAllongement: plat=%2.2f, proj=%2.2f, r=%2.1f%%",

		info.allongement, info.allongementProj, info.allongementProj/info.allongement*100.0f);

	printf("\n\tCorde mini = %2.3fm, maxi = %2.3fm", info.cordeMin, info.cordeMax);

	printf("\n\tLarg caisson mini = %2.3fm, maxi = %2.3fm", info.largMin, info.largMax);

}



/*****************/

/* EcritLigneDXF */

/*****************/

void EcritLigneDXF(FILE *fid, char *nom, int en3d,

				   float x1, float y1, float z1,

				   float x2, float y2, float z2)

{

	fprintf(fid,"0\nLINE\n8\n%s\n10\n%3.3f\n20\n%3.3f\n", nom, x1*1000.0f, y1*1000.0f);

	if(en3d) fprintf(fid,"30\n%3.3f\n", z1*1000.0f);

	fprintf(fid,"11\n%3.3f\n21\n%3.3f\n", x2*1000.0f, y2*1000.0f);

	if(en3d) fprintf(fid,"31\n%3.3f\n", z2*1000.0f);

}
/*****************/

/* EcritPolyBeginDXF */

/*****************/

void EcritPolyBeginDXF(FILE *fid, char *nom)

{

	fprintf(fid,"0\nPOLYLINE\n8\n%s\n66\n1\n70\n1\n", nom);

}
/*****************/

/* EcritPolyEndDXF */

/*****************/

void EcritPolyEndDXF(FILE *fid, char *nom)

{

	fprintf(fid,"0\nSEQEND\n", nom);

}


/*****************/

/* EcritPolyVertexDXF */

/*****************/

void EcritPolyVertexDXF(FILE *fid, char *nom, int en3d,

				   float x1, float y1, float z1)

{

	fprintf(fid,"0\nVERTEX\n8\n%s\n10\n%3.3f\n20\n%3.3f\n", nom, x1*1000.0f, y1*1000.0f);

	if(en3d) fprintf(fid,"30\n%3.3f\n", z1*1000.0f);

}

/**********************/
/* EcritureFichierFGen */
/**********************/
void EcritureFichierFGen(char *NomFichier, Forme *foil)
{
	ofstream out(NomFichier,ios::out);
	if ( out.good() )
	{
		out << dec << fixed << setfill('0');

		out 
			<< "scale = 1.0" << endl
			<< "foil = FROM_WINDWAX" << endl
			<< "(" << endl
			<< "\tinclude_profile ( USER \"USER.profile.fgen\" )" << endl
			;

		for( int i = 0; i < foil->m_nbProfils; i++ )
		{
			Profil *current = foil->m_pProfils[i];
			float height = current->m_fWidth;
			if ( height == 0.0f )
				height = 0.1f;
			height /= 100.0f;

			float a = current->m_fInclin;

			out 
				<< "\trib" << endl
				<< "\t(" << endl
				<< "\t\talias = Rib_" << setw(2) << i+1 << endl
				<< "\t\tmodel = USER" << endl
				<< "\t\tposition = ( (" 
						<<			current->m_fNezX * 1000
						<< ", " <<	-current->m_fNezZ * 1000
						<< ", " <<	current->m_fNezY * 1000
					<< "), ("
						<<			current->m_fNezX * 1000
						<< ", " <<	(-current->m_fNezZ + current->m_fLength )* 1000
						<< ", " <<	current->m_fNezY * 1000
					<< ") )" << endl
				<< "\t\torientation = ( (" 
						<<			current->m_fNezX * 1000
						<< ", " <<	-current->m_fNezZ * 1000
						<< ", " <<	current->m_fNezY * 1000
					<< "), ("
						<<			(current->m_fNezX + height*cos(a))* 1000 // TOCHECK
						<< ", " <<	(-current->m_fNezZ)* 1000
						<< ", " <<	(current->m_fNezY + height*sin(a)) * 1000 // TOCHECK
					<< ") )" << endl
				<< "\t\taoa = " << current->m_fWash*RAD2DEG << endl
				<< "\t\theight = " << height << endl
				<< "\t)" << endl
				;
		}

		out 
			<< ")" << endl // foil
			;

	}
	out.close();
}


/**********************/

/* EcritureFichierDXF */

/**********************/

void EcritureFichierDXF(char *NomFichier, TAxe *axe)

{

	int i, j, cpt;

	FILE *fid;

	Courbe *courb;

	TMesh *mesh;

	char text[50];



	/**** message ****/

	printf("\nEcriture fichier de DXF: '%s'",NomFichier);



	/**** ouverture fichier en ecriture ****/

	if( (fid = fopen( NomFichier, "wt" )) == NULL )

	{

		printf( "\nErreur ouverture fichier '%s'", NomFichier);

		exit(0);

	}

	

	/**** ecriture DXF *****/



	/* entete DXF minimal !!! */

    fprintf(fid,"0\nSECTION\n2\nENTITIES\n");



	/*ecriture courbes*/

	courb=axe->Courb; cpt=0;

	while (courb!=NULL)

	{

		printf ("\n courb [%s]", courb->name.c_str());
/*
		if ((courb->name.find("Marge") != string::npos) || (courb->name.find("Coin") != string::npos)) {
			sprintf(text, "Ext");
		} else 

		if (courb->name.find("suspentes") != string::npos) {
			sprintf(text, "X");
		} else 
		if ((courb->name.find("AR") != string::npos) || (courb->name.find("Avant") != string::npos)  || (courb->name.find("Patron") != string::npos)) {
			sprintf(text, "Int");
		} else 
		if (courb->name.length() == 0) {
			sprintf(text, "Text");
		} else {
			sprintf(text, "%s", "Other");
		} 
		*/

		cpt++; //sprintf(text, "COURB%03d", cpt);
		//sprintf(text, "COURB%03d%s", cpt, courb->name.c_str());

		sprintf(text, "%s", courb->name.c_str());

		//text = layerName;
		for(i=0; i<courb->pts->GetLignes()-1; i++)

		{

			EcritLigneDXF(fid, text, axe->axe3d,

				courb->pts->Element(i,0), courb->pts->Element(i,1), courb->pts->Element(i,2),

				courb->pts->Element(i+1,0), courb->pts->Element(i+1,1), courb->pts->Element(i+1,2));

		}

		courb=courb->CourbSuiv;

	}



	/*ecriture meshs*/

	mesh=axe->Mesh; cpt=0;

	while (mesh!=NULL)

	{

		cpt++; sprintf(text, "MESH%03d", cpt);

		for(i=0; i<mesh->x->GetLignes()-1; i++)

			for(j=0; j<mesh->x->GetColonnes()-1; j++)

			{

				EcritLigneDXF(fid, text, axe->axe3d,

					mesh->x->Element(i,j), mesh->y->Element(i,j), mesh->z->Element(i,j),

					mesh->x->Element(i,j+1), mesh->y->Element(i,j+1), mesh->z->Element(i,j+1));

				EcritLigneDXF(fid, text, axe->axe3d,

					mesh->x->Element(i,j), mesh->y->Element(i,j), mesh->z->Element(i,j),

					mesh->x->Element(i+1,j), mesh->y->Element(i+1,j), mesh->z->Element(i+1,j));

				if (i==mesh->x->GetLignes()-2)

					EcritLigneDXF(fid, text, axe->axe3d,

					mesh->x->Element(i+1,j), mesh->y->Element(i+1,j), mesh->z->Element(i+1,j),

					mesh->x->Element(i+1,j+1), mesh->y->Element(i+1,j+1), mesh->z->Element(i+1,j+1));

				if (j==mesh->x->GetColonnes()-2)

					EcritLigneDXF(fid, text, axe->axe3d,

					mesh->x->Element(i,j+1), mesh->y->Element(i,j+1), mesh->z->Element(i,j+1),

					mesh->x->Element(i+1,j+1), mesh->y->Element(i+1,j+1), mesh->z->Element(i+1,j+1));

			}

			mesh=mesh->MeshSuiv;

	}

	



	/*fin DXF*/

    fprintf(fid,"0\nENDSEC\n0\nEOF\n");





	/**** fermeture fichier ****/

	if(fclose(fid))

	{

		printf("\nProbleme � la fermeture du fichier");

		exit(0);

	}

}

//EcritureManyFichierPolyDXF(PtrNomFichier, AxePatronDXF, AxeMarginDXF, 1, AxeRepDXF, 0, AxeCercleDXF, Numerotation, AxePatronTextDXF, W, H);

void EcritureManyFichierPolyDXF(char *NomFichier, int n, TAxe **axe, TAxe **axe2, int rep, TAxe **axeR, int vent, TAxe **axeC, int num, TAxe **axeT, float* W, float* H)
{
    	FILE *fid;
	printf("\nEcriture POLY fichier de DXF: '%s'",NomFichier);
	if( (fid = fopen( NomFichier, "wt" )) == NULL ) {
		printf( "\nErreur ouverture fichier '%s'", NomFichier);
		exit(0);
	}
        fprintf(fid,"0\nSECTION\n2\nENTITIES\n");
        float dx=0.0f,dy=0.0f,dz=0.0f;
        for (int i=0; i<n;i++) {
            EcritureFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], vent, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*1.1f;
        }
       fprintf(fid,"0\nENDSEC\n0\nEOF\n");
	if(fclose(fid))	{
		printf("\nProbleme � la fermeture du fichier");
		exit(0);
	}
}

/**********************/

/* EcritureFichierPolyDXF */

/**********************/

void EcritureFichierPolyDXF(char *NomFichier, TAxe *axe, TAxe *axe2, int rep, TAxe *axeR, int vent, TAxe *axeC, int num, TAxe *axeT )
{
	int i, j, cpt;
	float fx, fy, fz;
	FILE *fid;
	Courbe *courb;
        TMesh *mesh;
	char text[50];
	/**** message ****/
	printf("\nEcriture POLY fichier de DXF: '%s'",NomFichier);
	/**** ouverture fichier en ecriture ****/
	if( (fid = fopen( NomFichier, "wt" )) == NULL )
	{
		printf( "\nErreur ouverture fichier '%s'", NomFichier);
		exit(0);
	}
	/**** ecriture DXF *****/
	/* entete DXF minimal !!! */
        fprintf(fid,"0\nSECTION\n2\nENTITIES\n");
	float curx, cury, x1, y1;
	/*ecriture courbes*/
	for (int _i=0; _i < 2; _i++){
		curx = -1000000.0f;
		cury = -2000000.0f;
		x1 = -3000000.0f;
		y1 = -4000000.0f;
		if (_i == 0) { 
			courb=axe->Courb;
			sprintf(text, "Patron");
		}
		if (_i == 1) {
			courb=axe2->Courb;
			sprintf(text, "Marge");
		}
		cpt=0;
		EcritPolyBeginDXF(fid, text);
		int q = 0;
		while (courb != NULL)
		{
//			printf ("\n while in ecritureFichierPolyDXF");
//
//			printf ("\n courb [%s]", courb->name.c_str());

			if (q == 0) {
				x1 = courb->pts->Element(0,0);
				y1 = courb->pts->Element(0,1);
				//printf ("XY1 (%f, %f)", x1, y1);
			}
	
/*			if ((courb->name.find("Marge") != string::npos) || (courb->name.find("Coin") != string::npos)) {
				sprintf(text, "Ext");
			} else 

			if (courb->name.find("suspentes") != string::npos) {
				sprintf(text, "X");
			} else 
			if ((courb->name.find("AR") != string::npos) || (courb->name.find("Avant") != string::npos)  || (courb->name.find("Patron") != string::npos)) {
				sprintf(text, "Int");
			} else 
			if (courb->name.length() == 0) {
				sprintf(text, "Text");
			} else {
				sprintf(text, "%s", courb->name.c_str());
			} */
			cpt++; 
                 	sprintf(text, "%s", courb->name.c_str());
			//text = layerName;
			for(i = 0; i < courb->pts->GetLignes(); i++)
			{
				fx = courb->pts->Element(i,0);
				fy = courb->pts->Element(i,1);
				fz = courb->pts->Element(i,2);
				//printf ("\n%s : f(%f, %f) cur (%f, %f) %d", courb->name.c_str(), fx, fy, curx, cury, q);
				if ( (q==0) || !((fx == curx) && (fy == cury)) ) {
					if ((q==0) || !((fx == x1) && (fy == y1)) ) {
						EcritPolyVertexDXF(fid, text, axe->axe3d, fx, fy, fz);
						//printf (" !print!");
						curx = fx;
						cury = fy;
					}
				}
				q++;
			}
			courb=courb->CourbSuiv;
		}
		EcritPolyEndDXF(fid, text);
	}

	if (vent == 1) {
		courb=axeC->Courb;
		sprintf(text, "Vent");
		while (courb != NULL)
		{
			EcritPolyBeginDXF(fid, text);
			for(i = 0; i < courb->pts->GetLignes(); i++)
			{
				fx = courb->pts->Element(i,0);
				fy = courb->pts->Element(i,1);
				fz = courb->pts->Element(i,2);
				//printf ("\n%s : f(%f, %f) cur (%f, %f) %d", courb->name.c_str(), fx, fy, curx, cury, q);
				EcritPolyVertexDXF(fid, text, axe->axe3d, fx, fy, fz);
			}
			EcritPolyEndDXF(fid, text);
			courb=courb->CourbSuiv;
		}
	} 



	if (rep) {
		courb=axeR->Courb; cpt=0;
		while (courb != NULL)
		{
			//printf ("\n courb [%s]", courb->name.c_str());
			cpt++; 
			sprintf(text, "Rep");
			for(i=0; i<courb->pts->GetLignes()-1; i++)
			{
				EcritLigneDXF(fid, text, axe->axe3d,
					courb->pts->Element(i,0), courb->pts->Element(i,1), courb->pts->Element(i,2),
					courb->pts->Element(i+1,0), courb->pts->Element(i+1,1), courb->pts->Element(i+1,2));
			}
			courb=courb->CourbSuiv;
		}
	}
	if (num == 1) {
		courb=axeT->Courb; cpt=0;
		while (courb != NULL)
		{
			//printf ("\n courb [%s]", courb->name.c_str());
			cpt++; 
			sprintf(text, "Text");
			for(i=0; i<courb->pts->GetLignes()-1; i++)
			{
				EcritLigneDXF(fid, text, axe->axe3d,
					courb->pts->Element(i,0), courb->pts->Element(i,1), courb->pts->Element(i,2),
					courb->pts->Element(i+1,0), courb->pts->Element(i+1,1), courb->pts->Element(i+1,2));
			}
			courb=courb->CourbSuiv;
		}
	}

	/*ecriture meshs*/
	mesh=axe->Mesh; cpt=0;
	while (mesh!=NULL)
	{
		cpt++; sprintf(text, "MESH%03d", cpt);
		for(i=0; i<mesh->x->GetLignes()-1; i++)
			for(j=0; j<mesh->x->GetColonnes()-1; j++)
			{
				EcritLigneDXF(fid, text, axe->axe3d,
					mesh->x->Element(i,j), mesh->y->Element(i,j), mesh->z->Element(i,j),
					mesh->x->Element(i,j+1), mesh->y->Element(i,j+1), mesh->z->Element(i,j+1));
				EcritLigneDXF(fid, text, axe->axe3d,
					mesh->x->Element(i,j), mesh->y->Element(i,j), mesh->z->Element(i,j),
					mesh->x->Element(i+1,j), mesh->y->Element(i+1,j), mesh->z->Element(i+1,j));
				if (i==mesh->x->GetLignes()-2)
					EcritLigneDXF(fid, text, axe->axe3d,
					mesh->x->Element(i+1,j), mesh->y->Element(i+1,j), mesh->z->Element(i+1,j),
					mesh->x->Element(i+1,j+1), mesh->y->Element(i+1,j+1), mesh->z->Element(i+1,j+1));
				if (j==mesh->x->GetColonnes()-2)
					EcritLigneDXF(fid, text, axe->axe3d,
					mesh->x->Element(i,j+1), mesh->y->Element(i,j+1), mesh->z->Element(i,j+1),
					mesh->x->Element(i+1,j+1), mesh->y->Element(i+1,j+1), mesh->z->Element(i+1,j+1));
			}
			mesh=mesh->MeshSuiv;
	}
	/*fin DXF*/
    fprintf(fid,"0\nENDSEC\n0\nEOF\n");
	/**** fermeture fichier ****/
	if(fclose(fid))
	{
		printf("\nProbleme � la fermeture du fichier");
		exit(0);
	}
}


/**********************/

/* EcritureFichierPolyDXFDelta */

/**********************/

void EcritureFichierPolyDXFDelta(FILE *fid, TAxe *axe, TAxe *axe2, int rep, TAxe *axeR, int vent, TAxe *axeC, int num, TAxe *axeT, float dx, float dy, float dz, int n )
{
	int i, j, cpt;
	float fx, fy, fz;
	Courbe *courb;
        TMesh *mesh;
	char text[50];
	float curx, cury, x1, y1;
	for (int _i=0; _i < 2; _i++){
		curx = -1000000.0f;
		cury = -2000000.0f;
		x1 = -3000000.0f;
		y1 = -4000000.0f;
		if (_i == 0) {
			courb=axe->Courb;
			sprintf(text, "Patron%d", n);
		}
		if (_i == 1) {
			courb=axe2->Courb;
			sprintf(text, "Marge%d", n);
		}
		cpt=0;
		EcritPolyBeginDXF(fid, text);
		int q = 0;
		while (courb != NULL)
		{
			if (q == 0) {
				x1 = courb->pts->Element(0,0);
				y1 = courb->pts->Element(0,1);
				//printf ("XY1 (%f, %f)", x1, y1);
			}
			cpt++;
                 	sprintf(text, "%s", courb->name.c_str());
			//text = layerName;
			for(i = 0; i < courb->pts->GetLignes(); i++)
			{
				fx = courb->pts->Element(i,0);
				fy = courb->pts->Element(i,1);
				fz = courb->pts->Element(i,2);
				//printf ("\n%s : f(%f, %f) cur (%f, %f) %d", courb->name.c_str(), fx, fy, curx, cury, q);
				if ( (q==0) || !((fx == curx) && (fy == cury)) ) {
					if ((q==0) || !((fx == x1) && (fy == y1)) ) {
						EcritPolyVertexDXF(fid, text, axe->axe3d, dx+fx, dy + fy, fz + dz);
						//printf (" !print!");
						curx = fx;
						cury = fy;
					}
				}
				q++;
			}
			courb=courb->CourbSuiv;
		}
		EcritPolyEndDXF(fid, text);
	}

	if (vent == 1) {
		courb=axeC->Courb;
		sprintf(text, "Vent%d", n);
		while (courb != NULL)
		{
			EcritPolyBeginDXF(fid, text);
			for(i = 0; i < courb->pts->GetLignes(); i++)
			{
				fx = courb->pts->Element(i,0);
				fy = courb->pts->Element(i,1);
				fz = courb->pts->Element(i,2);
				EcritPolyVertexDXF(fid, text, axe->axe3d, dx+fx , dy+fy, dz+fz);
			}
			EcritPolyEndDXF(fid, text);
			courb=courb->CourbSuiv;
		}
	}

	if (rep) {
		courb=axeR->Courb; cpt=0;
		while (courb != NULL)
		{
			cpt++;
			sprintf(text, "Rep%d", n);
			for(i=0; i<courb->pts->GetLignes()-1; i++)
			{
				EcritLigneDXF(fid, text, axe->axe3d,
					dx+courb->pts->Element(i,0), dy+courb->pts->Element(i,1), dz+courb->pts->Element(i,2),
					dx+courb->pts->Element(i+1,0), dy+courb->pts->Element(i+1,1), dz+courb->pts->Element(i+1,2));
			}
			courb=courb->CourbSuiv;
		}
	}
	if (num == 1) {
		courb=axeT->Courb; cpt=0;
		while (courb != NULL)
		{
			cpt++;
			sprintf(text, "Text%d", n);
			for(i=0; i<courb->pts->GetLignes()-1; i++)
			{
				EcritLigneDXF(fid, text, axe->axe3d,
					dx+courb->pts->Element(i,0), dy+courb->pts->Element(i,1), dz+courb->pts->Element(i,2),
					dx+courb->pts->Element(i+1,0), dy+courb->pts->Element(i+1,1), dz+courb->pts->Element(i+1,2));
			}
			courb=courb->CourbSuiv;
		}
	}

	/*ecriture meshs*/
	mesh=axe->Mesh; cpt=0;
	while (mesh!=NULL)
	{
		cpt++; sprintf(text, "MESH%03d", cpt);
		for(i=0; i<mesh->x->GetLignes()-1; i++)
			for(j=0; j<mesh->x->GetColonnes()-1; j++)
			{
				EcritLigneDXF(fid, text, axe->axe3d,
					dx+mesh->x->Element(i,j), dy+mesh->y->Element(i,j), dz+mesh->z->Element(i,j),
					dx+mesh->x->Element(i,j+1), dy+mesh->y->Element(i,j+1), dz+mesh->z->Element(i,j+1));
				EcritLigneDXF(fid, text, axe->axe3d,
					dx+mesh->x->Element(i,j), dy+mesh->y->Element(i,j), dz+mesh->z->Element(i,j),
					dx+mesh->x->Element(i+1,j), dy+mesh->y->Element(i+1,j), dz+mesh->z->Element(i+1,j));
				if (i==mesh->x->GetLignes()-2)
					EcritLigneDXF(fid, text, axe->axe3d,
					dx+mesh->x->Element(i+1,j), dy+mesh->y->Element(i+1,j), dz+mesh->z->Element(i+1,j),
					dx+mesh->x->Element(i+1,j+1), dy+mesh->y->Element(i+1,j+1),dz+ mesh->z->Element(i+1,j+1));
				if (j==mesh->x->GetColonnes()-2)
					EcritLigneDXF(fid, text, axe->axe3d,
					dx+mesh->x->Element(i,j+1), dy+mesh->y->Element(i,j+1), dz+mesh->z->Element(i,j+1),
					dx+mesh->x->Element(i+1,j+1), dy+mesh->y->Element(i+1,j+1), dz+mesh->z->Element(i+1,j+1));
			}
			mesh=mesh->MeshSuiv;
	}
}