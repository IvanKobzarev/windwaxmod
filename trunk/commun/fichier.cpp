/**************************************
* gestion chargement divers fichiers
* forme ,profil, patrons ...
***************************************/

//#include <afx.h>		//class CString
//#include <afxdlgs.h>	//class CFileDialog

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <conio.h>
#include <math.h>
#include <fstream>
#include <iomanip>

//using namespace std;
#define sqr(f1) ((f1)*(f1))
#define CHISLOPI	3.141592675f

#ifndef DEBUG
#define DEBUG false
#endif

#include "C:\windwax\QFLR5\src\Design\AFoil.h"
#include "C:\windwax\QFLR5\src\Objects\Wing.h"
#include "C:\windwax\QFLR5\src\Objects\WOpp.h"
#include "C:\windwax\QFLR5\src\Objects\Sf.h"
#include "C:\windwax\QFLR5\src\Objects\Pf.h"

#include "fichier.h"
#include "plot.h"
#include "matrice.h"
#include "profil.h"
#include "geom.h"
#include "rasklad.h"
#include "patternsproject.h"


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
	CoeffProgGeom = 0.93f;
        CoeffExp = 0.0f;
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

void Forme::SerializeWing (QDataStream &ar) {
	printf ("\nin Forme::SerializeWing()");

	CWing* wing = new CWing ();

	wing -> m_WingName = "Kite";
	char foilNameL[50];
	char foilNameR[50];

	double LongNerv, EpaiRel, xp,yp, xo,yo,zo, a,v,m;



/*

    Matrice *ExtProfCentN, *ExtProfBoutN;

	double EpaiRelProfCent, EpaiRelProfBout;
	double coeffx, coeffyCent, coeffyBout;
	int i,j;


	EpaiRelProfCent = EpaisseurRelative(ExtProfCent, IntProfCent);
	EpaiRelProfBout = EpaisseurRelative(ExtProfBout, IntProfBout);


	for (i=0; i<forme->m_nbProfils; i++)
	{
		//longueur nervure courante
		LongNerv = forme->m_pProfils[i]->m_fLength;
		//epaisseur relative
		EpaiRel = forme->m_pProfils[i]->m_fWidth;
		//position xo,yo,zo du nez
		xo = forme->m_pProfils[i]->m_fNezX;
		yo = forme->m_pProfils[i]->m_fNezY;
		zo = forme->m_pProfils[i]->m_fNezZ;
                //printf ("\n %d long=%f width=%f (%f, %f, %f)", i, forme->m_pProfils[i]->m_fLength, forme->m_pProfils[i]->m_fWidth, xo, yo, zo);
		//inclinaison de la nervure par rapport a l'horizontale
		a = forme->m_pProfils[i]->m_fInclin;
        if ((isCenterPanel) && (i == 0)) a = forme->m_pProfils[0]->m_fInclin - 0.33333333f * fabs(forme->m_pProfils[0]->m_fInclin - forme->m_pProfils[1]->m_fInclin);
                
		//angle vrillage
		v = forme->m_pProfils[i]->m_fWash;
		//coeff morphing
		m = forme->m_pProfils[i]->m_fMorph;
		coeffx = LongNerv/100.0f;
		coeffyCent = LongNerv*EpaiRel/(EpaiRelProfCent*100.0f);
		coeffyBout = LongNerv*EpaiRel/(EpaiRelProfBout*100.0f);

		for (j=0; j<IntProfCent->GetLignes(); j++)
		{
			xp = IntProfCent->Element(j,0)*coeffx;
			yp = IntProfCent->Element(j,1)*coeffyCent*m
				+ IntProfBout->Element(j,1)*coeffyBout*(1.0f-m);
			(*XInt)->SetElement(i,j, xo+(yp-xp*(double)sin(v))*(double)cos(a));
			(*YInt)->SetElement(i,j, yo+(yp-xp*(double)sin(v))*(double)sin(a));
			(*ZInt)->SetElement(i,j, zo-(xp*(double)cos(v)));

		}

		for (j=0; j<ExtProfCent->GetLignes(); j++)
		{
            xp = ExtProfCent->Element(j,0)*coeffx;
            yp = ExtProfCent->Element(j,1)*coeffyCent*m
                    + ExtProfBout->Element(j,1)*coeffyBout*(1.0f-m);
			(*XExt)->SetElement(i,j, xo+(yp-xp*(double)sin(v))*(double)cos(a));
			(*YExt)->SetElement(i,j, yo+(yp-xp*(double)sin(v))*(double)sin(a));
			(*ZExt)->SetElement(i,j, zo-(xp*(double)cos(v)));
		}

	}


*/

        int isCenterPanel = (1 & NbCaiss);
		double positionSum = 0;
		double	xprev = m_pProfils[0]->m_fNezX;
		double	yprev = m_pProfils[0]->m_fNezY;

		wing->m_NPanel = isCenterPanel + m_nbProfils-1;	
		if (isCenterPanel) {
			LongNerv = m_pProfils[0]->m_fLength;
			EpaiRel = m_pProfils[0]->m_fWidth;

			xo = m_pProfils[0]->m_fNezX;
			yo = m_pProfils[0]->m_fNezY;
			zo = m_pProfils[0]->m_fNezZ;

			wing->m_RFoil[0]= "P000";
			wing->m_LFoil[0]= "P000";

			wing->m_TDihedral[0]=0;
			wing->m_TChord[0] = LongNerv;
			wing->m_TPos[0] = 0.0;

			wing->m_TOffset[0]=0.0;
			wing->m_TTwist[0]=0.0;

			wing->m_NXPanels[0]=2;
			wing->m_NYPanels[0]=2;
			wing->m_XPanelDist[0]=1;
			wing->m_YPanelDist[0]=0;

			xprev = 0;
		}

		double x1, y1, x2, y2;

		for (int i = 0; i < m_nbProfils; i++) {
			
			LongNerv = m_pProfils[i]->m_fLength;
			EpaiRel = m_pProfils[i]->m_fWidth;
			
			xo = m_pProfils[i]->m_fNezX;
			yo = m_pProfils[i]->m_fNezY;
			zo = m_pProfils[i]->m_fNezZ;
			positionSum += dist2d(xo, yo, xprev, yprev);
			xprev = xo;
			yprev = yo;

			wing->m_TDihedral[i+isCenterPanel]=0;
			v = m_pProfils[i]->m_fWash;
			printf ("\n vrillage=%f", (v*180/CHISLOPI));
			if (i < m_nbProfils-1) {
				x1 = m_pProfils[i]->m_fNezX;
				y1 = m_pProfils[i]->m_fNezY;

				x2 = m_pProfils[i+1]->m_fNezX;
				y2 = m_pProfils[i+1]->m_fNezY;
				wing->m_TDihedral[i+isCenterPanel] = -180/CHISLOPI*atan(abs(y1-y2)/abs(x1-x2));
			} 

			printf ("\ni=%d panel %d", i, (i+isCenterPanel));
			sprintf(foilNameL, "P%03d", i);
			sprintf(foilNameR, "P%03d", i);
			wing->m_RFoil[i+isCenterPanel]= foilNameL;
			wing->m_LFoil[i+isCenterPanel]= foilNameR;

			/*
			if (i < m_nbProfils-1) {
				wing->m_RFoil[i+isCenterPanel]= foilNameL;
				wing->m_LFoil[i+isCenterPanel]= foilNameR;
			} else {
				sprintf(foilNameL, "P%03d", i);
				sprintf(foilNameR, "P%03d", i);
				wing->m_RFoil[i+isCenterPanel]= foilNameL;
				wing->m_LFoil[i+isCenterPanel]= foilNameR;
			} 
			*/
			

			wing->m_TChord[i+isCenterPanel] = LongNerv;
			wing->m_TPos[i+isCenterPanel] = positionSum;
			//wing->m_Twist[i+isCenterPanel] = ;

			wing->m_TOffset[i+isCenterPanel]=m_pProfils[0]->m_fNezZ - zo;
			wing->m_TTwist[i+isCenterPanel]=v*180/CHISLOPI;

			wing->m_NXPanels[i+isCenterPanel]=2;
			wing->m_NYPanels[i+isCenterPanel]=2;
			wing->m_XPanelDist[i+isCenterPanel]=1;
			wing->m_YPanelDist[i+isCenterPanel]=0;
		}

		printf ("\n ...for");

		wing->m_bVLMAutoMesh = 1;
		wing->m_bSymetric = 1;
		printf ("\n in Forme::SerializeWing() ComputeGeometry()");
		//wing->ComputeGeometry();
		
		wing -> SerializeWing(ar, true, 5);
 		printf ("\n...SerializeWing()");
}



void Forme::SerializeToWpa(QDataStream &ar) 
{
	printf ("\nSerializeToWpa()");
	int m_LengthUnit;
	int m_AreaUnit;
	int m_WeightUnit;
	int m_SpeedUnit;
	int m_ForceUnit;
	int m_MomentUnit;
	
	m_LengthUnit  = 3;
	m_AreaUnit    = 3;
	m_WeightUnit  = 1;
	m_SpeedUnit   = 0;
	m_ForceUnit   = 0;
	m_MomentUnit  = 0;

	int m_Type=0;
	int m_UnitType=2;//1= International, 2= English
	int m_AnalysisType=1; //0=LLT;1=VLM;2=Panel
	int m_RefAreaType=1; //1=classic or 2=projected on x-yplane

	bool m_bVLM1=true; //true if Classic, false if Quendez
	bool m_bThinSurfaces=true;//true if Plane Panel calculation on middle surface, false if on top & bottom
	bool m_bWakeRollUp=true;//true if wake roll up is to be taken into account in calculation
	bool m_bTiltedGeom=true;//true if calculation is performed on the tilted geometry, at alpha=0.0
	bool m_bViscous=true;
    bool m_bGround=true;

	double m_QInf = 9.0, m_Weight = 1.0, m_Alpha = 2, m_XCmRef = 0.25;

	double m_CoG_x=0.25;
	double m_CoG_y=0;
	double m_CoG_z=0;

	double m_Beta = 0;
	double m_Density = 1.225, m_Viscosity = 1.5e-5;
	double m_WingLoad = 1;
	double m_Height = 0;

    int wingSize = 1;
    int polarSize = 0;
    int bodySize = 0;
    int planeSize = 0;

	ar << 100013;
	ar << m_LengthUnit;
	ar << m_AreaUnit;
	ar << m_WeightUnit;
	ar << m_SpeedUnit;
	ar << m_ForceUnit;
	ar << m_MomentUnit;

	ar << m_Type;
	ar << (float)m_Weight;
	ar << (float)m_QInf;
	ar << (float)m_CoG_x;
	ar << (float)m_CoG_y;
	ar << (float)m_CoG_z;

	ar << (float)m_Density;
	ar << (float)m_Viscosity;
	ar << (float)m_Alpha;
	ar << (float)m_Beta;
	ar << m_AnalysisType;

	if (m_bVLM1) ar << 1; else ar << 0;
	ar << 1;
	if (m_bTiltedGeom) ar << 1; else ar << 0;
	if (m_bWakeRollUp) ar << 1; else ar << 0;

	ar << 1;
    SerializeWing(ar);

	// now store all the WPolars
	ar << (int)polarSize;
	/*for (i=0; i<m_oaWPolar.GetSize();i++)
	{
		pWPolar = (CWPolar*)m_oaWPolar.GetAt(i);
		pWPolar->m_pParent = this;
		pWPolar->SerializeWPlr(ar);
	}*/

	// next store all the WOpps
	/*if(m_bSaveWOpps)
	{
		ar << (int)m_oaWOpp.GetSize();
		for (i=0; i<m_oaWOpp.GetSize();i++)
		{
			pWOpp = (CWOpp*)m_oaWOpp.GetAt(i);
			pWOpp->SerializeWOpp(ar);
		}
	}
	else */
    ar << 0;

	/*ar << 1;
	CWOpp* wopp = new CWOpp();
	wopp->m_WingName="Kite";
	wopp->SerializeWOpp(ar);*/

	// then the foils,  polars and Opps
		
	WritePolars(ar);

	// next the bodies
	ar << (int)bodySize;
	/*for (i=0; i<bodySize;i++)
	{
		pBody = (CBody*)m_oaBody.GetAt(i);
		pBody->SerializeBody(ar);
	}*/

	// last write the planes...
	ar << (int)planeSize;
	/*for (i=0; i<m_oaPlane.GetSize();i++)
	{
		pPlane = (CPlane*)m_oaPlane.GetAt(i);
		pPlane->SerializePlane(ar);
	}*/

	/* if(m_bSaveWOpps)
	{
		// not forgetting their POpps
		ar << (int)m_oaPOpp.GetSize();
		for (i=0; i<m_oaPOpp.GetSize();i++)
		{
			pPOpp = (CPOpp*)m_oaPOpp.GetAt(i);
			pPOpp->SerializePOpp(ar);
		}
	}
	else */

    ar << 0;

	CSF *m_pSF;
	CPF *m_pPF;

	m_pSF = new CSF();
	m_pSF->m_bModified = false;
	m_pSF->InitSplineFoil();

	m_pPF = new CPF();
	m_pPF->m_bModified = false;
	m_pPF->InitSplinedFoil();

	m_pSF->Serialize(ar, true);
	m_pPF->Serialize(ar, true);

	printf ("\n...SerializeToWpa()");
}

void EcritureFichierPolyDXFDelta(FILE *fid, TAxe *axe, TAxe *axe2, int rep, TAxe *axeR, int vent, TAxe *axeC, int num, TAxe *axeT, double dx, double dy, double dz, int n );


enum trim_type {LEFT=1, RIGHT, LEFT_AND_RIGHT};

void trim(char *str, char c=' ', trim_type t=LEFT_AND_RIGHT)
{
     int len = strlen(str);
     int beg=0;
     if(t & RIGHT)
        while(str[len-1]==c) str[--len]=0;
    if(t & LEFT)
        while(str[beg]==c) str[beg++]=0;
    if(beg) memmove(str, str+beg, len-beg+1);
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
                // printf ("\nf1=%f1 f2=%f2", f1, f2);
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
                //printf ("\n extr f1=%f1 f2=%f2", f1, f2);
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
                //printf ("\n intr f1=%f1 f2=%f2", f1, f2);
			}

		}

		

		//fermeture fichier

		if(fclose(fid))

		{

			printf("\nProbleme ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ la fermeture du fichier '%s'", NomProf);

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
//        if (DEBUG) printf ("\n...LectureFichierProfil()");

}
int* LectureFichierVentHoles(char* NomFic, int* quant, int* central) {
    FILE *fid;
    int* m = new int[1000];
    for (int i = 0; i < 1000; i++) m[i] = 0;
    int n;
    int nerv;
    (*central) = 0;
    if( (fid = fopen( NomFic, "rt" )) == NULL )
    {
            printf( "\nErreur ouverture fichier '%s'", NomFic);
    } else {
            fscanf(fid,"%d", &n);
            //printf ("\n n=%d", n);
            if (n != 0) {
                for (int i = 0; i < n; i++) {
                    fscanf (fid, "%d", &nerv);
                    //printf (" nerv %d", nerv);
                    if (nerv == -1) (*central) = 1; else m[nerv] = 1;
                }
            }
            (*quant) = n;
    }
    return m;
}

/***********************/
/* LectureFichierReperPoints */
/***********************/

Matrice** LectureFichierReperPoints(char* NomFic) {
    FILE *fid;
    Matrice** m= new Matrice*[3];
    int n;
    float p;
    if( (fid = fopen( NomFic, "rt" )) == NULL )
    {
            printf( "\nErreur ouverture fichier '%s'", NomFic);
    } else {
        for (int i=1;i<3;i++) {
            fscanf(fid,"%d", &n);
            //printf ("\n n=%d", n);
            m[i] = new Matrice(n, 1);
            for (int j = 0; j < n; j++) {
                fscanf (fid, "%f", &p);
                //printf ("\n p=%f", p);
                m[i]->SetElement(j,0,p);
            }
        }
    }
    return m;
}

int* LectureFichierDiagNervs(char* NomFic, int* quant) {
    FILE *fid;
    int* res;
    int n;
    if( (fid = fopen( NomFic, "rt" )) == NULL )
    {
            printf( "\nErreur ouverture fichier '%s'", NomFic);
    } else {

            fscanf(fid,"%d", &n);
            //printf ("\n n=%d", n);
            (*quant) = n;
            res = new int[n];
            for (int i = 0; i < n; i++) {
                fscanf (fid, "%d", &(res[i]));
                //printf (" dn=%d", res[i]);
            }

    }
    return res;
}
/*
 PROJECT f17v2

MARGE_DEB 1.0
MARGE_FIN 1.0
MARGE1 1.0
MARGE2 1.0

PINCES 1
PINCES_POS_BA 14.0
PINCES_AMP_BA 2.5
PINCES_POS_BF 25.0
PINCES_AMP_BF 2.0

PINCE_TYPE_NEZ 0
PINCE_RADIUS_ALG_K_NEZ 0.985
PINCE_SET_EQUAL_AMP_TO_FUI 1
PINCE_POWER_NEZ 0
PINCE_POWER_NEZ_VALUE 1.3
PINCE_ARCTAN_NEZ 0
PINCE_ARCTAN_NEZ_K1 0
PINCE_ARCTAN_NEZ_K2 0
PINCE_ARCTAN_NEZ_K3 0
 
PINCE_TYPE_FUI 1
PINCE_RADIUS_ALG_K_FUI 0.985
PINCE_SET_EQUAL_AMP_TO_NEZ 1
PINCE_POWER_FUI 1
PINCE_POWER_FUI_VALUE 1.3
PINCE_ARCTAN_FUI 0
PINCE_ARCTAN_FUI_K1 0
PINCE_ARCTAN_FUI_K2 0
PINCE_ARCTAN_FUI_K3 0



DIAG_NERVURES 1
DIAG_NERVURES_FILE diagnervs.txt
DIAG_NERVURES_POS_INT_A 8.0
DIAG_NERVURES_POS_INT_F 71.0
DIAG_NERVURES_POS_EXT_A 2.0
DIAG_NERVURES_POS_EXT_F 80.0

VENT_HOLES 1
VENT_HOLES_FILE ventholes.txt
VENT_HOLES_DOUBLE 1
VENT_HOLES_DEB 2.0
VENT_HOLES_FIN 6.0

VENT_HOLES_KLAPANS 1
VENT_HOLES_KLAPANS_FIN 16.0

REP_POINTS 1
REP_POINTS_FROM_FILE 1
REP_POINTS_FILE rep_points.txt

LAYOUT_VENTILATION 1
LAYOUT_MARGE_EXT_FIN 4.0
LAYOUT_ACCURACY 0.7
 */
WindPatternsProject* LectureWindPatternsProject(char* NomFic) {
    printf ("\nLectureWindPatternsProject(%s)", NomFic);
    FILE *fid;
    WindPatternsProject *wpp = NULL;
    char* motsClef[54] =
                    {"PROJECT", // 1
				"MARGE_DEB", "MARGE_FIN", "MARGE1", "MARGE2", //4
				"PINCES","PINCES_POS_BA","PINCES_AMP_BA","PINCES_POS_BF","PINCES_AMP_BF", //5

                "PINCE_TYPE_NEZ","PINCE_RADIUS_ALG_K_NEZ","PINCE_EQUAL_AMP_TO_FUI", // 9
                "PINCE_POWER_NEZ", "PINCE_POWER_NEZ_VALUE",
                "PINCE_ARCTAN_NEZ", "PINCE_ARCTAN_NEZ_K1", "PINCE_ARCTAN_NEZ_K2", "PINCE_ARCTAN_NEZ_K3",

                "PINCE_TYPE_FUI", "PINCE_RADIUS_ALG_K_FUI", "PINCE_EQUAL_AMP_TO_NEZ", // 9
                "PINCE_POWER_FUI", "PINCE_POWER_FUI_VALUE",
                "PINCE_ARCTAN_FUI", "PINCE_ARCTAN_FUI_K1", "PINCE_ARCTAN_FUI_K2", "PINCE_ARCTAN_FUI_K3",

				"DIAG_NERVURES","DIAG_NERVURES_FILE", // 6
                "DIAG_NERVURES_POS_INT_A","DIAG_NERVURES_POS_INT_F", "DIAG_NERVURES_POS_EXT_A","DIAG_NERVURES_POS_EXT_F",

				"VENT_HOLES", "VENT_HOLES_FILE", "VENT_HOLES_DEB", "VENT_HOLES_FIN", // 7
                "VENT_HOLES_KLAPANS", "VENT_HOLES_KLAPANS_DOUBLE", "VENT_HOLES_KLAPANS_FIN",

                "REP_POINTS", "REP_POINTS_FROM_FILE", "REP_POINTS_FILE", // 3
                "LAYOUT_VENTILATION",
                "LAYOUT_MARGE_EXT_FIN","LAYOUT_MARGE_INT_FIN", "LAYOUT_MARGE_NERV_FIN", "LAYOUT_MARGE_DIAG_NERV_FIN",
                "LAYOUT_ACCURACY",//3

                "XMASHTAB", "TEXTX", "TEXTY", "FORME"}; //3

    if( (fid = fopen( NomFic, "rt" )) == NULL )
    {
            printf( "\nErreur ouverture fichier '%s'", NomFic);
    }
    else
    {
        wpp = new WindPatternsProject();
        strcpy(wpp->fromFile, NomFic);
        char s[255];
        for (int i = 0; i < 54; i++) {
            //printf ("\n i=%d", i);
            if(!TrouveMotDansFichierTexte(fid, motsClef[i])) {
                printf("\n\tWarning: mot clef '%s' not in file!", motsClef[i]);
            } else {
                switch(i)
                {
                    case 0: {  fgets (s, 255, fid);  trim(s); trim(s, '\r'); trim(s, '\n'); strcpy(wpp->name, s); break;}
                    case 1: fscanf(fid,"%f", &(wpp->MargeDeb)); break;
                    case 2: fscanf(fid,"%f", &(wpp->MargeFin));  break;
                    case 3: fscanf(fid,"%f", &(wpp->Marge[0]));  break;
                    case 4: fscanf(fid,"%f", &(wpp->Marge[1])); break;

                    case 5: fscanf(fid,"%d", &(wpp->isPinces)); break;
                    case 6: fscanf(fid,"%f", &(wpp->PosPinceBA[0])); break;
                    case 7: fscanf(fid,"%f", &(wpp->AmpPinceBA[0])); break;
                    case 8: fscanf(fid,"%f", &(wpp->PosPinceBF[0])); break;
                    case 9: fscanf(fid,"%f", &(wpp->AmpPinceBF[0])); break;

                    case 10: fscanf(fid,"%d", &(wpp->PinceNosRadio)); break;
                    case 11: fscanf(fid,"%f", &(wpp->PinceRadiusAlgKNos)); break;
                    case 12: fscanf(fid,"%d", &(wpp->PinceNosEqualAmp)); break;
                    case 13: fscanf(fid,"%d", &(wpp->PincePowerA)); break;
                    case 14: fscanf(fid,"%f", &(wpp->PincePowerValueA)); break;
                    case 15: fscanf(fid,"%d", &(wpp->PinceArctanA)); break;
                    case 16: fscanf(fid,"%f", &(wpp->PinceArctanK1ValueA)); break;
                    case 17: fscanf(fid,"%f", &(wpp->PinceArctanK2ValueA)); break;
                    case 18: fscanf(fid,"%f", &(wpp->PinceArctanK3ValueA)); break;

                    case 19: fscanf(fid,"%d", &(wpp->PinceHvostRadio)); break;
                    case 20: fscanf(fid,"%f", &(wpp->PinceRadiusAlgKHvost)); break;
                    case 21: fscanf(fid,"%d", &(wpp->PinceHvostEqualAmp)); break;
                    case 22: fscanf(fid,"%d", &(wpp->PincePowerF)); break;
                    case 23: fscanf(fid,"%f", &(wpp->PincePowerValueF)); break;
                    case 24: fscanf(fid,"%d", &(wpp->PinceArctanF)); break;
                    case 25: fscanf(fid,"%f", &(wpp->PinceArctanK1ValueF)); break;
                    case 26: fscanf(fid,"%f", &(wpp->PinceArctanK2ValueF)); break;
                    case 27: fscanf(fid,"%f", &(wpp->PinceArctanK3ValueF)); break;

                    case 28: fscanf(fid,"%d", &(wpp->DiagNervs)); break;
                    case 29: {  fgets (s, 255, fid);  trim(s); trim(s, '\r'); trim(s, '\n'); strcpy(wpp->NomFichierDiagNerv, s);  break;}
                    case 30: fscanf(fid,"%f", &(wpp->PosDiagNerv2A)); break;
                    case 31: fscanf(fid,"%f", &(wpp->PosDiagNerv2F)); break;
                    case 32: fscanf(fid,"%f", &(wpp->PosDiagNerv1A)); break;
                    case 33: fscanf(fid,"%f", &(wpp->PosDiagNerv1F)); break;

                    case 34: fscanf(fid,"%d", &(wpp->VentHoles)); break;
                    case 35: {  fgets (s, 255, fid); trim(s); trim(s, '\r'); trim(s, '\n'); strcpy(wpp->NomFichierVentHoles, s); break;}
                    case 36: fscanf(fid,"%f", &(wpp->VentHolesDeb)); break;
                    case 37: fscanf(fid,"%f", &(wpp->VentHolesFin)); break;

                    case 38: fscanf(fid,"%d", &(wpp->RaskladKlapans)); break;
                    case 39: fscanf(fid,"%d", &(wpp->VentHolesDouble)); break;
                    case 40: fscanf(fid,"%f", &(wpp->PosKlapanFin)); break;

                    case 41: fscanf(fid,"%d", &(wpp->RepPoints)); break;
                    case 42: fscanf(fid,"%d", &(wpp->ReperPointsFromFile)); break;
                    case 43: {  fgets (s, 255, fid); trim(s); trim(s, '\r'); trim(s, '\n'); strcpy(wpp->NomFichierRepPoints, s);  break;}

                    case 44: fscanf(fid,"%d", &(wpp->VentilationLayout)); break;

                    case 45: fscanf(fid,"%f", &(wpp->margeFinExt)); break;
                    case 46: fscanf(fid,"%f", &(wpp->margeFinInt)); break;
                    case 47: fscanf(fid,"%f", &(wpp->margeFinNerv)); break;
                    case 48: fscanf(fid,"%f", &(wpp->margeFinDiagNerv)); break;


                    case 49: fscanf(fid,"%f", &(wpp->tochnostRasklad2)); break;
                    case 50: fscanf(fid,"%f", &(wpp->XMashtab)); break;
                    case 51: fscanf(fid,"%f", &(wpp->textX)); break;
                    case 52: fscanf(fid,"%f", &(wpp->textY)); break;
                    case 53: {  fgets (s, 255, fid); trim(s); trim(s, '\r'); trim(s, '\n'); strcpy(wpp->NomFichierForme, s);  break;}
                    default:;
                }
               
            }
        }

	}
    if(fclose(fid))
    {
        printf("\nProbleme  la fermeture du fichier");
    }
    //printf ("\n..LectureWindPatternsProject(%s)", NomFic);
    return wpp;
}

void EcritureWindPatternsProject(char *NomFichier, WindPatternsProject *wpp) {
	FILE *fid;
	/**** message ****/
	printf("\nEcriture fichier WindPatternsProject: '%s'",NomFichier);
	/**** ouverture fichier en ecriture ****/
	if( (fid = fopen( NomFichier, "wt" )) == NULL )
	{
		printf( "\nErreur ouverture fichier '%s'", NomFichier);
		exit(0);
	}
    char* motsClef[54] =
                    {"PROJECT", // 1

		"MARGE_DEB", "MARGE_FIN", "MARGE1", "MARGE2", //4

		"PINCES","PINCES_POS_BA","PINCES_AMP_BA","PINCES_POS_BF","PINCES_AMP_BF", //5

                "PINCE_TYPE_NEZ","PINCE_RADIUS_ALG_K_NEZ","PINCE_EQUAL_AMP_TO_FUI", // 9
                "PINCE_POWER_NEZ", "PINCE_POWER_NEZ_VALUE",
                "PINCE_ARCTAN_NEZ", "PINCE_ARCTAN_NEZ_K1", "PINCE_ARCTAN_NEZ_K2", "PINCE_ARCTAN_NEZ_K3",

                "PINCE_TYPE_FUI", "PINCE_RADIUS_ALG_K_FUI", "PINCE_EQUAL_AMP_TO_NEZ", // 9
                "PINCE_POWER_FUI", "PINCE_POWER_FUI_VALUE",
                "PINCE_ARCTAN_FUI", "PINCE_ARCTAN_FUI_K1", "PINCE_ARCTAN_FUI_K2", "PINCE_ARCTAN_FUI_K3",

		"DIAG_NERVURES","DIAG_NERVURES_FILE", // 6
                "DIAG_NERVURES_POS_INT_A","DIAG_NERVURES_POS_INT_F", "DIAG_NERVURES_POS_EXT_A","DIAG_NERVURES_POS_EXT_F",

		"VENT_HOLES", "VENT_HOLES_FILE", "VENT_HOLES_DEB", "VENT_HOLES_FIN", // 7
                "VENT_HOLES_KLAPANS", "VENT_HOLES_KLAPANS_DOUBLE", "VENT_HOLES_KLAPANS_FIN",

                "REP_POINTS", "REP_POINTS_FROM_FILE", "REP_POINTS_FILE", // 3
                "LAYOUT_VENTILATION",
                "LAYOUT_MARGE_EXT_FIN","LAYOUT_MARGE_INT_FIN", "LAYOUT_MARGE_NERV_FIN", "LAYOUT_MARGE_DIAG_NERV_FIN",
                "LAYOUT_ACCURACY",//3
                "XMASHTAB", "TEXTX", "TEXTY", "FORME"}; //3

                    fprintf(fid,"\n%s %s", motsClef[0], wpp->name);
                    fprintf(fid,"\n\n%s %f", motsClef[1],wpp->MargeDeb);
                    fprintf(fid,"\n%s %f", motsClef[2],wpp->MargeFin);
                    fprintf(fid,"\n%s %f", motsClef[3],wpp->Marge[0]);
                    fprintf(fid,"\n%s %f", motsClef[4],wpp->Marge[1]);

                    fprintf(fid,"\n\n%s %d", motsClef[5],wpp->isPinces);
                    fprintf(fid,"\n%s %f", motsClef[6],wpp->PosPinceBA[0]);
                    fprintf(fid,"\n%s %f", motsClef[7],wpp->AmpPinceBA[0]);
                    fprintf(fid,"\n%s %f", motsClef[8],wpp->PosPinceBF[0]);
                    fprintf(fid,"\n%s %f", motsClef[9],wpp->AmpPinceBF[0]);

                    fprintf(fid,"\n\n%s %d", motsClef[10],wpp->PinceNosRadio);
                    fprintf(fid,"\n%s %f", motsClef[11],wpp->PinceRadiusAlgKNos);
                    fprintf(fid,"\n%s %d", motsClef[12],wpp->PinceNosEqualAmp);
                    fprintf(fid,"\n%s %d", motsClef[13],wpp->PincePowerA);
                    fprintf(fid,"\n%s %f", motsClef[14],wpp->PincePowerValueA);
                    fprintf(fid,"\n%s %d", motsClef[15],wpp->PinceArctanA);
                    fprintf(fid,"\n%s %f", motsClef[16],wpp->PinceArctanK1ValueA);
                    fprintf(fid,"\n%s %f", motsClef[17],wpp->PinceArctanK2ValueA);
                    fprintf(fid,"\n%s %f", motsClef[18],wpp->PinceArctanK3ValueA);

                    fprintf(fid,"\n\n%s %d", motsClef[19],wpp->PinceHvostRadio);
                    fprintf(fid,"\n%s %f", motsClef[20],wpp->PinceRadiusAlgKHvost);
                    fprintf(fid,"\n%s %d", motsClef[21],wpp->PinceHvostEqualAmp);
                    fprintf(fid,"\n%s %d", motsClef[22],wpp->PincePowerF);
                    fprintf(fid,"\n%s %f", motsClef[23],wpp->PincePowerValueF);
                    fprintf(fid,"\n%s %d", motsClef[24],wpp->PinceArctanF);
                    fprintf(fid,"\n%s %f", motsClef[25],wpp->PinceArctanK1ValueF);
                    fprintf(fid,"\n%s %f", motsClef[26],wpp->PinceArctanK2ValueF);
                    fprintf(fid,"\n%s %f", motsClef[27],wpp->PinceArctanK3ValueF);

                    fprintf(fid,"\n\n%s %d", motsClef[28],wpp->DiagNervs);
                    fprintf(fid,"\n%s %s", motsClef[29],wpp->NomFichierDiagNerv);
                    fprintf(fid,"\n%s %f", motsClef[30],wpp->PosDiagNerv2A);
                    fprintf(fid,"\n%s %f", motsClef[31],wpp->PosDiagNerv2F);
                    fprintf(fid,"\n%s %f", motsClef[32],wpp->PosDiagNerv1A);
                    fprintf(fid,"\n%s %f", motsClef[33],wpp->PosDiagNerv1F);

                    fprintf(fid,"\n\n%s %d", motsClef[34],wpp->VentHoles);
                    fprintf(fid,"\n%s %s", motsClef[35],wpp->NomFichierVentHoles);
                    fprintf(fid,"\n%s %f", motsClef[36],wpp->VentHolesDeb);
                    fprintf(fid,"\n%s %f", motsClef[37],wpp->VentHolesFin);

                    fprintf(fid,"\n\n%s %d", motsClef[38],wpp->RaskladKlapans);
                    fprintf(fid,"\n%s %d", motsClef[39],wpp->VentHolesDouble);
                    fprintf(fid,"\n%s %f", motsClef[40],wpp->PosKlapanFin);

                    fprintf(fid,"\n\n%s %d", motsClef[41],wpp->RepPoints);
                    fprintf(fid,"\n%s %d", motsClef[42],wpp->ReperPointsFromFile);
                    fprintf(fid,"\n%s %s", motsClef[43],wpp->NomFichierRepPoints);

                    fprintf(fid,"\n\n%s %d", motsClef[44],wpp->VentilationLayout);

                    fprintf(fid,"\n%s %f", motsClef[45],wpp->margeFinExt);
                    fprintf(fid,"\n%s %f", motsClef[46],wpp->margeFinInt);
                    fprintf(fid,"\n%s %f", motsClef[47],wpp->margeFinNerv);
                    fprintf(fid,"\n%s %f", motsClef[48],wpp->margeFinDiagNerv);

                    fprintf(fid,"\n%s %f", motsClef[49],wpp->tochnostRasklad2);
                    fprintf(fid,"\n%s %f", motsClef[50],wpp->XMashtab);
                    fprintf(fid,"\n%s %f", motsClef[51],wpp->textX);
                    fprintf(fid,"\n%s %f", motsClef[52],wpp->textY);
                    fprintf(fid,"\n%s %s", motsClef[53],wpp->NomFichierForme);

	if(fclose(fid))
	{
		printf("\nProbleme  la fermeture du fichier WindPatternsProject");
		exit(0);
	}
    
}

Forme* LectureFichierForme(char* NomFic)
{
//    if (DEBUG) printf ("\n LectureFichierForme");
	int i,j,n;
	FILE *fid;
	Matrice *m=NULL;
	double noVer=0.0f;
	Forme *f=NULL;
	//liste des mots clefs du fichier
	char* motsClef[19] = {"VERSION",
		"NB_ALVEOLES", "COEFF_PROG_GEOM", "COEFF_EXP", "EPAI_REL_CENTRE",
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
		for (int i = 0; i < 19; i++)
		{
                    //if (DEBUG) printf ("\n\n\n i=%d", i);
                    bool go = true;
                    //recherche mot clef
                    if(!TrouveMotDansFichierTexte(fid, motsClef[i])) {
                        go=false;
                        printf("\n\tWarning: mot clef '%s' not in file!", motsClef[i]);
                        if ( i == 11) {
//                            if (DEBUG) printf (" i==11");
                            f->mCtrlE = CloneMat(f->mCtrlD);
							f->mCtrlE->SetElement(0, 1, f->mCtrlD->Element(0, 1)-0.2f);
							f->mCtrlE->SetElement(1, 1, f->mCtrlD->Element(0, 1)-0.2f);
                			f->mCtrlE->SetElement(2, 1, f->mCtrlD->Element(2, 1)-0.2f);
							f->mCtrlE->SetElement(3, 1, f->mCtrlD->Element(3, 1)-0.01f);
                            continue;
                        }
                        if (i == 3) {
//                            if (DEBUG) printf (" i==3");
                            f->CoeffExp = 0.0f;
                            go = false;
                        }
                    }
                    if (go)
                    {
//                        if (DEBUG) printf (" go!");
                            //action en fonction du mot clef
                            switch(i)
                            {
                            case 0: fscanf(fid,"%lf", &noVer); break;
                            case 1: fscanf(fid,"%d", &(f->NbCaiss)); break;
                            case 2: fscanf(fid,"%lf", &(f->CoeffProgGeom));  break;
                            case 3: fscanf(fid,"%lf", &(f->CoeffExp));  break;
                            case 4: fscanf(fid,"%lf", &(f->EpaiRelCent)); break;
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
                            case 15: {
                                    m = new Matrice(4,2);
                                    double f1, f2, f3, f4, f5, f6, f7, f8;
                                    fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf", &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8 );
//                                    if (DEBUG) printf (" %lf %lf %lf %lf %lf %lf %lf %lf", f1, f2, f3, f4, f5, f6, f7, f8);
                                    m->SetElement(0,0,f1); m->SetElement(0,1,f2);
                                    m->SetElement(1,0,f3); m->SetElement(1,1,f4);
                                    m->SetElement(2,0,f5); m->SetElement(2,1,f6);
                                    m->SetElement(3,0,f7); m->SetElement(3,1,f8);
                                    switch(i)
                                    {
                                    case 5: f->mCtrlNez = m; break;
                                    case 6: f->mCtrlFui = m; break;
                                    case 7: f->mCtrlA = m; break;
                                    case 8: f->mCtrlB = m; break;
                                    case 9: f->mCtrlC = m; break;
                                    case 10: f->mCtrlD = m; break;
                                    case 11: f->mCtrlE = m; break;
                                    case 12: f->mCtrlDiedre = m; break;
                                    case 13: f->mCtrlMorphing = m; break;
                                    case 14: f->mCtrlVrillage = m; break;
                                    case 15: f->mCtrlEpaiRel = m; break;
                                    default:;
                                    }
                                    break;
                            }
                            case 16: {
                                    char temp[1024];
                                                    fscanf(fid,"%s", temp);
                                                    f->m_strNomProfilCent = temp;
                                                    break;
                                             }

                            case 17: {
                                                    char temp[1024];
                                                    fscanf(fid,"%s", temp);
                                                    f->m_strNomProfilBout = temp;
                                                    break;
                                             }

                            case 18: {
                                    f->DeleteProfils();
                                    fscanf(fid,"\n%d",&n);
                                //fscanf(fid,"%d",&n);
                                f->AllocateProfils(n);
                                char* s = new char[255];
                                fgets (s, 255, fid);
                                for (j = 0; j<f->m_nbProfils; j++)
                                {
                                        double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13 = 0.0f;
                                        fgets (s, 255, fid);
                                        sscanf(s,"\n%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                                                                  &n, //no Nervure
                                                                  &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11, &f12, &f13 );
//                                        if (DEBUG) printf ("\n\n__ %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",n, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13);
                                        if (f13 == 0.0f) f13=98.0f;
                                        /*fscanf(fid,"\n%d %f %f %f %f %f %f %f %f %f %f %f %f %f",
                                                                  &n, //no Nervure
                                                                  &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11, &f12, &f13 );*/
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
                                        f->m_pProfils[j]->m_fPosD	= f12;	//pos relative ligne D*/
                                        f->m_pProfils[j]->m_fPosE	= f13;	//pos relative ligne E
                                }//boucle lecture tableau de forme

                                break;
                            }
                                    default:;

                            }//switch mot clef

                    }//test trouve mot clef

		}//boucle

		printf ("\n f->m_strNomProfilCent=%s", f->m_strNomProfilCent.c_str());
		printf ("\n f->m_strNomProfilBout=%s", f->m_strNomProfilBout.c_str());

		/**** fermeture fichier ****/
		if(fclose(fid))
		{
			printf("\nProbleme  la fermeture du fichier");
		}

	}
	//test ouverture fichier ok
	// if (DEBUG) printf ("\n pg=%f ce=%f er=%f", f->CoeffProgGeom, f->CoeffExp, f->EpaiRelCent);
    // if (DEBUG) printf ("\n ...LectureFichierForme");
	return f;
}

/***********************/

/* LectureFichierForme2 */

/***********************/

Forme* LectureFichierForme2(char* NomFic)

{
//    if (DEBUG) printf ("\n LectureFichierForme2()");
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

	char* motsClef[30] = {"VERSION",

		"NB_ALVEOLES", "COEFF_PROG_GEOM","COEFF_EXP", "EPAI_REL_CENTRE",

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

				printf("\n\tWarning: mot clef '%s' non trouvÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ !", motsClef[i]);

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

						//boucle ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ partir de la 1ere nervure du centre vers l'extrÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½mitÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½

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

			printf("\nProbleme ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ la fermeture du fichier");

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
	Matrice* m[11];
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
	fprintf(fid,"\nCOEFF_EXP %1.5f", f->CoeffExp);
	fprintf(fid,"\nEPAI_REL_CENTRE %2.3f", f->EpaiRelCent);

	

	/**** Ecriture des points de controle ****/

	

	/*init pointeur des matrices de coordonnÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½es des points de controle*/
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
	/*boucle ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ partir de la 1ere nervure du centre vers l'extrÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½mitÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½*/
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

		printf("\nProbleme ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ la fermeture du fichier");

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

	

	/*init pointeur des matrices de coordonnÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½es des points de controle*/

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

	/*boucle ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ partir de la 1ere nervure du centre vers l'extrÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½mitÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½*/

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

		printf("\nProbleme ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ la fermeture du fichier");

		exit(0);

	}

}


/********************/

/* CalculInfoForme */

/********************/

void CalculInfoForme( Forme* F, TInfoForme* info )

{
	double xNerv, xNervPrec, lNerv, lNervPrec, larg;
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
			+ (double)sqrt(
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



	/*envergure projetÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½e*/

	info->envergureProj =  2.0f * F->m_pProfils[nNerv-1]->m_fNezX;

	/*surface projetÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½e*/

	info->surfaceProj = 0.0f; xNervPrec = 0.0f;

	lNervPrec = F->m_pProfils[0]->m_fLength;
	for(i=0; i<nNerv; i++)
	{
		lNerv = F->m_pProfils[i]->m_fLength;
		xNerv = F->m_pProfils[i]->m_fNezX;
		info->surfaceProj += (lNerv + lNervPrec) * (xNerv - xNervPrec);
		lNervPrec = lNerv; xNervPrec = xNerv;
	}

	/*allongement projetÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½*/

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
				   double x1, double y1, double z1,
				   double x2, double y2, double z2)
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
				   double x1, double y1, double z1)
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
			double height = current->m_fWidth;
			if ( height == 0.0f )
				height = 0.1f;
			height /= 100.0f;

			double a = current->m_fInclin;

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




void Forme::SerializeFoil (QDataStream &ar, Matrice* ExtProf, Matrice* IntProf) { 
	CFoil* foil = new CFoil();
	foil->InitFromWindWax(ExtProf, IntProf);
	foil->m_FoilName="DummyFoil";
	foil->Serialize(ar, true);
}

/*
void CalculForme3D(Forme *forme, int isPercent, double percent,
				   Matrice *ExtProfCent, Matrice *IntProfCent,
				   Matrice *ExtProfBout, Matrice *IntProfBout,
				   Matrice **XExt, Matrice **YExt, Matrice **ZExt,
				   Matrice **XInt, Matrice **YInt, Matrice **ZInt)

{
    Matrice *ExtProfCentN, *ExtProfBoutN;
	double LongNerv, EpaiRel, xp,yp, xo,yo,zo, a,v,m;
	double EpaiRelProfCent, EpaiRelProfBout;
	double coeffx, coeffyCent, coeffyBout;
	int i,j;
    bool isCenterPanel = (1 & forme->NbCaiss);

	EpaiRelProfCent = EpaisseurRelative(ExtProfCent, IntProfCent);
	EpaiRelProfBout = EpaisseurRelative(ExtProfBout, IntProfBout);
	
	*XExt = Zeros(forme->m_nbProfils, ExtProfCent->GetLignes());
	*YExt = Zeros(forme->m_nbProfils, ExtProfCent->GetLignes());
	*ZExt = Zeros(forme->m_nbProfils, ExtProfCent->GetLignes());

	*XInt = Zeros(forme->m_nbProfils, IntProfCent->GetLignes());
	*YInt = Zeros(forme->m_nbProfils, IntProfCent->GetLignes());
	*ZInt = Zeros(forme->m_nbProfils, IntProfCent->GetLignes());

	for (i=0; i<forme->m_nbProfils; i++)
	{
		LongNerv = forme->m_pProfils[i]->m_fLength;
		EpaiRel = forme->m_pProfils[i]->m_fWidth;
		xo = forme->m_pProfils[i]->m_fNezX;
		yo = forme->m_pProfils[i]->m_fNezY;
		zo = forme->m_pProfils[i]->m_fNezZ;
		a = forme->m_pProfils[i]->m_fInclin;
        if ((isCenterPanel) && (i == 0)) a = forme->m_pProfils[0]->m_fInclin - 0.33333333f * fabs(forme->m_pProfils[0]->m_fInclin - forme->m_pProfils[1]->m_fInclin);
                
		v = forme->m_pProfils[i]->m_fWash;
		m = forme->m_pProfils[i]->m_fMorph;

		coeffx = LongNerv/100.0f;
		coeffyCent = LongNerv*EpaiRel/(EpaiRelProfCent*100.0f);
		coeffyBout = LongNerv*EpaiRel/(EpaiRelProfBout*100.0f);

		for (j=0; j<IntProfCent->GetLignes(); j++)
		{
			xp = IntProfCent->Element(j,0)*coeffx;
			yp = IntProfCent->Element(j,1)*coeffyCent*m
				+ IntProfBout->Element(j,1)*coeffyBout*(1.0f-m);
			(*XInt)->SetElement(i,j, xo+(yp-xp*(double)sin(v))*(double)cos(a));
			(*YInt)->SetElement(i,j, yo+(yp-xp*(double)sin(v))*(double)sin(a));
			(*ZInt)->SetElement(i,j, zo-(xp*(double)cos(v)));

		}

		for (j=0; j<ExtProfCent->GetLignes(); j++)
		{
            xp = ExtProfCent->Element(j,0)*coeffx;
            yp = ExtProfCent->Element(j,1)*coeffyCent*m
		                        + ExtProfBout->Element(j,1)*coeffyBout*(1.0f-m);

			(*XExt)->SetElement(i,j, xo+(yp-xp*(double)sin(v))*(double)cos(a));
			(*YExt)->SetElement(i,j, yo+(yp-xp*(double)sin(v))*(double)sin(a));
			(*ZExt)->SetElement(i,j, zo-(xp*(double)cos(v)));
  		}
	}
}
*/

void Forme::WritePolars (QDataStream &ar) {
		printf ("\n WritePolars..");
	    ar << 100003;
        int m_oaFoilSize=m_nbProfils;
        ar << (int)m_oaFoilSize;

		//SerializeFoil(ar, ExtProfCent, IntProfCent);
        //for (int i=0; i<m_oaFoilSize; i++){
                //pFoil = (CFoil*)m_oaFoil.GetAt(i);
                //pFoil->Serialize(ar);
                //SerializeFoil(ar);
		//}
		double LongNerv, EpaiRel, xp,yp, xo,yo,zo, a,v,m, w, kw;
		double EpaiRelProfCent, EpaiRelProfBout;
		double coeffx, coeffyCent, coeffyBout;
		int i, j;
        //bool isCenterPanel = (1 & forme->NbCaiss);

		EpaiRelProfCent = EpaisseurRelative(ExtProfCent, IntProfCent);
		EpaiRelProfBout = EpaisseurRelative(ExtProfBout, IntProfBout);
		char foilName[50];
		printf ("\n Forme->m_nbProfils=%d [0 -> %d]", m_nbProfils, (m_nbProfils-1));
		double w0 = m_pProfils[0]->m_fWidth;
		for (i=0; i<m_nbProfils; i++)
		{
			Matrice *ExtProf, *IntProf;
			ExtProf = Zeros(ExtProfCent->GetLignes(), 2);
			IntProf = Zeros(IntProfCent->GetLignes(), 2);

			//LongNerv = forme->m_pProfils[i]->m_fLength;
			w = m_pProfils[i]->m_fWidth;
			
			kw = w / w0;
			if (i == (m_nbProfils-1))  kw=0.1;
			//xo = forme->m_pProfils[i]->m_fNezX;
			//yo = forme->m_pProfils[i]->m_fNezY;
			//zo = forme->m_pProfils[i]->m_fNezZ;
			//a = forme->m_pProfils[i]->m_fInclin;
			//if ((isCenterPanel) && (i == 0)) a = forme->m_pProfils[0]->m_fInclin - 0.33333333f * fabs(forme->m_pProfils[0]->m_fInclin - forme->m_pProfils[1]->m_fInclin);
	                
			//v = forme->m_pProfils[i]->m_fWash;
			m = m_pProfils[i]->m_fMorph;
			//if (DEBUG) printf ("\n%3d -> a=%f v=%f m=%f", i, a * 180.0f/pi, v, m);
			coeffx = 1.0f;//LongNerv/100.0f;
			coeffyCent = 1.0f;// LongNerv*EpaiRel/(EpaiRelProfCent*100.0f);
			coeffyBout = 1.0f;// LongNerv*EpaiRel/(EpaiRelProfBout*100.0f);
			//printf ("\n i=%d kw=%f", i , kw);
			for (j=0; j<IntProfCent->GetLignes(); j++)
			{
				xp = IntProfCent->Element(j,0)*coeffx;
				yp = IntProfCent->Element(j,1)*coeffyCent*m
					+ IntProfBout->Element(j,1)*coeffyBout*(1.0f-m);
				if (kw != 1) yp = yp * kw;
				if (i == (m_nbProfils-1)) printf (" %f", yp);
				IntProf->SetElement(j,0,xp);
				IntProf->SetElement(j,1,yp);
			}

			for (j=0; j<ExtProfCent->GetLignes(); j++)
			{
				xp = ExtProfCent->Element(j,0)*coeffx;
				yp = ExtProfCent->Element(j,1)*coeffyCent*m
						+ ExtProfBout->Element(j,1)*coeffyBout*(1.0f-m);
				if (kw != 1) yp = yp * kw;
				ExtProf->SetElement(j,0,xp);
				ExtProf->SetElement(j,1,yp);
			}

			CFoil* foil = new CFoil();
			foil->InitFromWindWax(ExtProf, IntProf);
			sprintf(foilName, "P%03d", i);
			foil->m_FoilName = foilName;
			printf ("\n Serialize: %s", foilName);
			foil->Serialize(ar, true, 5);

			delete (foil);
			delete (ExtProf);
			delete (IntProf);
		}

        //then write polars
        ar << 0;
        ar << 0;
		printf ("\n ...WritePolars..");
}



bool EcritureFichierWpa(char *NomFichier, Forme *forme)
{
	//CWaitCursor wait;
	CFileException fe;
    string nomFichierStr = NomFichier;
	QFile fp(NomFichier);

	if (!fp.open(QIODevice::WriteOnly))
	{
		printf ("\nCould not open the file for writing");
		return false;
	}

	QDataStream ar(&fp);//, CArchive::store);
	ar.setByteOrder(QDataStream::LittleEndian);
	printf ("\n EcritureFichierWpa");
	forme -> SerializeToWpa(ar);
	printf ("\n...EcritureFichierWpa");
	//ar.Close();
	fp.close();
	return true;
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
                //	printf ("\n courb [%s]", courb->name.c_str());
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
		printf("\nProbleme ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ la fermeture du fichier");
		exit(0);
	}
}

//EcritureManyFichierPolyDXF(PtrNomFichier, AxePatronDXF, AxeMarginDXF, 1, AxeRepDXF, 0, AxeCercleDXF, Numerotation, AxePatronTextDXF, W, H);

void EcritureManyFichierPolyDXF(char *NomFichier, int np, int n, TAxe **axe, TAxe **axe2, int rep, TAxe **axeR, int vent, TAxe **axeC, int num, TAxe **axeT, double* W, double* H)
{
    	FILE *fid;
	printf("\nEcritureMANY POLY fichier de DXF: '%s'",NomFichier);
	if( (fid = fopen( NomFichier, "wt" )) == NULL ) {
		printf( "\nErreur ouverture fichier '%s'", NomFichier);
		exit(0);
	}
        fprintf(fid,"0\nSECTION\n2\nENTITIES\n");
        double dx=0.0f,dy=0.0f,dz=0.0f;

        int povorot1=0, povorot2=0, len1=0, len2=0;
        if (np&1) {
            povorot1=(2*(np-1));
            len1=2*(np-1);
            povorot2=povorot1 + ((np-1));
            len2=np-1;
        } else {
            povorot1=(2*(np-1)+1);
            len1=2*(np-1)+1;
            povorot2=povorot1 + ((np));
            len2=np;
        }
        double maxW1=-1000000.0f,maxW2=-1000000.0f;
        /*
        for (int i=0; i<n;i++) {

            if (i==povorot1) { dx=maxW1*2.0; dy=0.0f; maxW2=-1000000.0f;}
            if (i==povorot2) { dx=maxW1*2.0 + maxW2*2.0; dy=0.0f;}
            printf ("\nMFP i=%d",i);
            EcritureFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], vent, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
            if (W[i] > maxW1) maxW1 = W[i];
            if (W[i] > maxW2) maxW2 = W[i];
        }*/

        for (int i=0; i<povorot1;i++) {
            EcritureFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], vent, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
            if (W[i] > maxW1) maxW1 = W[i];
        }
        dx=maxW1*2.0; dy=0.0f; maxW2=-1000000.0f;
        /*int puk=0;
        if (np&1) {
            puk = povorot1+1;
        } else {
            puk = povorot1;
        }*/
        for (int i=povorot2-1;i>=povorot1;i--) {
            EcritureFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], vent, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
            if (W[i] > maxW2) maxW2 = W[i];
        }
        for (int i=povorot1+(np&1);i<povorot2;i++) {
            EcritureFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], vent, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
            if (W[i] > maxW2) maxW2 = W[i];
        }
        dx=maxW1*2.0 + maxW2*4.0; dy=0.0f;
        for (int i=povorot2; i<n;i++) {
            EcritureFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], vent, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
        }

       fprintf(fid,"0\nENDSEC\n0\nEOF\n");
       //printf ("\n...fclose(fid)");
	if(fclose(fid))	{
		printf("\nProbleme ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ la fermeture du fichier");
		exit(0);
	}
       //printf ("\n end of EcrFichManyF...");

}
void EcritureManyFichierPolyDXF2(char *NomFichier, int np, int n, TAxe **axe, TAxe **axe2, int rep, TAxe **axeR, int vent, TAxe **axeC, int num, TAxe **axeT, double* W, double* H, int* numncon)
{
    	FILE *fid;
	printf("\nEcritureMANY POLY fichier de DXF: '%s'",NomFichier);
	if( (fid = fopen( NomFichier, "wt" )) == NULL ) {
		printf( "\nErreur ouverture fichier '%s'", NomFichier);
		exit(0);
	}
        fprintf(fid,"0\nSECTION\n2\nENTITIES\n");
        double dx=0.0f,dy=0.0f,dz=0.0f;

        int povorot1=0, povorot2=0,povorot2d=0, len1=0, len2=0;
        povorot1=numncon[1];
        povorot2=numncon[2];
        povorot2d=numncon[3];
/*        if (np&1) {
            povorot1=(2*(np-1));
            len1=2*(np-1);
            povorot2=povorot1 + ((np-1));
            len2=np-1;
        } else {
            povorot1=(2*(np-1)+1);
            len1=2*(np-1)+1;
            povorot2=povorot1 + ((np));
            len2=np;
        }*/
        
        double maxW1=-1000000.0f,maxW2=-1000000.0f,maxW3=-1000000.0f;

        for (int i=0; i<povorot1;i++) {
            EcritureFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], 0, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
            if (W[i] > maxW1) maxW1 = W[i];
        }
        dx=maxW1*2.0; dy=0.0f; maxW2=-1000000.0f; maxW3=-1000000.0f;
        /* +(np&1) */
        for (int i=povorot1;i<povorot2;i++) {
            EcritureFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], vent, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
			//printf ("\n%d dy+%f",i, H[i]);
            if (W[i] > maxW2) maxW2 = W[i];
        }
        dx=maxW1*2.0 + maxW2*4.0; dy=0.0f;
        for (int i=povorot2; i<povorot2d;i++) {
            EcritureFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], vent, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
			//printf ("\n%d dy+%f",i, H[i]);
            if (W[i] > maxW3) maxW3 = W[i];
        }

        dx=maxW1*2.0 + maxW2*4.0 + maxW3*4.0; dy=0.0f;
        for (int i=povorot2d; i<n;i++) {
            EcritureFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], 0, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
        }

       fprintf(fid,"0\nENDSEC\n0\nEOF\n");
       //printf ("\n...fclose(fid)");
	if(fclose(fid))	{
		printf("\nProbleme  la fermeture du fichier");
		exit(0);
	}

}

/**********************/
/* EcritureFichierPolyDXF */
/**********************/

void EcritureFichierPolyDXF(char *NomFichier, TAxe *axe, TAxe *axe2, int rep, TAxe *axeR, int vent, TAxe *axeC, int num, TAxe *axeT )
{

	int i, j, cpt;
	double fx, fy, fz;
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
	double curx, cury, x1, y1;
	/*ecriture courbes*/
	for (int _i=0; _i < 2; _i++){
		curx = -1000000.0f;
		cury = -2000000.0f;
		x1 = -3000000.0f;
		y1 = -4000000.0f;
//                printf ("\n _i=%d", _i);
		if (_i == 0) { 
			courb=axe->Courb;
			sprintf(text, "Patron");
		}
		if (_i == 1) {
			courb=axe2->Courb;
			sprintf(text, "Marge");
		}
		cpt=0;
//              printf ("\n before polybegin");
		EcritPolyBeginDXF(fid, text);
//              printf ("\n after polybegin");
		int q = 0;
		while (courb != NULL)
		{
//			printf ("\n while in ecritureFichierPolyDXF");
//			printf ("\n courb next");

			if (q == 0) {
  //                          printf ("\n q==0");
        			x1 = courb->pts->Element(0,0);
				y1 = courb->pts->Element(0,1);
				//printf ("XY1 (%f, %f)", x1, y1);
//                            printf ("\n... q==0");
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
  //                      printf ("\nbefore sprintf");
                 	sprintf(text, "%s", courb->name.c_str());
    //                    printf ("\n... sprintf");
			//text = layerName;
			for(i = 0; i < courb->pts->GetLignes(); i++)
			{
//                            printf ("\n %d/%d",i, courb->pts->GetLignes());
//                            printf ("\n getx");
				fx = courb->pts->Element(i,0);
//                            printf ("\n gety");
				fy = courb->pts->Element(i,1);
//                            printf ("\n getz");
				//fz = courb->pts->Element(i,2);
                                fz = 0;
//                            printf ("\n ....after get");
//				printf ("\n%s : f(%f, %f) cur (%f, %f) %d", courb->name.c_str(), fx, fy, curx, cury, q);
				if ( (q==0) || !((fx == curx) && (fy == cury)) ) {
					if ((q==0) || !((fx == x1) && (fy == y1)) ) {
						EcritPolyVertexDXF(fid, text, axe->axe3d, fx, fy, fz);
//						printf (" !print!");
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
        //printf ("\n before vent");
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
				//fz = courb->pts->Element(i,2);
                                fz = 0;
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
		printf("\nProbleme ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½ la fermeture du fichier");
		exit(0);
	}
}


/**********************/
/* EcritureFichierPolyDXFDelta */
/**********************/

void EcritureFichierPolyDXFDelta(FILE *fid, TAxe *axe, TAxe *axe2, int rep, TAxe *axeR, int vent, TAxe *axeC, int num, TAxe *axeT, double dx, double dy, double dz, int n )
{
    //printf ("\n EcritureFichPolyDxfDelta");
	int i, j, cpt;
	double fx, fy, fz;
	Courbe *courb;
        TMesh *mesh;
	char text[50];
	double curx, cury, x1, y1;
	for (int _i=0; _i < 2; _i++){
  //          printf ("\n _i=%d", _i);
		curx = -1000000.0f;
		cury = -2000000.0f;
		x1 = -3000000.0f;
		y1 = -4000000.0f;
		if (_i == 0) {
			courb=axe->Courb;
			//sprintf(text, "Patron%d", n);
                        sprintf(text, "Patron");
		}
		if (_i == 1) {
			courb=axe2->Courb;
			//sprintf(text, "Marge%d", n);
                        sprintf(text, "Marge");
		}
		cpt=0;
//                printf ("\n 1 ");
//                printf ("\n text=[%s] ", text);
		EcritPolyBeginDXF(fid, text);
//              printf ("after EcritPolyBegin...");
		int q = 0;
		while (courb != NULL)
		{
			if (q == 0) {
    //                        printf ("\nq==0");
				x1 = courb->pts->Element(0,0);
				y1 = courb->pts->Element(0,1);
				//printf ("XY1 (%f, %f)", x1, y1);
      //                      printf ("\n...q==0");
			}
			cpt++;
                 	//sprintf(text, "%s", courb->name.c_str());
			//text = layerName;
			for(i = 0; i < courb->pts->GetLignes(); i++)
			{
				fx = courb->pts->Element(i,0);
				fy = courb->pts->Element(i,1);
				fz = courb->pts->Element(i,2);
//				printf ("\n%s : f(%f, %f) cur (%f, %f) %d", courb->name.c_str(), fx, fy, curx, cury, q);
				if ( (q==0) || !((fx == curx) && (fy == cury)) ) {
					if ((q==0) || !((fx == x1) && (fy == y1)) ) {
						EcritPolyVertexDXF(fid, text, axe->axe3d, dx+fx, dy + fy, fz + dz);
//						printf (" !print!");
						curx = fx;
						cury = fy;
					}
				}
				q++;
			}
			courb=courb->CourbSuiv;
		}
		EcritPolyEndDXF(fid, text);
//                printf ("...end...");
	}

	if (vent == 1) {
//            printf ("\n Vent");
                courb=axeC->Courb;
		//sprintf(text, "Vent%d", n);
                sprintf(text, "Vent");
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
//            printf ("\n in Rep");
		courb=axeR->Courb; cpt=0;
		while (courb != NULL)
		{
			cpt++;
			//sprintf(text, "Rep%d", n);
                        sprintf(text, courb->name.c_str());
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
//                printf ("\n Num=1");
		courb=axeT->Courb; cpt=0;
		while (courb != NULL)
		{
			cpt++;
			//sprintf(text, "Text%d", n);
                        sprintf(text, "Text");
			for(i=0; i<courb->pts->GetLignes()-1; i++)
			{
				EcritLigneDXF(fid, text, axe->axe3d,
					dx+courb->pts->Element(i,0), dy+courb->pts->Element(i,1), dz+courb->pts->Element(i,2),
					dx+courb->pts->Element(i+1,0), dy+courb->pts->Element(i+1,1), dz+courb->pts->Element(i+1,2));
			}
			courb=courb->CourbSuiv;
		}
	}
//        printf ("\n before meshes");
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
        //printf ("...poly dxf delta");
}


void GenerateCourbe(WindPatternsProject* gfd, Matrice *Xd1,
					Matrice *Yd1, Matrice *P1,
					int nerv1, double deb1, int faceDeb1, double fin1, int faceFin1,
					Matrice *Xd2, Matrice *Yd2, Matrice *P2,
					int nerv2, double deb2, int faceDeb2, double fin2, int faceFin2, char *text,
			        TAxe **AxePatronP, TAxe **AxePatronDXFP, TAxe **AxePatronTextDXFP, TAxe **AxeMarginDXFP, TAxe **AxeCercleDXFP, TAxe **AxeRepDXFP, int Ventilation,
                    double marge1, double marge2, double margeDeb, double margeFin,bool makeRep, bool debug,
					bool isPince, Matrice *Xd01, Matrice *Yd01,double coeff1, Matrice *Xd02, Matrice *Yd02,double coeff2)
{
    if (debug) printf ("\n GenerateCourbe()");
    double Marge[2]={marge1, marge2};

    *AxePatronP = CreerAxe(gfd->FenetrePatron);
    *AxePatronDXFP = CreerAxe(gfd->FenetrePatron);
    *AxeRepDXFP = CreerAxe(gfd->FenetrePatron);
    *AxePatronTextDXFP = CreerAxe(gfd->FenetrePatron);
    *AxeMarginDXFP = CreerAxe(gfd->FenetrePatron);
    *AxeCercleDXFP = CreerAxe(gfd->FenetrePatron);
    int ajCoAr = 0, ajCoMaAr = 0, ajCoAv = 0, ajCoMaAv = 0;
    int ajCo1[2] = {0, 0};
    int ajCo2[2] = {0, 0};
    int nerv[2] = {nerv1, nerv2};
    if (nerv1 == -1) nerv[0]=0;
    int face[2] = {faceFin1, faceFin2};
    int FaceDeb[2]= {faceDeb1, faceDeb2};
    int FaceFin[2]= {faceFin1, faceFin2};
    double Deb[2] = {deb1, deb2};
    double Fin[2] = {fin1, fin2};
    Matrice *distance;
    int n;
    double xText, yText;
    char texte[100], texteExtInt[3] = "EI";
    Matrice *interpXSuspente, *interpYSuspente;
    double posSuspente, xSuspente[2][205], ySuspente[2][205];

    Matrice * Xd[2] = {Xd1, Xd2};
	Matrice * Xd0[2] = {Xd01, Xd02};
    Matrice * Yd[2] = {Yd1, Yd2};
	Matrice * Yd0[2] = {Yd01, Yd02};
    Matrice * P[2] = {P1, P2};
	double coeff[2]= {coeff1, coeff2};
    int i = 0, j = 0;
    Courbe * CourbPatron[2], *CourbPatronDXF[2], *CourbPatronBack[2], *CourbMarge[2], *CourbMargeBack[2], *CourbMargeDXF[2];
    Courbe *CourbAv, *CourbAvBack, *CourbAr, *CourbArDXF, *CourbMargeAv, *CourbMargeAvBack, *CourbMargeAr, *CourbCoin, *CourbMargeArDXF, *CourbCoin1[2], *CourbCoin1Back[2], *CourbCoin2[2], *CourbCoin2Back[2], *CourbRep, *CourbRepDXF;
    Courbe *CourbCercle, *CourbCercleDXF;
    for (i = 0; i < 2; i++) {
      if (debug) printf ("\n GenerateCourbe 1for: %d", i);
        CourbPatron[i] = new Courbe("Patron");
        CourbPatronDXF[i] = new Courbe("Patron");
        CourbPatronBack[i] = new Courbe("PatronBack");
        CourbPatron[i]->points = OFF;
        CourbPatronDXF[i]->points = OFF;
        CourbPatronBack[i]->points = OFF;
        CourbPatron[i]->symX = OFF;
        CourbPatronDXF[i]->symX = OFF;
        CourbPatronBack[i]->symX = OFF;
        CourbPatron[i]->pts = new Matrice(Xd[i]->GetLignes(), 2);
        CourbPatronDXF[i]->pts = new Matrice(Xd[i]->GetLignes(), 2);
        CourbPatronBack[i]->pts = new Matrice(Xd[i]->GetLignes(), 2);
        for (j = 0; j < Xd[i]->GetLignes(); j++) {
            CourbPatron[i]->pts->SetElement(j, 0, Xd[i]->Element(j, 0));
            CourbPatronDXF[i]->pts->SetElement(j, 0, Xd[i]->Element(j, 0));
            CourbPatronBack[i]->pts->SetElement(Xd[i]->GetLignes() - j - 1, 0, Xd[i]->Element(j, 0));

            CourbPatron[i]->pts->SetElement(j, 1, Yd[i]->Element(j, 0));
            CourbPatronDXF[i]->pts->SetElement(j, 1, Yd[i]->Element(j, 0));
            CourbPatronBack[i]->pts->SetElement(Xd[i]->GetLignes() - j - 1, 1, Yd[i]->Element(j, 0));
        }
        AjoutCourbe(*AxePatronP, CourbPatron[i]);
        distance = Ones(CourbPatron[i]->pts->GetLignes(), 1);
        for (j = 0; j < distance->GetLignes(); j++) distance->MultiplyElement(j, 0, Marge[i] / 100.0f);
        if (debug) printf ("\n CC Marge/Marge/MargeBack");
        CourbMarge[i] = new Courbe("Marge");
        CourbMargeDXF[i] = new Courbe("Marge");
        CourbMargeBack[i] = new Courbe("MargeBack");
        if (i == 0) {
            CourbMarge[i]->pts = CalculContour(CourbPatron[i]->pts, distance, -1);
            CourbMargeDXF[i]->pts = CalculContour(CourbPatron[i]->pts, distance, -1);
            CourbMargeBack[i]->pts = CalculContour(CourbPatron[i]->pts, distance, -1);
        } else {
            CourbMarge[i]->pts = CalculContour(CourbPatron[i]->pts, distance, +1);
            CourbMargeDXF[i]->pts = CalculContour(CourbPatron[i]->pts, distance, +1);
            CourbMargeBack[i]->pts = CalculContour(CourbPatron[i]->pts, distance, +1);
        }
        CourbMarge[i]->points = OFF;
        CourbMarge[i]->symX = OFF;
        CourbMargeDXF[i]->points = OFF;
        CourbMargeDXF[i]->symX = OFF;
        CourbMargeBack[i]->points = OFF;
        CourbMargeBack[i]->symX = OFF;

        for (j = 0; j < CourbMarge[i]->pts->GetLignes(); j++) {
            CourbMargeBack[i]->pts->SetElement(CourbMarge[i]->pts->GetLignes() - j - 1, 0, CourbMarge[i]->pts->Element(j, 0));
            CourbMargeBack[i]->pts->SetElement(CourbMarge[i]->pts->GetLignes() - j - 1, 1, CourbMarge[i]->pts->Element(j, 1));
        }
        AjoutCourbe(*AxePatronP, CourbMarge[i]);
        delete(distance);
    }
    AjoutCourbe(*AxePatronDXFP, CourbPatronDXF[0]);
    if ((Xd[0]->Element(0, 0) != Xd[1]->Element(0, 0)) //test points Av cote 1&2 confondus
            || (Yd[0]->Element(0, 0) != Yd[1]->Element(0, 0))) {
        if (debug) printf ("\n GenerateCourbe Avant");
        CourbAv = new Courbe("Avant");
        CourbAvBack = new Courbe("AvantBack");
        CourbAv->pts = new Matrice(2, 2);
        CourbAvBack->pts = new Matrice(2, 2);

        CourbAv->pts->SetElement(0, 0, Xd[0]->Element(0, 0));
        CourbAvBack->pts->SetElement(1, 0, Xd[0]->Element(0, 0));

        CourbAv->pts->SetElement(1, 0, Xd[1]->Element(0, 0));
        CourbAvBack->pts->SetElement(0, 0, Xd[1]->Element(0, 0));

        CourbAv->pts->SetElement(0, 1, Yd[0]->Element(0, 0));
        CourbAvBack->pts->SetElement(1, 1, Yd[0]->Element(0, 0));

        CourbAv->pts->SetElement(1, 1, Yd[1]->Element(0, 0));
        CourbAvBack->pts->SetElement(0, 1, Yd[1]->Element(0, 0));

        CourbAv->points = OFF;
        CourbAv->symX = OFF;
        CourbAvBack->points = OFF;
        CourbAv->symX = OFF;

        AjoutCourbe(*AxePatronP, CourbAv);
        ajCoAv = 1;

        //marge avant
        distance = Ones(2, 1);
        for (j = 0; j < distance->GetLignes(); j++)
            distance->MultiplyElement(j, 0, margeDeb / 100.0f);

        CourbMargeAv = new Courbe("MargeAV");
        CourbMargeAvBack = new Courbe("MargeAvBack");
        if (debug) printf ("\n CC MargeAv/MargeAvBack");
        CourbMargeAv->pts = CalculContour(CourbAv->pts, distance, +1);
        CourbMargeAvBack->pts = CalculContour(CourbAv->pts, distance, +1);
        for (j = 0; j < CourbMargeAv->pts->GetLignes(); j++) {
            CourbMargeAvBack->pts->SetElement(CourbMargeAv->pts->GetLignes() - j - 1, 0, CourbMargeAv->pts->Element(j, 0));
            CourbMargeAvBack->pts->SetElement(CourbMargeAv->pts->GetLignes() - j - 1, 1, CourbMargeAv->pts->Element(j, 1));
        }

        CourbMargeAv->points = OFF;
        CourbMargeAv->symX = OFF;
        CourbMargeAvBack->points = OFF;
        CourbMargeAvBack->symX = OFF;

        AjoutCourbe(*AxePatronP, CourbMargeAv);
        ajCoMaAv = 1;

        delete(distance);
        for (i = 0; i < 2; i++) {

            CourbCoin1[i] = new Courbe("Coin1");
            CourbCoin = new Courbe("Coin");
            CourbCoin1Back[i] = new Courbe("Coin1Back");

            CourbCoin1[i]->points = OFF;
            CourbCoin->points = OFF;
            CourbCoin1Back[i]->points = OFF;

            CourbCoin1[i]->pts = Zeros(3, 2);
            CourbCoin->pts = Zeros(3, 2);
            CourbCoin1Back[i]->pts = Zeros(3, 2);

            CourbCoin1[i]->pts->SetElement(0, 0, CourbMarge[i]->pts->Element(0, 0));
            CourbCoin->pts->SetElement(0, 0, CourbMarge[i]->pts->Element(0, 0));

            CourbCoin1[i]->pts->SetElement(0, 1, CourbMarge[i]->pts->Element(0, 1));
            CourbCoin->pts->SetElement(0, 1, CourbMarge[i]->pts->Element(0, 1));

            //			CourbCoin1[i]->pts->SetElement(2,0, CourbMargeAv->pts->Element(0,0));
            //			CourbCoin->pts->SetElement(2,0, CourbMargeAv->pts->Element(0,0));
            CourbCoin1[i]->pts->SetElement(2, 0, CourbMargeAv->pts->Element(i, 0));
            CourbCoin->pts->SetElement(2, 0, CourbMargeAv->pts->Element(i, 0));

            //			CourbCoin1[i]->pts->SetElement(2,1, CourbMargeAv->pts->Element(0,1));
            //			CourbCoin->pts->SetElement(2,1, CourbMargeAv->pts->Element(0,1));
            CourbCoin1[i]->pts->SetElement(2, 1, CourbMargeAv->pts->Element(i, 1));
            CourbCoin->pts->SetElement(2, 1, CourbMargeAv->pts->Element(i, 1));

            double x, y;
            Inter2Vecteurs(
                    CourbMarge[i]->pts->Element(1, 0), CourbMarge[i]->pts->Element(1, 1),
                    CourbMarge[i]->pts->Element(0, 0), CourbMarge[i]->pts->Element(0, 1),
                    CourbMargeAv->pts->Element(1, 0), CourbMargeAv->pts->Element(1, 1),
                    CourbMargeAv->pts->Element(0, 0), CourbMargeAv->pts->Element(0, 1),
                    &x, &y);

            CourbCoin1[i]->pts->SetElement(1, 0, x);
            CourbCoin->pts->SetElement(1, 0, x);

            CourbCoin1[i]->pts->SetElement(1, 1, y);
            CourbCoin->pts->SetElement(1, 1, y);

            for (int _i = 0; _i < 3; _i++) {
                CourbCoin1Back[i]->pts -> SetElement(2 - _i, 0, CourbCoin1[i]->pts->Element(_i, 0));
                CourbCoin1Back[i]->pts -> SetElement(2 - _i, 1, CourbCoin1[i]->pts->Element(_i, 1));
            }
            AjoutCourbe(*AxePatronP, CourbCoin);
            ajCo1[i] = 1;
        }
    }
    /*else //relie les marges cote 1&2 par un segment
    {
        if (debug) printf ("\n GenerateCourbe AV");
        CourbAv = new Courbe("AV");
        CourbAv->points = OFF;
        CourbAvBack = new Courbe("AVBack");
        CourbAvBack->points = OFF;
        CourbAv->pts = Zeros(2, 2);
        CourbAvBack->pts = Zeros(2, 2);
        CourbAv->pts->SetElement(0, 0, CourbMarge[0]->pts->Element(0, 0));
        CourbAvBack->pts->SetElement(1, 0, CourbMarge[0]->pts->Element(0, 0));

        CourbAv->pts->SetElement(0, 1, CourbMarge[0]->pts->Element(0, 1));
        CourbAvBack->pts->SetElement(1, 1, CourbMarge[0]->pts->Element(0, 1));

        CourbAv->pts->SetElement(1, 0, CourbMarge[1]->pts->Element(0, 0));
        CourbAvBack->pts->SetElement(0, 0, CourbMarge[1]->pts->Element(0, 0));

        CourbAv->pts->SetElement(1, 1, CourbMarge[1]->pts->Element(0, 1));
        CourbAvBack->pts->SetElement(0, 1, CourbMarge[1]->pts->Element(0, 1));

        AjoutCourbe(*AxePatronP, CourbAv);
        // 1
        ajCoAv = 0;
    }*/
    int n0 = Xd[0]->GetLignes() - 1;
    int n1 = Xd[1]->GetLignes() - 1;
    if ((Xd[0]->Element(n0, 0) != Xd[1]->Element(n1, 0)) //test points Ar cote 1&2 confondus
            || (Yd[0]->Element(n0, 0) != Yd[1]->Element(n1, 0))) {
        if (debug) printf ("\n GenerateCourbe AR");
        CourbAr = new Courbe("AR");
        CourbArDXF = new Courbe("AR");

        CourbAr->pts = new Matrice(2, 2);
        CourbArDXF->pts = new Matrice(2, 2);

        CourbAr->pts->SetElement(0, 0, Xd[0]->Element(n0, 0));
        CourbArDXF->pts->SetElement(0, 0, Xd[0]->Element(n0, 0));

        CourbAr->pts->SetElement(1, 0, Xd[1]->Element(n1, 0));
        CourbArDXF->pts->SetElement(1, 0, Xd[1]->Element(n1, 0));

        CourbAr->pts->SetElement(0, 1, Yd[0]->Element(n0, 0));
        CourbArDXF->pts->SetElement(0, 1, Yd[0]->Element(n0, 0));

        CourbAr->pts->SetElement(1, 1, Yd[1]->Element(n1, 0));
        CourbArDXF->pts->SetElement(1, 1, Yd[1]->Element(n1, 0));

        CourbAr->points = OFF;
        CourbAr->symX = OFF;
        CourbArDXF->points = OFF;
        CourbAr->symX = OFF;

        AjoutCourbe(*AxePatronP, CourbAr);
        ajCoAr = 1;

        distance = Ones(2, 1);

        for (j = 0; j < distance->GetLignes(); j++)
            distance->MultiplyElement(j, 0, margeFin / 100.0f);

        CourbMargeAr = new Courbe("MargeAR");
        CourbMargeArDXF = new Courbe("MargeAR");
        if (debug) printf ("\n CC MargeAr/MargeArDXF");
        CourbMargeAr->pts = CalculContour(CourbAr->pts, distance, -1);
        CourbMargeArDXF->pts = CalculContour(CourbAr->pts, distance, -1);

        CourbMargeAr->points = OFF;
        CourbMargeAr->symX = OFF;
        CourbMargeArDXF->points = OFF;
        CourbMargeArDXF->symX = OFF;

        AjoutCourbe(*AxePatronP, CourbMargeAr);
        ajCoMaAr = 1;

        delete(distance);
        for (i = 0; i < 2; i++) {
            CourbCoin2[i] = new Courbe("Coin2");
            CourbCoin = new Courbe("Coin2");
            CourbCoin2Back[i] = new Courbe("Coin2Back");

            CourbCoin2[i]->points = OFF;
            CourbCoin->points = OFF;
            CourbCoin2Back[i]->points = OFF;

            CourbCoin2[i]->pts = Zeros(3, 2);
            CourbCoin->pts = Zeros(3, 2);
            CourbCoin2Back[i]->pts = Zeros(3, 2);
            if (debug) printf ("\n GenerateCourbe CourbMarge(%d", i);
            n = CourbMarge[i]->pts->GetLignes() - 1;

            CourbCoin2[i]->pts->SetElement(0, 0, CourbMarge[i]->pts->Element(n, 0));
            CourbCoin->pts->SetElement(0, 0, CourbMarge[i]->pts->Element(n, 0));

            CourbCoin2[i]->pts->SetElement(0, 1, CourbMarge[i]->pts->Element(n, 1));
            CourbCoin->pts->SetElement(0, 1, CourbMarge[i]->pts->Element(n, 1));

            //CourbCoin2[i]->pts->SetElement(2,0, CourbMargeAr->pts->Element(0,0));
            //CourbCoin->pts->SetElement(2,0, CourbMargeAr->pts->Element(0,0));
            CourbCoin2[i]->pts->SetElement(2, 0, CourbMargeAr->pts->Element(i, 0));
            CourbCoin->pts->SetElement(2, 0, CourbMargeAr->pts->Element(i, 0));

            //CourbCoin2[i]->pts->SetElement(2,1, CourbMargeAr->pts->Element(0,1));
            //CourbCoin->pts->SetElement(2,1, CourbMargeAr->pts->Element(0,1));
            CourbCoin2[i]->pts->SetElement(2, 1, CourbMargeAr->pts->Element(i, 1));
            CourbCoin->pts->SetElement(2, 1, CourbMargeAr->pts->Element(i, 1));

            double x, y;
            Inter2Vecteurs(
                    CourbMarge[i]->pts->Element(n, 0), CourbMarge[i]->pts->Element(n, 1),
                    CourbMarge[i]->pts->Element(n - 1, 0), CourbMarge[i]->pts->Element(n - 1, 1),
                    CourbMargeAr->pts->Element(1, 0), CourbMargeAr->pts->Element(1, 1),
                    CourbMargeAr->pts->Element(0, 0), CourbMargeAr->pts->Element(0, 1),
                    &x, &y);

            if ( x < CourbMarge[i]->pts->Element(n, 0)) {
                if (debug) printf ("\n\n x < ... %f\n\n", 1000.0f* (CourbMarge[i]->pts->Element(n, 0)-x));
                CourbMarge[i]->pts->SetElement(n, 0, x);
                CourbMarge[i]->pts->SetElement(n, 1, y);
                CourbMargeDXF[i]->pts->SetElement(n, 0, x);
                CourbMargeDXF[i]->pts->SetElement(n, 1, y);
                CourbMargeBack[i]->pts->SetElement(0, 0, x);
                CourbMargeBack[i]->pts->SetElement(0, 1, y);

                CourbCoin2[i]->pts->SetElement(0, 0, CourbMarge[i]->pts->Element(n, 0));
                CourbCoin2[i]->pts->SetElement(0, 1, CourbMarge[i]->pts->Element(n, 1));
                CourbCoin->pts->SetElement(0, 0, CourbMarge[i]->pts->Element(n, 0));
                CourbCoin->pts->SetElement(0, 1, CourbMarge[i]->pts->Element(n, 1));
            }


            CourbCoin2[i]->pts->SetElement(1, 0, x);
            CourbCoin->pts->SetElement(1, 0, x);

            CourbCoin2[i]->pts->SetElement(1, 1, y);
            CourbCoin->pts->SetElement(1, 1, y);

            for (int _i = 0; _i < 3; _i++) {
                CourbCoin2Back[i]->pts -> SetElement(2 - _i, 0, CourbCoin2[i]->pts->Element(_i, 0));
                CourbCoin2Back[i]->pts -> SetElement(2 - _i, 1, CourbCoin2[i]->pts->Element(_i, 1));
            }

            AjoutCourbe(*AxePatronP, CourbCoin);
            ajCo2[i] = 1;
        }
    } /*else //relie les marges cote 1&2 par un segment
    {
        CourbAr = new Courbe("AR");
        CourbAr->points = OFF;
        CourbArDXF = new Courbe("AR");
        CourbArDXF->points = OFF;

        CourbAr->pts = Zeros(2, 2);
        CourbArDXF->pts = Zeros(2, 2);

        CourbAr->pts->SetElement(0, 0, CourbMarge[0]->pts->Element(n0, 0));
        CourbArDXF->pts->SetElement(0, 0, CourbMarge[0]->pts->Element(n0, 0));

        CourbAr->pts->SetElement(0, 1, CourbMarge[0]->pts->Element(n0, 1));
        CourbArDXF->pts->SetElement(0, 1, CourbMarge[0]->pts->Element(n0, 1));

        CourbAr->pts->SetElement(1, 0, CourbMarge[1]->pts->Element(n1, 0));
        CourbArDXF->pts->SetElement(1, 0, CourbMarge[1]->pts->Element(n1, 0));

        CourbAr->pts->SetElement(1, 1, CourbMarge[1]->pts->Element(n1, 1));
        CourbArDXF->pts->SetElement(1, 1, CourbMarge[1]->pts->Element(n1, 1));

        AjoutCourbe(*AxePatronP, CourbAr);
        // 1
        ajCoAr = 0;
    }*/
    if (ajCoAr) AjoutCourbe(*AxePatronDXFP, CourbArDXF);

    AjoutCourbe(*AxePatronDXFP, CourbPatronBack[1]);
    if (ajCoAv) AjoutCourbe(*AxePatronDXFP, CourbAvBack);


    AjoutCourbe(*AxeMarginDXFP, CourbMargeDXF[0]);
    if (ajCo2[0]) AjoutCourbe(*AxeMarginDXFP, CourbCoin2[0]);

    if (ajCoMaAr) AjoutCourbe(*AxeMarginDXFP, CourbMargeArDXF);
    if (ajCo2[1]) AjoutCourbe(*AxeMarginDXFP, CourbCoin2Back[1]);


    AjoutCourbe(*AxeMarginDXFP, CourbMargeBack[1]);
    if (ajCo1[1]) AjoutCourbe(*AxeMarginDXFP, CourbCoin1[1]);
    if (ajCoMaAv) AjoutCourbe(*AxeMarginDXFP, CourbMargeAvBack);
    if (ajCo1[0]) AjoutCourbe(*AxeMarginDXFP, CourbCoin1Back[0]);

    if (debug) printf ("\n GC go in Reper()");
    int ventisave=0;
    double xk0, xk1, xk2, yk0, yk1, yk2, xkf, ykf;
    if (makeRep)
    for (i = 0; i < 2; i++) {
        if (debug) printf ("\n i=%d\n", i);
            for (j = 0;j < 5; j++) {
                if (debug) printf (" set j=%d ", j);
                xSuspente[i][j] = -100000.0f;
                ySuspente[i][j] = -100000.0f;
            }

			Matrice* m = 0;

			if (gfd->CorrectRepPoints && isPince && (coeff[i] != -1)) {
				printf ("\n coeff=%f", coeff[i]);
				m = getReperPointsPince (gfd, Xd0[i], Yd0[i], coeff[i], 
										Xd[i], Yd[i], P[i], nerv[i], Deb[i], FaceDeb[i], Fin[i], FaceFin[i], (nerv1==nerv2));
			} else {
				m = getReperPoints(gfd, Xd[i], Yd[i], P[i], nerv[i], Deb[i], FaceDeb[i], Fin[i], FaceFin[i], (nerv1==nerv2));
				if (coeff[i] == -1) 
					printf ("\nlength XY=%6.2f", 1000*calculCourbeLength(Xd[i], Yd[i]));
			} 
				

            if (debug) for (int _i=0; _i<m->GetLignes(); _i++) {
                printf ("\n grp %d (%f, %f, %f)", _i, m->Element(_i,0), m->Element(_i,1), m->Element(_i,2));
            }
            //ajout croix au graphe
            for (j = 0;j<m->GetLignes();j++) {
                if (debug)    printf ("\nj3=%d", j);
                    int mark = m->Element(j,2);
                    double xSus = m->Element(j,0);
                    double ySus = m->Element(j,1);

                    if (((mark == REP_TRIANGLE) || (mark == REP_V)) && (m->Element(j, 3) > 0.0f)) {
                        ventisave = (int) m->Element(j, 3) - 1;
                        xSuspente[i][ventisave] = xSus;
                        ySuspente[i][ventisave] = ySus;
                    }
                    if ((i==0) && (nerv1==nerv2)) {
                        if (mark == REP_0) {
                            xk0 = xSus;
                            yk0 = ySus;
                            //printf ("\n REP_0 (%f, %f)", xk0, yk0);
                        }
                        if (mark == REP_VENTDEB) {
                            xk1 = xSus;
                            yk1 = ySus;
                            //printf ("\n REP_VENTDEB (%f, %f)", xk1, yk1);
                        }
                        if (mark == REP_VENTFIN) {
                            xk2 = xSus;
                            yk2 = ySus;
                            //printf ("\n REP_VENTFIN (%f, %f)", xk2, yk2);
                        }
                        if (mark == REP_KLAPAN_FIN) {
                            xkf = xSus;
                            ykf = ySus;
                            //printf ("\n REP_KLAPAN_FIN (%f, %f)", xkf, ykf);
                        }
                    }


                    if ((mark == REP_CROSS) && (gfd->ReperesProfile[i])) {
                            CourbRep = new Courbe("Repers points");
                            CourbRep->points = OFF;
                            CourbRepDXF = new Courbe("Repers points");
                            CourbRepDXF->points = OFF;
                            CourbRep->pts = Zeros(2, 2);
                            CourbRepDXF->pts = Zeros(2, 2);
                            CourbRep->pts->SetElement(0, 0, xSus - gfd->XMashtab * 0.01f);
                            CourbRep->pts->SetElement(0, 1, ySus - gfd->XMashtab * 0.01f);
                            CourbRep->pts->SetElement(1, 0, xSus + gfd->XMashtab * 0.01f);
                            CourbRep->pts->SetElement(1, 1, ySus + gfd->XMashtab * 0.01f);
                            int _i, _j;
                            for (_i = 0; _i < 2; _i++)
                                for (_j = 0; _j < 2; _j++)
                                    CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                            AjoutCourbe(*AxeRepDXFP, CourbRepDXF);
                            CourbRep = new Courbe("Repers points");
                            CourbRep->points = OFF;
                            CourbRepDXF = new Courbe("Repers points");
                            CourbRepDXF->points = OFF;
                            CourbRep->pts = Zeros(2, 2);
                            CourbRepDXF->pts = Zeros(2, 2);
                            CourbRep->pts->SetElement(0, 0, xSus - gfd->XMashtab * 0.01f);
                            CourbRep->pts->SetElement(0, 1, ySus + gfd->XMashtab * 0.01f);
                            CourbRep->pts->SetElement(1, 0, xSus + gfd->XMashtab * 0.01f);
                            CourbRep->pts->SetElement(1, 1, ySus - gfd->XMashtab * 0.01f);
                            for (_i = 0; _i < 2; _i++)
                                for (_j = 0; _j < 2; _j++)
                                    CourbRepDXF->pts->SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                            AjoutCourbe(*AxeRepDXFP, CourbRepDXF);
                        } 
                        if ((mark==REP_LINE) || ((gfd->ReperesProfile[i]) && (mark == REP_MIDDLE_LINE))) {
                                char courbName[255];
                                if (mark==REP_MIDDLE_LINE) 
                                    strcpy(courbName, "Repers points");
                                else
                                    strcpy(courbName, "Repers diag");
                                //perpendiular
                                int direction = -1;
                                if (i == 1) direction = 1;
                                double _x0 = xSus;
                                double _y0 = ySus;
                                int ix, icol;
                                Ind(_x0, Xd[i], &ix, &icol);
                                //Ind(ySuspente[i][j], interpYSuspente, &iy, &icol);
                                double xc, yc;
                                double _l = 3 * gfd->XMashtab * 0.01f;
                                if (mark == REP_MIDDLE_LINE) {
                                    _l = 0.4 * _l;
                                }
                                if ((Xd[i]->Element(ix, 0) > _x0) && (ix != 0)) ix--;
                                if ((Xd[i]->Element(ix,0) == _x0) && (ix != 0)) ix--;
                                if (ix != 0)
                                    CalculVecteurNormal(Xd[i]->Element(ix,0), Yd[i]->Element(ix, 0),_x0, _y0, &xc, &yc, _l, direction);
                                else
                                    CalculVecteurNormal(_x0, _y0,Xd[i]->Element(1,0), Yd[i]->Element(1, 0), &xc, &yc, _l, direction);
                                double _x1 = _x0;
                                double _y1 = _y0;
                                double _x2 = xc;
                                double _y2 = yc;
                                if (ix == 0) {
                                    _x2 = _x0 + (xc - Xd[i]->Element(1,0));
                                    _y2 = _y0 + (yc - Yd[i]->Element(1,0));
                                }
                                int _i, _j;
                                CourbRep = new Courbe(courbName);
                                CourbRep->points = OFF;
                                CourbRepDXF = new Courbe(courbName);
                                CourbRepDXF->points = OFF;
                                CourbRep->pts = Zeros(2, 2);
                                CourbRepDXF->pts = Zeros(2, 2);
                                CourbRep->pts->SetElement(0, 0, _x1);
                                CourbRep->pts->SetElement(0, 1, _y1);
                                CourbRep->pts->SetElement(1, 0, _x2);
                                CourbRep->pts->SetElement(1, 1, _y2);
                                for (_i = 0; _i < 2; _i++)
                                    for (_j = 0; _j < 2; _j++)
                                        CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                                AjoutCourbe(*AxeRepDXFP, CourbRepDXF);
                            }
                            if ((gfd -> ReperesSuspentes[i]) && ((mark == REP_TRIANGLE) || (mark==REP_V)) ) {
                                int direction = -1;
                                if (i == 1) direction = 1;
                                double _x0 = xSus;
                                double _y0 = ySus;
                                double _l = 2*gfd->XMashtab * 0.01f;
                                double _x1 = _x0;
                                double _y1 = _y0;

                                int ix, icol;
                                Ind(xSus, Xd[i], &ix, &icol);
                                double xc, yc;
                                if ((Xd[i]->Element(ix,0) > _x0) && (ix != 0 )) ix--;
                                if ((Xd[i]->Element(ix,0) == _x0) && (ix != 0) ) ix--;
                                if  (ix != 0)
                                    CalculVecteurNormal(Xd[i]->Element(ix,0), Yd[i]->Element(ix, 0), _x0, _y0, &xc, &yc, _l, direction);
                                else
                                    CalculVecteurNormal( _x0, _y0, Xd[i]->Element(1,0), Yd[i]->Element(1, 0), &xc, &yc, _l, direction);
                                if (ix == 0) ix=1;

                                double nx = (Xd[i]->Element(ix, 0)-_x0);
                                double ny = (Yd[i]->Element(ix, 0)-_y0);
                                double rasst = sqrt (sqr(nx) + sqr(ny));
                                nx=_l/2.0f*nx/rasst;
                                ny=_l/2.0f*ny/rasst;

                                double _x2 = xc + nx;
                                double _y2 = yc + ny;

                                double _x3 = xc - nx;
                                double _y3 = yc - ny;
                                int _i, _j;
                                CourbRep = new Courbe("Repers suspente");
                                CourbRep->points = OFF;
                                CourbRepDXF = new Courbe("Repers suspente");
                                CourbRepDXF->points = OFF;
                                CourbRep->pts = Zeros(2, 2);
                                CourbRepDXF->pts = Zeros(2, 2);
                                CourbRep->pts->SetElement(0, 0, _x1);
                                CourbRep->pts->SetElement(0, 1, _y1);
                                CourbRep->pts->SetElement(1, 0, _x2);
                                CourbRep->pts->SetElement(1, 1, _y2);
                                for (_i = 0; _i < 2; _i++)
                                    for (_j = 0; _j < 2; _j++)
                                        CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));

                                AjoutCourbe(*AxeRepDXFP, CourbRepDXF);

                                if (mark != REP_V) {
                                    CourbRep = new Courbe("Repers suspente");
                                    CourbRep->points = OFF;
                                    CourbRepDXF = new Courbe("Repers suspente");
                                    CourbRepDXF->points = OFF;
                                    CourbRep->pts = Zeros(2, 2);
                                    CourbRepDXF->pts = Zeros(2, 2);
                                    CourbRep->pts->SetElement(0, 0, _x2);
                                    CourbRep->pts->SetElement(0, 1, _y2);
                                    CourbRep->pts->SetElement(1, 0, _x3);
                                    CourbRep->pts->SetElement(1, 1, _y3);
                                    for (_i = 0; _i < 2; _i++)
                                        for (_j = 0; _j < 2; _j++)
                                            CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                                    AjoutCourbe(*AxeRepDXFP, CourbRepDXF);
                                }
                                CourbRep = new Courbe("Repers suspente");
                                CourbRep->points = OFF;
                                CourbRepDXF = new Courbe("Repers suspente");
                                CourbRepDXF->points = OFF;
                                CourbRep->pts = Zeros(2, 2);
                                CourbRepDXF->pts = Zeros(2, 2);
                                CourbRep->pts->SetElement(0, 0, _x3);
                                CourbRep->pts->SetElement(0, 1, _y3);
                                CourbRep->pts->SetElement(1, 0, _x1);
                                CourbRep->pts->SetElement(1, 1, _y1);
                                for (_i = 0; _i < 2; _i++)
                                    for (_j = 0; _j < 2; _j++)
                                        CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                                AjoutCourbe(*AxeRepDXFP, CourbRepDXF);
                            }
            }

            if ((i==0) && (nerv1==nerv2) && (gfd->RaskladKlapans) && (gfd->VentHoles)) {
                //printf ("\n go make shov!");
                bool match = false;
				int _i, _j;
                if ((nerv1==-1) || (nerv1==0)){
                        if (gfd->VentCentralNerv) match = true;
                }
                else {
                        if ((gfd->noNervVH[nerv1] == 1) || (gfd->noNervVH[nerv1 - 1]==1)) match = true;
                }

                if (match) {
                    //printf (" \n MATCH true!");
                    if (gfd->VentHolesDouble) {
                        // xk0, yk0 -> xkfin, ykfin
                        CourbRep = new Courbe("Klapan line");
                        CourbRep->points = OFF;
                        CourbRepDXF = new Courbe("Klapan line");
                        CourbRepDXF->points = OFF;
                        CourbRep->pts = Zeros(2, 2);
                        CourbRepDXF->pts = Zeros(2, 2);
                        CourbRep->pts->SetElement(0, 0, xk0);
                        CourbRep->pts->SetElement(0, 1, yk0);
                        CourbRep->pts->SetElement(1, 0, xkf);
                        CourbRep->pts->SetElement(1, 1, 0.0f);
                        for (_i = 0; _i < 2; _i++)
                            for (_j = 0; _j < 2; _j++)
                                CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                        AjoutCourbe(*AxeRepDXFP, CourbRepDXF);
                    }
                    // xk1, yk1 -> xkfin, ykfin
                    CourbRep = new Courbe("Klapan line");
                    CourbRep->points = OFF;
                    CourbRepDXF = new Courbe("Klapan line");
                    CourbRepDXF->points = OFF;
                    CourbRep->pts = Zeros(2, 2);
                    CourbRepDXF->pts = Zeros(2, 2);
                    CourbRep->pts->SetElement(0, 0, xk1);
                    CourbRep->pts->SetElement(0, 1, yk1);
                    CourbRep->pts->SetElement(1, 0, xkf);
                    CourbRep->pts->SetElement(1, 1, 0.0f);
                    for (_i = 0; _i < 2; _i++)
                        for (_j = 0; _j < 2; _j++)
                            CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                    AjoutCourbe(*AxeRepDXFP, CourbRepDXF);
                    // xk2, yk2 -> xkfin, ykfin
                    CourbRep = new Courbe("Klapan line");
                    CourbRep->points = OFF;
                    CourbRepDXF = new Courbe("Klapan line");
                    CourbRepDXF->points = OFF;
                    CourbRep->pts = Zeros(2, 2);
                    CourbRepDXF->pts = Zeros(2, 2);
                    CourbRep->pts->SetElement(0, 0, xk2);
                    CourbRep->pts->SetElement(0, 1, yk2);
                    CourbRep->pts->SetElement(1, 0, xkf);
                    CourbRep->pts->SetElement(1, 1, 0.0f);
                    for (_i = 0; _i < 2; _i++)
                        for (_j = 0; _j < 2; _j++)
                            CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                    AjoutCourbe(*AxeRepDXFP, CourbRepDXF);
                } else {
                    printf (" \n NOOOOT MATCH!");
                }

            }

    }

    if (debug) printf ("\n... in Reper()");

    if (Ventilation == 1) {
      if (debug)   printf ("\nif Ventil");
        for (i = 1; i < 4; i++) {
            if ((xSuspente[0][i] != -100000.0f)
                    && (xSuspente[1][i] != -100000.0f)
                    && (xSuspente[0][i - 1] != -100000.0f)
                    && (xSuspente[1][i - 1] != -100000.0f)) {

                CourbCercle = new Courbe("Cercle");
                CourbCercleDXF = new Courbe("Cercle");

                CourbCercle->points = OFF;
                CourbCercleDXF->points = OFF;

                CourbCercle->pts = Cercle(
                        (xSuspente[0][i] + xSuspente[1][i] + xSuspente[0][i - 1] + xSuspente[1][i - 1]) / 4.0,
                        (ySuspente[0][i] + ySuspente[1][i] + ySuspente[0][i - 1] + ySuspente[1][i - 1]) / 4.0,
                        (ySuspente[0][i] + ySuspente[0][i - 1] - ySuspente[1][i] - ySuspente[1][i - 1]) / 8.0,
                        36);
                CourbCercleDXF->pts = Cercle(
                        (xSuspente[0][i] + xSuspente[1][i] + xSuspente[0][i - 1] + xSuspente[1][i - 1]) / 4.0,
                        (ySuspente[0][i] + ySuspente[1][i] + ySuspente[0][i - 1] + ySuspente[1][i - 1]) / 4.0,
                        (ySuspente[0][i] + ySuspente[0][i - 1] - ySuspente[1][i] - ySuspente[1][i - 1]) / 8.0,
                        36);

                AjoutCourbe(*AxePatronP, CourbCercle);
                AjoutCourbe(*AxeCercleDXFP, CourbCercleDXF);
            }

        }
      if (debug)   printf ("\n...if Ventil");
    }

    //if (Numerotation == 1) {
    xText = (Xd[0]->Element(0, 0) + Xd[0]->Element(Xd[0]->GetLignes() - 1, 0)
            + Xd[1]->Element(0, 0) + Xd[1]->Element(Xd[1]->GetLignes() - 1, 0)) / 4.0f;

    double x0len = Xd[0]->Element(Xd[0]->GetLignes() - 1, 0) - Xd[0]->Element(0, 0);
    double x1len = Xd[1]->Element(Xd[1]->GetLignes() - 1, 0) - Xd[1]->Element(0, 0);
    double xlen = (x0len + x1len) / 2.0f;

    //		printf ("\nXd[0]->Element(0,0) %f",Xd[0]->Element(0,0));
    //		printf ("\nXd[0]->Element(Xd[0]->GetLignes()-1,0) %f", Xd[0]->Element(Xd[0]->GetLignes()-1,0));
    //		printf ("\nXd[1]->Element(0,0) %f", Xd[1]->Element(0,0));
    //		printf ("\nXd[1]->Element(Xd[1]->GetLignes()-1,0) %f", Xd[1]->Element(Xd[1]->GetLignes()-1,0));

    xText = gfd->textX * xlen;
    yText = (Yd[0]->Element(0, 0) + Yd[0]->Element(Yd[0]->GetLignes() - 1, 0)
            + Yd[1]->Element(0, 0) + Yd[1]->Element(Yd[1]->GetLignes() - 1, 0)) / 4.0f;

    double y0len = Yd[1]->Element(0, 0) - Yd[0]->Element(0, 0);
    double y1len = Yd[1]->Element(Yd[1]->GetLignes() - 1, 0) - Yd[0]->Element(Yd[0]->GetLignes() - 1, 0);
    double ylen = (y0len + y1len) / 2.0f;

    yText = Yd[1]->Element(0, 0) - gfd->textY*ylen;

    //		printf ("\nYd[0]->Element(0,0) %f",Yd[0]->Element(0,0));
    //		printf ("\nYd[0]->Element(Yd[0]->GetLignes()-1,0) %f", Yd[0]->Element(Yd[0]->GetLignes()-1,0));
    //		printf ("\nYd[1]->Element(0,0) %f", Yd[1]->Element(0,0));
    //		printf ("\nYd[1]->Element(Yd[1]->GetLignes()-1,0) %f", Yd[1]->Element(Yd[1]->GetLignes()-1,0));
    //sprintf(texte, "%d%c%03.1f%c%03.1f_%d%c%03.1f%c%03.1f",
    //	NoNerv[0],texteExtInt[FaceDeb[0]-1],Deb[0],texteExtInt[FaceFin[0]-1],Fin[0],
    //	NoNerv[1],texteExtInt[FaceDeb[1]-1],Deb[1],texteExtInt[FaceFin[1]-1],Fin[1]);
    sprintf(texte, "%s", text);
    AjoutTexte(*AxePatronP, texte, 0.02f, 0.0f, xText, yText);
    sprintf(texte, "%s", text);
    AjoutTexte(*AxePatronTextDXFP, texte, 0.02f, 0.0f, xText, yText);
    if (debug) printf ("\n...GenerateCourbe()");
    //}
    //    printf ("\n ...GenerateCourbe()");
}
