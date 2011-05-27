/**************************************
* gestion chargement divers fichiers
* forme ,profil, patrons ...
***************************************/

#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <string>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

#include "fichier.h"
#include "profil.h"
#include "geom.h"
#include "design.h"

#define sqr(f1) ((f1)*(f1))
#define CHISLOPI	3.141592675f

#ifndef DEBUG
#define DEBUG false
#endif

void writeFichierPolyDXFDelta(FILE *fid, TAxe *axe, TAxe *axe2, int rep, TAxe *axeR, int vent, TAxe *axeC, int num, TAxe *axeT, double dx, double dy, double dz, int n );


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
/* readFichierProfil */
/************************/

KiteDesign* readKiteDesignFromFile(const char* FilePath) {
	cout << "readKiteDesignFromFile: " << FilePath << endl;
	KiteDesign* kd = new KiteDesign();
	ifstream in(FilePath);
	int n=0;
	string word;
	in >> n;
	//cout <<"n:"<< n << endl;
	kd->n_elements=n;
	kd->kiteDesignElements = new KiteDesignElement*[n];

	for (int i = 0; i < n; i++) {
		in >> word;
		//cout << "word:" << word << endl;
		if (word.compare("LINE") == 0) {
			//cout << "It's LINE" << endl;
			Line* line = new Line(in);
			//line->print();
			kd->kiteDesignElements[i]=line;
		}
	}

	kd->colorTable = new ColorTable(in);
	kd->colorTable->print();
	return kd;
}

void readFichierProfil(const char* NomProf, Matrix** extrados, Matrix** intrados)
{
	FILE *fid;
	char line[255];
	int NbInt, NbExt, i;
	float x,y,coeff;
	//radde pour format Selig ...

	int nbData, iNez;
	Matrix *data = NULL;

	/* ouverture fichier profil */
	printf("\nread fichier profil: '%s'", NomProf);
	if( (fid = fopen( NomProf, "rt" )) == NULL )
	{
		printf( "\nErreur ouverture fichier: '%s'", NomProf );
		char ex[100];
		sprintf (ex, " could not open `%s`", NomProf);
		throw ex;
	}
	else
	{
		/* read entete du fichier */
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
		
			//1ere read pour determiner le nombre de data
			nbData = 0;
			while( fscanf(fid,"%f %f", &x, &y) != EOF ) nbData++;
			//read des datas
			rewind(fid);
			fgets( line, 255, fid ); //premiere ligne
			data = new Matrix(nbData, 2);
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
			*extrados=new Matrix(NbExt,2);

			for(i=0; i<NbExt; i++)
			{
				(*extrados)->SetElement(NbExt-i-1,0,data->Element(i,0));
				(*extrados)->SetElement(NbExt-i-1,1,data->Element(i,1));
			}
			//recopie data dans le tableau d'intrados
			delete(*intrados); NbInt = nbData-iNez;
			*intrados=new Matrix(NbInt,2);
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
			// read nombre de points en intrados et extrados
			rewind(fid);
			fgets( line, 255, fid ); //premiere ligne
			fscanf(fid,"%f %f", &x, &y);
			NbExt=(int)floor(x);	NbInt=(int)floor(y);		
			// read pts Extrados
			delete(*extrados);
			*extrados=new Matrix(NbExt,2);
			for(i=0; i<NbExt; i++)
        		{
				float f1, f2;
				fscanf(fid,"%f %f", &f1, &f2 );
				(*extrados)->SetElement(i,0,f1);
				(*extrados)->SetElement(i,1,f2);
                //printf ("\n extr f1=%f1 f2=%f2", f1, f2);
			}
			// read pts Intrados
			delete(*intrados);
			*intrados=new Matrix(NbInt,2);
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
			printf("\nProbleme fermeture du fichier '%s'", NomProf);
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

int* readFichierVentHoles(char* NomFic, int* quant, int* central) {
    FILE *fid;
    int* m = new int[1000];
    for (int i = 0; i < 1000; i++) m[i] = 0;
    int n;
    int nerv;
    (*central) = 0;
    if( (fid = fopen( NomFic, "rt" )) == NULL )
    {
            printf( "\nErreur ouverture fichier '%s'", NomFic);
			char ex[100];
			sprintf (ex, " could not open `%s`", NomFic);
			throw (ex);
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
/* readFichierReperPoints */
/***********************/

Matrix** readFichierReperPoints(char* NomFic) {
    FILE *fid;
    Matrix** m= new Matrix*[3];
    int n;
    float p;
    if( (fid = fopen( NomFic, "rt" )) == NULL )
    {
            printf( "\nErreur ouverture fichier '%s'", NomFic);
			char ex[100];
			sprintf (ex, " could not open `%s`", NomFic);
			throw (ex);
    } else {
        for (int i=1;i<3;i++) {
            fscanf(fid,"%d", &n);
            //printf ("\n n=%d", n);
            m[i] = new Matrix(n, 1);
            for (int j = 0; j < n; j++) {
                fscanf (fid, "%f", &p);
                //printf ("\n p=%f", p);
                m[i]->SetElement(j,0,p);
            }
        }
    }
    return m;
}

int* readFichierDiagNervs(char* NomFic, int* quant) {
    FILE *fid;
    int* res;
    int n;
    if( (fid = fopen( NomFic, "rt" )) == NULL )
    {
            printf( "\nErreur ouverture fichier '%s'", NomFic);
			char ex[100];
			sprintf (ex, " could not open `%s`", NomFic);
			throw (ex);
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

WindPatternsProject* readWindPatternsProject(char* NomFic) {
    FILE *fid;
    WindPatternsProject *wpp = NULL;
    if( (fid = fopen( NomFic, "rt" )) == NULL )
    {
            printf( "\nErreur ouverture fichier '%s'", NomFic);
    }
    else
    {
    char* motsClef[55] =
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

                "XMASHTAB", "TEXTX", "TEXTY", "FORME", "BALLONEMENT"}; //3


        wpp = new WindPatternsProject();
        strcpy(wpp->fromFile, NomFic);
        char s[255];
        for (int i = 0; i < 55; i++) {
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
                    case 29: {  fgets (s, 255, fid);  trim(s); trim(s, '\r'); trim(s, '\n'); strcpy(wpp->fileNameDiagNerv, s);  break;}
                    case 30: fscanf(fid,"%f", &(wpp->PosDiagNerv2A)); break;
                    case 31: fscanf(fid,"%f", &(wpp->PosDiagNerv2F)); break;
                    case 32: fscanf(fid,"%f", &(wpp->PosDiagNerv1A)); break;
                    case 33: fscanf(fid,"%f", &(wpp->PosDiagNerv1F)); break;

                    case 34: fscanf(fid,"%d", &(wpp->VentHoles)); break;
                    case 35: {  fgets (s, 255, fid); trim(s); trim(s, '\r'); trim(s, '\n'); strcpy(wpp->fileNameVentHoles, s); break;}
                    case 36: fscanf(fid,"%f", &(wpp->VentHolesDeb)); break;
                    case 37: fscanf(fid,"%f", &(wpp->VentHolesFin)); break;

                    case 38: fscanf(fid,"%d", &(wpp->LayoutKlapans)); break;
                    case 39: fscanf(fid,"%d", &(wpp->VentHolesDouble)); break;
                    case 40: fscanf(fid,"%f", &(wpp->PosKlapanFin)); break;

                    case 41: fscanf(fid,"%d", &(wpp->RepPoints)); break;
                    case 42: fscanf(fid,"%d", &(wpp->ReperPointsFromFile)); break;
                    case 43: {  fgets (s, 255, fid); trim(s); trim(s, '\r'); trim(s, '\n'); strcpy(wpp->fileNameRepPoints, s);  break;}

                    case 44: fscanf(fid,"%d", &(wpp->VentilationLayout)); break;

                    case 45: fscanf(fid,"%f", &(wpp->margeFinExt)); break;
                    case 46: fscanf(fid,"%f", &(wpp->margeFinInt)); break;
                    case 47: fscanf(fid,"%f", &(wpp->margeFinNerv)); break;
                    case 48: fscanf(fid,"%f", &(wpp->margeFinDiagNerv)); break;


                    case 49: fscanf(fid,"%f", &(wpp->tochnostLayout2)); break;
                    case 50: fscanf(fid,"%f", &(wpp->XMashtab)); break;
                    case 51: fscanf(fid,"%f", &(wpp->textX)); break;
                    case 52: fscanf(fid,"%f", &(wpp->textY)); break;
                    case 53: {  fgets (s, 255, fid); trim(s); trim(s, '\r'); trim(s, '\n'); strcpy(wpp->fileNameForm, s);  break;}
					case 54: {  fgets (s, 255, fid); trim(s); trim(s, '\r'); trim(s, '\n'); strcpy(wpp->ballonementPath, s);  break;}
                    default:;
                }
               
            }
        }

	}

    if(fclose(fid))
    {
        printf("\nProbleme  la fermeture du fichier");
    }
    return wpp;
}

void writeWindPatternsProject(char *fileName, WindPatternsProject *wpp) {
	FILE *fid;
	/**** message ****/
	printf("\nwrite fichier WindPatternsProject: '%s'",fileName);
	/**** ouverture fichier en write ****/
	if( (fid = fopen( fileName, "wt" )) == NULL )
	{
		printf( "\nErreur ouverture fichier '%s'", fileName);
		exit(0);
	}
    char* motsClef[55] =
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

                "XMASHTAB", "TEXTX", "TEXTY", "FORME", "BALLONEMENT"}; //3


					fprintf(fid,"\n%s %s", motsClef[0], fileName);
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
                    fprintf(fid,"\n%s %s", motsClef[29],wpp->fileNameDiagNerv);
                    fprintf(fid,"\n%s %f", motsClef[30],wpp->PosDiagNerv2A);
                    fprintf(fid,"\n%s %f", motsClef[31],wpp->PosDiagNerv2F);
                    fprintf(fid,"\n%s %f", motsClef[32],wpp->PosDiagNerv1A);
                    fprintf(fid,"\n%s %f", motsClef[33],wpp->PosDiagNerv1F);

                    fprintf(fid,"\n\n%s %d", motsClef[34],wpp->VentHoles);
                    fprintf(fid,"\n%s %s", motsClef[35],wpp->fileNameVentHoles);
                    fprintf(fid,"\n%s %f", motsClef[36],wpp->VentHolesDeb);
                    fprintf(fid,"\n%s %f", motsClef[37],wpp->VentHolesFin);

                    fprintf(fid,"\n\n%s %d", motsClef[38],wpp->LayoutKlapans);
                    fprintf(fid,"\n%s %d", motsClef[39],wpp->VentHolesDouble);
                    fprintf(fid,"\n%s %f", motsClef[40],wpp->PosKlapanFin);

                    fprintf(fid,"\n\n%s %d", motsClef[41],wpp->RepPoints);
                    fprintf(fid,"\n%s %d", motsClef[42],wpp->ReperPointsFromFile);
                    fprintf(fid,"\n%s %s", motsClef[43],wpp->fileNameRepPoints);

                    fprintf(fid,"\n\n%s %d", motsClef[44],wpp->VentilationLayout);

                    fprintf(fid,"\n%s %f", motsClef[45],wpp->margeFinExt);
                    fprintf(fid,"\n%s %f", motsClef[46],wpp->margeFinInt);
                    fprintf(fid,"\n%s %f", motsClef[47],wpp->margeFinNerv);
                    fprintf(fid,"\n%s %f", motsClef[48],wpp->margeFinDiagNerv);

                    fprintf(fid,"\n%s %f", motsClef[49],wpp->tochnostLayout2);
                    fprintf(fid,"\n%s %f", motsClef[50],wpp->XMashtab);
                    fprintf(fid,"\n%s %f", motsClef[51],wpp->textX);
                    fprintf(fid,"\n%s %f", motsClef[52],wpp->textY);
                    fprintf(fid,"\n%s %s", motsClef[53],wpp->fileNameForm);
                    fprintf(fid,"\n%s %s", motsClef[54],wpp->ballonementPath);

	if(fclose(fid))
	{
		printf("\nProbleme  la fermeture du fichier WindPatternsProject");
		exit(0);
	}
    
}


Ballonement* readBallonementFromFile(char* NomFic) {
	FILE *fid;
	Matrix *m = NULL;
	Ballonement* bal = new Ballonement();
	//kChord, kMf, wN, dyw
	Matrix *kChord, *kMf, *wN, *dyw, *powerTail;
	if( (fid = fopen( NomFic, "rt" )) == NULL )
	{
		printf( "\nErreur ouverture fichier '%s'", NomFic);
		char ex[100];
		sprintf (ex, " could not open `%s`", NomFic);
		throw (ex);
	}
	else
	{
		// readBallonement
		int n=-1, ind=0;
		fscanf(fid, "%d", &n);
		double _kChord=0, _kMf=0, _wN=0, _dyw=0, _powerTail=0;
		kChord = new Matrix(n, 1);
		kMf = new Matrix(n, 1);
		wN = new Matrix(n, 1);
		dyw = new Matrix(n, 1);
		powerTail = new Matrix(n, 1);

		for (int i = 0; i < n; i ++) {
			fscanf(fid,"%d %lf %lf %lf %lf %lf", &ind, &_kChord, &_kMf, &_wN, &_dyw, &_powerTail );
			kChord->SetElement(i, 0, _kChord);
			kMf->SetElement(i, 0, _kMf);
			wN->SetElement(i, 0, _wN);
			dyw->SetElement(i, 0, _dyw);
			powerTail->SetElement(i, 0, _powerTail);
		}
		bal->kChord = kChord;
		bal->kMf = kMf;
		bal->wN = wN;
		bal->dyw = dyw;
		bal->powerTail = powerTail;
	}
	if(fclose(fid))
	{
		printf("\nProbleme  reading ballonement");
		exit(0);
	}
	return bal;
}

Form* readFichierForm(char* NomFic)
{
//    if (DEBUG) printf ("\n readFichierForm");
	int i,j,n;
	FILE *fid;
	Matrix *m=NULL;
	double noVer=0.0f;
	Form *f=NULL;
	//liste des mots clefs du fichier
	char* motsClef[19] = {"VERSION",
		"NB_ALVEOLES", "COEFF_PROG_GEOM", "COEFF_EXP", "EPAI_REL_CENTRE",
		"BORD_ATTAQUE_FORME","BORD_DE_FUITE_FORME","SUSPENTAGE_LIGNE_A",
		"SUSPENTAGE_LIGNE_B","SUSPENTAGE_LIGNE_C","SUSPENTAGE_LIGNE_D", "SUSPENTAGE_LIGNE_E",
		"DIEDRE","MORPHING","VRILLAGE","EPAISSEUR_RELATIVE",
		"PROFIL_CENTRE", "PROFIL_BOUT", "TABLEAU_FORME" };
	/**** message ****/
	printf("\nread fichier de Form: '%s'",NomFic);
	/**** ouverture fichier en read ****/
	if( (fid = fopen( NomFic, "rt" )) == NULL )
	{
		printf( "\nErreur ouverture fichier '%s'", NomFic);
		char ex[100];
		sprintf (ex, " could not open `%s`", NomFic);
		throw (ex);
	}
	else
	{
		/**** allocation memoire ****/
		f = new Form();
		/**** read des points de controle ****/
		/*Boucle de read*/
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
                                    m = new Matrix(4,2);
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
                                }//boucle read tableau de forme

                                break;
                            }
                                    default:;
                            }//switch mot clef
                    }//test trouve mot clef
		}
		if(fclose(fid))
		{
			printf("\nProbleme  la fermeture du fichier");
		}
	}
	return f;
}


/***********************/
/* readFichierForm2 */
/***********************/

Form* readFichierForm2(char* NomFic)

{
//    if (DEBUG) printf ("\n readFichierForm2()");
	int i,j,n;

	FILE *fid;

	Matrix *m=NULL;

	float noVer;

	Form *f=NULL;

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

	printf("\nread fichier de Form: '%s'",NomFic);

	

	/**** ouverture fichier en read ****/

	if( (fid = fopen( NomFic, "rt" )) == NULL )

	{

		printf( "\nErreur ouverture fichier '%s'", NomFic);

	}

	else

	{

		/**** allocation memoire ****/

		f = new Form();

		f -> courbInput = true;

		int NbNerv = 0;

		/**** read des points de controle ****/

		/*Boucle de read*/

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

					m = new Matrix(NN,2);

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
						}//boucle read tableau de forme

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
/* writeFichierForm */
/************************/

void writeFichierForm(char *fileName, Form *f)

{
	int i;
	FILE *fid;
	char* texte[11]={"BORD_ATTAQUE_FORME","BORD_DE_FUITE_FORME","SUSPENTAGE_LIGNE_A",
		"SUSPENTAGE_LIGNE_B","SUSPENTAGE_LIGNE_C","SUSPENTAGE_LIGNE_D", "SUSPENTAGE_LIGNE_E",
		"DIEDRE","MORPHING","VRILLAGE","EPAISSEUR_RELATIVE"};
	Matrix* m[11];
	/**** message ****/
	printf("\nwrite fichier de Form: '%s'",fileName);
	/**** ouverture fichier en write ****/
	if( (fid = fopen( fileName, "wt" )) == NULL )
	{
		printf( "\nErreur ouverture fichier '%s'", fileName);
		exit(0);
	}

	/**** write No de version ****/
	fprintf(fid,"\nVERSION %1.2f", NO_VERSION);

	/**** write des parametres ****/
	fprintf(fid,"\n\nNB_ALVEOLES %d", f->NbCaiss);
	fprintf(fid,"\nCOEFF_PROG_GEOM %1.5f", f->CoeffProgGeom);
	fprintf(fid,"\nCOEFF_EXP %1.5f", f->CoeffExp);
	fprintf(fid,"\nEPAI_REL_CENTRE %2.3f", f->EpaiRelCent);

	/**** write des points de controle ****/
	/*init pointeur des matrices de coordonnÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¯ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â¿ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â½es des points de controle*/
	m [0]=f->mCtrlNez; m [1]=f->mCtrlFui;
	m [2]=f->mCtrlA; m [3]=f->mCtrlB; m [4]=f->mCtrlC; m [5]=f->mCtrlD;
	m [6]=f->mCtrlE;
	m [7]=f->mCtrlDiedre; m [8]=f->mCtrlMorphing;
	m [9]=f->mCtrlVrillage; m [10]=f->mCtrlEpaiRel;

//	m [6]=f->mCtrlDiedre; m [7]=f->mCtrlMorphing;
//	m [8]=f->mCtrlVrillage; m [9]=f->mCtrlEpaiRel;
	/*boucle write des points de controle*/

	for (i=0; i<11; i++)
		fprintf(fid,"\n%s %1.3f %1.3f  %1.3f %1.3f  %1.3f %1.3f  %1.3f %1.3f",
		texte[i],
		m [i]->Element(0,0), m [i]->Element(0,1),
		m [i]->Element(1,0), m [i]->Element(1,1),
		m [i]->Element(2,0), m [i]->Element(2,1),
		m [i]->Element(3,0), m [i]->Element(3,1));


	/**** write noms des profils ****/
	fprintf(fid,"\n\nPROFIL_CENTRE %s", f->m_strNomProfilCent.c_str());
	fprintf(fid,"\nPROFIL_BOUT %s", f->m_strNomProfilBout.c_str());

	/**** write du tableau de forme ****/
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
/* writeFichierForm2 */
/************************/

void writeFichierForm2(char *fileName, Form *f)
{
	printf ("writeFichierForm2\n");
	int i;
	FILE *fid;
	char* texte[11]={"BORD_ATTAQUE_FORME_COURBE","BORD_DE_FUITE_FORME_COURBE","SUSPENTAGE_LIGNE_A_COURBE",
		"SUSPENTAGE_LIGNE_B_COURBE","SUSPENTAGE_LIGNE_C_COURBE","SUSPENTAGE_LIGNE_D_COURBE", "SUSPENTAGE_LIGNE_E_COURBE",
		"DIEDRE_COURBE","MORPHING_COURBE","VRILLAGE_COURBE","EPAISSEUR_RELATIVE_COURBE"};
	char* texte2[11]={"BORD_ATTAQUE_FORME_CTRL","BORD_DE_FUITE_FORME_CTRL","SUSPENTAGE_LIGNE_A_CTRL",
		"SUSPENTAGE_LIGNE_B_CTRL","SUSPENTAGE_LIGNE_C_CTRL","SUSPENTAGE_LIGNE_D_CTRL", "SUSPENTAGE_LIGNE_E_CTRL",
		"DIEDRE_CTRL","MORPHING_CTRL","VRILLAGE_CTRL","EPAISSEUR_RELATIVE_CTRL"};
	Matrix* m[10];
	Matrix* m2[10];

	/**** message ****/
	printf("\nwrite fichier de Form: '%s'",fileName);
	/**** ouverture fichier en write ****/
	if( (fid = fopen( fileName, "wt" )) == NULL )
	{
		printf( "\nErreur ouverture fichier '%s'", fileName);
		exit(0);
	}
	/**** write No de version ****/
	fprintf(fid,"\nVERSION %1.2f", NO_VERSION);
	/**** write des parametres ****/
	fprintf(fid,"\n\nNB_ALVEOLES %d", f->NbCaiss);
	fprintf(fid,"\nCOEFF_PROG_GEOM %1.5f", f->CoeffProgGeom);
	fprintf(fid,"\nEPAI_REL_CENTRE %2.3f", f->EpaiRelCent);
	/**** write des points de controle ****/
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
	/*boucle write des points de controle*/
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
	
	/**** write noms des profils ****/
	fprintf(fid,"\n\nPROFIL_CENTRE %s", f->m_strNomProfilCent.c_str());
	fprintf(fid,"\nPROFIL_BOUT %s", f->m_strNomProfilBout.c_str());
	/**** write du tableau de forme ****/
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
/* calcInfoForm */
/********************/

void calcInfoForm( Form* F, TInfoForm* info )

{
	double xNerv, xNervPrec, lNerv, lNervPrec, larg;
	int i, nNerv;
	Matrix *tabXNerv = NULL;

	/*calc position nervures a plat*/
	nNerv = F->m_nbProfils;
	tabXNerv = new Matrix(nNerv,1);

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
/* AfficheInfoForm */
/********************/

void AfficheInfoForm( TInfoForm info )
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



/**********************/
/* writeFichierFGen */
/**********************/
void writeFichierFGen(char *fileName, Form *foil)
{
	ofstream out(fileName,ios::out);
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
