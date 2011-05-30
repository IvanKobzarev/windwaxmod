#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <string>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

#include "dxf.h"
#include "profil.h"
#include "geom.h"
#include "design.h"
#include "form.h"

#define sqr(f1) ((f1)*(f1))
#define CHISLOPI	3.141592675f

#ifndef DEBUG
#define DEBUG false
#endif

void EcritLigneDXF(FILE *fid, char *nom, int en3d,
				   double x1, double y1, double z1,
				   double x2, double y2, double z2)
{
	fprintf(fid,"0\nLINE\n8\n%s\n10\n%3.3f\n20\n%3.3f\n", nom, x1*1000.0f, y1*1000.0f);
	if(en3d) fprintf(fid,"30\n%3.3f\n", z1*1000.0f);
	fprintf(fid,"11\n%3.3f\n21\n%3.3f\n", x2*1000.0f, y2*1000.0f);
	if(en3d) fprintf(fid,"31\n%3.3f\n", z2*1000.0f);
}

void EcritPolyBeginDXF(FILE *fid, char *nom)
{
	fprintf(fid,"0\nPOLYLINE\n8\n%s\n66\n1\n70\n1\n", nom);
}

void EcritPolyEndDXF(FILE *fid, char *nom)
{
	fprintf(fid,"0\nSEQEND\n", nom);
}


void EcritPolyVertexDXF(FILE *fid, char *nom, int en3d,
				   double x1, double y1, double z1)
{
	fprintf(fid,"0\nVERTEX\n8\n%s\n10\n%3.3f\n20\n%3.3f\n", nom, x1*1000.0f, y1*1000.0f);
	if(en3d) fprintf(fid,"30\n%3.3f\n", z1*1000.0f);
}


void writeFichierDXF(char *fileName, TAxe *axe)

{
	int i, j, cpt;
	FILE *fid;
	Courbe *courb;
	TMesh *mesh;
	char text[50];
	/**** message ****/
	printf("\nwrite fichier de DXF: '%s'",fileName);
	/**** ouverture fichier en write ****/
	if( (fid = fopen( fileName, "wt" )) == NULL )
	{
		printf( "\nErreur ouverture fichier '%s'", fileName);
		exit(0);
        }
	/**** write DXF *****/
	/* entete DXF minimal !!! */
        fprintf(fid,"0\nSECTION\n2\nENTITIES\n");
	/*write courbes*/
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
	/*write meshs*/
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
		printf("\nProbleme ÃƒÆ’Ã‚Â¯Ãƒâ€šÃ‚Â¿Ãƒâ€šÃ‚Â½ la fermeture du fichier");
		exit(0);
	}
}

//writeManyFichierPolyDXF(PtrfileName, AxePatronDXF, AxeMarginDXF, 1, AxeRepDXF, 0, AxeCercleDXF, Numerotation, AxePatronTextDXF, W, H);

void writeManyFichierPolyDXF(char *fileName, int np, int n, TAxe **axe, TAxe **axe2, int rep, TAxe **axeR, int vent, TAxe **axeC, int num, TAxe **axeT, double* W, double* H)
{
    	FILE *fid;
	printf("\nwriteMANY POLY fichier de DXF: '%s'",fileName);
	if( (fid = fopen( fileName, "wt" )) == NULL ) {
		printf( "\nErreur ouverture fichier '%s'", fileName);
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
            writeFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], vent, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
            if (W[i] > maxW1) maxW1 = W[i];
            if (W[i] > maxW2) maxW2 = W[i];
        }*/

        for (int i=0; i<povorot1;i++) {
            writeFichierPolyDXFDelta(fid, axe[i], axe2[i],
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
            writeFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], vent, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
            if (W[i] > maxW2) maxW2 = W[i];
        }
        for (int i=povorot1+(np&1);i<povorot2;i++) {
            writeFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], vent, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
            if (W[i] > maxW2) maxW2 = W[i];
        }
        dx=maxW1*2.0 + maxW2*4.0; dy=0.0f;
        for (int i=povorot2; i<n;i++) {
            writeFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], vent, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
        }

       fprintf(fid,"0\nENDSEC\n0\nEOF\n");
       //printf ("\n...fclose(fid)");
	if(fclose(fid))	{
		printf("\nProbleme close file");
		exit(0);
	}
}

//writeManyFichierPolyDXF2(PtrfileName, n, q, AxePD, AxeMD, 1, AxeRepD, gfd->VentilationLayout, AxeCD, 1, AxePTD, W, H, numncol);
//writeManyFichierPolyDXF(char *fileName, int np, int n, TAxe **axe, TAxe **axe2, int rep, TAxe **axeR, int vent, TAxe **axeC, int num, TAxe **axeT, double* W, double* H)
void writeLayoutToDXF(char *fileName, WindPatternsProject* gfd, Layout* layout) {
    FILE *fid;
	printf("\nwrite MANY POLY fichier de DXF: '%s'",fileName);
	if( (fid = fopen( fileName, "wt" )) == NULL ) {
		printf( "\nError write file '%s'", fileName);
		exit(0);
	}
    fprintf(fid,"0\nSECTION\n2\nENTITIES\n");
    double maxW1 = -1000000.0f, maxW2 = -1000000.0f, maxW3 = -1000000.0f;
    //panelsExt
    double dx = 0.0f, dy = 0.0f, dz = 0.0f;
    for (int i = 0; i < layout->panelsExt.size(); i++) {
        LayoutElementExport* lee = layout->panelsExt[i]->leexport;
        writeFichierPolyDXFDelta(fid, lee->AxePD, lee->AxeMD,
            1, lee->AxeRepD, gfd->VentilationLayout, lee->AxeCD, 1, lee->AxePTD, dx, dy, dz, 0);
        dy = dy + lee->H * 2.0f;
        if (lee->W > maxW1) maxW1 = lee->W;
    }

    //profs
    dx = maxW1*2.0; dy = 0.0f; maxW2 = -1000000.0f; maxW3 = -1000000.0f;
    for (int i = 0; i < layout->profs.size(); i++) {
        LayoutElementExport* lee = layout->profs[i]->leexport;
        writeFichierPolyDXFDelta(fid, lee->AxePD, lee->AxeMD,
            1, lee->AxeRepD, gfd->VentilationLayout, lee->AxeCD, 1, lee->AxePTD, dx, dy, dz, 0);
        dy = dy + lee->H * 2.0f;
        if (lee->W > maxW2) maxW2 = lee->W;
    }

    //diagNervs
    dx = maxW1 * 2.0 + maxW2 * 4.0; dy = 0.0f;
    for (int i = 0; i < layout->diagNervs.size(); i++) {
        LayoutElementExport* lee = layout->diagNervs[i]->leexport;
        writeFichierPolyDXFDelta(fid, lee->AxePD, lee->AxeMD,
            1, lee->AxeRepD, gfd->VentilationLayout, lee->AxeCD, 1, lee->AxePTD, dx, dy, dz, 0);
        dy = dy + lee->H * 2.0f;
        if (lee->W > maxW3) maxW3 = lee->W;
    }

    //panelsInt
    dx = maxW1 * 2.0 + maxW2 * 4.0 + maxW3 * 4.0; dy = 0.0f;
    for (int i = 0; i < layout->panelsInt.size(); i++) {
        LayoutElementExport* lee = layout->panelsInt[i]->leexport;
        writeFichierPolyDXFDelta(fid, lee->AxePD, lee->AxeMD,
            1, lee->AxeRepD, gfd->VentilationLayout, lee->AxeCD, 1, lee->AxePTD, dx, dy, dz, 0);
        dy = dy + lee->H * 2.0f;
    }

    fprintf(fid,"0\nENDSEC\n0\nEOF\n");
    if(fclose(fid))	{
	    printf("\nError close file");
	    exit(0);
    }
}

void writeManyFichierPolyDXF2(char *fileName, int np, int n, TAxe **axe, TAxe **axe2, int rep, TAxe **axeR, int vent, TAxe **axeC, int num, TAxe **axeT, double* W, double* H, int* numncon)
{
    FILE *fid;
	printf("\nwriteMANY POLY fichier de DXF: '%s'",fileName);
	if( (fid = fopen( fileName, "wt" )) == NULL ) {
		printf( "\nErreur ouverture fichier '%s'", fileName);
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
            writeFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], 0, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
            if (W[i] > maxW1) maxW1 = W[i];
        }
        dx=maxW1*2.0; dy=0.0f; maxW2=-1000000.0f; maxW3=-1000000.0f;
        /* +(np&1) */
        for (int i=povorot1;i<povorot2;i++) {
            writeFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], vent, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
			//printf ("\n%d dy+%f",i, H[i]);
            if (W[i] > maxW2) maxW2 = W[i];
        }
        dx=maxW1*2.0 + maxW2*4.0; dy=0.0f;
        for (int i=povorot2; i<povorot2d;i++) {
            writeFichierPolyDXFDelta(fid, axe[i], axe2[i],
                                            rep, axeR[i], vent, axeC[i], num, axeT[i], dx, dy, dz, i);
            dy=dy + H[i]*2.0f;
			//printf ("\n%d dy+%f",i, H[i]);
            if (W[i] > maxW3) maxW3 = W[i];
        }

        dx=maxW1*2.0 + maxW2*4.0 + maxW3*4.0; dy=0.0f;
        for (int i=povorot2d; i<n;i++) {
            writeFichierPolyDXFDelta(fid, axe[i], axe2[i],
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
/* writeFichierPolyDXF */
/**********************/

void writeFichierPolyDXF(char *fileName, TAxe *axe, TAxe *axe2, int rep, TAxe *axeR, int vent, TAxe *axeC, int num, TAxe *axeT )
{

	int i, j, cpt;
	double fx, fy, fz;
	FILE *fid;
	Courbe *courb;
        TMesh *mesh;
	char text[50];
	/**** message ****/
	printf("\nwrite POLY fichier de DXF: '%s'",fileName);
	/**** ouverture fichier en write ****/
	if( (fid = fopen( fileName, "wt" )) == NULL )
	{
		printf( "\nErreur ouverture fichier '%s'", fileName);
		exit(0);
	}
	/**** write DXF *****/
	/* entete DXF minimal !!! */
        fprintf(fid,"0\nSECTION\n2\nENTITIES\n");
	double curx, cury, x1, y1;
	/*write courbes*/
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
//			printf ("\n while in writeFichierPolyDXF");
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

	/*write meshs*/
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
		printf("\nProbleme ÃƒÆ’Ã‚Â¯Ãƒâ€šÃ‚Â¿Ãƒâ€šÃ‚Â½ la fermeture du fichier");
		exit(0);
	}
}


/**********************/
/* writeFichierPolyDXFDelta */
/**********************/

void writeFichierPolyDXFDelta(FILE *fid, TAxe *axe, TAxe *axe2, int rep, TAxe *axeR, int vent, TAxe *axeC, int num, TAxe *axeT, double dx, double dy, double dz, int n )
{
    //printf ("\n writeFichPolyDxfDelta");
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
	/*write meshs*/
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
