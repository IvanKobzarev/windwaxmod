// ProfileBallone.cpp : Defines the entry point for the console application.
//
#pragma warning(disable:4786)
#pragma warning(disable:4514)

//#include "stdafx.h"
#include "afx.h"		//class CString
#include "afxdlgs.h"	//class CFileDialog

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "../commun/plot.h"
#include "../commun/matrice.h"
#include "../commun/geom.h"
#include "../commun/fichier.h"
#include "../commun/profil.h"
#include "../commun/pince.h"
#include "../commun/rasklad.h"
#include "../commun/logger.h"

#include "../GL/glui.h"
#include "../GL/glut.h"

GLUI *glui;
GLUI_StaticText *FicProject;
TAxe *AxeProfiles;

Forme *F;
char ProjectName[255];

/*divers tableaux*/
Matrice *IntProfCent, *ExtProfCent, *IntProfBout, *ExtProfBout;
Matrice** ReperPoints;
char NomFichierForme[255];
char NomFichierProject[255];
char NomFichierRepPoints[255];
char NomFichierVentHoles[255];
char NomFichierDiagNerv[255];

WindPatternsProject* gfd = new WindPatternsProject();

void display(void) {
    //CalculVue3dEtPatron();
    VisuAxe(AxeProfiles);
    glutSwapBuffers();
}


void motion(int x, int y) {
}

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    display();
}

void keyboard(unsigned char key, int x, int y) {
}

void BoutonSouris(int button, int state, int x, int y) {
}

void InitFenetre(void) {
    glEnable(GL_LINE_STIPPLE);
    glutDisplayFunc(&display);
    glutReshapeFunc(&reshape);
    glutKeyboardFunc(&keyboard);
    glutMotionFunc(&motion);
    glutMouseFunc(&BoutonSouris);
}

void LoadFromWindPatternsProject(WindPatternsProject* gfd) {
	Forme* tmpF=0;
	Matrice *tmpIntProfCent=0, *tmpExtProfCent=0, *tmpIntProfBout=0, *tmpExtProfBout=0;
	Matrice** tmpReperPoints=0;
	int tmpquantDiag=0;
	int tmpquantVH=0;
	int *tmpnoNervD=0;
	int *tmpnoNervVH=0;
	int tmpVentCentralNerv=0;

	try {
		tmpF = LectureFichierForme(gfd->NomFichierForme);
		try {
		tmpF->Validate();
		} catch (char* msg) {
			AfxMessageBox(msg);
			return;
		}

		LectureFichierProfil(tmpF->m_strNomProfilCent.c_str(), &tmpExtProfCent, &tmpIntProfCent);
		LectureFichierProfil(tmpF->m_strNomProfilBout.c_str(), &tmpExtProfBout, &tmpIntProfBout);
	} catch (char* sexception) {
		printf ("\n Problem loading project: %s", sexception);
		char msg[100];
		sprintf (msg, "Problem loading project: %s, check project file %s, forme file %s", sexception, gfd-> name, gfd->NomFichierForme);
		AfxMessageBox(msg);
		return;
	}

	/* -------- seems that loading file success ------------ */
	F = tmpF;
	ExtProfCent = tmpExtProfCent;
	IntProfCent = tmpIntProfCent;
	ExtProfBout = tmpExtProfBout;
	IntProfBout = tmpIntProfBout;
	/* --------                                 ------------ */
}

void ChargerFichierProject(int /*control*/) {
    CString NomFichier;
    char* PtrNomFichier;
    WindPatternsProject* oldWpp;
    CFileDialog DlgOpen(TRUE, NULL, "*.wpp", OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        NomFichier = DlgOpen.GetPathName();
        PtrNomFichier = NomFichier.GetBuffer(1);
        strcpy(NomFichierProject, PtrNomFichier);
        oldWpp = gfd;
        gfd = LectureWindPatternsProject(NomFichierProject);
        LoadFromWindPatternsProject(gfd);
        delete(oldWpp);
        FicProject->set_text(NomFichierProject);
    }
    //Appliquer(0);
    display();
    glui->sync_live();
}

int main(int argc, char** argv)
{
	printf ("\nProfile ballone");

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(100, 100);

	int FenetreProfiles = glutCreateWindow("Profile view");
    InitFenetre(); 

    AxeProfiles = CreerAxe(FenetreProfiles);

	int ws, hs;
    ws = glutGet(GLUT_SCREEN_WIDTH);
    hs = glutGet(GLUT_SCREEN_HEIGHT);

    glutSetWindow(FenetreProfiles);
    glutPositionWindow((int) (0.0 * (double) ws), (int) (0.50 * (double) hs));
    glutReshapeWindow((int) (1.0 * (double) ws), (int) (0.45 * (double) hs));

    glui = GLUI_Master.create_glui("Profile ballone", 0, 0, 0);

    GLUI_Button *btnLoadForme =
            glui->add_button("Load project", 0, &ChargerFichierProject);
    btnLoadForme->set_w(10);

    glui->add_column(true);

	strcpy(NomFichierProject, "f17project.wpp");

    FicProject = glui->add_statictext("???");
    FicProject->set_text(NomFichierProject);

    GLUI_Panel *panel = glui->add_panel("");

    GLUI_Master.set_glutIdleFunc(NULL);

    glui->sync_live();

    display();
    glutMainLoop();
	return 0;
}

