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
TAxe *AxeSel;
Courbe * CourbProfile;
Forme *F;
char ProjectName[255];
int FenetreProfiles;
/*divers tableaux*/
Matrice *IntProfCent, *ExtProfCent, *IntProfBout, *ExtProfBout;
Matrice** ReperPoints;
char NomFichierForme[255];
char NomFichierProject[255];
char NomFichierRepPoints[255];
char NomFichierVentHoles[255];
char NomFichierDiagNerv[255];

int XYGrid = 0;

WindPatternsProject* gfd = new WindPatternsProject();
GLUI_Spinner *SpinNoNerv[2];
int NoNerv[2] = {2, 3};
float dyw=2.0f, wN=25.0f, kChord=1.05f, kMf=1.0f;
int xSouris, ySouris;
int zoomIN = false, debZoomIN = false, finZoomIN = false, zoomOUT = false, quitZoom = false;

Courbe *CourbZoom;

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
    int Fenetre;
    xSouris = x;
    ySouris = y;
    AxeSel = NULL;

    Fenetre = glutGetWindow();

	if (Fenetre == FenetreProfiles) {
        AxeSel = AxeProfiles;
    }

    /*test debut/fin, rotation/translation/zoom pour Axe3d*/
    if (AxeSel == AxeProfiles) {
        if (state == GLUT_DOWN) {
            switch (button) {
                case GLUT_LEFT_BUTTON: //zoomIN
                    if (zoomOUT) {
                        quitZoom = true;
                        break;
                    }
                    //test si debut de zoom
                    if (!zoomIN) {
                        zoomIN = true, debZoomIN = true;
                    }
                    CourbZoom->segments = ON;
                    AxeProfiles->XAuto = OFF;
                    AxeProfiles->YAuto = OFF;
                    AxeProfiles->ZAuto = OFF;
                    break;
                case GLUT_RIGHT_BUTTON: //zoomOUT x2

                    if (zoomIN) {
                        quitZoom = true;
                        break;
                    }

                    zoomOUT = true;

                    AxeProfiles->XAuto = OFF;
                    AxeProfiles->YAuto = OFF;
                    AxeProfiles->ZAuto = OFF;

                    break;

                case GLUT_MIDDLE_BUTTON:

                    quitZoom = true;

                    break;

                default: printf("\n type de bouton souris inconnu !!!");

            }

        } else /*state = GLUT_UP*/ {
            if (zoomIN) {
                zoomIN = false;
                finZoomIN = true;
                CourbZoom->segments = OFF;
                motion(x, y); //appel a motion pour maj axe

            } else if (zoomOUT) {
                motion(x, y); //appel a motion pour maj axe
            } else if (quitZoom) {
                quitZoom = false;
                zoomIN = false;
                debZoomIN = false;
                finZoomIN = false;
                zoomOUT = false;

                AxeProfiles->XAuto = ON;
                AxeProfiles->YAuto = ON;
                AxeProfiles->ZAuto = ON;
            }

            glutPostRedisplay();

        }

    }

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
	gfd->ExtProfCent = ExtProfCent;
	gfd->ExtProfBout = ExtProfBout;
	gfd->IntProfCent = IntProfCent;
	gfd->IntProfBout = IntProfBout;
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



void Apply(int /*control*/) {
    AxeProfiles->XAuto = ON;
    AxeProfiles->YAuto = ON;
    AxeProfiles->ZAuto = ON;
    AxeProfiles->xmin = -2.0;
    AxeProfiles->xmax = 2.0;
    AxeProfiles->ymin = -2.0;
    AxeProfiles->ymax = 2.0;

	AxeProfiles->XGrid = XYGrid;
	AxeProfiles->YGrid = XYGrid;
	AxeProfiles->ZGrid = OFF;

    LibererCourbesAxe(AxeProfiles);
	
	ProfilGeom* pg1 = getProfile(gfd, F, NoNerv[0]);
	pg1->print();
	Courbe* cext1, *cint1;
	getCourbeFromProfilGeom(pg1, &cext1, &cint1);
	AjoutCourbe(AxeProfiles, cext1);
	AjoutCourbe(AxeProfiles, cint1);

/*	ProfilGeom* pg2 = getProfile(gfd, F, NoNerv[1]);
	Courbe* cext2, *cint2;
	getCourbeFromProfilGeom(pg2, &cext2, &cint2);
	AjoutCourbe(AxeProfiles, cext2);
	AjoutCourbe(AxeProfiles, cint2);*/

	ProfilGeom* pg1b = getBalloneProfilGeom(pg1, kChord, kMf, F->m_pProfils[NoNerv[0]]->m_fWidth, wN, dyw);
	Courbe* cext3, *cint3;
	getCourbeFromProfilGeom(pg1b, &cext3, &cint3);
	cext3->CouleurSegments[0]=1.0f;
	cext3->CouleurSegments[1]=0.0f;
	cext3->CouleurSegments[2]=0.0f;
	cint3->CouleurSegments[0]=1.0f;
	cint3->CouleurSegments[1]=0.0f;
	cint3->CouleurSegments[2]=0.0f;

	AjoutCourbe(AxeProfiles, cext3);
	AjoutCourbe(AxeProfiles, cint3);


	display();
}

void ModifXYGrid(int /*control*/) {
	printf ("\n ModifXYGrid()");
	AxeProfiles->XGrid = XYGrid;
	AxeProfiles->YGrid = XYGrid;

	display();
}

int main(int argc, char** argv)
{
	printf ("Profile ballone...");

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(100, 100);

	FenetreProfiles = glutCreateWindow("Profile view");
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

	GLUI_Panel *panel1 = glui->add_panel("");

	SpinNoNerv[0] = glui->add_spinner_to_panel(panel1, "No Nerv 1", GLUI_SPINNER_INT, &(NoNerv[0]));
    //SpinNoNerv[0] -> set_int_limits(-1, F->m_nbProfils - 1);

	SpinNoNerv[1] = glui->add_spinner_to_panel(panel1, "No Nerv 2", GLUI_SPINNER_INT, &(NoNerv[1]));
    //SpinNoNerv[1] -> set_int_limits(-1, F->m_nbProfils - 1);


	glui->add_checkbox("XY grid", &XYGrid, 0,  &ModifXYGrid);

    GLUI_Button *btnApply =
            glui->add_button("Apply", 0, &Apply);
    btnApply->set_w(10);

    glui->add_column(true);

	strcpy(NomFichierProject, "f17project.wpp");

	gfd = LectureWindPatternsProject(NomFichierProject);
    LoadFromWindPatternsProject(gfd);

    FicProject = glui->add_statictext("???");
    FicProject->set_text(NomFichierProject);

	GLUI_Panel *panel2 = glui->add_panel("");
	
    GLUI_Spinner *SpinKchord = glui->add_spinner_to_panel(
            panel2, "K chord", GLUI_SPINNER_FLOAT, &(kChord));
    SpinKchord -> set_float_limits(0.8, 1.5);

    GLUI_Spinner *SpinKmf = glui->add_spinner_to_panel(
            panel2, "K move forward", GLUI_SPINNER_FLOAT, &(kMf));
    SpinKmf -> set_float_limits(0.0, 1.0);

    GLUI_Spinner *SpinWn = glui->add_spinner_to_panel(
            panel2, "W new", GLUI_SPINNER_FLOAT, &(wN));
    SpinWn -> set_float_limits(0.0, 100.0);

    GLUI_Spinner *SpinDyw = glui->add_spinner_to_panel(
            panel2, "D yw", GLUI_SPINNER_FLOAT, &(dyw));
    SpinDyw -> set_float_limits(0.0, 100.0);

    //GLUI_Panel *panel = glui->add_panel("");

    GLUI_Master.set_glutIdleFunc(NULL);
	Apply(0);
    glui->sync_live();

    display();
    glutMainLoop();
	return 0;
}

