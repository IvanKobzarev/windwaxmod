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
#include "../commun/design.h"

#include "GL/glui.h"
#include "GL/glut.h"

GLUI *glui;
// static text on dialog form, that shows name of loaded project
GLUI_StaticText *FicProject;
// static text on dialog form, that shows name of loaded design 
GLUI_StaticText *FicDesign;
TAxe *AxeProjections;
TAxe *AxeSel;
TAxe *Axe3d;

int VisuFace = 0; 
int VisuSymetrique = 0;

Courbe * CourbProfile;
Forme *F;
char ProjectName[255];
int Fenetre2d, Fenetre3d;
/*divers tableaux*/
Matrice *IntProfCent, *ExtProfCent, *IntProfBout, *ExtProfBout;
Matrice** ReperPoints;
char NomFichierForme[255];

char NomFichierProject[255];
char NomFichierDesign[255];

char NomFichierRepPoints[255];
char NomFichierVentHoles[255];
char NomFichierDiagNerv[255];
int ProjOrthoPers = 0;
int XYGrid = 0;

int ExtPrjNoseUp = 0;
int IntPrjNoseUp = 1;

WindPatternsProject* gfd = new WindPatternsProject();
GLUI_Spinner *SpinNoNerv[2];
int NoNerv[2] = {0, 1};
float dyw = 0.5f, wN=22.0f, kChord=1.05f, kMf=0.25f, power=1.1f;
int xSouris, ySouris;
int zoomIN = false, debZoomIN = false, finZoomIN = false, zoomOUT = false, quitZoom = false;
int Rotation3d = OFF;
int Translation3d = OFF;
int ZoomTwist3d = OFF;

Courbe *CourbZoom;

void Apply(int /*control*/);

void ModifVisu3d(int /*control*/) {
    VisuAxe(Axe3d);
    glutSwapBuffers();
}

void ModifVisuSymetrique(int /*control*/) {
	printf ("\n ModifVisySymetrique()");
	Apply(0);
    VisuAxe(Axe3d);
    glutSwapBuffers();
}

void display(void) {
    VisuAxe(Axe3d);
    glutSwapBuffers();
    //CalculVue3dEtPatron();
    VisuAxe(AxeProjections);
    glutSwapBuffers();
}

/*****************/
/* reshape       */
/*****************/

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    display();
}

/*****************/
/* motion        */
/*****************/

void motion(int x, int y) {
    //pour gestion zoom dans fenetre patron
    double xs, ys;
    double xmin, xmax, ymin, ymax;
    int h, w; /*taille de la fenetre*/
    double dx, dy; //pour zoom OUT

    /*message pour test*/
    //printf("\n<%d,%d>",x,y);

    /*axe 3d -> rotation*/

    if (Rotation3d) {
        Axe3d->azimuth = Axe3d->azimuth + (x - xSouris);
        Axe3d->incidence = Axe3d->incidence + (y - ySouris);
        //Axe3dBal->azimuth = Axe3dBal->azimuth + (x - xSouris);
        //Axe3dBal->incidence = Axe3dBal->incidence + (y - ySouris);

        xSouris = x;
        ySouris = y;
    }/*axe 3d -> translation*/
    else if (Translation3d) {
        Axe3d->xcam = Axe3d->xcam + (double) (x - xSouris) / 100.0;
        Axe3d->ycam = Axe3d->ycam - (double) (y - ySouris) / 100.0;
        //Axe3dBal->xcam = Axe3dBal->xcam + (double) (x - xSouris) / 100.0;
        //Axe3dBal->ycam = Axe3dBal->ycam - (double) (y - ySouris) / 100.0;

        xSouris = x;
        ySouris = y;

    }/*axe 3d -> zoom*/
    else if (ZoomTwist3d) {
        Axe3d->twist = Axe3d->twist - (double) (x - xSouris);
        Axe3d->zcam = Axe3d->zcam + (double) (y - ySouris) / 100.0;
        //Axe3dBal->twist = Axe3dBal->twist - (double) (x - xSouris);
        //Axe3dBal->zcam = Axe3dBal->zcam + (double) (y - ySouris) / 100.0;

        xSouris = x;
        ySouris = y;
    }/*axe patron  -> zoom IN*/
    else if (zoomIN) {
        w = glutGet(GLUT_WINDOW_WIDTH);
        h = glutGet(GLUT_WINDOW_HEIGHT);

        xmin = AxeSel->xmin;
        xmax = AxeSel->xmax;

        ymin = AxeSel->ymin;
        ymax = AxeSel->ymax;

        xs = xmin + (double) x * (xmax - xmin) / (double) w;
        ys = ymin + ((double) h - (double) y)*(ymax - ymin) / h;
        
        if (debZoomIN) {
            debZoomIN = false;

            CourbZoom->pts->SetElement(0, 0, xs);
            CourbZoom->pts->SetElement(0, 1, ys);

            CourbZoom->pts->SetElement(4, 0, xs);
            CourbZoom->pts->SetElement(4, 1, ys);

            CourbZoom->pts->SetElement(3, 0, xs);
            CourbZoom->pts->SetElement(1, 1, ys);

        }

        CourbZoom->pts->SetElement(1, 0, xs);
        CourbZoom->pts->SetElement(3, 1, ys);

        CourbZoom->pts->SetElement(2, 0, xs);
        CourbZoom->pts->SetElement(2, 1, ys);
    }

    /*else if (finZoomIN) {
        finZoomIN = false;

        //test points non confondus

        if ((CourbZoom->pts->Element(0, 0) != CourbZoom->pts->Element(2, 0))
                && (CourbZoom->pts->Element(0, 1) != CourbZoom->pts->Element(2, 1))) {
            //assure xmin<xmax
            if (CourbZoom->pts->Element(0, 0) < CourbZoom->pts->Element(2, 0)) {
                AxePatron->xmin = CourbZoom->pts->Element(0, 0);
                AxePatron->xmax = CourbZoom->pts->Element(2, 0);
            } else {
                AxePatron->xmax = CourbZoom->pts->Element(0, 0);
                AxePatron->xmin = CourbZoom->pts->Element(2, 0);
            }
            //assure ymin<ymax
            if (CourbZoom->pts->Element(0, 1) < CourbZoom->pts->Element(2, 1)) {
                AxePatron->ymin = CourbZoom->pts->Element(0, 1);
                AxePatron->ymax = CourbZoom->pts->Element(2, 1);
            } else {
                AxePatron->ymax = CourbZoom->pts->Element(0, 1);
                AxePatron->ymin = CourbZoom->pts->Element(2, 1);
            }
        }
    } else if (zoomOUT) {
        w = glutGet(GLUT_WINDOW_WIDTH);
        h = glutGet(GLUT_WINDOW_HEIGHT);

        xmin = AxeSel->xmin;
        xmax = AxeSel->xmax;

        ymin = AxeSel->ymin;
        ymax = AxeSel->ymax;

        xs = xmin + (double) x * (xmax - xmin) / (double) w;
        ys = ymin + ((double) h - (double) y)*(ymax - ymin) / h;

        zoomOUT = false;

        dx = AxePatron->xmax - AxePatron->xmin;
        dy = AxePatron->ymax - AxePatron->ymin;

        AxePatron->xmin = xs - dx;
        AxePatron->xmax = xs + dx;

        AxePatron->ymin = ys - dy;
        AxePatron->ymax = ys + dy;
    } */
    glutPostRedisplay();

}

void keyboard(unsigned char key, int x, int y) {
}

void BoutonSouris(int button, int state, int x, int y) {
    int Fenetre;
    xSouris = x;
    ySouris = y;
    AxeSel = NULL;

    Fenetre = glutGetWindow();
	if (Fenetre == Fenetre2d) {
        AxeSel = AxeProjections;
    }
    if (Fenetre == Fenetre3d) {
        AxeSel = Axe3d;
    }

    /*test debut/fin, rotation/translation/zoom pour Axe3d*/
    if (AxeSel == Axe3d) {
        if (state == GLUT_DOWN) {
            switch (button) {
                case GLUT_LEFT_BUTTON:
                    if (Translation3d == ON) {
                        Translation3d = OFF;
                        ZoomTwist3d = ON;
                    } else Rotation3d = ON;
                    break;

                case GLUT_RIGHT_BUTTON:

                    if (Rotation3d == ON) {
                        Rotation3d = OFF;
                        ZoomTwist3d = ON;
                    } else Translation3d = ON;
                    break;

                case GLUT_MIDDLE_BUTTON:
                    ZoomTwist3d = ON;
                    break;

                default: printf("\n type de bouton souris inconnu !!!");
            }
        } else /*state = GLUT_UP*/ {
            finZoomIN = true;
            Rotation3d = OFF;
            Translation3d = OFF;
            ZoomTwist3d = OFF;
        }
    }


    /*test debut/fin, rotation/translation/zoom pour Axe3d*/
    if (AxeSel == AxeProjections) {
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
                    AxeProjections->XAuto = OFF;
                    AxeProjections->YAuto = OFF;
                    AxeProjections->ZAuto = OFF;
                    break;
                case GLUT_RIGHT_BUTTON: //zoomOUT x2

                    if (zoomIN) {
                        quitZoom = true;
                        break;
                    }

                    zoomOUT = true;

                    AxeProjections->XAuto = OFF;
                    AxeProjections->YAuto = OFF;
                    AxeProjections->ZAuto = OFF;

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

                AxeProjections->XAuto = ON;
                AxeProjections->YAuto = ON;
                AxeProjections->ZAuto = ON;
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

void InitLumiere(void) {
    float mat_specular[] = {1.0, 1.0, 1.0, 1.0};
    float mat_shininess[] = {50.0};
    float light_position[] = {0.0, 0.0, 1.0, 0.0};
    float light_color[] = {1.0, 1.0, 0.0, 1.0};
    glShadeModel(GL_SMOOTH);
    /*glClearColor(0.0, 0.0, 0.0, 0.0);*/
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_color);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_color);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
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

//TODO: design load calling
void ChargerFichierDesign(int /*control*/) {
    CString NomFichier;
    char* PtrNomFichier;
    //WindDesign* oldWpp;
    CFileDialog DlgOpen(TRUE, NULL, "*.wdn", OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        NomFichier = DlgOpen.GetPathName();
        PtrNomFichier = NomFichier.GetBuffer(1);
        strcpy(NomFichierDesign, PtrNomFichier);
        //oldWpp = gfd;
        //gfd = LectureWindPatternsProject(NomFichierProject);
        //LoadFromWindPatternsProject(gfd);
        //delete(oldWpp);
        FicDesign->set_text(NomFichierDesign);
    }
    //Appliquer(0);
    display();
    glui->sync_live();
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

void Apply2d(int /*control*/) {
}

void Apply3d(int /*control*/) {
}


void Apply(int /*control*/) {
	// 3D
	// forme visualisation
    Matrice *XExt, *YExt, *ZExt, *YExt0;
    Matrice *XInt, *YInt, *ZInt, *YInt0;

    CalculForme3D(F, 0, 0.0f,
            ExtProfCent, IntProfCent, ExtProfBout, IntProfBout,
            &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);

	printf ("\n WindDesigner.test()");
	Forme3D* f3d = getForme3D(F, 0, 0.0f,
            ExtProfCent, IntProfCent, ExtProfBout, IntProfBout);

    //AjoutForme3D(Axe3d, XExt, YExt, ZExt, XInt, YInt, ZInt, VisuFace, VisuSymetrique);
	//YExt0 = Zeros(f3d->YExt->GetLignes(), f3d->YExt->GetColonnes());
	//YInt0 = Zeros(f3d->YInt->GetLignes(), f3d->YInt->GetColonnes());

	AjoutForme3D(Axe3d, f3d->XExt, f3d->YExt, f3d->ZExt, f3d->XInt, f3d->YInt, f3d->ZInt, VisuFace, VisuSymetrique);




	// 2D
	//------------------------ Projections Fenetre calculation -------------------------------------
    AxeProjections->XAuto = ON;
    AxeProjections->YAuto = ON;
    AxeProjections->ZAuto = ON;
    AxeProjections->xmin = -100.0;
    AxeProjections->xmax = 100.0;
    AxeProjections->ymin = -100.0;
    AxeProjections->ymax = 100.0;
	AxeProjections->XGrid = XYGrid;
	AxeProjections->YGrid = XYGrid;
	AxeProjections->ZGrid = OFF;
    LibererCourbesAxe(AxeProjections);
	
	//AjoutCourbe(AxeProjections, cint4);
	printf ("\n FormeProjection* fp = getFormeProjection(f3d);");
	FormeProjection* fp = getFormeProjection(f3d);
	printf ("\n ajoutFormeProjectionCourbesToAxe(fp, AxeProjections);");
	ajoutFormeProjectionCourbesToAxe(AxeProjections, fp, VisuSymetrique, 0.0, IntPrjNoseUp);
	ajoutFormeProjectionCourbesToAxe(AxeProjections, fp, VisuSymetrique, 2.0, ExtPrjNoseUp);
	//ajoutFormeProjectionCourbesToAxe(FormeProjection* fp, TAxe* axe)
	printf ("\n display()");

	// Zoom
    CourbZoom = new Courbe("Zoom");
    CourbZoom->pts = Zeros(5, 2);
    CourbZoom->points = OFF;
    CourbZoom->segments = OFF;
    CourbZoom->symX = OFF;
    CourbZoom->CouleurSegments[0] = 0.0f;
    CourbZoom->CouleurSegments[2] = 0.0f;
    AjoutCourbe(AxeProjections, CourbZoom);
	//-------------------------
	display();
}

void ModifXYGrid(int /*control*/) {
	printf ("\n ModifXYGrid()");
	AxeProjections->XGrid = XYGrid;
	AxeProjections->YGrid = XYGrid;
	display();
}

void ModifExtPrjNoseUp(int /*control*/) {
	Apply(0);
	display();
}

void ModifIntPrjNoseUp(int /*control*/) {
	Apply(0);
	display();
}

void ModifProjection3d(int /*control*/) {
    /*test type de proj.*/
	if (ProjOrthoPers == 0) {
        Axe3d->proj = PROJ_ORTHOGONALE;
        //Axe3dBal->proj = PROJ_ORTHOGONALE;
	}
	else {
        Axe3d->proj = PROJ_PERSPECTIVE;
        //Axe3dBal->proj = PROJ_PERSPECTIVE;
	}
    /*retrace uniquement la vue 3d*/
    VisuAxe(Axe3d);
    //VisuAxe(Axe3dBal);
    glutSwapBuffers();
}


int main(int argc, char** argv)
{
	printf ("\nWind designer");
	// glut initialisation
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(100, 100);

	// 2d view for projections
	Fenetre2d = glutCreateWindow("2d design");
    InitFenetre(); 
    AxeProjections = CreerAxe(Fenetre2d);
	int ws, hs;
    ws = glutGet(GLUT_SCREEN_WIDTH);
    hs = glutGet(GLUT_SCREEN_HEIGHT);
    glutSetWindow(Fenetre2d);
    glutPositionWindow((int) (0.0 * (double) ws), (int) (0.50 * (double) hs));
    glutReshapeWindow((int) (1.0 * (double) ws), (int) (0.45 * (double) hs));

	Fenetre3d = glutCreateWindow("3d design");
    InitFenetre();
    InitLumiere();
    Axe3d = CreerAxe(Fenetre3d);
    Axe3d->axe3d = ON;
    glutSetWindow(Fenetre3d);
    glutPositionWindow((int) (0.25 * (double) ws), (int) (0.0 * (double) hs));
    glutReshapeWindow((int) (0.75 * (double) ws), (int) (0.48 * (double) hs));

	// main form whith GUI
    glui = GLUI_Master.create_glui("Wind designer", 0, 0, 0);
	GLUI_Panel *panelProject = glui->add_panel("");

    GLUI_Button *btnLoadForme = glui->add_button_to_panel(panelProject, "Load project", 0, &ChargerFichierProject);
    btnLoadForme->set_w(10);

    glui->add_column_to_panel(panelProject, false);
	// loading Default project
	strcpy(NomFichierProject, "f17project.wpp");
	gfd = LectureWindPatternsProject(NomFichierProject);
    LoadFromWindPatternsProject(gfd);
    FicProject = glui->add_statictext_to_panel(panelProject, "???");
    FicProject->set_text(NomFichierProject);

	GLUI_Panel *panelDesign = glui->add_panel("");
    GLUI_Button *btnLoadDesign = glui->add_button_to_panel(panelDesign, "Load design", 0, &ChargerFichierDesign);
    btnLoadDesign->set_w(10);
    glui->add_column_to_panel(panelDesign, false);
	// load default design
	strcpy(NomFichierDesign, "f17design.wdn");
    FicDesign = glui->add_statictext_to_panel(panelDesign, "???");
    FicDesign->set_text(NomFichierDesign);
	//FicDesign->set_h(10);

	GLUI_Panel *panel1 = glui->add_panel("");
	GLUI_Button *btn2d = glui->add_button_to_panel(panel1, "2d", 0, &Apply2d);
	btn2d->set_w(5);
	glui->add_column_to_panel(panel1, false);
	GLUI_Button *btn3d = glui->add_button_to_panel(panel1, "3d", 0, &Apply3d);
	btn3d->set_w(5);

	//SpinNoNerv[0] = glui->add_spinner_to_panel(panel1, "No Nerv 1", GLUI_SPINNER_INT, &(NoNerv[0]));
    ////SpinNoNerv[0] -> set_int_limits(-1, F->m_nbProfils - 1);

	//SpinNoNerv[1] = glui->add_spinner_to_panel(panel1, "No Nerv 2", GLUI_SPINNER_INT, &(NoNerv[1]));
    ////SpinNoNerv[1] -> set_int_limits(-1, F->m_nbProfils - 1);
	GLUI_Panel *panel2dOptions = glui->add_panel("2D");
	glui->add_checkbox_to_panel(panel2dOptions, "XY grid", &XYGrid, 0,  &ModifXYGrid);

	glui->add_checkbox_to_panel(panel2dOptions, "Ext Projection nose up", &ExtPrjNoseUp, 0,  &ModifExtPrjNoseUp);
	glui->add_checkbox_to_panel(panel2dOptions, "Int Projection nose up", &IntPrjNoseUp, 0,  &ModifIntPrjNoseUp);

    //GLUI_Button *btnApply = glui->add_button("Apply", 0, &Apply);
    //btnApply->set_w(10);

	/*GLUI_Panel *panel2 = glui->add_panel("");
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

    GLUI_Spinner *SpinPower = glui->add_spinner_to_panel(
            panel2, "power", GLUI_SPINNER_FLOAT, &(power));
    SpinPower -> set_float_limits(0.05, 3.0);*/

    GLUI_Panel *panelViewOptions = glui->add_panel("");
    glui->add_checkbox_to_panel(panelViewOptions, "symetrique", &VisuSymetrique, 0, &ModifVisuSymetrique);
    GLUI_Panel *panel_proj = glui->add_panel_to_panel(panelViewOptions, "Projection");
    GLUI_RadioGroup *radio_proj =
            glui->add_radiogroup_to_panel(panel_proj, &ProjOrthoPers, 0, &ModifProjection3d);
    glui->add_radiobutton_to_group(radio_proj, "Orthogonale");
    glui->add_radiobutton_to_group(radio_proj, "Perspective");

    GLUI_Master.set_glutIdleFunc(NULL);
	glui->sync_live();
	Apply(0);
	AxeSel = Axe3d;
    display();
    glutMainLoop();
	return 0;
}

