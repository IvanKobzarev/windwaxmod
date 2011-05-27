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
#include "../commun/layout.h"
#include "../commun/logger.h"
#include "../commun/design.h"

#include "GL/glui.h"
#include "GL/glut.h"

GLUI *glui;
// static text on dialog form, that shows name of loaded project
GLUI_StaticText *FicProject;
// static text on dialog form, that shows name of loaded design 
GLUI_StaticText *FicDesignExt;
GLUI_StaticText *FicDesignInt;
TAxe *AxeProjections;
TAxe *AxeSel;
TAxe *Axe3d;

int ViewFace = 0; 
int ViewSymetrique = 0;

Courbe * CourbProfile;
Form *F;
char ProjectName[255];
int window2d, window3d;
/*divers tableaux*/
Matrix *IntProfCent, *ExtProfCent, *IntProfBout, *ExtProfBout;
Matrix** ReperPoints;
char fileNameForm[255];

char fileNameProject[255];

char fileNameDesignExt[255];
char fileNameDesignInt[255];

char fileNameRepPoints[255];
char fileNameVentHoles[255];
char fileNameDiagNerv[255];
int ProjOrthoPers = 0;
int XYGrid = 0;

int ExtPrjNoseUp = 0;
int IntPrjNoseUp = 1;

WindPatternsProject* gfd = new WindPatternsProject();
KiteDesign* kiteDesignExt = 0;
KiteDesign* kiteDesignInt = 0;
GLUI_Spinner *SpinNoNerv[2];
int NoNerv[2] = {0, 1};
float dyw = 0.5f, wN=22.0f, kChord=1.05f, kMf=0.25f, power=1.1f;

float opacExt = 0.8f, opacInt = 0.6f;

int xSouris, ySouris;
int zoomIN = false, debZoomIN = false, finZoomIN = false, zoomOUT = false, quitZoom = false;
int Rotation3d = OFF;
int Translation3d = OFF;
int ZoomTwist3d = OFF;

Courbe *CourbZoom;

void apply(int /*control*/);

void modifView3d(int /*control*/) {
    ViewAxe(Axe3d);
    glutSwapBuffers();
}

void modifViewSymetrique(int /*control*/) {
	//printf ("\n modifVisySymetrique()");
	apply(0);
    ViewAxe(Axe3d);
    glutSwapBuffers();
}

void display(void) {
    ViewAxe(Axe3d);
    glutSwapBuffers();
    //calcVue3dEtPatron();
    ViewAxe(AxeProjections);
    glutSwapBuffers();
}

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    display();
}

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
    int window;
    xSouris = x;
    ySouris = y;
    AxeSel = NULL;

    window = glutGetWindow();
	if (window == window2d) {
        AxeSel = AxeProjections;
    }
    if (window == window3d) {
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

void InitWindow(void) {
    glEnable(GL_LINE_STIPPLE);
    glutDisplayFunc(&display);
    glutReshapeFunc(&reshape);
    glutKeyboardFunc(&keyboard);
    glutMotionFunc(&motion);
    glutMouseFunc(&BoutonSouris);
}

void InitLight(void) {
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
	glEnable(GL_COLOR_MATERIAL);
}


void LoadFromWindPatternsProject(WindPatternsProject* gfd) {
	Form* tmpF=0;
	Matrix *tmpIntProfCent=0, *tmpExtProfCent=0, *tmpIntProfBout=0, *tmpExtProfBout=0;
	Matrix** tmpReperPoints=0;
	int tmpquantDiag=0;
	int tmpquantVH=0;
	int *tmpnoNervD=0;
	int *tmpnoNervVH=0;
	int tmpVentCentralNerv=0;

	try {
		tmpF = readFichierForm(gfd->fileNameForm);
		try {
		tmpF->Validate();
		} catch (char* msg) {
			AfxMessageBox(msg);
			return;
		}

		readFichierProfil(tmpF->m_strNomProfilCent.c_str(), &tmpExtProfCent, &tmpIntProfCent);
		readFichierProfil(tmpF->m_strNomProfilBout.c_str(), &tmpExtProfBout, &tmpIntProfBout);


	} catch (char* sexception) {
		printf ("\n Problem loading project: %s", sexception);
		char msg[100];
		sprintf (msg, "Problem loading project: %s, check project file %s, forme file %s", sexception, gfd-> name, gfd->fileNameForm);
		AfxMessageBox(msg);
		return;
	}

	/* -------- seems that loading file success ------------ */
	F = tmpF;
	ExtProfCent = tmpExtProfCent;
	IntProfCent = tmpIntProfCent;
	ExtProfBout = tmpExtProfBout;
	IntProfBout = tmpIntProfBout;

	F->ExtProfCent = ExtProfCent;
	F->IntProfCent = IntProfCent;
	F->ExtProfBout = ExtProfBout;
	F->IntProfBout = IntProfBout;

	F->EpaiRelProfCent = EpaisseurRelative(ExtProfCent, IntProfCent);
	F->EpaiRelProfBout = EpaisseurRelative(ExtProfBout, IntProfBout);


	gfd->ExtProfCent = ExtProfCent;
	gfd->ExtProfBout = ExtProfBout;
	gfd->IntProfCent = IntProfCent;
	gfd->IntProfBout = IntProfBout;
	/* --------                                 ------------ */
}

//TODO: design load calling
void readFileDesignExt(int /*control*/) {
    CString fileName;
    char* PtrfileName;
    KiteDesign* oldKd;

    CFileDialog DlgOpen(TRUE, NULL, "*.wdn", OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        fileName = DlgOpen.GetPathName();
        PtrfileName = fileName.GetBuffer(1);
        strcpy(fileNameDesignExt, PtrfileName);
        oldKd = kiteDesignExt;
        kiteDesignExt = readKiteDesignFromFile(fileNameDesignExt);
		delete(oldKd);
        FicDesignExt->set_text(fileNameDesignExt);
    }
    apply(0);
    display();
    glui->sync_live();
}

void readFileDesignInt(int /*control*/) {
    CString fileName;
    char* PtrfileName;
    KiteDesign* oldKd;

    CFileDialog DlgOpen(TRUE, NULL, "*.wdn", OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        fileName = DlgOpen.GetPathName();
        PtrfileName = fileName.GetBuffer(1);
        strcpy(fileNameDesignInt, PtrfileName);
        oldKd = kiteDesignInt;
        kiteDesignInt = readKiteDesignFromFile(fileNameDesignInt);
		delete(oldKd);
        FicDesignInt->set_text(fileNameDesignInt);
    }
    apply(0);
    display();
    glui->sync_live();
}


void readProject(int /*control*/) {
    CString fileName;
    char* PtrfileName;
    WindPatternsProject* oldWpp;
    CFileDialog DlgOpen(TRUE, NULL, "*.wpp", OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        fileName = DlgOpen.GetPathName();
        PtrfileName = fileName.GetBuffer(1);
        strcpy(fileNameProject, PtrfileName);
        oldWpp = gfd;
        gfd = readWindPatternsProject(fileNameProject);
        LoadFromWindPatternsProject(gfd);
        delete(oldWpp);
        FicProject->set_text(fileNameProject);
    }
    //apply(0);
    display();
    glui->sync_live();
}


void modifSpinnerExt(int /*control*/) {
	TMesh* MeshCour;
	MeshCour=Axe3d->Mesh;
	while (MeshCour != NULL)
	{	
		if (MeshCour->side == EXT_SIDE) {
			MeshCour->colorFaces[3] = opacExt;
		}
		MeshCour=MeshCour->MeshSuiv;
	}
	display();
}

void modifSpinnerInt(int /*control*/) {
	TMesh* MeshCour;
	MeshCour=Axe3d->Mesh;
	while (MeshCour != NULL)
	{	
		if (MeshCour->side == INT_SIDE) {
			MeshCour->colorFaces[3] = opacInt;
		}
		MeshCour=MeshCour->MeshSuiv;
	}
	display();
}


void initAxe2d() {
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
}

void apply(int /*control*/) {
	printf ("\nclearMeshsAxe(Axe3d);");
    clearMeshsAxe(Axe3d);

    Matrix *XExt, *YExt, *ZExt, *YExt0;
    Matrix *XInt, *YInt, *ZInt, *YInt0;
    calcForm3D(F, 0, 0.0f,
            ExtProfCent, IntProfCent, ExtProfBout, IntProfBout,
            &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);
	Form3D* f3d = getForm3D(F, 0, 0.0f);
	FormProjection* fp = getFormProjection(f3d);

	printf ("\nOpacity:");
	printf ("\nopacExt=%f opacInt=%f", opacExt, opacInt);

	addForm3D(Axe3d, f3d->XExt, f3d->YExt, f3d->ZExt, f3d->XInt, f3d->YInt, f3d->ZInt, ViewFace, ViewSymetrique);
	
	// 2D
	//------------------------ Projections window calcation -------------------------------------
	initAxe2d();

	clearCourbesAxe(AxeProjections);
    clearMeshsAxe(AxeProjections);
	
	addForm3d2dKiteDesign( Axe3d, AxeProjections, f3d, fp, kiteDesignExt, ExtPrjNoseUp, opacExt, kiteDesignInt, IntPrjNoseUp, opacInt, ViewFace, ViewSymetrique);

	addFormProjectionCourbesToAxe( AxeProjections, fp, kiteDesignInt, ViewSymetrique, 0.0, IntPrjNoseUp );
	addFormProjectionCourbesToAxe( AxeProjections, fp, kiteDesignExt, ViewSymetrique, 2.0, ExtPrjNoseUp );

	// Zoom
    CourbZoom = new Courbe("Zoom");
    CourbZoom->pts = Zeros(5, 2);
    CourbZoom->points = OFF;
    CourbZoom->segments = OFF;
    CourbZoom->symX = OFF;
    CourbZoom->colorSegments[0] = 0.0f;
    CourbZoom->colorSegments[2] = 0.0f;
    addCourbe(AxeProjections, CourbZoom);

	display();
}

void modifXYGrid(int /*control*/) {
	//printf ("\n modifXYGrid()");
	AxeProjections->XGrid = XYGrid;
	AxeProjections->YGrid = XYGrid;
	display();
}

void modifExtPrjNoseUp(int /*control*/) {
	apply(0);
	display();
}

void modifIntPrjNoseUp(int /*control*/) {
	apply(0);
	display();
}

void modifProjection3d(int /*control*/) {
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
    ViewAxe(Axe3d);
    //ViewAxe(Axe3dBal);
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
	window2d = glutCreateWindow("2d design");
    InitWindow(); 
    AxeProjections = createAxe(window2d);
	int ws, hs;
    ws = glutGet(GLUT_SCREEN_WIDTH);
    hs = glutGet(GLUT_SCREEN_HEIGHT);
    glutSetWindow(window2d);
    glutPositionWindow((int) (0.0 * (double) ws), (int) (0.50 * (double) hs));
    glutReshapeWindow((int) (1.0 * (double) ws), (int) (0.45 * (double) hs));

	window3d = glutCreateWindow("3d design");
    InitWindow();
    InitLight();

    Axe3d = createAxe(window3d);
    Axe3d->axe3d = ON;

    glutSetWindow(window3d);
    glutPositionWindow((int) (0.25 * (double) ws), (int) (0.0 * (double) hs));
    glutReshapeWindow((int) (0.75 * (double) ws), (int) (0.48 * (double) hs));

	// main form whith GUI
    glui = GLUI_Master.create_glui("Wind designer", 0, 0, 0);
	GLUI_Panel *panelProject = glui->add_panel("");

    GLUI_Button *btnLoadForm = glui->add_button_to_panel(panelProject, "Load project", 0, &readProject);
    btnLoadForm->set_w(10);

    glui->add_column_to_panel(panelProject, false);
	// loading Default project
	strcpy(fileNameProject, "f17project.wpp");
	gfd = readWindPatternsProject(fileNameProject);
    LoadFromWindPatternsProject(gfd);
    FicProject = glui->add_statictext_to_panel(panelProject, "???");
    FicProject->set_text(fileNameProject);

	GLUI_Panel *panelDesign = glui->add_panel("");
    GLUI_Button *btnLoadDesignExt = glui->add_button_to_panel(panelDesign, "Load design Ext", 0, &readFileDesignExt);
    btnLoadDesignExt->set_w(10);

	GLUI_Button *btnLoadDesignInt = glui->add_button_to_panel(panelDesign, "Load design Int", 0, &readFileDesignInt);
    btnLoadDesignInt->set_w(10);

    glui->add_column_to_panel(panelDesign, false);
	// load default design
	strcpy(fileNameDesignExt, "f17ext.wdn");
    FicDesignExt = glui->add_statictext_to_panel(panelDesign, "???");
    FicDesignExt->set_text(fileNameDesignExt);
	kiteDesignExt = readKiteDesignFromFile(fileNameDesignExt);

	strcpy(fileNameDesignInt, "f17int.wdn");
    FicDesignInt = glui->add_statictext_to_panel(panelDesign, "???");
    FicDesignInt->set_text(fileNameDesignInt);
	kiteDesignInt = readKiteDesignFromFile(fileNameDesignInt);

	//FicDesign->set_h(10);
	GLUI_Panel *panel2dOptions = glui->add_panel("2D");
	glui->add_checkbox_to_panel(panel2dOptions, "XY grid", &XYGrid, 0,  &modifXYGrid);

	glui->add_checkbox_to_panel(panel2dOptions, "Ext Projection nose up", &ExtPrjNoseUp, 0,  &modifExtPrjNoseUp);
	glui->add_checkbox_to_panel(panel2dOptions, "Int Projection nose up", &IntPrjNoseUp, 0,  &modifIntPrjNoseUp);

    GLUI_Panel *panelViewOptions = glui->add_panel("");
    glui->add_checkbox_to_panel(panelViewOptions, "symetrique", &ViewSymetrique, 0, &modifViewSymetrique);
    GLUI_Panel *panel_proj = glui->add_panel_to_panel(panelViewOptions, "Projection");
    GLUI_RadioGroup *radio_proj =
            glui->add_radiogroup_to_panel(panel_proj, &ProjOrthoPers, 0, &modifProjection3d);
    glui->add_radiobutton_to_group(radio_proj, "Orthogonale");
    glui->add_radiobutton_to_group(radio_proj, "Perspective");

    GLUI_Panel *panelOpacity = glui->add_panel("Opacity");
	GLUI_Spinner *extTranspSpinner = glui->add_spinner_to_panel(panelOpacity, "Ext", GLUI_SPINNER_FLOAT, &(opacExt), 0, &modifSpinnerExt);
	extTranspSpinner->set_float_limits(0.0f,1.0f);

	GLUI_Spinner *intTranspSpinner = glui->add_spinner_to_panel(panelOpacity, "Int", GLUI_SPINNER_FLOAT, &(opacInt), 0, &modifSpinnerInt);
	
	intTranspSpinner->set_float_limits(0.0f,1.0f);
    glui->add_button("apply", 0, &apply);
	// ----- end of constructing GUI ----- 

    GLUI_Master.set_glutIdleFunc(NULL);
	glui->sync_live();

	apply(0);
	AxeSel = Axe3d;
    display();
    glutMainLoop();
	return 0;
}

