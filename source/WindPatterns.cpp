#pragma warning(disable:4786)
#pragma warning(disable:4514)

/***********/
/* include */
/***********/
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


#include "GL/glui.h"
#include "GL/glut.h"

/**********/
/* define */
/**********/
#define AUTEUR_DATE "(Thierry Pebayle - 18/08/2002)\n[Export FGen + Refactoring Olivier REMAUD 2007]\n[ Roman Lybimcev, Ivan Kobzarev 2009]"
#define sqr(f1) ((f1)*(f1))
#define pi	3.141592675f
//#define DEBUG 0
#define TEST_BUTTON true

/**********************/
/* variables globales */
/**********************/

/*live-variable GLUI*/
int ViewFace = 0; /*choix visu nervures(0) ou surfaces(1)*/
int ProjOrthoPers = 0; /*choix projection orthogonale ou perspective*/
int PinceType = 0;

int PinceNosRadio = 1;
int PinceHvostRadio = 1;

int ReperPointsFromFile = 1;
int ReperPointsPlotterFormat = 1;
int VentHoles = 1;
int DiagNervs = 1;
int VentHolesDouble = 0;
int isKlapan = 0;
int VentCentralNerv = 0;
float VentHolesDeb = 2.0f;
float VentHolesFin = 6.0f;

int ViewSymetrique = 0;
int LayoutSymetrique = 0;
int LayoutKlapans = 0;
int ViewPtsSuspentes = 1; /*visualise pts de suspentage*/
int ReperesCoutures = 1;
int ReperesSuspentes[2] = {1, 1};
int ReperesProfile[2] = {1, 1};

int Numerotation = 0;
int Ventilation = 0;
int RepPoints = 1;
int VentilationLayout = 1;
int layoutWithDesign = 0;
int CorrectRepPoints = 1;

int PincePowerA = 0;
int PincePowerF = 1;

int PinceArctanA = 1;
int PinceArctanF = 0;

int PinceNosEqualAmp = 0;
int PinceHvostEqualAmp = 0;
//int PinceNos0Amp = 0;
//int PinceHvost0Amp = 0;

int GoPince = 0;

float XMashtab = 0.5f;
float coeffMult = 1.0f;
//pour calc patrons
int NoNerv[2] = {2, 3};
float Deb[2] = {0.0f, 20.0f}, Fin[2] = {100.0f, 100.0f};

float PincePowerValueA=1.3f, PinceArctanK1ValueA=0.5f,  PinceArctanK2ValueA=3.0f,  PinceArctanK3ValueA=0.5f;
float PincePowerValueF=1.3f, PinceArctanK1ValueF=0.5f,  PinceArctanK2ValueF=3.0f, PinceArctanK3ValueF=0.5f;
float PosDiagNerv2A=8.0f, PosDiagNerv2F=71.0f, PosDiagNerv1A=2.0f, PosDiagNerv1F=80.0f;
float PosKlapanInt=8.0f, PosKlapanFin=16.0f;

float PinceRadiusAlgKHvost=0.985f, PinceRadiusAlgKNos=0.985f;

int FaceDeb[2], FaceFin[2];
int isPosNerv[2]={0, 0};
float PosNerv[2]={0.0f, 100.0f};
float Marge[2] = {1.0f, 1.0f};
//double PosRep[2];
int FaceRep[2];
float MargeDeb = 1.0f, MargeFin = 0.0f;

float PosPinceBA[2] = {14.0, 14.0};
float PosPinceBF[2] = {25.0, 25.0};
float AmpPinceBA[2] = {2.5f, 2.5f};
float AmpPinceBF[2] = {2.0f, 2.0f};

float tochnostLayout2 = 0.7f;
float margeFinExt = 4.0f;
float margeFinInt = 1.0f;
float margeFinNerv = 1.0f;
float margeFinDiagNerv = 1.0f;

float textX = 0.45f;
float textY = 0.3f;

double *pincePercentage1, *pincePercentage2;
bool *pinceCalc;

int zoomIN = false, debZoomIN = false, finZoomIN = false, zoomOUT = false, quitZoom = false;

Courbe *CourbZoom;

/*global pour control rotation, ...*/
int xSouris, ySouris;
int Rotation3d = OFF; /*flag rotation 3d -> bouton gauche souris*/
int Translation3d = OFF; /*flag translation 3d -> bouton droit souris*/
int ZoomTwist3d = OFF; /*flag twist 3d -> bouton gauche et droit souris*/

/*fenetres OpenGL*/
int window3d; /*visu 3D*/
//int window3dBal; /*visu 3D*/
int windowPatron; /*visu patron*/

double *pinceLAAmp1;
double *pinceLFAmp1;
double *pinceRAAmp1;
double *pinceRFAmp1;

Matrix** funcL1;
Matrix** funcL2;
Matrix** funcR1;
Matrix** funcR2;

double *lenp1;

double *pinceLAAmp2;
double *pinceLFAmp2;
double *pinceRAAmp2;
double *pinceRFAmp2;
double *lenp2;
double *coeffn;
double *coeffd;
int quantDiag = 0;
int quantVH = 0;
int *noNervD;
int *noNervVH;
int *noNervD1;
int *noNervD2;

/* axes graphiques*/
TAxe *Axe3d;
//TAxe *Axe3dBal;
TAxe *AxeSel;
TAxe *AxePatron;
TAxe *AxePatronDXF, *AxePatronTextDXF, *AxeMarginDXF, *AxeCercleDXF, *AxeRepDXF;
//Courbe *CourbeIntDXF;

/*forme*/
Form *F;
char ProjectName[255];

/*divers tableaux*/
Matrix *IntProfCent, *ExtProfCent, *IntProfBout, *ExtProfBout;
Matrix** ReperPoints;
Ballonement* ballonement;
char ballonementPath[255];
char fileNameForm[255];
char fileNameProject[255];
char fileNameRepPoints[255];
char fileNameVentHoles[255];
char fileNameDiagNerv[255];

char fileNameDesignExt[255];
char fileNameDesignInt[255];


GLUI *glui;
GLUI_StaticText *FicForm, *FicRepPoints, *FicVentHoles, *FicDiagNerv, *BlankST, *BlankST2, *FicProject;
GLUI_StaticText *FicDesignExt;
GLUI_StaticText *FicDesignInt;
GLUI_EditText *NumText, *DiagNervText, *VentHolesNervsText;
GLUI_Spinner *SpinNoNerv[2];

KiteDesign* kiteDesignExt = 0;
KiteDesign* kiteDesignInt = 0;

/***********************************************/
/*liste des procedures definies dans ce fichier*/
/***********************************************/

void writeManyFichierPolyDXF(char *fileName, int np, int n, TAxe **axe, TAxe **axe2, int rep, TAxe **axeR, int vent, TAxe **axeC, int num, TAxe **axeT, double* W, double* H);
void writeManyFichierPolyDXF2(char *fileName, int np, int n, TAxe **axe, TAxe **axe2, int rep, TAxe **axeR, int vent, TAxe **axeC, int num, TAxe **axeT, double* W, double* H, int* numncon);

void display(void);
void reshape(int w, int h);
void motion(int x, int y);
void BoutonSouris(int button, int state, int x, int y);
void Positionnewindows(void);
void keyboard(unsigned char key, int x, int y);
void InitWindow(void);
void InitLight(void);
void calcVue3dEtPatron(void);
void calcMagicPinceR(int noNerv1, double perc, int face, double* percNew);
void modifProjection3d(int control);
void modifView3d(int control);
void modifViewSymetrique(int control);
void apply(int control);
void applyMagic(int control);

void quit(int control);
void Info(int control);
void readForm(int control);

Matrix** readFichierReperPoints(char* NomFic);

void initValueDialogue(void);
void adder(int control);
void getPointByPos (Matrix *Xd, Matrix *Yd, Matrix *P, double Pos, double *xr, double *yr);
void calcMagicPincesLayout();

WindPatternsProject* gfd = new WindPatternsProject();

WindPatternsProject* getWindPatternsProject() {
    strcpy(gfd->name,ProjectName);
    strcpy(gfd->fileNameForm, fileNameForm);

    strcpy(gfd->ballonementPath, ballonementPath);
	gfd->ballonement=ballonement;
    gfd->ReperPointsPlotterFormat = ReperPointsPlotterFormat;
    gfd->PincePowerA = PincePowerA;
    gfd->PincePowerF = PincePowerF;
    gfd->PinceArctanA =PinceArctanA;
    gfd->PinceArctanF = PinceArctanF;
    gfd->PinceArctanK1ValueA = PinceArctanK1ValueA;
    gfd->PinceArctanK1ValueF = PinceArctanK1ValueF;
    gfd->PinceArctanK2ValueA = PinceArctanK2ValueA;
    gfd->PinceArctanK2ValueF = PinceArctanK2ValueF;
    gfd->PinceArctanK3ValueA = PinceArctanK3ValueA;
    gfd->PinceArctanK3ValueF = PinceArctanK3ValueF;
    gfd->PincePowerValueA = PincePowerValueA;
    gfd->PincePowerValueF = PincePowerValueF;

    gfd->PinceNosRadio = PinceNosRadio;
    gfd->PinceHvostRadio = PinceHvostRadio;
    gfd->PosPinceBA[0]=PosPinceBA[0];
    gfd->PosPinceBA[1]=PosPinceBA[1];
    gfd->AmpPinceBA[0]=AmpPinceBA[0];
    gfd->AmpPinceBA[1]=AmpPinceBA[1];

    gfd->PosPinceBF[0]=PosPinceBF[0];
    gfd->PosPinceBF[1]=PosPinceBF[1];
    gfd->AmpPinceBF[0]=AmpPinceBF[0];
    gfd->AmpPinceBF[1]=AmpPinceBF[1];

    gfd->PinceRadiusAlgKNos=PinceRadiusAlgKNos;
    gfd->PinceRadiusAlgKHvost=PinceRadiusAlgKHvost;
    gfd->Form=F;
    gfd->LayoutSymetrique=LayoutSymetrique;
    gfd->PosDiagNerv1A=PosDiagNerv1A;
    gfd->PosDiagNerv1F=PosDiagNerv1F;
    gfd->PosDiagNerv2A=PosDiagNerv2A;
    gfd->PosDiagNerv2F=PosDiagNerv2F;
    gfd->windowPatron=windowPatron;
    gfd->MargeDeb = MargeDeb;
    gfd->MargeFin = MargeFin;
    gfd->Marge[0] = Marge[0];
    gfd->Marge[1] = Marge[1];

    gfd->ReperPointsFromFile=ReperPointsFromFile;
    gfd->ReperPoints=ReperPoints;
    gfd->ReperesSuspentes[0]=ReperesSuspentes[0];
    gfd->ReperesSuspentes[1]=ReperesSuspentes[1];
    gfd->ReperesProfile[0]=ReperesProfile[0];
    gfd->ReperesProfile[1]=ReperesProfile[1];
    gfd->XMashtab=XMashtab;

	gfd->CorrectRepPoints = CorrectRepPoints;
    gfd->Ventilation=Ventilation;
    gfd->VentilationLayout=VentilationLayout;
    gfd->textX=textX;
    gfd->textY=textY;

    gfd->tochnostLayout2=tochnostLayout2;

    gfd->PinceNosEqualAmp=PinceNosEqualAmp;
    gfd->PinceHvostEqualAmp=PinceHvostEqualAmp;

    gfd->Deb[0]=Deb[0];
    gfd->Deb[1]=Deb[1];
    gfd->Fin[0]=Fin[0];
    gfd->Fin[1]=Fin[1];

    gfd->VentHoles = VentHoles;
    gfd->noNervVH = noNervVH;
    gfd->VentHolesDeb = VentHolesDeb;
    gfd->VentHolesFin = VentHolesFin;
    gfd->VentHolesDouble = VentHolesDouble;
    gfd->LayoutKlapans = LayoutKlapans;

    gfd->quantVH = quantVH;
    gfd->VentCentralNerv=VentCentralNerv;
    //ParseDiagNervs();
	gfd->DiagNervs=DiagNervs;
    gfd->quantDiag=quantDiag;
    gfd->noNervD = noNervD;

    gfd->margeFinExt = margeFinExt;
    gfd->margeFinInt = margeFinInt;
    gfd->margeFinNerv = margeFinNerv;
    gfd->margeFinDiagNerv = margeFinDiagNerv;

    gfd->PosKlapanFin = PosKlapanFin;

    gfd->isPosNerv = isPosNerv;
    gfd->PosNerv = PosNerv;

    gfd->ExtProfCent=ExtProfCent;
    gfd->ExtProfBout=ExtProfBout;
    gfd->IntProfCent=IntProfCent;
    gfd->IntProfBout=IntProfBout;
    strcpy(gfd->fileNameDiagNerv, fileNameDiagNerv);
    strcpy(gfd->fileNameVentHoles, fileNameVentHoles);
    strcpy(gfd->fileNameRepPoints, fileNameRepPoints);

    FicForm->set_text(fileNameForm);
    FicVentHoles->set_text(fileNameVentHoles);
    FicRepPoints->set_text(fileNameRepPoints);
    FicDiagNerv->set_text(fileNameDiagNerv);
    return gfd;
}

void LoadFromWindPatternsProject(WindPatternsProject* gfd) {
	Form* tmpF=0;
	Matrix *tmpIntProfCent=0, *tmpExtProfCent=0, *tmpIntProfBout=0, *tmpExtProfBout=0;
	Matrix** tmpReperPoints=0;
	Ballonement* tmpBallonement;
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
		tmpnoNervVH = readFichierVentHoles(gfd->fileNameVentHoles, &tmpquantVH, &tmpVentCentralNerv);
		tmpnoNervD = readFichierDiagNervs(gfd->fileNameDiagNerv, &tmpquantDiag);
		tmpReperPoints = readFichierReperPoints(gfd->fileNameRepPoints);
		tmpBallonement = readBallonementFromFile(gfd->ballonementPath);
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
	noNervVH = tmpnoNervVH;
	quantVH = tmpquantVH;
	VentCentralNerv = tmpVentCentralNerv;
	noNervD = tmpnoNervD;
	quantDiag = tmpquantDiag;
	ReperPoints = tmpReperPoints;
	ballonement = tmpBallonement;
	/* --------                                 ------------ */
 
    strcpy(fileNameRepPoints, gfd->fileNameRepPoints);
    strcpy(fileNameDiagNerv, gfd->fileNameDiagNerv);
    strcpy(fileNameVentHoles, gfd->fileNameVentHoles);
	strcpy(fileNameForm, gfd->fileNameForm);
	strcpy(ballonementPath, gfd->ballonementPath);
    strcpy(ProjectName, gfd -> name);	
	PincePowerA = gfd->PincePowerA;
    PincePowerF = gfd->PincePowerF;
    PinceArctanA = gfd->PinceArctanA;
    PinceArctanF = gfd->PinceArctanF;
    PinceArctanK1ValueA = gfd->PinceArctanK1ValueA;
    PinceArctanK1ValueF = gfd->PinceArctanK1ValueF;
    PinceArctanK2ValueA = gfd->PinceArctanK2ValueA;
    PinceArctanK2ValueF = gfd->PinceArctanK2ValueF;
    PinceArctanK3ValueA = gfd->PinceArctanK3ValueA;
    PinceArctanK3ValueF = gfd->PinceArctanK3ValueF;
    PincePowerValueA = gfd->PincePowerValueA;
    PincePowerValueF = gfd->PincePowerValueF;
    
    PinceNosRadio = gfd->PinceNosRadio;
    PinceHvostRadio = gfd->PinceHvostRadio;
    PosPinceBA[0]=gfd->PosPinceBA[0];
    PosPinceBA[1]=gfd->PosPinceBA[0];
    AmpPinceBA[0]=gfd->AmpPinceBA[0];
    AmpPinceBA[1]=gfd->AmpPinceBA[0];

    PosPinceBF[0]=gfd->PosPinceBF[0];
    PosPinceBF[1]=gfd->PosPinceBF[0];
    AmpPinceBF[0]=gfd->AmpPinceBF[0];
    AmpPinceBF[1]=gfd->AmpPinceBF[0];

    PinceRadiusAlgKNos=gfd->PinceRadiusAlgKNos;
    PinceRadiusAlgKHvost=gfd->PinceRadiusAlgKHvost;
    //LayoutSymetrique=gfd->LayoutSymetrique;
    PosDiagNerv1A=gfd->PosDiagNerv1A;
    PosDiagNerv1F=gfd->PosDiagNerv1F;
    PosDiagNerv2A=gfd->PosDiagNerv2A;
    PosDiagNerv2F=gfd->PosDiagNerv2F;

    MargeDeb = gfd->MargeDeb;
    MargeFin = gfd->MargeFin;
    
    margeFinExt=gfd->margeFinExt;
    margeFinInt=gfd->margeFinInt;
    margeFinNerv=gfd->margeFinNerv;
    margeFinDiagNerv=gfd->margeFinDiagNerv;
    
    Marge[0] = gfd->Marge[0];
    Marge[1] = gfd->Marge[1];
    
    ReperPointsFromFile = gfd->ReperPointsFromFile;
    RepPoints = gfd -> RepPoints;

    ReperesProfile[0]=gfd->RepPoints;
    ReperesProfile[1]=gfd->RepPoints;
    ReperesSuspentes[0]=gfd->RepPoints;
    ReperesSuspentes[1]=gfd->RepPoints;

    XMashtab=gfd->XMashtab;
    Ventilation=gfd->Ventilation;
    VentilationLayout=gfd->VentilationLayout;
    textX=gfd->textX;
    textY=gfd->textY;

    tochnostLayout2=gfd->tochnostLayout2;
    PinceNosEqualAmp=gfd->PinceNosEqualAmp;
    PinceHvostEqualAmp=gfd->PinceHvostEqualAmp;
    VentHoles = gfd->VentHoles;
    VentHolesDeb = gfd->VentHolesDeb;
    VentHolesFin = gfd->VentHolesFin;
    VentHolesDouble = gfd->VentHolesDouble;
    LayoutKlapans = gfd->LayoutKlapans;
    PosKlapanFin = gfd->PosKlapanFin;
}


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

void initValueDialogue(void) {
    NoNerv[0] = 7;
    NoNerv[1] = 7;

    Deb[0] = 0.0f;
    Deb[1] = 0.0f;

    Fin[0] = 100.0f;
    Fin[1] = 100.0f;

    FaceDeb[0] = 2;
    FaceDeb[1] = 1;

    FaceFin[0] = 2;
    FaceFin[1] = 1;

    SpinNoNerv[0] -> set_int_limits(-1, F->m_nbProfils - 1);
    SpinNoNerv[1] -> set_int_limits(-1, F->m_nbProfils - 1);
}



void readRepPoints(int /*control*/) {
    CString fileName;
    char* PtrfileName;
    Matrix** oldRP;
    CFileDialog DlgOpen(TRUE, NULL, "*.txt", OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        fileName = DlgOpen.GetPathName();
        PtrfileName = fileName.GetBuffer(1);
        strcpy(fileNameRepPoints, PtrfileName);
        oldRP = ReperPoints;
        ReperPoints = readFichierReperPoints(fileNameRepPoints);
        delete(oldRP);
        FicRepPoints->set_text(fileNameRepPoints);
    }
    if (ReperPointsFromFile==1) apply(0);
}


void readDiagNervs(int /*control*/) {
    CString fileName;
    char* PtrfileName;
    int* oldNoNervD;
    CFileDialog DlgOpen(TRUE, NULL, "*.txt", OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        fileName = DlgOpen.GetPathName();
        PtrfileName = fileName.GetBuffer(1);
        strcpy(fileNameDiagNerv, PtrfileName);
        oldNoNervD = noNervD;
        
        noNervD = readFichierDiagNervs(fileNameDiagNerv, &quantDiag);
        delete(oldNoNervD);
        FicDiagNerv->set_text(fileNameDiagNerv);
    }
}


void readVentHoles(int /*control*/) {
    //printf ("\nCharger Fichier Vent Holes()");
    CString fileName;
    char* PtrfileName;
    int* oldVH;
    //ouverture boite de dialogue
    CFileDialog DlgOpen(TRUE, NULL, "*.txt", OFN_OVERWRITEPROMPT, NULL, NULL);

    if (DlgOpen.DoModal() == IDOK) {
        //recupere nom de fichier
        fileName = DlgOpen.GetPathName();
        PtrfileName = fileName.GetBuffer(1);
        strcpy(fileNameVentHoles, PtrfileName);
        oldVH = noNervVH;
        VentCentralNerv = 0;
        noNervVH = readFichierVentHoles(fileNameVentHoles, &quantVH, &VentCentralNerv);
        delete(oldVH);
        FicVentHoles->set_text(fileNameVentHoles);
    }
}


void readProject(int /*control*/) {
    CString fileName;
    char* PtrfileName;
    WindPatternsProject* oldWpp;
    //ouverture boite de dialogue
    CFileDialog DlgOpen(TRUE, NULL, "*.wpp", OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        //recupere nom de fichier
        fileName = DlgOpen.GetPathName();
        PtrfileName = fileName.GetBuffer(1);
        strcpy(fileNameProject, PtrfileName);
        oldWpp = gfd;
        gfd = readWindPatternsProject(fileNameProject);
        LoadFromWindPatternsProject(gfd);
        delete(oldWpp);
        FicProject->set_text(fileNameProject);
    }
    apply(0);
    display();
    glui->sync_live();
}

void readForm(int /*control*/) {
    CString fileName;
    char* PtrfileName;
    Form* oldF;
    //ouverture boite de dialogue
    CFileDialog DlgOpen(TRUE, NULL, "*.txt", OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        //recupere nom de fichier
        fileName = DlgOpen.GetPathName();
        PtrfileName = fileName.GetBuffer(1);
        strcpy(fileNameForm, PtrfileName);
        //read fichier de forme
        //clearForm(F);
        oldF = F;
        F = readFichierForm(fileNameForm);
       //libere memoire ancienne forme
        delete(oldF);
        //maj interface GLUI
        FicForm->set_text(fileNameForm);

        //chargement profils
        readFichierProfil(F->m_strNomProfilCent.c_str(), &ExtProfCent, &IntProfCent);
        readFichierProfil(F->m_strNomProfilBout.c_str(), &ExtProfBout, &IntProfBout);
        InterpoleProfilBout(&ExtProfBout, ExtProfCent);
        InterpoleProfilBout(&IntProfBout, IntProfCent);
        

        //par defaut pas de zoom, et axe automatique
        quitZoom = false;
        zoomIN = false;
        debZoomIN = false;
        finZoomIN = false;
        zoomOUT = false;

        AxePatron->XAuto = ON;
        AxePatron->YAuto = ON;
        AxePatron->ZAuto = ON;
        //mise a jour limite/valeur boite de dialogue
        initValueDialogue();
    }
    //mise a jour affichage
    //calcVue3dEtPatron();
    apply(0);
    display();

}

void readForm2(int /*control*/) {
    CString fileName;
    char* PtrfileName;
    Form* oldF;
    //ouverture boite de dialogue
    CFileDialog DlgOpen(TRUE, NULL, "*.txt", OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        //recupere nom de fichier
        fileName = DlgOpen.GetPathName();
        PtrfileName = fileName.GetBuffer(1);
        strcpy(fileNameForm, PtrfileName);
        //read fichier de forme
        //clearForm(F);
        oldF = F;
        F = readFichierForm2(fileNameForm);
        //libere memoire ancienne forme
        delete(oldF);
        //maj interface GLUI
        FicForm->set_text(fileNameForm);
        //chargement profils
        readFichierProfil(F->m_strNomProfilCent.c_str(), &ExtProfCent, &IntProfCent);
        readFichierProfil(F->m_strNomProfilBout.c_str(), &ExtProfBout, &IntProfBout);
        InterpoleProfilBout(&ExtProfBout, ExtProfCent);
        InterpoleProfilBout(&IntProfBout, IntProfCent);
        //par defaut pas de zoom, et axe automatique
        quitZoom = false;
        zoomIN = false;
        debZoomIN = false;
        finZoomIN = false;
        zoomOUT = false;

        AxePatron->XAuto = ON;
        AxePatron->YAuto = ON;
        AxePatron->ZAuto = ON;

        //mise a jour limite/valeur boite de dialogue
        initValueDialogue();
    }
    //mise a jour affichage
    calcVue3dEtPatron();
    display();
}

void modifPinceNosRadio(int /*control*/) {
    calcVue3dEtPatron();
    display();
}

void modifRepPts(int /*control*/) {
    calcVue3dEtPatron();
    display();
}

void modifVentilation(int /*control*/) {
    apply(0);
}

void modifPinceHvostRadio(int /*control*/) {
    calcVue3dEtPatron();
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

/***************/
/* modifView3d */
/***************/
void modifView3d(int /*control*/) {
    /*recalc et retrace la vue 3d*/
    calcVue3dEtPatron();
    ViewAxe(Axe3d);
    glutSwapBuffers();
    //ViewAxe(Axe3dBal);
    //glutSwapBuffers();

}

/***********************/
/* modifViewSymetrique */
/***********************/
void modifViewSymetrique(int /*control*/) {
    /*recalc et retrace la vue 3d*/
    calcVue3dEtPatron();
    ViewAxe(Axe3d);
    glutSwapBuffers();
    //ViewAxe(Axe3dBal);
    //glutSwapBuffers();

}

void modifGoPince(int /*control*/) {
    calcVue3dEtPatron();
    display();
}

void modifLayoutSymetrique(int /*control*/) {

}

void modifVentilationLayout(int /*control*/) {

}

void modifLayoutWithDesign(int /*control*/) {

}

void adder(int /*control*/) {

}

void apply(int /*control*/) {
    int i; //,j;
    Matrix *xe, *xi;
    Matrix *extProf, *intProf;
    xe = new Matrix(ExtProfCent->GetLignes(), 1);
    xi = new Matrix(IntProfCent->GetLignes(), 1);
    for (i = 0; i < ExtProfCent->GetLignes(); i++)        xe->SetElement(i, 0, ExtProfCent->Element(i, 0));
    for (i = 0; i < IntProfCent->GetLignes(); i++)        xi->SetElement(i, 0, IntProfCent->Element(i, 0));
    for (i = 0; i < 2; i++) {
        if ((FaceDeb[i] == 1) && (!ValeurPresente(Deb[i], xe)))            addeValeurCroissant(Deb[i], &xe);
        if ((FaceDeb[i] == 2) && (!ValeurPresente(Deb[i], xi)))            addeValeurCroissant(Deb[i], &xi);
        if ((FaceFin[i] == 1) && (!ValeurPresente(Fin[i], xe)))            addeValeurCroissant(Fin[i], &xe);
        if ((FaceFin[i] == 2) && (!ValeurPresente(Fin[i], xi)))            addeValeurCroissant(Fin[i], &xi);
    }

    if ((!ValeurPresente(PosDiagNerv2A, xi))) addeValeurCroissant(PosDiagNerv2A, &xi);
    if ((!ValeurPresente(PosDiagNerv2F, xi))) addeValeurCroissant(PosDiagNerv2F, &xi);
    if ((!ValeurPresente(PosDiagNerv1A, xe))) addeValeurCroissant(PosDiagNerv1A, &xe);
    if ((!ValeurPresente(PosDiagNerv1F, xe))) addeValeurCroissant(PosDiagNerv1F, &xe);

    if ((!ValeurPresente(PosPinceBA[0], xi))) addeValeurCroissant(PosPinceBA[0], &xi);
    if ((!ValeurPresente(100.0f-PosPinceBF[0], xi))) addeValeurCroissant(100.0f-PosPinceBF[0], &xi);
    if ((!ValeurPresente(PosPinceBA[0], xe))) addeValeurCroissant(PosPinceBA[0], &xe);
    if ((!ValeurPresente(100.0f-PosPinceBF[0], xe))) addeValeurCroissant(100.0f-PosPinceBF[0], &xe);

    if (VentHoles) {
        if ((!ValeurPresente(VentHolesDeb, xi))) addeValeurCroissant(VentHolesDeb, &xi);
        if ((!ValeurPresente(VentHolesFin, xi))) addeValeurCroissant(VentHolesFin, &xi);
    }
    if (isKlapan) {
        if ((!ValeurPresente(PosKlapanInt, xi))) addeValeurCroissant(PosKlapanInt, &xi);
        if ((!ValeurPresente(PosKlapanFin, xe))) addeValeurCroissant(PosKlapanFin, &xe);
    }

    extProf = Zeros(xe->GetLignes(), 2);
    intProf = Zeros(xi->GetLignes(), 2);
    for (i = 0; i < extProf->GetLignes(); i++)
        extProf->SetElement(i, 0, xe->Element(i, 0));

    for (i = 0; i < intProf->GetLignes(); i++)
        intProf->SetElement(i, 0, xi->Element(i, 0));

    InterpoleProfilBout(&ExtProfBout, extProf);
    InterpoleProfilBout(&IntProfBout, intProf);
    InterpoleProfilBout(&ExtProfCent, extProf);
    InterpoleProfilBout(&IntProfCent, intProf);

    //liberation memoire

    delete(xe);
    delete(xi);
    delete(extProf);
    delete(intProf);

    //par defaut pas de zoom, et axe automatique pour AxePatron

    quitZoom = false;
    zoomIN = false;
    debZoomIN = false;
    finZoomIN = false;
    zoomOUT = false;

    AxePatron->XAuto = ON;
    AxePatron->YAuto = ON;
    AxePatron->ZAuto = ON;
    //ParseVentHolesNervs();
    //recalc forme 3D et patron avec "nouveau" profil
    calcVue3dEtPatron();
    display();

}


void applyMagic2(int /*control*/) {
    int i; //,j;
    Matrix *xe, *xi;
    Matrix *extProf, *intProf;
    xe = new Matrix(ExtProfCent->GetLignes(), 1);
    xi = new Matrix(IntProfCent->GetLignes(), 1);
    for (i = 0; i < ExtProfCent->GetLignes(); i++)
        xe->SetElement(i, 0, ExtProfCent->Element(i, 0));
    for (i = 0; i < IntProfCent->GetLignes(); i++)
        xi->SetElement(i, 0, IntProfCent->Element(i, 0));
    for (i = 0; i < 2; i++) {
        if ((FaceDeb[i] == 1) && (!ValeurPresente(Deb[i], xe))) addeValeurCroissant(Deb[i], &xe);
        if ((FaceDeb[i] == 2) && (!ValeurPresente(Deb[i], xi))) addeValeurCroissant(Deb[i], &xi);
        if ((FaceFin[i] == 1) && (!ValeurPresente(Fin[i], xe))) addeValeurCroissant(Fin[i], &xe);
        if ((FaceFin[i] == 2) && (!ValeurPresente(Fin[i], xi))) addeValeurCroissant(Fin[i], &xi);
    }

    if ((!ValeurPresente(PosDiagNerv2A, xi))) addeValeurCroissant(PosDiagNerv2A, &xi);
    if ((!ValeurPresente(PosDiagNerv2F, xi))) addeValeurCroissant(PosDiagNerv2F, &xi);
    if ((!ValeurPresente(PosDiagNerv1A, xe))) addeValeurCroissant(PosDiagNerv1A, &xe);
    if ((!ValeurPresente(PosDiagNerv1F, xe))) addeValeurCroissant(PosDiagNerv1F, &xe);

    if ((!ValeurPresente(PosPinceBA[0], xi))) addeValeurCroissant(PosPinceBA[0], &xi);
    if ((!ValeurPresente(100.0f - PosPinceBF[0], xi))) addeValeurCroissant(100.0f - PosPinceBF[0], &xi);
    if ((!ValeurPresente(PosPinceBA[0], xe))) addeValeurCroissant(PosPinceBA[0], &xe);
    if ((!ValeurPresente(100.0f - PosPinceBF[0], xe))) addeValeurCroissant(100.0f - PosPinceBF[0], &xe);

    if (VentHoles) {
        if ((!ValeurPresente(VentHolesDeb, xi))) addeValeurCroissant(VentHolesDeb, &xi);
        if ((!ValeurPresente(VentHolesFin, xi))) addeValeurCroissant(VentHolesFin, &xi);
    }


    extProf = Zeros(xe->GetLignes(), 2);
    intProf = Zeros(xi->GetLignes(), 2);
    for (i = 0; i < extProf->GetLignes(); i++)
        extProf->SetElement(i, 0, xe->Element(i, 0));
    for (i = 0; i < intProf->GetLignes(); i++)
        intProf->SetElement(i, 0, xi->Element(i, 0));
    InterpoleProfilBout(&ExtProfBout, extProf);
    InterpoleProfilBout(&IntProfBout, intProf);
    InterpoleProfilBout(&ExtProfCent, extProf);
    InterpoleProfilBout(&IntProfCent, intProf);
    delete(xe);
    delete(xi);
    delete(extProf);
    delete(intProf);
    //printf ("gocalcPinceLayout()");
	Layout* Layout = calcIndepPinceLayout(getWindPatternsProject(), F);
    //printf ("...gocalcPinceLayout()");
	//printf ("goSaveLayout2()");
    SaveLayout2(getWindPatternsProject(), Layout);
	//printf ("...goSaveLayout2()");
}

void saveFichier3dDXF( int ) //
{
	CString fileName;
	LPTSTR PtrfileName;
	//ouverture boite de dialogue
	CFileDialog DlgOpen(FALSE, NULL, "*.dxf", OFN_OVERWRITEPROMPT, NULL, NULL);
	if (DlgOpen.DoModal() == IDOK)
	{
		//recupere nom de fichier
		fileName=DlgOpen.GetPathName();
		PtrfileName = fileName.GetBuffer(1);
       	//write fichier DXF
		writeFichierDXF(PtrfileName, Axe3d);
	}
}

void ExportWpa( int ) //
{
/*	CString fileName;
	LPTSTR PtrfileName;
	//ouverture boite de dialogue
	CFileDialog DlgOpen(FALSE, NULL, "*.wpa", OFN_OVERWRITEPROMPT, NULL, NULL);
	if (DlgOpen.DoModal() == IDOK)
	{
		//recupere nom de fichier
		fileName=DlgOpen.GetPathName();
		PtrfileName = fileName.GetBuffer(1);
		
		F->ExtProfCent = ExtProfCent;
		F->ExtProfBout = ExtProfBout;
		F->IntProfCent = IntProfCent;
		F->IntProfBout = IntProfBout;

		if (writeFichierWpa(PtrfileName, F))
			printf ("\nForm exported in WPA sucessfully");
	}*/
}


/*************/
/* save    */
/*************/

void save(int dxf) {
    //boite de dialogue fichier
    CString fileName;
    LPTSTR PtrfileName;
    //parcours liste des courbes
    Courbe *Cour;
    FILE *fid;
    char ext[10];
    int i;
    //ouverture boite de dialogue
    if (dxf > 0) strcpy(ext, "*.dxf");
    else strcpy(ext, "*.txt");
    CFileDialog DlgOpen(FALSE, NULL, ext, OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        //recupere nom de fichier
        fileName = DlgOpen.GetPathName();
        PtrfileName = fileName.GetBuffer(1);
        //write fichier patron
        if (dxf > 0) //format DXF
        {
            if (dxf == 1)
                writeFichierDXF(PtrfileName, AxePatron);
            else {
                writeFichierPolyDXF(PtrfileName, AxePatronDXF, AxeMarginDXF, (ReperesSuspentes[0] || ReperesSuspentes[1]),
                        AxeRepDXF, Ventilation, AxeCercleDXF, Numerotation, AxePatronTextDXF);
                
                /*writeFichierPolyDXFDelta(fid, AxePatronDXF, AxeMarginDXF,
                                            (ReperesSuspentes[0] || ReperesSuspentes[1]), AxeRepDXF, Ventilation, AxeCercleDXF,
                        Numerotation, AxePatronTextDXF, 0.0f, 0.0f, 0.0f, 0);*/
            }
        } else //format TXT
        {
            fid = fopen(PtrfileName, "wt");
            Cour = AxePatron->Courb;
            while (Cour != NULL) {
                fprintf(fid, "\n%d", Cour->pts->GetLignes());
                for (i = 0; i < Cour->pts->GetLignes(); i++)
                    fprintf(fid, "\n%1.5f %1.5f", Cour->pts->Element(i, 0), Cour->pts->Element(i, 1));

                Cour = Cour->CourbSuiv;

            }

        }

    }

}

void saveProject(int /*control*/) {
    //boite de dialogue fichier
    CString fileName;
    LPTSTR PtrfileName;
    //parcours liste des courbes
    Courbe *Cour;
    FILE *fid;
    char ext[10];
    int i;
    //ouverture boite de dialogue
    strcpy(ext, "*.wpp");
    CFileDialog DlgOpen(FALSE, NULL, ext, OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        //recupere nom de fichier
        fileName = DlgOpen.GetPathName();
        PtrfileName = fileName.GetBuffer(1);
        writeWindPatternsProject(PtrfileName, getWindPatternsProject());
    }

}


/***********/
/* quit */
/***********/

void quit(int control) {
    if (MessageBox(NULL, "Voulez-vous vraiment quit l'application ?", "Confirmation", MB_YESNO) == IDYES) {
        /*liberation espace memoire*/
        clearAxe(Axe3d);
        //clearAxe(Axe3dBal);
        clearAxe(AxePatron);
        clearAxe(AxePatronDXF);
        clearAxe(AxeRepDXF);
        clearAxe(AxePatronTextDXF);
        clearAxe(AxeMarginDXF);
        clearAxe(AxeCercleDXF);
        delete(F);
        delete(IntProfCent);
        delete(ExtProfCent);
        delete(IntProfBout);
        delete(ExtProfBout);
        exit(0);
    }
}

/********/
/* Info */


void Info(int /*control*/) {

}

void calcMagicPincesLayout() {
    int n = F->m_nbProfils;
    pincePercentage1 = new double[n];
    pincePercentage2 = new double[n];

    for (int _i = 1; _i < n - 1; _i++) {
        pincePercentage1[_i] = 0.0f;
        pincePercentage2[_i] = 0.0f;
    }
    if (n & 1) pincePercentage1[0] = 1.0f;
    else calcMagicPinceR(0, 1.0f, 1, &pincePercentage1[0]);
    for (int _i = 1; _i < n - 1; _i++) {
        calcMagicPinceR(_i, 2.0f - pincePercentage1[_i - 1], 1, &pincePercentage1[_i]);
    }
     printf("\n");
    if (n & 1) pincePercentage2[0] = 1.0f;
    else calcMagicPinceR(0, 1.0f, 2, &pincePercentage2[0]);
    for (int _i = 1; _i < n - 1; _i++) {
        calcMagicPinceR(_i, 2.0f - pincePercentage2[_i - 1], 2, &pincePercentage2[_i]);
    }
}

void GetFormProfile (int nerv, int face, Matrix** XProf, Matrix** YProf) {
    Matrix *XExt, *YExt, *ZExt, *XInt, *YInt, *ZInt;
    calcForm3D(F, 0, 0.0f,
        ExtProfCent, IntProfCent, ExtProfBout, IntProfBout,
        &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);
    int i = nerv;
    double LongNerv = F->m_pProfils[i]->m_fLength;
    double EpaiRel = F->m_pProfils[i]->m_fWidth;
    double m = F->m_pProfils[i]->m_fMorph;
    double coeffx = LongNerv/100.0f;
    double EpaiRelProfCent = EpaisseurRelative(ExtProfCent, IntProfCent);
    double EpaiRelProfBout = EpaisseurRelative(ExtProfBout, IntProfBout);
    double coeffyCent = LongNerv*EpaiRel/(EpaiRelProfCent*100.0f);
    double coeffyBout = LongNerv*EpaiRel/(EpaiRelProfBout*100.0f);
    double xp, yp;
    int n = 0, j = 0;
    if (face == 1) n = ExtProfCent->GetLignes(); else n = IntProfCent->GetLignes();
    *XProf = new Matrix(n, 1);
    *YProf = new Matrix(n, 1);

    if (face == 1) {
        for (j = 0; j < ExtProfCent -> GetLignes(); j++)
        {
                xp = ExtProfCent->Element(j,0)*coeffx;
                yp = ExtProfCent->Element(j,1)*coeffyCent*m
                        + ExtProfBout->Element(j,1)*coeffyBout*(1.0f-m);
                (*XProf)->SetElement(j, 0, xp);
                (*YProf)->SetElement(j, 0, yp);
        }
    } else {
        for (j = 0; j < IntProfCent -> GetLignes(); j++)
        {
                xp = IntProfCent->Element(j,0) * coeffx;
                yp = IntProfCent->Element(j,1) * coeffyCent * m
                        + IntProfBout->Element(j,1) * coeffyBout*(1.0f - m);
                (*XProf)->SetElement(j, 0, xp);
                (*YProf)->SetElement(j, 0, yp);
        }
    }
}

void Test(int control) {
    // printf ("\n TEST()");

	Matrix* Xd0 = new Matrix(4, 0);
	Matrix* Yd0 = new Matrix(4, 0);

	Xd0->SetElement (0, 0, 0.0);
	Yd0->SetElement (0, 0, 0.0);

	Xd0->SetElement (1, 0, 2.0);
	Yd0->SetElement (1, 0, 2.0);

	Xd0->SetElement (2, 0, 2.0);
	Yd0->SetElement (2, 0, 0.0);

	Xd0->SetElement (3, 0, 5.0);
	Yd0->SetElement (3, 0, 0.0);

	printf ("\n 0.0, 0.0: %f", calcPolyLength(Xd0, Yd0, 0.0, 0.0));
	printf ("\n 2.0, 2.0: %f", calcPolyLength(Xd0, Yd0, 2.0, 2.0));
	printf ("\n 2.0, 0.0: %f", calcPolyLength(Xd0, Yd0, 2.0, 0.0));
	printf ("\n 4.0, 0.0: %f", calcPolyLength(Xd0, Yd0, 4.0, 0.0));
	printf ("\n 5.0, 0.0: %f", calcPolyLength(Xd0, Yd0, 5.0, 0.0));

	// int res = calcRepPointByLength(Xd, Yd, l*coeff, &newx, &newy);
	int res=0;
	double x=0, y=0;
	res = calcRepPointByLength(Xd0, Yd0, 0.0, &x, &y);
	printf ("\n 0.0: res=%d x=%f y=%f", res, x, y);

	res = calcRepPointByLength(Xd0, Yd0, 2.82, &x, &y);
	printf ("\n 2.82: res=%d x=%f y=%f", res, x, y);

	res = calcRepPointByLength(Xd0, Yd0, 4.828427, &x, &y);
	printf ("\n 4.828427: res=%d x=%f y=%f", res, x, y);

	res = calcRepPointByLength(Xd0, Yd0, 6.82, &x, &y);
	printf ("\n 6.82: res=%d x=%f y=%f", res, x, y);

	res = calcRepPointByLength(Xd0, Yd0, 6.828427, &x, &y);
	printf ("\n 6.828427: res=%d x=%f y=%f", res, x, y);

/*    char filename[255];
    strcpy(filename, "bu.dxf");
    strcat (filename,".log");
    printf ("\n logfilename=[%s]", filename);
    CopyFile("ConfigDefaut.txt", filename, true); */

/*     int n = F -> m_nbProfils;

    Pince** PinceFunction1;
    Pince** PinceFunction2;
    PinceFunction1 = new Pince*[n];
    PinceFunction2 = new Pince*[n];*/


    //Matrix *FF1,*FF2, *RF;
//    double ampFfA, ampFfF, ampRfA, ampRfF;
    /*getPinceFunctions(2, 3, 1,
                        &FF1,&FF2, &ampFfA, &ampFfF,
                        &RF, &ampRfA, &ampRfF);*/
/*    Matrix* m = new Matrix (5, 2);
    for (int i = 0; i < 5; i++) {
        m -> SetElement(i, 0, i);
        m -> SetElement(i, 1, 4*i);
    }
    Matrix* p = new Matrix (2, 1);
    p->SetElement(0, 0, 1.5f);
    p->SetElement(1, 0, 3.5f);

    Matrix* mn;
    AddPointsToCourb(m, p, &mn);
    for (i = 0; i<mn->GetLignes();i++) {
        printf ("\n %d -> (%f, %f)", i, mn->Element(i,0),  mn->Element(i,1));
    }*/

    /*Matrix* addRepPts = getAddRepPts();
    for  (int _i = 0; _i < addRepPts->GetLignes();_i++){
        printf ("\n arp %d -> %f", _i, addRepPts->Element(_i, 0));
    }*/
    //writeWindPatternsProject("testproject.txt", getWindPatternsProject());
}


void SaveLayout() {
    int startNoNerv = 0, noNerv = 0, face = 1;
    Matrix * Xd[2], *Yd[2];//, *Xdp[2], *Ydp[2];
    Matrix * X[2], *Y[2], *Z[2], *P[2];//, *newP[2];
    int n = F->m_nbProfils;
    /*
     * int n = F->m_nbProfils;
     * if (n&1)  pincePercentage1[0]=1.0f; else calcMagicPinceR(0, 1.0f, 1, &pincePercentage1[0]);
        for (_i = 1; _i < n - 1; _i++) {
            calcMagicPinceR(_i, 2.0f - pincePercentage1[_i - 1], 1, &pincePercentage1[_i]);
        }
        printf("\n");
        if (n&1)  pincePercentage2[0]=1.0f; else calcMagicPinceR(0, 1.0f, 2, &pincePercentage2[0]);
        for (_i = 1; _i < n - 1; _i++) {
            calcMagicPinceR(_i, 2.0f - pincePercentage2[_i - 1], 2, &pincePercentage2[_i]);
        }
        calcPince(Xd[0], Yd[0], Xd[1], Yd[1],
                X[0], Y[0], Z[0], P[0],
                X[1], Y[1], Z[1], P[1],
                PosPinceBA[0], AmpPinceBA[0] * 1.0f, PosPinceBF[0], AmpPinceBF[0] * 1.0f, 1.5f,
                &Xdp[0], &Ydp[0], &newP[0], &Xdp[1], &Ydp[1], &newP[1]);
        delete(P[0]);
        delete(P[1]);
        P[0]=newP[0];
        P[1]=newP[1];
     */
    int size = 6 * n + 2 + 10;
    double *H, *W;
    TAxe **AxeP, **AxePD, **AxePTD, **AxeMD, **AxeCD, **AxeRepD;
    H = new double[size];
    W = new double[size];
    AxeP = new TAxe*[size];
    AxePD = new TAxe*[size];
    AxePTD = new TAxe*[size];
    AxeMD = new TAxe*[size];
    AxeCD = new TAxe*[size];
    AxeRepD = new TAxe*[size];
    face = 1;
    int *n1, *n2, *f1, *f2;
    int *isPince;
    bool *s1, *s2;
    double *p;
    f1 = new int[size];
    f2 = new int[size];
    n1 = new int[size];
    n2 = new int[size];
    s1 = new bool[size];
    s2 = new bool[size];
    p = new double[size];
    isPince = new int[size];
    //p2 = new double[size];
    int isave = 0, i = 0;
    for (face = 1; face <= 2; face++) {
        for (i = n - 1; i > 0; i--) {
            n1[isave] = i;
            n2[isave] = i - 1;
            s1[isave] = true;
            s2[isave] = true;
            f1[isave] = face;
            f2[isave] = face;
            if (face == 1) {
                p[isave] = 1.0f - pincePercentage1[i - 1];
            } else {
                p[isave] = 1.0f - pincePercentage2[i - 1];
            }
            isPince[isave] = 1;
            isave++;
        }
        if (!(n & 1)) {
            n1[isave] = -1;
            n2[isave] = 0;
            s1[isave] = false;
            s2[isave] = false;
            p[isave] = 1.0f;
            f1[isave] = face;
            f2[isave] = face;
            isPince[isave] = 1;
            isave++;
        }
        for (i = 0; i < n - 1; i++) {
            n1[isave] = i;
            n2[isave] = i + 1;
            s1[isave] = false;
            s2[isave] = false;
            if (face == 1) {
                p[isave] = pincePercentage1[i];
            } else {
                p[isave] = pincePercentage2[i];
            }
            f1[isave] = face;
            f2[isave] = face;
            isPince[isave] = 1;
            isave++;
        }
        if (face == 1) {
            // profile
            /*            for (i = n - 2; i > 0; i--) {
                            n1[isave] = i;
                            n2[isave] = i;
                            s1[isave] = true;
                            s2[isave] = true;
                            f1[isave] = 2;
                            f2[isave] = 1;
                            isPince[isave]=0;
                            isave++;
                        }*/
            /* if (!(n & 1)) {
                 n1[isave] = 0;
                 n2[isave] = 0;
                 s1[isave] = true;
                 s2[isave] = true;
                 p[isave] = 1.0f;
                 f1[isave] = 2;
                 f2[isave] = 1;
                 isPince[isave]=0;
                 isave++;
             }*/
            for (i = 0; i < n - 1; i++) {
                n1[isave] = i;
                n2[isave] = i;
                s1[isave] = false;
                s2[isave] = false;
                isPince[isave] = 0;
                f1[isave] = 2;
                f2[isave] = 1;
                isave++;
            }
        }
    }
    int q = isave;
    char text[100];
    //return;
    for (int t = 0; t < q; t++) {
        // printf ("\n t=%d ", t);
        //  printf ("\n n%d %d", n1[t], f1[t]);
        //		if (s1[t]) printf ("true"); else printf ("false");
        //  printf ("\n n%d %d", n2[t], f2[t]);
        //		if (s2[t]) printf ("true"); else printf ("false");
        calcPatron(getWindPatternsProject(),n1[t], s1[t], f1[t], f1[t], 0.0f, 100.0f,
                n2[t], s2[t], f2[t], f2[t], 0.0f, 100.0f,
                &Xd[0], &Yd[0], &Xd[1], &Yd[1],
                &X[0], &Y[0], &Z[0], &P[0],
                &X[1], &Y[1], &Z[1], &P[1]);
        //	printf ("\n after calcPatron");
        if (isPince[t]) {
            Pince* pince=new Pince();

/*            calcPince(getWindPatternsProject(),Xd[0], Yd[0], Xd[1], Yd[1],
                    X[0], Y[0], Z[0], P[0],
                    X[1], Y[1], Z[1], P[1],
                    PosPinceBA[0], AmpPinceBA[0] * p[t], PosPinceBF[0], AmpPinceBF[0] * (1.0f - p[t]), modePinces,
                    &Xdp[0], &Ydp[0], &newP[0], &Xdp[1], &Ydp[1], &newP[1]);*/
            //        printf ("\n after calcPince");
            delete(P[0]);
            delete(P[1]);
//            P[0] = newP[0];
//            P[1] = newP[1];
        }
        sprintf(text, "N%dF%dN%dF%d", n1[t], f1[t], n2[t], f2[t]);
/*
        if (isPince[t])
            GenerateCourbe(getWindPatternsProject(),Xdp[0], Ydp[0], P[0], n1[t],  f1[t],
                Xdp[1], Ydp[1], P[1], n2[t], f2[t], text,
                &AxeP[t], &AxePD[t], &AxePTD[t], &AxeMD[t], &AxeCD[t], &AxeRepD[t]);
        else
            GenerateCourbe(getWindPatternsProject(), Xd[0], Yd[0], P[0], n1[t], f1[t],
                Xd[1], Yd[1], P[1], n2[t], f2[t], text,
                &AxeP[t], &AxePD[t], &AxePTD[t], &AxeMD[t], &AxeCD[t], &AxeRepD[t]);
*/
//        calcMaxWH(Xdp[0], Ydp[0], Xdp[1], Ydp[1], &W[t], &H[t]);


        /*delete (X[0]);
        delete (Y[0]);
        delete (Z[0]);
        delete (P[0]);
        delete (X[1]);
        delete (Y[1]);
        delete (Z[1]);
        delete (P[1]);
        delete (Xd[0]);
        delete (Yd[0]);
        delete (Xd[1]);
        delete (Yd[1]);*/
    }

    CString fileName;
    LPTSTR PtrfileName;
    char ext[10];
    strcpy(ext, "*.dxf");
    CFileDialog DlgOpen(FALSE, NULL, ext, OFN_OVERWRITEPROMPT, NULL, NULL);
    if (DlgOpen.DoModal() == IDOK) {
        fileName = DlgOpen.GetPathName();
        PtrfileName = fileName.GetBuffer(1);
        writeManyFichierPolyDXF(PtrfileName, n, q, AxePD, AxeMD, 1, AxeRepD, 0, AxeCD, 1, AxePTD, W, H);
    }
    //TAxe *AxeP[n], *AxePD[n], *AxePTD[n], *AxeMD[n], *AxeCD[n], *AxeRepD[n];
    for (i = 0; i < q; i++) {
        //printf("\n cleare i=%d", i);
        clearCourbesAxe(AxeP[i]);
        clearCourbesAxe(AxePD[i]);
        clearCourbesAxe(AxePTD[i]);
        clearCourbesAxe(AxeMD[i]);
        clearCourbesAxe(AxeCD[i]);
        clearCourbesAxe(AxeRepD[i]);
    }
    //printf("\n end of clear...");
}


void calcMagicPinceR(int noNerv1, double perc, int face, double* percNew) {
    Matrix * Xd1[2], *Yd1[2], *Xd1p[2], *Yd1p[2];
    Matrix * Xd2[2], *Yd2[2], *Xd2p[2], *Yd2p[2];
    Matrix * X1[2], *Y1[2], *Z1[2], *P1[2], *newP1[2];
    Matrix * X2[2], *Y2[2], *Z2[2], *P2[2], *newP2[2];
    bool showPoints = false;
    //if ( DEBUG ) printf("\ncalcMagicPinceR(%d, %f, %d)", noNerv1, perc, face);
    int i = 0;
    calcPatron(getWindPatternsProject(), noNerv1 - 1, false, face, face, Deb[0], Fin[0],
            noNerv1, false, face, face, Deb[1], Fin[1],
            &Xd1[0], &Yd1[0], &Xd1[1], &Yd1[1],
            &X1[0], &Y1[0], &Z1[0], &P1[0],
            &X1[1], &Y1[1], &Z1[1], &P1[1]);
    Matrix *lenm1 = Longueur(Xd1[1], Yd1[1]);
    double len1 = lenm1->Element(lenm1->GetLignes() - 1, 0);
    double width, w1, w2;
    calcWidth(Xd1[0], Yd1[0], Xd1[1], Yd1[1], &w1);
    int jt = 0;
    if (showPoints) {
        for (jt = 0; jt < Xd1[0]->GetLignes(); jt = jt + 6) {
            printf("\n0 %d(%7.4f, %7.4f)", jt, Xd1[0]->Element(jt, 0), Yd1[0]->Element(jt, 0));
        }
        for (jt = 0; jt < Xd1[0]->GetLignes(); jt = jt + 6) {
            printf("\n1 %d(%7.4f, %7.4f)", jt, Xd1[1]->Element(jt, 0), Yd1[1]->Element(jt, 0));
        }
    }
    calcWidth(Xd1[0], Yd1[0], Xd1[1], Yd1[1], &width);
    calcPatron(getWindPatternsProject(),noNerv1, false, face, face, Deb[0], Fin[0],
            noNerv1 + 1, false, face, face, Deb[1], Fin[1],
            &Xd2[0], &Yd2[0], &Xd2[1], &Yd2[1],
            &X2[0], &Y2[0], &Z2[0], &P2[0],
            &X2[1], &Y2[1], &Z2[1], &P2[1]);
    Matrix *lenm2 = Longueur(Xd2[0], Yd2[0]);
    double len2 = lenm2->Element(lenm2->GetLignes() - 1, 0);
    double minDelta = 1000000.0f, minapAI1 = 0.0f, minapFI1 = 0.0f, minp2 = 1000000.0f;
    double ppA = PosPinceBA[1];
    double apA = AmpPinceBA[1];
    double ppF = PosPinceBF[1];
    double apF = AmpPinceBF[1];
    double startPerc = 0.0, endPerc = 2.0, hPerc = 0.001, l1p1, l2p0, delta;
    Matrix *long1p1, *long2p0;
    double mode1 = 1.5;
    double mode2 = 1.5;
    calcWidth(Xd2[0], Yd2[0], Xd2[1], Yd2[1], &w2);
    if (showPoints) {
        for (jt = 0; jt < Xd2[0]->GetLignes(); jt = jt + 6) {
            printf("\n0 %d(%7.4f, %7.4f)", jt, Xd2[0]->Element(jt, 0), Yd2[0]->Element(jt, 0));
        }
        for (jt = 0; jt < Xd2[0]->GetLignes(); jt = jt + 6) {
            printf("\n1 %d(%7.4f, %7.4f)", jt, Xd2[1]->Element(jt, 0), Yd2[1]->Element(jt, 0));
        }
    }

    calcPinceAlone(getWindPatternsProject(),Xd1[1], Yd1[1], Xd1[0], Yd1[0], X1[1], Y1[1], Z1[1], P1[1], X1[0], Y1[0], Z1[0], P1[0],
            ppA, apA*perc, ppF, apF*perc, mode1, &Xd1p[1], &Yd1p[1], &newP1[1]);
    long1p1 = Longueur(Xd1p[1], Yd1p[1]);
    l1p1 = long1p1->Element(long1p1->GetLignes() - 1, 0);

    for (double apAI1 = startPerc; apAI1 <= endPerc; apAI1 += hPerc) {
        calcPinceAlone(getWindPatternsProject(),Xd2[0], Yd2[0], Xd2[1], Yd2[1], X2[0], Y2[0], Z2[0], P2[0], X2[1], Y2[1], Z2[1], P2[1],
                ppA, apA*apAI1, ppF, apF*apAI1, mode2, &Xd2p[0], &Yd2p[0], &newP2[0]);
        long2p0 = Longueur(Xd2p[0], Yd2p[0]);
        l2p0 = long2p0->Element(long2p0->GetLignes() - 1, 0);
        delta = fabs(l1p1 - l2p0);
        if (delta < minDelta) {
            minDelta = delta;
            minapAI1 = apAI1;
            minp2 = l2p0;
        }
        delete (long2p0);
        delete (Xd2p[0]);
        delete (Yd2p[0]);
        delete (newP2[0]);
    }
    delete (long1p1);
    delete (Xd1p[1]);
    delete (Yd1p[1]);
    printf("\n%2df%d %6.4f (%5.3f) %6.4f(%7.4f)[%4.3f] M%7.4f(%7.4f) w%5.3f/%5.3f",
            noNerv1, face, (minDelta * 1000), minapAI1, l1p1, len1, perc, minp2, len2, w1, w2);
    *percNew = minapAI1;
    if (minDelta < (tochnostLayout2 * 0.001)) {
        printf(" OK:)");
    } else {
        printf(" BAD!");
    }
    delete (newP1[1]);
    delete (X1[0]);
    delete (Y1[0]);
    delete (Z1[0]);
    delete (P1[0]);
    delete (X1[1]);
    delete (Y1[1]);
    delete (Z1[1]);
    delete (P1[1]);
    delete (X2[0]);
    delete (Y2[0]);
    delete (Z2[0]);
    delete (P2[0]);
    delete (X2[1]);
    delete (Y2[1]);
    delete (Z2[1]);
    delete (P2[1]);
    delete (Xd1[0]);
    delete (Yd1[0]);
    delete (Xd2[0]);
    delete (Yd2[0]);
    delete (Xd1[1]);
    delete (Yd1[1]);
    delete (Xd2[1]);
    delete (Yd2[1]);
}

/***********************/
/* calcPinces        */
/***********************/

void calcPinces(Matrix *Xd1, Matrix *Yd1,
        double PosPinceBA1, double AmpPinceBA1, double PosPinceBF1, double AmpPinceBF1,
        Matrix *Xd2, Matrix *Yd2,
        double PosPinceBA2, double AmpPinceBA2, double PosPinceBF2, double AmpPinceBF2, double modePincesT,
        Matrix **Xdp1, Matrix **Ydp1,
        Matrix **Xdp2, Matrix **Ydp2) {
    double PosPinceBAT[2] = {PosPinceBA1, PosPinceBA2};
    double PosPinceBFT[2] = {PosPinceBF1, PosPinceBF2};
    double AmpPinceBAT[2] = {AmpPinceBA1, AmpPinceBA2};
    double AmpPinceBFT[2] = {AmpPinceBF1, AmpPinceBF2};

    int n, i;
    double longPince, ampPince, amp, larg, longCoteMax;
    Matrix * longCote[2], *newXd[2], *newYd[2];
    int iPince, iBidon, i0, i1;

    Matrix * Xd[2], *Yd[2];
    Xd[0] = Xd1;
    Yd[0] = Yd1;
    Xd[1] = Xd2;
    Yd[1] = Yd2;

    longCote[0] = Longueur(Xd[0], Yd[0]);
    longCote[1] = Longueur(Xd[1], Yd[1]);

    newXd[0] = new Matrix(Xd[0]->GetLignes(), 1);
    newYd[0] = new Matrix(Xd[0]->GetLignes(), 1);
    newXd[1] = new Matrix(Xd[1]->GetLignes(), 1);
    newYd[1] = new Matrix(Xd[1]->GetLignes(), 1);

    for (i = 0; i < Xd[0]->GetLignes(); i++) {
        newXd[0]->SetElement(i, 0, Xd[0]->Element(i, 0));
        newXd[1]->SetElement(i, 0, Xd[1]->Element(i, 0));
        newYd[0]->SetElement(i, 0, Yd[0]->Element(i, 0));
        newYd[1]->SetElement(i, 0, Yd[1]->Element(i, 0));
    }

    for (n = 0; n < 2; n++) {
        if (n == 0) {
            i0 = 0;
            i1 = 1;
        } else {
            i0 = 1;
            i1 = 0;
        }
        if (PosPinceBAT[i0] > 0.0f) {
            longCoteMax = longCote[i0]->Element(longCote[i0]->GetLignes() - 1, 0);
            longPince = PosPinceBAT[i0] * longCoteMax / 100.0f;
            Ind(longPince, longCote[i0], &iPince, &iBidon);
            larg = sqrt(sqr(Xd[i0]->Element(0, 0) - Xd[i1]->Element(0, 0)) + sqr(Yd[i0]->Element(0, 0) - Yd[i1]->Element(0, 0)));
            ampPince = AmpPinceBAT[i0] * larg / 100.0f;
            for (i = 0; i < iPince; i++) {
                amp = ampPince * pow(1.0f - longCote[i0]->Element(i, 0) / longPince, modePincesT);
                larg = sqrt(sqr(Xd[i0]->Element(i, 0) - Xd[i1]->Element(i, 0)) + sqr(Yd[i0]->Element(i, 0) - Yd[i1]->Element(i, 0)));
                newXd[i0]->SetElement(i, 0, Xd[i0]->Element(i, 0) + (Xd[i1]->Element(i, 0) - Xd[i0]->Element(i, 0)) * amp / larg);
                newYd[i0]->SetElement(i, 0, Yd[i0]->Element(i, 0) + (Yd[i1]->Element(i, 0) - Yd[i0]->Element(i, 0)) * amp / larg);
            }
        }
    }
    for (n = 0; n < 2; n++) {
        if (n == 0) {
            i0 = 0;
            i1 = 1;
        } else {
            i0 = 1;
            i1 = 0;
        }
        if (PosPinceBFT[i0] > 0.0f) {
            longCoteMax = longCote[i0]->Element(longCote[i0]->GetLignes() - 1, 0);
            longPince = PosPinceBFT[i0] * longCoteMax / 100.0f;
            Ind(longCoteMax - longPince, longCote[i0], &iPince, &iBidon);
            larg = sqrt(sqr(Xd[i0]->Element(0, 0) - Xd[i1]->Element(0, 0)) + sqr(Yd[i0]->Element(0, 0) - Yd[i1]->Element(0, 0)));
            ampPince = AmpPinceBFT[i0] * larg / 100.0f;
            for (i = iPince; i < Xd[i0]->GetLignes(); i++) {
                amp = ampPince * pow(1.0f - (longCoteMax - longCote[i0]->Element(i, 0)) / longPince, modePincesT);
                larg = sqrt(sqr(Xd[i0]->Element(i, 0) - Xd[i1]->Element(i, 0)) + sqr(Yd[i0]->Element(i, 0) - Yd[i1]->Element(i, 0)));
                newXd[i0]->SetElement(i, 0, Xd[i0]->Element(i, 0) + (Xd[i1]->Element(i, 0) - Xd[i0]->Element(i, 0)) * amp / larg);
                newYd[i0]->SetElement(i, 0, Yd[i0]->Element(i, 0) + (Yd[i1]->Element(i, 0) - Yd[i0]->Element(i, 0)) * amp / larg);
            }
        }
    }
    *Xdp1 = CloneMat(newXd[0]);
    *Ydp1 = CloneMat(newYd[0]);
    *Xdp2 = CloneMat(newXd[1]);
    *Ydp2 = CloneMat(newYd[1]);
}

/***********************/
/* calcVue3dEtPatron */
/***********************/

void calcVue3dEtPatron(void)
{
    //printf ("\n calcVue3dEtPatron");
    Matrix *XExt, *YExt, *ZExt;
    Matrix *XInt, *YInt, *ZInt;
    Matrix *XExtBal, *YExtBal, *ZExtBal;
    Matrix *XIntBal, *YIntBal, *ZIntBal;

    Matrix *PtsSuspentes; //pour recuperer la position 3D des pts de suspentage
    
    int ajCoAr = 0, ajCoMaAr = 0, ajCoAv = 0, ajCoMaAv = 0;
    int ajCo1[2] = {0, 0};
    int ajCo2[2] = {0, 0};
    //pour developpe patron ...
    int i, j;
    bool symetrique[2] = {false, false};
    int nerv[2];
    Matrix * X[2], *Y[2], *Z[2];
    Matrix * P[2], *rP[2]; // % de la corde du profil
    Matrix * Xd[2], *Yd[2], *_Xd[2],  *_Yd[2], *rXd[2], *rYd[2];
    int iDeb, iFin, iCol;
    Matrix *xExtProf, *xIntProf;
    Matrix *iAv, *iAp, *Res;
    Courbe * CourbPatron[2], *CourbPatronDXF[2], *CourbPatronBack[2], *CourbMarge[2], *CourbMargeBack[2], *CourbMargeDXF[2];
    Courbe *CourbAv, *CourbAvBack, *CourbAr, *CourbArDXF, *CourbMargeAv, *CourbMargeAvBack, *CourbMargeAr, *CourbCoin, *CourbMargeArDXF, *CourbCoin1[2], *CourbCoin1Back[2], *CourbCoin2[2], *CourbCoin2Back[2], *CourbRep, *CourbRepDXF;
    Courbe *CourbPin[2];
    Matrix *distance;
    int n; //, iRep;
    double xText, yText;
    char texte[100], texteExtInt[3] = "EI";
    //pour calc position suspentes
    //Matrix *interpXSuspente, *interpYSuspente, *newP[2];
    //double posSuspente, 
    double xSuspente[2][205], ySuspente[2][205];
    TMesh *MeshPatron;
    Matrix *newXd[2],  *newYd[2], *rnewXd[2],  *rnewYd[2];//, *Xdp, *Ydp;
    Courbe *CourbCercle, *CourbCercleDXF;
    clearCourbesAxe(Axe3d);
    clearMeshsAxe(Axe3d);
    //clearCourbesAxe(Axe3dBal);
    //clearMeshsAxe(Axe3dBal);

    calcForm3D(F, 0, 0.0f,
            ExtProfCent, IntProfCent, ExtProfBout, IntProfBout,
            &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);
    //addForm3D(Axe3d, XExt, YExt, ZExt, XInt, YInt, ZInt, ViewFace, ViewSymetrique);
	int showBal = 1;
	if (showBal) {
		//printf ("\n if showBal");
		//printf ("\n readBallonementFile(...");
		//Ballonement *bal = readBallonementFromFile("bal.txt");
		//printf ("\n calcForm3DBallonement(...");
		calcForm3DBallonement(getWindPatternsProject(), F,  0, 0.0f,
				ExtProfCent, IntProfCent, ExtProfBout, IntProfBout,
				&XExtBal, &YExtBal, &ZExtBal, &XIntBal, &YIntBal, &ZIntBal);
		//printf ("\n addForm3D(...");
		addForm3D(Axe3d, XExtBal, YExtBal, ZExtBal, XIntBal, YIntBal, ZIntBal, ViewFace, ViewSymetrique);
		//printf ("\n !!!...");
	}

    /*adde points de suspentage*/
    if (ViewPtsSuspentes) {
        addPtsSuspentage(Axe3d, F, IntProfCent, XExt, YExt, ZExt, XInt, YInt, ZInt,
                ViewSymetrique, &PtsSuspentes);
    }
    clearCourbesAxe(AxePatron);
    clearCourbesAxe(AxeRepDXF);
    clearCourbesAxe(AxePatronDXF);
    clearCourbesAxe(AxePatronTextDXF);
    clearCourbesAxe(AxeMarginDXF);
    clearCourbesAxe(AxeCercleDXF);

    xExtProf = new Matrix(ExtProfCent->GetLignes(), 1);
    xIntProf = new Matrix(IntProfCent->GetLignes(), 1);

    for (i = 0; i < ExtProfCent->GetLignes(); i++) {
        xExtProf->SetElement(i, 0, ExtProfCent->Element(i, 0));
    //    printf ("\n ExtProf[%d]=%f %f", i, ExtProfCent->Element(i, 0), ExtProfCent->Element(i,1));
    }
    for (i = 0; i < IntProfCent->GetLignes(); i++) {
        xIntProf->SetElement(i, 0, IntProfCent->Element(i, 0));
      //  printf ("\n IntProf[%d]=%f %f", i, xIntProf->Element(i, 0),  IntProfCent->Element(i,1));
    }
    
    double myDeb = Deb[0];
    double myFin = Fin[0];

    if (GoPince) {
        if (FaceDeb[0] == FaceFin[0]) {
            if (Deb[0] < Deb[1]) myDeb = Deb[0]; else myDeb = Deb[1];
        } else {
            if (Deb[0] > Deb[1]) myDeb = Deb[0]; else myDeb = Deb[1];
        }
    }

    for (i = 0; i < 2; i++) {
        nerv[i] = NoNerv[i];
        if (NoNerv[i] == -1) {
            nerv[i] = 0;
            symetrique[i] = true;
        }
    }
    if (isKlapan) {
        calcPatronKlapan(getWindPatternsProject(), nerv[0], symetrique[0], nerv[1], symetrique[1], PosKlapanInt, PosKlapanFin,
                    &Xd[0], &Yd[0], &Xd[1], &Yd[1],
                    &X[0], &Y[0], &Z[0], &P[0],
                    &X[1], &Y[1], &Z[1], &P[1]);
        
    } else
        if (GoPince)
            calcPatron(getWindPatternsProject(), nerv[0], symetrique[0], FaceDeb[0], FaceFin[0], myDeb, myFin,
                        nerv[1], symetrique[1], FaceDeb[1], FaceFin[1], myDeb, myFin,
                        &Xd[0], &Yd[0], &Xd[1], &Yd[1],
                        &X[0], &Y[0], &Z[0], &P[0],
                        &X[1], &Y[1], &Z[1], &P[1]);
        else {
            if (isPosNerv[0] || isPosNerv[1]) {
                if (!isPosNerv[0])
                    if (FaceDeb[0]==1) PosNerv[0]=100.0f; else PosNerv[0]=0.0f;
                if (!isPosNerv[1])
                    if (FaceDeb[1]==1) PosNerv[1]=100.0f; else PosNerv[1]=0.0f;

                calcPatronPosNerv(getWindPatternsProject(), nerv[0], symetrique[0], FaceDeb[0], FaceFin[0], Deb[0], Fin[0],
                            nerv[1], symetrique[1], FaceDeb[1], FaceFin[1], Deb[1], Fin[1],
                            &Xd[0], &Yd[0], &Xd[1], &Yd[1],
                            &X[0], &Y[0], &Z[0], &P[0],
                            &X[1], &Y[1], &Z[1], &P[1]);
            }
            else
            {
            calcPatron(getWindPatternsProject(), nerv[0], symetrique[0], FaceDeb[0], FaceFin[0], Deb[0], Fin[0],
                        nerv[1], symetrique[1], FaceDeb[1], FaceFin[1], Deb[1], Fin[1],
                        &Xd[0], &Yd[0], &Xd[1], &Yd[1],
                        &X[0], &Y[0], &Z[0], &P[0],
                        &X[1], &Y[1], &Z[1], &P[1]);
            }
        }


    /*
        //test face int/ext
        if (FaceDeb[i] == 1 && FaceFin[i] == 1) //debut et fin en extrados
        {
            if (Deb[i] > Fin[i]) //inversion Deb<->Fin
            {
                Ind(Fin[i], xExtProf, &iDeb, &iCol);
                Ind(Deb[i], xExtProf, &iFin, &iCol);
            } else {
                Ind(Deb[i], xExtProf, &iDeb, &iCol);
                Ind(Fin[i], xExtProf, &iFin, &iCol);
            }
            X[i] = new Matrix(iFin - iDeb + 1, 1);
            Y[i] = new Matrix(iFin - iDeb + 1, 1);
            Z[i] = new Matrix(iFin - iDeb + 1, 1);
            P[i] = new Matrix(iFin - iDeb + 1, 1);
            for (j = 0; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Y[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Z[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], iDeb + j));
            for (j = 0; j < P[i]->GetLignes(); j++) {
                P[i]->SetElement(j, 0, xExtProf->Element(iDeb + j, 0));
            }
        } else if (FaceDeb[i] == 1) //debut en extrados et fin en intrados
        {
            Ind(Deb[i], xExtProf, &iDeb, &iCol);
            Ind(Fin[i], xIntProf, &iFin, &iCol);
            X[i] = new Matrix(iFin + iDeb - 1, 1);
            Y[i] = new Matrix(iFin + iDeb - 1, 1);
            Z[i] = new Matrix(iFin + iDeb - 1, 1);
            P[i] = new Matrix(iFin + iDeb - 1, 1);
            for (j = 0; j <= iDeb; j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) P[i]->SetElement(j, 0, xExtProf->Element(iDeb - j, 0));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xIntProf->Element(j - iDeb, 0));
        } else if (FaceFin[i] == 1) //debut en intrados et fin en extrados
        {
            Ind(Deb[i], xIntProf, &iDeb, &iCol);
            Ind(Fin[i], xExtProf, &iFin, &iCol);
            X[i] = new Matrix(iFin + iDeb + 1, 1);
            Y[i] = new Matrix(iFin + iDeb + 1, 1);
            Z[i] = new Matrix(iFin + iDeb + 1, 1);
            P[i] = new Matrix(iFin + iDeb + 1, 1);
            for (j = 0; j <= iDeb; j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], iDeb - j));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZExt->Element(nerv[i], j - iDeb));
            for (j = 0; j <= iDeb; j++) P[i]->SetElement(j, 0, xIntProf->Element(iDeb - j, 0));
            for (j = iDeb + 1; j < X[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xExtProf->Element(j - iDeb, 0));
            //                printf ("\n r2 3");
        } else //debut et fin en intrados
        {
            if (Deb[i] > Fin[i]) //inversion Deb<->Fin
            {
                Ind(Fin[i], xIntProf, &iDeb, &iCol);
                Ind(Deb[i], xIntProf, &iFin, &iCol);
            } else {
                Ind(Deb[i], xIntProf, &iDeb, &iCol);
                Ind(Fin[i], xIntProf, &iFin, &iCol);
            }
            X[i] = new Matrix(iFin - iDeb + 1, 1);
            Y[i] = new Matrix(iFin - iDeb + 1, 1);
            Z[i] = new Matrix(iFin - iDeb + 1, 1);
            P[i] = new Matrix(iFin - iDeb + 1, 1);
            for (j = 0; j < X[i]->GetLignes(); j++) X[i]->SetElement(j, 0, XInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Y[i]->GetLignes(); j++) Y[i]->SetElement(j, 0, YInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < Z[i]->GetLignes(); j++) Z[i]->SetElement(j, 0, ZInt->Element(nerv[i], iDeb + j));
            for (j = 0; j < P[i]->GetLignes(); j++) P[i]->SetElement(j, 0, xIntProf->Element(iDeb + j, 0));
        }

    }
*/
/*
    if (X[0]->GetLignes() > X[1]->GetLignes()) {
        iAv = LinSpace(0.0f, (double) X[1]->GetLignes() - 1, X[1]->GetLignes());
        iAp = LinSpace(0.0f, (double) X[1]->GetLignes() - 1, X[0]->GetLignes());
        Res = new Matrix(X[0]->GetLignes(), 1);
        InterpLinMat(iAv, X[1], iAp, Res);
        delete(X[1]);
        X[1] = Res;
        Res = new Matrix(X[0]->GetLignes(), 1);
        InterpLinMat(iAv, Y[1], iAp, Res);
        delete(Y[1]);
        Y[1] = Res;
        Res = new Matrix(X[0]->GetLignes(), 1);
        InterpLinMat(iAv, Z[1], iAp, Res);
        delete(Z[1]);
        Z[1] = Res;
        Res = new Matrix(X[0]->GetLignes(), 1);
        InterpLinMat(iAv, P[1], iAp, Res);
        delete(P[1]);
        P[1] = Res;
        delete(iAv);
        delete(iAp);
    } else if (X[1]->GetLignes() > X[0]->GetLignes()) {
        iAv = LinSpace(0.0f, (double) X[0]->GetLignes() - 1, X[0]->GetLignes());
        iAp = LinSpace(0.0f, (double) X[0]->GetLignes() - 1, X[1]->GetLignes());
        Res = new Matrix(X[1]->GetLignes(), 1);
        InterpLinMat(iAv, X[0], iAp, Res);
        delete(X[0]);
        X[0] = Res;
        Res = new Matrix(X[1]->GetLignes(), 1);
        InterpLinMat(iAv, Y[0], iAp, Res);
        delete(Y[0]);
        Y[0] = Res;
        Res = new Matrix(X[1]->GetLignes(), 1);
        InterpLinMat(iAv, Z[0], iAp, Res);
        delete(Z[0]);
        Z[0] = Res;
        Res = new Matrix(X[1]->GetLignes(), 1);
        InterpLinMat(iAv, P[0], iAp, Res);
        delete(P[0]);
        P[0] = Res;
        delete(iAv);
        delete(iAp);
    }

    for (i = 0; i < 2; i++) {
        if (symetrique[i] == true) {
            for (j = 0; j < X[i]->GetLignes(); j++)
                X[i]->MultiplyElement(j, 0, -1.0f);
        }
    }
    
    if (DEBUG) printf ("\n gocalcDevelope()");
    calcDeveloppe(
            X[0], Y[0], Z[0], X[1], Y[1], Z[1],
            &Xd[0], &Yd[0], &Xd[1], &Yd[1]); */
    
    //printf ("\n P[0]->GetLignes()=%d", P[0]->GetLignes());
    //for (i=0; i<P[0]->GetLignes();i++) printf ("\nP[0]->Element(%d,0)=%f",i,P[0]->Element(i, 0));
    //printf ("\n P[1]->GetLignes()=%d", P[1]->GetLignes());
    //for (i=0; i<P[1]->GetLignes();i++) printf ("\nP[1]->Element(%d,0)=%f",i,P[1]->Element(i, 0));


        //printf ("\nXd[0] %d", Xd[0]->GetLignes());
        //Xd[0]->print(0);
        //printf ("\nYd[0] %d", Yd[0]->GetLignes());
        //Yd[0]->print(0);
        //printf ("\nXd[1] %d", Xd[1]->GetLignes());
        //Xd[1]->print(0);
        //printf ("\nYd[1] %d", Yd[1]->GetLignes());
        //Yd[1]->print(0);
  

  //init tableau pour calc des pinces !!!
    /*Matrix* long1 = Longueur(Xd[0], Yd[0]);
    double l1 = long1->Element(long1->GetLignes() - 1, 0);
    Matrix* long2 = Longueur(Xd[1], Yd[1]);
    double l2 = long2->Element(long2->GetLignes() - 1, 0);
    printf("\n(%4.2f)L1=%9.7f L2=%9.7f", AmpPinceBA[0], l1, l2);
    delete(long1);
    delete(long2);*/
    //printf ("\n before CloneMat...()");
    newXd[0] = CloneMat(Xd[0]);
    newYd[0] = CloneMat(Yd[0]);
    newXd[1] = CloneMat(Xd[1]);
    newYd[1] = CloneMat(Yd[1]);
    //printf ("\n begore goPince()");
    bool wasPince = false;
    if ((GoPince == 1)
            && ((nerv[0] != nerv[1]) || ((nerv[0] == nerv[1])&&(symetrique[0] != symetrique[1])  ))
                && (FaceFin[0] == FaceFin[1])  && (FaceDeb[0] == FaceDeb[1])) {
                //printf ("\n go go go Pince ()");
                wasPince = true;

                //printf ("\nbefore GetPincePlus() myDeb=%f", myDeb);
                Pince* pince = getPincePlus (getWindPatternsProject(), NoNerv[0], NoNerv[1], myDeb, myFin, myDeb, myFin, FaceDeb[0], FaceFin[0]);
                //printf ("\nafter GetPincePlus()");

                //printf ("\n function1 L=%d", pince->function1->GetLignes());
                /*for (int i = 0;i<pince->function1->GetLignes();i++) {
                    printf ("\nfunction1 %d -> %f", i, pince->function1->Element(i,0));
                }*/
                //printf ("\n function2 L=%d", pince->function2->GetLignes());
                /*for (i = 0;i<pince->function2->GetLignes();i++) {
                    printf ("\nfunction2 %d -> %f", i, pince->function2->Element(i,0));
                }*/
                    
                pince -> debug = false;
                pince -> SetDiffAmps(pince->AmpA, pince->AmpF, pince->AmpA, pince->AmpF);
                if (pince -> glue) pince -> SetDiffAmps0(pince->AmpA0, pince->AmpA0);
                
                //printf ("\nbefore calcPincePlusNew");
                calcPincePlusNew(Xd[0], Yd[0], Xd[1], Yd[1],  pince, &newXd[0], &newYd[0], &newXd[1], &newYd[1]);
                //printf ("\nafter calcPincePlusNew");

                if ((myDeb != Deb[0])||(myDeb != Deb[1])) {
                    // rezem Xd, newXd
                    rXd[0] = GetFunctionSrezDeb(P[0], Xd[0], Deb[0]);
                    rYd[0] = GetFunctionSrezDeb(P[0], Yd[0], Deb[0]);
                    delete(Xd[0]);
                    Xd[0]=rXd[0];
                    delete(Yd[0]);
                    Yd[0]=rYd[0];
                    rnewXd[0] = GetFunctionSrezDeb(P[0], newXd[0], Deb[0]);
                    rnewYd[0] = GetFunctionSrezDeb(P[0], newYd[0], Deb[0]);
                    delete(newXd[0]);
                    newXd[0]=rnewXd[0];
                    delete(newYd[0]);
                    newYd[0]=rnewYd[0];

                    rXd[1] = GetFunctionSrezDeb(P[1], Xd[1], Deb[1]);
                    rYd[1] = GetFunctionSrezDeb(P[1], Yd[1], Deb[1]);
                    delete(Xd[1]);
                    Xd[1]=rXd[1];
                    delete(Yd[1]);
                    Yd[1]=rYd[1];

                    rnewXd[1] = GetFunctionSrezDeb(P[1], newXd[1], Deb[1]);
                    rnewYd[1] = GetFunctionSrezDeb(P[1], newYd[1], Deb[1]);
                    delete(newXd[1]);
                    newXd[1]=rnewXd[1];
                    delete(newYd[1]);
                    newYd[1]=rnewYd[1];

                    rP[0] = GetFunctionSrezDeb(P[0], P[0], Deb[0]);
                    delete(P[0]);
                    P[0]=rP[0];
                    rP[1] = GetFunctionSrezDeb(P[1], P[1], Deb[1]);
                    delete(P[1]);
                    P[1]=rP[1];
                }

                _Xd[0]=Xd[0];
                _Xd[1]=Xd[1];
                _Yd[0]=Yd[0];
                _Yd[1]=Yd[1];

                //delete(Xd[0]);
                //delete(Xd[1]);
                //delete(Yd[0]);
                //delete(Yd[1]);
                //delete(P[0]);
                //delete(P[1]);
                //printf ("\nnewXd[0] %d", newXd[0]->GetLignes());
                //newXd[0]->print(0);
                //printf ("\nnewYd[0] %d", newYd[0]->GetLignes());
                //newYd[0]->print(0);

                //printf ("\nnewXd[1] %d", newXd[1]->GetLignes());
                //newXd[1]->print(0);
                //printf ("\nnewYd[1] %d", newYd[1]->GetLignes());
                //newYd[1]->print(0);

                Xd[0] = newXd[0];
                Xd[1] = newXd[1];
                Yd[0] = newYd[0];
                Yd[1] = newYd[1];
                //P[0] = newP[0];
                //P[1] = newP[1];
    }


    if (coeffMult != 1.0f) {
        calcPatronWithCoeff(Xd[0], Yd[0], Xd[1], Yd[1], coeffMult, &newXd[0], &newYd[0], &newXd[1], &newYd[1]);
        delete (Xd[0]);
        delete (Yd[0]);
        delete(Xd[1]);
        delete (Yd[1]);
        Xd[0] = newXd[0];
        Yd[0] = newYd[0];
        Xd[1] = newXd[1];
        Yd[1] = newYd[1];
    }

/*
    printf ("\nXd[0]");
    Xd[0]->print(0);
    printf ("\nYd[0]");
    Yd[0]->print(0);
    printf ("\nXd[1]");
    Xd[1]->print(0);
    printf ("\nYd[1]");
    Yd[1]->print(0);
*/
/*    long1 = Longueur(Xd[0], Yd[0]);
    l1 = long1->Element(long1->GetLignes() - 1, 0);
    long2 = Longueur(Xd[1], Yd[1]);
    l2 = long2->Element(long2->GetLignes() - 1, 0);
    printf(" LP1=%9.7f LP2=%9.7f", l1, l2);*/
    for (i = 0; i < 2; i++) {
        CourbPatron[i] = new Courbe("Patron");
        CourbPatronDXF[i] = new Courbe("Patron");
        CourbPatronBack[i] = new Courbe("PatronBack");
        CourbPatron[i]->points = OFF;
        CourbPatronDXF[i]->points = OFF;
        CourbPatronBack[i]->points = OFF;
        CourbPatron[i]->symX = OFF;
        CourbPatronDXF[i]->symX = OFF;
        CourbPatronBack[i]->symX = OFF;
        CourbPatron[i]->pts = new Matrix(Xd[i]->GetLignes(), 2);
        CourbPatronDXF[i]->pts = new Matrix(Xd[i]->GetLignes(), 2);
        CourbPatronBack[i]->pts = new Matrix(Xd[i]->GetLignes(), 2);
        for (j = 0; j < Xd[i]->GetLignes(); j++) {
            CourbPatron[i]->pts->SetElement(j, 0, Xd[i]->Element(j, 0));
            CourbPatronDXF[i]->pts->SetElement(j, 0, Xd[i]->Element(j, 0));
            CourbPatronBack[i]->pts->SetElement(Xd[i]->GetLignes() - j - 1, 0, Xd[i]->Element(j, 0));

            CourbPatron[i]->pts->SetElement(j, 1, Yd[i]->Element(j, 0));
            CourbPatronDXF[i]->pts->SetElement(j, 1, Yd[i]->Element(j, 0));
            CourbPatronBack[i]->pts->SetElement(Xd[i]->GetLignes() - j - 1, 1, Yd[i]->Element(j, 0));
        }

        /*for (int t = 0; t < CourbPatronBack[i]->pts->GetLignes(); t++) {
            printf ("\nPatronBack [%d] (%d) -> (%f, %f) ", i, t,
                    1000*CourbPatronBack[i]->pts->Element(t, 0),
                    1000*CourbPatronBack[i]->pts->Element(t, 1));
        }*/

        addCourbe(AxePatron, CourbPatron[i]);

        if (wasPince){
            //CourbPin
            CourbPin[i] = new Courbe("Pin");
            CourbPin[i]->points = OFF;
            CourbPin[i]->symX = OFF;
            CourbPin[i]->pts = new Matrix(_Xd[i]->GetLignes(), 2);
            for (j = 0; j < _Xd[i]->GetLignes(); j++) {
                CourbPin[i]->pts->SetElement(j, 0, _Xd[i]->Element(j, 0));
                CourbPin[i]->pts->SetElement(j, 1, _Yd[i]->Element(j, 0));
            }
            addCourbe(AxePatron, CourbPin[i]);
        }

        distance = Ones(CourbPatron[i]->pts->GetLignes(), 1);
        for (j = 0; j < distance->GetLignes(); j++) distance->MultiplyElement(j, 0, Marge[i] / 100.0f);
        CourbMarge[i] = new Courbe("Marge");
        CourbMargeDXF[i] = new Courbe("Marge");
        CourbMargeBack[i] = new Courbe("MargeBack");
        if (i == 0)  CourbMarge[i]->pts = calcContour(CourbPatron[i]->pts, distance, -1);
            else CourbMarge[i]->pts = calcContour(CourbPatron[i]->pts, distance, +1);

        CourbMargeDXF[i]->pts = CloneMat(CourbMarge[i]->pts);
        CourbMargeBack[i]->pts = CloneMat(CourbMarge[i]->pts);

        CourbMarge[i]->points = OFF;
        CourbMarge[i]->symX = OFF;
        CourbMargeDXF[i]->points = OFF;
        CourbMargeDXF[i]->symX = OFF;
        CourbMargeBack[i]->points = OFF;
        CourbMargeBack[i]->symX = OFF;

        for (j = 0; j < CourbMarge[i]->pts->GetLignes(); j++) {
            CourbMargeBack[i]->pts->SetElement(CourbMarge[i]->pts->GetLignes() - j - 1, 0, CourbMarge[i]->pts->Element(j, 0));
            CourbMargeBack[i]->pts->SetElement(CourbMarge[i]->pts->GetLignes() - j - 1, 1, CourbMarge[i]->pts->Element(j, 1));
            /*printf ("\nMarge[%d](%d)    -> (%f,%f) ", i, j,
                    1000*CourbMargeDXF[i]->pts->Element(j, 0),
                    1000*CourbMargeDXF[i]->pts->Element(j, 1));*/
        }

/*        for (t = 0; t < CourbMargeBack[i]->pts->GetLignes(); t++) {
            printf ("\nMargeBack [%d](%d) -> (%f, %f) ", i, t,
                    1000*CourbMargeBack[i]->pts->Element(t, 0),
                    1000*CourbMargeBack[i]->pts->Element(t, 1));
        }*/


        addCourbe(AxePatron, CourbMarge[i]);
        delete(distance);
    }

    addCourbe(AxePatronDXF, CourbPatronDXF[0]);
    if ((Xd[0]->Element(0, 0) != Xd[1]->Element(0, 0)) //test points Av cote 1&2 confondus
            || (Yd[0]->Element(0, 0) != Yd[1]->Element(0, 0))) {
        
        //left points not equal
        //printf ("\n left points not equal");
        CourbAv = new Courbe("Avant");
        CourbAvBack = new Courbe("AvantBack");
        CourbAv->pts = new Matrix(2, 2);
        CourbAvBack->pts = new Matrix(2, 2);

        CourbAv->pts->SetElement(0, 0, Xd[0]->Element(0, 0));
        CourbAv->pts->SetElement(0, 1, Yd[0]->Element(0, 0));
        CourbAv->pts->SetElement(1, 0, Xd[1]->Element(0, 0));
        CourbAv->pts->SetElement(1, 1, Yd[1]->Element(0, 0));
        
        CourbAvBack->pts->SetElement(0, 0, Xd[1]->Element(0, 0));
        CourbAvBack->pts->SetElement(0, 1, Yd[1]->Element(0, 0));
        CourbAvBack->pts->SetElement(1, 0, Xd[0]->Element(0, 0));
        CourbAvBack->pts->SetElement(1, 1, Yd[0]->Element(0, 0));


        CourbAv->points = OFF;
        CourbAv->symX = OFF;
        CourbAvBack->points = OFF;
        CourbAv->symX = OFF;
        addCourbe(AxePatron, CourbAv);
        ajCoAv = 1;
        //marge avant
        distance = Ones(2, 1);
        for (j = 0; j < distance->GetLignes(); j++)
            distance->MultiplyElement(j, 0, MargeDeb / 100.0f);
        CourbMargeAv = new Courbe("MargeAV");
        CourbMargeAvBack = new Courbe("MargeAvBack");
        CourbMargeAv->pts = calcContour(CourbAv->pts, distance, +1);
        CourbMargeAvBack->pts = calcContour(CourbAv->pts, distance, +1);
        for (j = 0; j < CourbMargeAv->pts->GetLignes(); j++) {
            CourbMargeAvBack->pts->SetElement(CourbMargeAv->pts->GetLignes() - j - 1, 0, CourbMargeAv->pts->Element(j, 0));
            CourbMargeAvBack->pts->SetElement(CourbMargeAv->pts->GetLignes() - j - 1, 1, CourbMargeAv->pts->Element(j, 1));
        }

        CourbMargeAv->points = OFF;
        CourbMargeAv->symX = OFF;
        CourbMargeAvBack->points = OFF;
        CourbMargeAvBack->symX = OFF;

        addCourbe(AxePatron, CourbMargeAv);
        ajCoMaAv = 1;

        delete(distance);

        if  ( 	( abs(CourbMargeAv->pts->Element(0, 0) -  CourbMargeAv->pts->Element(1, 0)) > 0.000001 ) 
				|| 
				( abs(CourbMargeAv->pts->Element(0, 1) - CourbMargeAv->pts->Element(1, 1)) > 0.000001 )	) {

			for (i = 0; i < 2; i++)
			{
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

				// CourbCoin1[i]->pts->SetElement(2,0, CourbMargeAv->pts->Element(0,0));
				// CourbCoin->pts->SetElement(2,0, CourbMargeAv->pts->Element(0,0));
				CourbCoin1[i]->pts->SetElement(2, 0, CourbMargeAv->pts->Element(i, 0));
				CourbCoin->pts->SetElement(2, 0, CourbMargeAv->pts->Element(i, 0));

				// CourbCoin1[i]->pts->SetElement(2,1, CourbMargeAv->pts->Element(0,1));
				// CourbCoin->pts->SetElement(2,1, CourbMargeAv->pts->Element(0,1));
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
				addCourbe(AxePatron, CourbCoin);
				ajCo1[i] = 1;

			}
		}

    }
    /*else //relie les marges cote 1&2 par un segment
    {
        printf ("\n left points equal");
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

        addCourbe(AxePatron, CourbAv);
        //
        ajCoAv = 0;
        //
    }*/
    int n0 = Xd[0]->GetLignes() - 1;
    int n1 = Xd[1]->GetLignes() - 1;
    if ((Xd[0]->Element(n0, 0) != Xd[1]->Element(n1, 0)) //test points Ar cote 1&2 confondus
            || (Yd[0]->Element(n0, 0) != Yd[1]->Element(n1, 0))) {
        CourbAr = new Courbe("AR");
        CourbArDXF = new Courbe("AR");

        CourbAr->pts = new Matrix(2, 2);
        CourbArDXF->pts = new Matrix(2, 2);

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

        addCourbe(AxePatron, CourbAr);
        ajCoAr = 1;

        distance = Ones(2, 1);
        for (j = 0; j < distance->GetLignes(); j++)
            distance->MultiplyElement(j, 0, MargeFin / 100.0f);

        CourbMargeAr = new Courbe("MargeAR");
        CourbMargeArDXF = new Courbe("MargeAR");
        CourbMargeAr->pts = calcContour(CourbAr->pts, distance, -1);
        CourbMargeArDXF->pts = calcContour(CourbAr->pts, distance, -1);

        CourbMargeAr->points = OFF;
        CourbMargeAr->symX = OFF;
        CourbMargeArDXF->points = OFF;
        CourbMargeArDXF->symX = OFF;

        addCourbe(AxePatron, CourbMargeAr);
        ajCoMaAr = 1;

        delete(distance);

        if  ( 	( abs(CourbMargeAr->pts->Element(0, 0) -  CourbMargeAr->pts->Element(1, 0)) > 0.000001 ) 
				|| 
				( abs(CourbMargeAr->pts->Element(0, 1) - CourbMargeAr->pts->Element(1, 1)) > 0.000001 )	) 
			{
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

				n = CourbMarge[i]->pts->GetLignes() - 1;

				CourbCoin2[i]->pts->SetElement(0, 0, CourbMarge[i]->pts->Element(n, 0));
				CourbCoin->pts->SetElement(0, 0, CourbMarge[i]->pts->Element(n, 0));

				CourbCoin2[i]->pts->SetElement(0, 1, CourbMarge[i]->pts->Element(n, 1));
				CourbCoin->pts->SetElement(0, 1, CourbMarge[i]->pts->Element(n, 1));

				CourbCoin2[i]->pts->SetElement(2, 0, CourbMargeAr->pts->Element(i, 0));
				CourbCoin->pts->SetElement(2, 0, CourbMargeAr->pts->Element(i, 0));

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
					CourbMarge[i]->pts->SetElement(n, 0, x);
					CourbMarge[i]->pts->SetElement(n, 1, y);
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

				addCourbe(AxePatron, CourbCoin);
				ajCo2[i] = 1;

			}
		}

    }
    /*else //relie les marges cote 1&2 par un segment
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

        addCourbe(AxePatron, CourbAr);
        ajCoAr = 0;
    }*/

    if (ajCoAr) {
        //printf ("\n ajCoAR!!!");
        addCourbe(AxePatronDXF, CourbArDXF);
    } else {
        //printf ("\n NOT ajCoAR!!!");
    }

    addCourbe(AxePatronDXF, CourbPatronBack[1]);
    if (ajCoAv) {
        //printf ("\n ajCoAvBack!!!");
        /*printf ("\n(%f, %f) -> (%f, %f)", 1000*CourbAvBack->pts->Element(0,0), 1000*CourbAvBack->pts->Element(0,1),
                                            1000*CourbAvBack->pts->Element(1,0), 1000*CourbAvBack->pts->Element(1,1) );*/
        addCourbe(AxePatronDXF, CourbAvBack);
    } else {
        //printf ("\n NOT ajCoAvBack!!!");
    }
    //printf ("\n");
    addCourbe(AxeMarginDXF, CourbMargeDXF[0]);
    if (ajCo2[0]) addCourbe(AxeMarginDXF, CourbCoin2[0]);


    if (ajCoMaAr) {
        addCourbe(AxeMarginDXF, CourbMargeArDXF);
        //printf ("\n MargeajCoMaAr!!!");
    } else {
        //printf ("\n NOT MargeajCoMaAr!!!");
    }
    if (ajCo2[1]) addCourbe(AxeMarginDXF, CourbCoin2Back[1]);

    addCourbe(AxeMarginDXF, CourbMargeBack[1]);
    if (ajCo1[1]) addCourbe(AxeMarginDXF, CourbCoin1[1]);
    if (ajCoMaAv) {
        //printf ("\n MargeajCoMaAvBack!!!");
        addCourbe(AxeMarginDXF, CourbMargeAvBack);
    } else {
        //printf ("\n NOT MargeajCoMaAvBack!!!");
    }
    if (ajCo1[0]) addCourbe(AxeMarginDXF, CourbCoin1Back[0]);

/*    i=0;
    Matrix* m = getReperPoints(getWindPatternsProject(), Xd[i], Yd[i], P[i], nerv[i], Deb[i], FaceDeb[i], Fin[i], FaceFin[i]);
    for (i = 0; i < m->GetLignes(); i++) printf ("\n%d -> %f %f %f %f", i, m->Element(i,0), m->Element(i,1), m->Element(i,2), m->Element(i,3));
    delete (m);*/

    int ventisave=0;
    double xk0, xk1, xk2, yk0, yk1, yk2, xkf, ykf;
    if (!isKlapan)
    for (i = 0; i < 2; i++) {
        //printf ("\n make ReperPoints i=%d", i);
        for (j = 0;j < 5; j++) {
            xSuspente[i][j] = -100000.0f;
            ySuspente[i][j] = -100000.0f;
        }
            Matrix* m = getReperPoints(getWindPatternsProject(), Xd[i], Yd[i], P[i], nerv[i], Deb[i], FaceDeb[i], Fin[i], FaceFin[i], (NoNerv[0]==NoNerv[1]));
            for (j = 0; j < m->GetLignes(); j++)
            {
                    int mark = m->Element(j,2);
                    double xSus = m->Element(j,0);
                    double ySus = m->Element(j,1);

                    if ((mark == REP_TRIANGLE) && (m->Element(j, 3) > 0.0f)) {
                        ventisave = (int) m->Element(j, 3) - 1;
                        xSuspente[i][ventisave] = xSus;
                        ySuspente[i][ventisave] = ySus;
                    }
                    if (i==0) {
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

                    if ((ReperesProfile[i]) && (mark == REP_CROSS)) {
                            CourbRep = new Courbe("Reperes suspentes 1");
                            CourbRep->points = OFF;
                            CourbRepDXF = new Courbe("Reperes suspentes 1");
                            CourbRepDXF->points = OFF;
                            CourbRep->pts = Zeros(2, 2);
                            CourbRepDXF->pts = Zeros(2, 2);
                            CourbRep->pts->SetElement(0, 0, xSus - XMashtab * 0.01f);
                            CourbRep->pts->SetElement(0, 1, ySus - XMashtab * 0.01f);
                            CourbRep->pts->SetElement(1, 0, xSus + XMashtab * 0.01f);
                            CourbRep->pts->SetElement(1, 1, ySus + XMashtab * 0.01f);
                            int _i, _j;
                            for (_i = 0; _i < 2; _i++)
                                for (_j = 0; _j < 2; _j++)
                                    CourbRepDXF -> pts -> SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                            addCourbe(AxeRepDXF, CourbRepDXF);
                            addCourbe(AxePatron, CourbRep);
                            CourbRep = new Courbe("Reperes suspentes 2");
                            CourbRep->points = OFF;
                            CourbRepDXF = new Courbe("Reperes suspentes 2");
                            CourbRepDXF->points = OFF;
                            CourbRep->pts = Zeros(2, 2);
                            CourbRepDXF->pts = Zeros(2, 2);
                            CourbRep->pts->SetElement(0, 0, xSus - XMashtab * 0.01f);
                            CourbRep->pts->SetElement(0, 1, ySus + XMashtab * 0.01f);
                            CourbRep->pts->SetElement(1, 0, xSus + XMashtab * 0.01f);
                            CourbRep->pts->SetElement(1, 1, ySus - XMashtab * 0.01f);
                            for (_i = 0; _i < 2; _i++)
                                for (_j = 0; _j < 2; _j++)
                                    CourbRepDXF->pts->SetElement(_i, _j, CourbRep->pts->Element(_i, _j));
                            addCourbe(AxeRepDXF, CourbRepDXF);
                            addCourbe(AxePatron, CourbRep);
                        } 
                        if (mark==REP_LINE || ((ReperesProfile[i]) && (mark == REP_MIDDLE_LINE))) {
                                //perpendiular
                                int direction = -1;
                                if (i == 1) direction = 1;
                                double _x0 = xSus;
                                double _y0 = ySus;
                                int ix, icol;
                                Ind(_x0, Xd[i], &ix, &icol);
                                //Ind(ySuspente[i][j], interpYSuspente, &iy, &icol);
                                double xc, yc;
                                double _l = 3 * XMashtab * 0.01f;
                                if (mark == REP_MIDDLE_LINE) {
                                    _l = 0.4 * _l;
                                }
                                if ((Xd[i]->Element(ix, 0) > _x0) && (ix != 0)) ix--;
                                if ((Xd[i]->Element(ix,0) == _x0) && (ix != 0)) ix--;
                                if (ix != 0)
                                    calcVecteurNormal(Xd[i]->Element(ix,0), Yd[i]->Element(ix, 0),_x0, _y0, &xc, &yc, _l, direction);
                                else
                                    calcVecteurNormal(_x0, _y0,Xd[i]->Element(1,0), Yd[i]->Element(1, 0), &xc, &yc, _l, direction);
                                double _x1 = _x0;
                                double _y1 = _y0;
                                double _x2 = xc;
                                double _y2 = yc;
                                if (ix == 0) {
                                    _x2 = _x0 + (xc - Xd[i]->Element(1,0));
                                    _y2 = _y0 + (yc - Yd[i]->Element(1,0));
                                }

                                int _i, _j;
                                CourbRep = new Courbe("RepP");
                                CourbRep->points = OFF;
                                CourbRepDXF = new Courbe("RepP");
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
                                addCourbe(AxeRepDXF, CourbRepDXF);
                                addCourbe(AxePatron, CourbRep);
                            }
                            if ((ReperesSuspentes[i] == 1) && ((mark==REP_TRIANGLE) || (mark==REP_V)))  {
                                    int direction = -1;
                                    if (i == 1) direction = 1;
                                    double _x0 = xSus;
                                    double _y0 = ySus;
                                    double _l = 2*XMashtab * 0.01f;
                                    double _x1 = _x0;
                                    double _y1 = _y0;

                                    int ix, icol;

                                    Ind(xSus, Xd[i], &ix, &icol);
                                    double xc, yc;
                                    if ((Xd[i]->Element(ix,0) > _x0) && (ix != 0 )) ix--;
                                    if ((Xd[i]->Element(ix,0) == _x0) && (ix != 0) ) ix--;
                                    if  (ix != 0)
                                        calcVecteurNormal(Xd[i]->Element(ix,0), Yd[i]->Element(ix, 0), _x0, _y0, &xc, &yc, _l, direction);
                                    else
                                        calcVecteurNormal( _x0, _y0, Xd[i]->Element(1,0), Yd[i]->Element(1, 0), &xc, &yc, _l, direction);
                                    if (ix == 0) ix = 1;
                                    double nx = (Xd[i]->Element(ix, 0)-_x0);
                                    double ny = (Yd[i]->Element(ix, 0)-_y0);
                                    double rasst = sqrt (sqr(nx) + sqr(ny));
                                    nx = _l/2.0f*nx/rasst;
                                    ny = _l/2.0f*ny/rasst;

                                    double _x2 = xc + nx;
                                    double _y2 = yc + ny;

                                    double _x3 = xc - nx;
                                    double _y3 = yc - ny;

                                    int _i, _j;
                                    CourbRep = new Courbe("RepP");
                                    CourbRep->points = OFF;
                                    CourbRepDXF = new Courbe("RepP");
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
                                    addCourbe(AxeRepDXF, CourbRepDXF);
                                    addCourbe(AxePatron, CourbRep);
                                    if (mark != REP_V) {
                                        CourbRep = new Courbe("RepP");
                                        CourbRep->points = OFF;
                                        CourbRepDXF = new Courbe("RepP");
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
                                        addCourbe(AxeRepDXF, CourbRepDXF);
                                        addCourbe(AxePatron, CourbRep);
                                    }

                                    CourbRep = new Courbe("RepP");
                                    CourbRep->points = OFF;
                                    CourbRepDXF = new Courbe("RepP");
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
                                    addCourbe(AxeRepDXF, CourbRepDXF);
                                    addCourbe(AxePatron, CourbRep);
                            }
                        
                    }

            if ((i==0) && (NoNerv[0]==NoNerv[1]) && (LayoutKlapans) && (VentHoles)) {
                bool match = false;
                int _i, _j;
                //printf ("\nnoNervVH: %d %d", noNervVH[NoNerv[0]], noNervVH[NoNerv[0]-1]);
                if ((NoNerv[0]==-1) || (NoNerv[0]==0)){
                  //  printf ("\n if 1");
                        if (VentCentralNerv) match = true;
                }
                    
                else {
                    //        printf ("\n if 2");
                        if ((noNervVH[NoNerv[0]] == 1) || (noNervVH[NoNerv[0] - 1]==1)) match = true;
                }

                if (match) {
                    //printf (" \n MATCH!");
                    if (VentHolesDouble) {
                        // xk0, yk0 -> xkfin, ykfin
                        CourbRep = new Courbe("RepP");
                        CourbRep->points = OFF;
                        CourbRepDXF = new Courbe("RepP");
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
                        addCourbe(AxeRepDXF, CourbRepDXF);
                        addCourbe(AxePatron, CourbRep);

                    }
                    // xk1, yk1 -> xkfin, ykfin
                    CourbRep = new Courbe("RepP");
                    CourbRep->points = OFF;
                    CourbRepDXF = new Courbe("RepP");
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
                    addCourbe(AxeRepDXF, CourbRepDXF);
                    addCourbe(AxePatron, CourbRep);

                    // xk2, yk2 -> xkfin, ykfin
                    CourbRep = new Courbe("RepP");
                    CourbRep->points = OFF;
                    CourbRepDXF = new Courbe("RepP");
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
                    addCourbe(AxeRepDXF, CourbRepDXF);
                    addCourbe(AxePatron, CourbRep);

                } else {
                    //printf (" \n NOOOOT MATCH!");
                }

            }

    }// endfor i

    if (Ventilation == 1) {
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
                addCourbe(AxePatron, CourbCercle);
                addCourbe(AxeCercleDXF, CourbCercleDXF);
            }
        }
    }
       // add texte !!!

    if (Numerotation == 1) {
        xText = (Xd[0]->Element(0, 0) + Xd[0]->Element(Xd[0]->GetLignes() - 1, 0)
                + Xd[1]->Element(0, 0) + Xd[1]->Element(Xd[1]->GetLignes() - 1, 0)) / 4.0f;
        double x0len = Xd[0]->Element(Xd[0]->GetLignes() - 1, 0) - Xd[0]->Element(0, 0);
        double x1len = Xd[1]->Element(Xd[1]->GetLignes() - 1, 0) - Xd[1]->Element(0, 0);
        double xlen = (x0len + x1len) / 2.0f;
        xText = textX * xlen;
        yText = (Yd[0]->Element(0, 0) + Yd[0]->Element(Yd[0]->GetLignes() - 1, 0)
                + Yd[1]->Element(0, 0) + Yd[1]->Element(Yd[1]->GetLignes() - 1, 0)) / 4.0f;
        double y0len = Yd[1]->Element(0, 0) - Yd[0]->Element(0, 0);
        double y1len = Yd[1]->Element(Yd[1]->GetLignes() - 1, 0) - Yd[0]->Element(Yd[0]->GetLignes() - 1, 0);
        double ylen = (y0len + y1len) / 2.0f;
        yText = Yd[1]->Element(0, 0) - textY*ylen;
        //sprintf(texte, "%d%c%03.1f%c%03.1f_%d%c%03.1f%c%03.1f",
        //	NoNerv[0],texteExtInt[FaceDeb[0]-1],Deb[0],texteExtInt[FaceFin[0]-1],Fin[0],
        //	NoNerv[1],texteExtInt[FaceDeb[1]-1],Deb[1],texteExtInt[FaceFin[1]-1],Fin[1]);
        sprintf(texte, "%s", NumText->get_text());
        addTexte(AxePatron, texte, 0.02f, 0.0f, xText, yText);
        sprintf(texte, "%s", NumText->get_text());
        addTexte(AxePatronTextDXF, texte, 0.02f, 0.0f, xText, yText);
    }

    //add d'une courbe bidon pour le zoom
    CourbZoom = new Courbe("Zoom");
    CourbZoom->pts = Zeros(5, 2);
    CourbZoom->points = OFF;
    CourbZoom->segments = OFF;
    CourbZoom->symX = OFF;
    CourbZoom->colorSegments[0] = 0.0f;
    CourbZoom->colorSegments[2] = 0.0f;
    addCourbe(AxePatron, CourbZoom);

    //creation Mesh pour visualisation 3d du patron
    MeshPatron = createMesh();
    MeshPatron->InvNormales = OFF;
    MeshPatron->segments = OFF;
    MeshPatron->faces = ON;
    MeshPatron->x = new Matrix(X[0]->GetLignes(), 2);
    MeshPatron->y = new Matrix(X[0]->GetLignes(), 2);
    MeshPatron->z = new Matrix(X[0]->GetLignes(), 2);
    for (i = 0; i < X[0]->GetLignes(); i++) {
        MeshPatron->x->SetElement(i, 0, X[0]->Element(i, 0));
        MeshPatron->x->SetElement(i, 1, X[1]->Element(i, 0));
        MeshPatron->y->SetElement(i, 0, Y[0]->Element(i, 0));
        MeshPatron->y->SetElement(i, 1, Y[1]->Element(i, 0));
        MeshPatron->z->SetElement(i, 0, Z[0]->Element(i, 0));
        MeshPatron->z->SetElement(i, 1, Z[1]->Element(i, 0));
    }
    addMesh(Axe3d, MeshPatron);
    //definition eclairage
    Axe3d->eclairage = ON;

    //printf ("Libere");
    for (i = 0; i < 2; i++) {
        delete(X[i]);
        delete(Y[i]);
        delete(Z[i]);
        delete(P[i]);
        delete(Xd[i]);
        delete(Yd[i]);
    } 
    delete(XExt);
    delete(YExt);
    delete(ZExt);

    delete(XInt);
    delete(YInt);
    delete(ZInt);
}

/****************/
/* display      */
/****************/

void display(void) {
    //calcVue3dEtPatron();
    ViewAxe(Axe3d);
    glutSwapBuffers();
	//ViewAxe(Axe3dBal);
    //glutSwapBuffers();
    ViewAxe(AxePatron);
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
        
	/*si debut du zoom, init cadre*/
        if (debZoomIN) {
            debZoomIN = false;

            CourbZoom->pts->SetElement(0, 0, xs);
            CourbZoom->pts->SetElement(0, 1, ys);

            CourbZoom->pts->SetElement(4, 0, xs);
            CourbZoom->pts->SetElement(4, 1, ys);

            CourbZoom->pts->SetElement(3, 0, xs);
            CourbZoom->pts->SetElement(1, 1, ys);

        }

        /*mise a jour position cadre*/
        CourbZoom->pts->SetElement(1, 0, xs);
        CourbZoom->pts->SetElement(3, 1, ys);

        CourbZoom->pts->SetElement(2, 0, xs);
        CourbZoom->pts->SetElement(2, 1, ys);
    }/*si fin du zoom, init position axe*/

    else if (finZoomIN) {
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
    }

    /*affichage*/
    glutPostRedisplay();

}

/****************/
/* BoutonSouris */
/****************/

void BoutonSouris(int button, int state, int x, int y) {
    int window;
    //double xmil, ymil; //pour zoom OUT
    /*memorise position de la souris dans la fenetre*/
    xSouris = x;
    ySouris = y;
    AxeSel = NULL;
    /*Rotation3d = OFF;
    Translation3d = OFF;
    ZoomTwist3d = OFF;*/

    window = glutGetWindow();
    if (window == window3d) {
        AxeSel = Axe3d;
    }
    //if (window == window3dBal) {
//        AxeSel = Axe3dBal;
    //}

    if (window == windowPatron) {
        AxeSel = AxePatron;
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

    /*test debut/fin zoomIN, zoomOUT, retour Axe Automatique pour Axe3d*/
    if (AxeSel == AxePatron) {
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
                    AxePatron->XAuto = OFF;
                    AxePatron->YAuto = OFF;
                    AxePatron->ZAuto = OFF;
                    break;
                case GLUT_RIGHT_BUTTON: //zoomOUT x2

                    if (zoomIN) {
                        quitZoom = true;
                        break;
                    }

                    zoomOUT = true;

                    AxePatron->XAuto = OFF;
                    AxePatron->YAuto = OFF;
                    AxePatron->ZAuto = OFF;

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

                AxePatron->XAuto = ON;
                AxePatron->YAuto = ON;
                AxePatron->ZAuto = ON;
            }

            glutPostRedisplay();

        }

    }

}

/*****************/
/* keyboard      */
/*****************/

void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case 27: quit(0);
            break;
        default:;
    }
}

/*****************/
/* InitWindow   */

/*****************/

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
}


int main(int argc, char** argv) {
    int ws, hs;
    double EpaiRelProfCent, EpaiRelProfBout;
    /*message*/
    printf("\nWindPatterns ");
    printf(AUTEUR_DATE);
    printf("\n");

    /*init glut*/
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(100, 100);

    /* creation fenetre/axe/courbes Axe3d*/
    window3d = glutCreateWindow("Viewalisation 3D");
    InitWindow();
    InitLight();
    Axe3d = createAxe(window3d);
    Axe3d->axe3d = ON;

    //window3dBal = glutCreateWindow("Viewalisation 3D Balone");
    //InitWindow();
    //InitLight();
    //Axe3dBal = createAxe(window3dBal);
    //Axe3dBal->axe3d = ON;


    /* creation fenetre/axe/courbes AxePatron*/
    windowPatron = glutCreateWindow("Patron");
    InitWindow(); //InitLight();
    AxePatron = createAxe(windowPatron);
    AxePatronDXF = createAxe(windowPatron);
    AxeRepDXF = createAxe(windowPatron);
    AxePatronTextDXF = createAxe(windowPatron);
    AxeMarginDXF = createAxe(windowPatron);
    AxeCercleDXF = createAxe(windowPatron);

    /*recupere taille de l'ecran en pixel*/
    ws = glutGet(GLUT_SCREEN_WIDTH);
    hs = glutGet(GLUT_SCREEN_HEIGHT);
    /*redefini taille et position des fenetre 3D*/
    glutSetWindow(window3d);
    glutPositionWindow((int) (0.70 * (double) ws), (int) (0.05 * (double) hs));
    glutReshapeWindow((int) (0.30 * (double) ws), (int) (0.40 * (double) hs));
    //glutSetWindow(window3dBal);
    //glutPositionWindow((int) (0.70 * (double) ws), (int) (0.05 * (double) hs));
    //glutReshapeWindow((int) (0.30 * (double) ws), (int) (0.40 * (double) hs));

    /*redefini taille et position des fenetre Patron*/
    glutSetWindow(windowPatron);
    glutPositionWindow((int) (0.0 * (double) ws), (int) (0.50 * (double) hs));
    glutReshapeWindow((int) (1.0 * (double) ws), (int) (0.45 * (double) hs));

    strcpy(fileNameProject, "f17project.wpp");
    gfd = readWindPatternsProject(fileNameProject);

    LoadFromWindPatternsProject(gfd);

    EpaiRelProfCent = EpaisseurRelative(ExtProfCent, IntProfCent);
    EpaiRelProfBout = EpaisseurRelative(ExtProfBout, IntProfBout);

    InterpoleProfilBout(&ExtProfBout, ExtProfCent);
    InterpoleProfilBout(&IntProfBout, IntProfCent);

	//printf("Ce programme utilise GLUI version: %3.2f\n", GLUI_Master.get_version());
    glui = GLUI_Master.create_glui("Boite de dialogue", 0, 0, 0);

    /*cote 1*/
    /*No Nerv 1*/
    SpinNoNerv[0] = glui->add_spinner("No Nerv 1", GLUI_SPINNER_INT, &(NoNerv[0]));
    SpinNoNerv[0] -> set_int_limits(-1, F->m_nbProfils - 1);
    /*% Deb 1*/

    GLUI_Panel *panel_Deb1 = glui->add_panel(""); //, GLUI_PANEL_NONE);
    GLUI_Spinner *SpinDeb1 = glui->add_spinner_to_panel(
            panel_Deb1, "% Deb 1", GLUI_SPINNER_FLOAT, &(Deb[0]));

    SpinDeb1 -> set_float_limits(0.0, 100.0);
    GLUI_Listbox *ListDeb1 = glui->add_listbox_to_panel(panel_Deb1, "Face", &FaceDeb[0]);
    ListDeb1->add_item(1, "Ext");
    ListDeb1->add_item(2, "Int");

    /*% Deb 1*/
    //GLUI_Panel *panel_Fin1 = glui->add_panel(""); //, GLUI_PANEL_NONE);
    GLUI_Spinner *SpinFin1 = glui->add_spinner_to_panel(
            panel_Deb1, "% Fin 1", GLUI_SPINNER_FLOAT, &(Fin[0]));
    SpinFin1 -> set_float_limits(0.0, 100.0);
    GLUI_Listbox *ListFin1 = glui->add_listbox_to_panel(panel_Deb1, "Face", &FaceFin[0]);
    ListFin1->add_item(1, "Ext");
    ListFin1->add_item(2, "Int");

    GLUI_Rollout *panel_PosNerv1 =glui->add_rollout("Pos Nerv 1",false);
    glui->add_checkbox_to_panel(panel_PosNerv1, "Pos Nerv 1", &isPosNerv[0]);
    GLUI_Spinner *SpinPosNerv1 = glui->add_spinner_to_panel(panel_PosNerv1, "Pos", GLUI_SPINNER_FLOAT, &(PosNerv[0]));
    SpinPosNerv1 -> set_float_limits(0.0, 100.0);

    /*Reperes 1*/
    /*
        GLUI_Panel *panel_Rep1 = glui->add_panel("");//, GLUI_PANEL_NONE);
        GLUI_Spinner *SpinPosRep1 = glui->add_spinner_to_panel(
                panel_Rep1,"% Rep", GLUI_SPINNER_FLOAT, &(PosRep[0]));
        SpinPosRep1 -> set_float_limits( 0.0, 100.0 );
        GLUI_Listbox *ListRep1 = glui->add_listbox_to_panel(panel_Rep1,"Face",&FaceRep[0]);
        ListRep1->add_item( 1, "Ext" ); ListRep1->add_item( 2, "Int" );
        glui->add_button_to_panel(panel_Rep1, "adder", 0, adder );
     */

    /*Marge de couture en cm*/
    GLUI_Spinner *SpinMarge1 = glui->add_spinner("Marge 1", GLUI_SPINNER_FLOAT, &(Marge[0]));
    SpinMarge1 -> set_float_limits(0.0, 5.0);

    /*reperes suspentes*/
    glui->add_checkbox("reperes suspentes", &(ReperesSuspentes[0]));
    glui->add_checkbox("reperes profile", &(ReperesProfile[0]));

    /*panel pinces!!!*/
    GLUI_Panel *panel_Pince1 = glui->add_panel("Pinces 1"); //, GLUI_PANEL_NONE);
    /*position debut pincement au BA en %*/
    GLUI_Spinner *SpinPosPinceBA1 = glui->add_spinner_to_panel(panel_Pince1,
            "Pos BA (%)", GLUI_SPINNER_FLOAT, &(PosPinceBA[0]));
    SpinPosPinceBA1 -> set_float_limits(0.0, 40.0);
    /*Amplitude pincement au BA en %*/
    GLUI_Spinner *SpinAmpPinceBA1 = glui->add_spinner_to_panel(panel_Pince1,
            "Amp BA (%)", GLUI_SPINNER_FLOAT, &(AmpPinceBA[0]));
    SpinAmpPinceBA1 -> set_float_limits(0.0, 40.0);
    /*position debut pincement au BF en %*/
    GLUI_Spinner *SpinPosPinceBF1 = glui->add_spinner_to_panel(panel_Pince1,
            "Pos BF (%)", GLUI_SPINNER_FLOAT, &(PosPinceBF[0]));
    SpinPosPinceBF1 -> set_float_limits(0.0, 40.0);
    /*Amplitude pincement au BF en %*/
    GLUI_Spinner *SpinAmpPinceBF1 = glui->add_spinner_to_panel(panel_Pince1,
            "Amp BF (%)", GLUI_SPINNER_FLOAT, &(AmpPinceBF[0]));
    SpinAmpPinceBF1 -> set_float_limits(0.0, 40.0);

    //pinceType
    GLUI_Rollout *panel_PinceTypeA=glui->add_rollout("Pince type A(Nos)", false);
    //GLUI_Panel *panel_PinceTypeA = glui->add_panel_to_panel(RolloutPinceTypeA, "Pince type A(Nos)");
    GLUI_Spinner *pinceRadiusAlgKNosSpinner = glui->add_spinner_to_panel(panel_PinceTypeA,
            "k (d.width)", GLUI_SPINNER_FLOAT, &(PinceRadiusAlgKNos));
    pinceRadiusAlgKNosSpinner -> set_float_limits(0.0001, 5.0);

    glui->add_checkbox_to_panel(panel_PinceTypeA, "Set = Amp to Hvost", &PinceNosEqualAmp, 0, &modifPinceNosRadio);
    //glui->add_checkbox_to_panel(panel_PinceTypeA, "Amp 0", &PinceNos0Amp, 0,  &modifPinceNosRadio);
    GLUI_RadioGroup *radio_pinceNos =
            glui->add_radiogroup_to_panel(panel_PinceTypeA, &PinceNosRadio, 0, &modifPinceNosRadio);

    glui->add_radiobutton_to_group(radio_pinceNos, "Radius auto alg");

    glui->add_radiobutton_to_group(radio_pinceNos, "Function");
    glui->add_checkbox_to_panel(panel_PinceTypeA, "Power", &PincePowerA);
    GLUI_Spinner *pincePowerSpinnerA = glui->add_spinner_to_panel(panel_PinceTypeA,
            "power", GLUI_SPINNER_FLOAT, &(PincePowerValueA));
    pincePowerSpinnerA -> set_float_limits(0.0, 5.0);
    glui->add_checkbox_to_panel(panel_PinceTypeA, "Arctan", &PinceArctanA);
    GLUI_Spinner *pinceArctanK1SpinnerA = glui->add_spinner_to_panel(panel_PinceTypeA,
            "k1", GLUI_SPINNER_FLOAT, &(PinceArctanK1ValueA));
    pinceArctanK1SpinnerA -> set_float_limits(0.0, 10.0);
    GLUI_Spinner *pinceArctanK2SpinnerA = glui->add_spinner_to_panel(panel_PinceTypeA,
            "k2", GLUI_SPINNER_FLOAT, &(PinceArctanK2ValueA));
    pinceArctanK2SpinnerA -> set_float_limits(0.0, 20.0);
    GLUI_Spinner *pinceArctanK3SpinnerA = glui->add_spinner_to_panel(panel_PinceTypeA,
            "k3", GLUI_SPINNER_FLOAT, &(PinceArctanK3ValueA));
    pinceArctanK3SpinnerA -> set_float_limits(0.0, 1.0);

    //DiagNervText = glui->add_edittext_to_panel(panel_DiagNerv, "Diag Nervures:", GLUI_EDITTEXT_TEXT, "3 5 7");
    //DiagNervText -> set_w(80);
    GLUI_Rollout *panel_DiagNerv = glui->add_rollout("Diagonal nervures",false);
    //GLUI_Panel *panel_DiagNerv = glui->add_panel_to_panel(RolloutDiagonal, "Diagonal nerv");
    //DiagNervText = glui->add_edittext_to_panel(panel_DiagNerv, "D. Nerv:", GLUI_EDITTEXT_TEXT, "2 5 8 11 14 17 20 23 26 29 32");
    //DiagNervText -> set_w(150);
    glui->add_checkbox_to_panel(panel_DiagNerv, "Diagonal nervures", &DiagNervs);

    GLUI_Spinner *posDiagNerv1A = glui->add_spinner_to_panel(panel_DiagNerv,
            "Pos Int A", GLUI_SPINNER_FLOAT, &(PosDiagNerv2A));
    posDiagNerv1A -> set_float_limits(0.0, 100.0);
    GLUI_Spinner *posDiagNerv1F = glui->add_spinner_to_panel(panel_DiagNerv,
            "Pos Int F", GLUI_SPINNER_FLOAT, &(PosDiagNerv2F));
    posDiagNerv1F -> set_float_limits(0.0, 100.0);

    GLUI_Spinner *posDiagNerv2A = glui->add_spinner_to_panel(panel_DiagNerv,
            "Pos Ext A", GLUI_SPINNER_FLOAT, &(PosDiagNerv1A));
    posDiagNerv2A -> set_float_limits(0.0, 100.0);
    GLUI_Spinner *posDiagNerv2F = glui->add_spinner_to_panel(panel_DiagNerv,
            "Pos Ext F", GLUI_SPINNER_FLOAT, &(PosDiagNerv1F));
    posDiagNerv2F -> set_float_limits(0.0, 100.0);

    FicDiagNerv = glui->add_statictext_to_panel(panel_DiagNerv, "???");
    FicDiagNerv->set_text(fileNameDiagNerv);

    GLUI_Button *boutonLoadDiagNervs =
            glui->add_button_to_panel(panel_DiagNerv, "Load", 0, &readDiagNervs);
    boutonLoadDiagNervs->set_w(10);

    GLUI_Rollout *panel_Klapan = glui->add_rollout("Klapan",false);
    glui->add_checkbox_to_panel(panel_Klapan, "Klapan", &isKlapan);
    GLUI_Spinner *SpinPosKlapanInt = glui->add_spinner_to_panel(panel_Klapan,
                                    "Pos Klapan Int", GLUI_SPINNER_FLOAT, &(PosKlapanInt));
	SpinPosKlapanInt->set_float_limits(0.0, 100.0);
    GLUI_Spinner *SpinPosKlapanFin = glui->add_spinner_to_panel(panel_Klapan,
                                    "Fin Klapan", GLUI_SPINNER_FLOAT, &(PosKlapanFin));

	SpinPosKlapanFin->set_float_limits(0.0, 100.0);
    /*cote 2*/
    glui->add_column(true);
    /*No Nerv 2*/
    SpinNoNerv[1] = glui->add_spinner("No Nerv 2", GLUI_SPINNER_INT, &(NoNerv[1]));
    SpinNoNerv[1] -> set_int_limits(-1, F->m_nbProfils - 1);
    /*% Deb 2*/
    GLUI_Panel *panel_Deb2 = glui->add_panel(""); //, GLUI_PANEL_NONE);
    GLUI_Spinner *SpinDeb2 = glui->add_spinner_to_panel(
            panel_Deb2, "% Deb 2", GLUI_SPINNER_FLOAT, &(Deb[1]));
    SpinDeb2 -> set_float_limits(0.0, 100.0);
    GLUI_Listbox *ListDeb2 = glui->add_listbox_to_panel(panel_Deb2, "Face", &FaceDeb[1]);
    ListDeb2->add_item(1, "Ext");
    ListDeb2->add_item(2, "Int");
    /*% Deb 2*/
    //GLUI_Panel *panel_Fin2 = glui->add_panel(""); //, GLUI_PANEL_NONE);
    GLUI_Spinner *SpinFin2 = glui->add_spinner_to_panel(
            panel_Deb2, "% Fin 2", GLUI_SPINNER_FLOAT, &(Fin[1]));
    SpinFin2 -> set_float_limits(0.0, 100.0);
    GLUI_Listbox *ListFin2 = glui->add_listbox_to_panel(panel_Deb2, "Face", &FaceFin[1]);
    ListFin2->add_item(1, "Ext");
    ListFin2->add_item(2, "Int");


    GLUI_Rollout *panel_PosNerv2 =glui->add_rollout("Pos Nerv 2", false);
    glui->add_checkbox_to_panel(panel_PosNerv2, "Pos Nerv 2", &isPosNerv[1]);
    GLUI_Spinner *SpinPosNerv2 = glui->add_spinner_to_panel(panel_PosNerv2, "Pos", GLUI_SPINNER_FLOAT, &(PosNerv[1]));
    SpinPosNerv2 -> set_float_limits(0.0, 100.0);

    //ListFin2->set_int_val(1);
    /*Reperes 2*/
    /*
        GLUI_Panel *panel_Rep2 = glui->add_panel("");//, GLUI_PANEL_NONE);
        GLUI_Spinner *SpinPosRep2 = glui->add_spinner_to_panel(
                panel_Rep2,"% Rep", GLUI_SPINNER_FLOAT, &(PosRep[1]));
        SpinPosRep2 -> set_float_limits( 0.0, 100.0 );
        GLUI_Listbox *ListRep2 = glui->add_listbox_to_panel(panel_Rep2,"Face",&FaceRep[1]);
        ListRep2->add_item( 1, "Ext" ); ListRep2->add_item( 2, "Int" );
        glui->add_button_to_panel(panel_Rep2, "adder", 1, adder );
     */

    /*Marge de couture en cm*/
    GLUI_Spinner *SpinMarge2 = glui->add_spinner(
            "Marge 2", GLUI_SPINNER_FLOAT, &(Marge[1]));
    SpinMarge2 -> set_float_limits(0.0, 5.0);

    /*reperes suspentes*/
    glui->add_checkbox("reperes suspentes", &(ReperesSuspentes[1]));
    glui->add_checkbox("reperes profile", &(ReperesProfile[1]));
    /*panel pinces!!!*/
    GLUI_Panel *panel_Pince2 = glui->add_panel("Pinces 2"); //, GLUI_PANEL_NONE);
    /*position debut pincement au BA en %*/
    GLUI_Spinner *SpinPosPinceBA2 = glui->add_spinner_to_panel(panel_Pince2,
            "Pos BA (%)", GLUI_SPINNER_FLOAT, &(PosPinceBA[1]));
    SpinPosPinceBA2 -> set_float_limits(0.0, 40.0);
    /*Amplitude pincement au BA en %*/
    GLUI_Spinner *SpinAmpPinceBA2 = glui->add_spinner_to_panel(panel_Pince2,
            "Amp BA (%)", GLUI_SPINNER_FLOAT, &(AmpPinceBA[1]));
    SpinAmpPinceBA2 -> set_float_limits(0.0, 40.0);
    /*position debut pincement au BF en %*/
    GLUI_Spinner *SpinPosPinceBF2 = glui->add_spinner_to_panel(panel_Pince2,
            "Pos BF (%)", GLUI_SPINNER_FLOAT, &(PosPinceBF[1]));
    SpinPosPinceBF2 -> set_float_limits(0.0, 40.0);
    /*Amplitude pincement au BF en %*/
    GLUI_Spinner *SpinAmpPinceBF2 = glui->add_spinner_to_panel(panel_Pince2,
            "Amp BF (%)", GLUI_SPINNER_FLOAT, &(AmpPinceBF[1]));
    SpinAmpPinceBF2 -> set_float_limits(0.0, 40.0);

    GLUI_Rollout *panel_PinceTypeF = glui->add_rollout("Pince type F(Hvost)", false);
    //GLUI_Panel *panel_PinceTypeF = glui->add_panel_to_panel(RolloutPinceTypeF, "Pince type F(Hvost)");
    GLUI_Spinner *pinceRadiusAlgKHvostSpinner = glui->add_spinner_to_panel(panel_PinceTypeF,
            "k (d.width)", GLUI_SPINNER_FLOAT, &(PinceRadiusAlgKHvost));
    pinceRadiusAlgKHvostSpinner -> set_float_limits(0.0001, 5.0);

    glui->add_checkbox_to_panel(panel_PinceTypeF, "Set = Amp to Nos", &PinceHvostEqualAmp, 0, &modifPinceHvostRadio);
    //glui->add_checkbox_to_panel(panel_PinceTypeF, "Amp 0", &PinceHvost0Amp, 0,  &modifPinceNosRadio);
    GLUI_RadioGroup *radio_pinceHvost =
            glui->add_radiogroup_to_panel(panel_PinceTypeF, &PinceHvostRadio, 0, &modifPinceHvostRadio);
    glui->add_radiobutton_to_group(radio_pinceHvost, "Radius auto alg");
    glui->add_radiobutton_to_group(radio_pinceHvost, "Function");
    glui->add_checkbox_to_panel(panel_PinceTypeF, "Power", &PincePowerF);
    GLUI_Spinner *pincePowerSpinnerF = glui->add_spinner_to_panel(panel_PinceTypeF,
            "power", GLUI_SPINNER_FLOAT, &(PincePowerValueF));
    pincePowerSpinnerF -> set_float_limits(0.0, 5.0);
    glui->add_checkbox_to_panel(panel_PinceTypeF, "Arctan", &PinceArctanF);
    GLUI_Spinner *pinceArctanK1SpinnerF = glui->add_spinner_to_panel(panel_PinceTypeF,
            "k1", GLUI_SPINNER_FLOAT, &(PinceArctanK1ValueF));
    pinceArctanK1SpinnerF -> set_float_limits(0.0, 10.0);
    GLUI_Spinner *pinceArctanK2SpinnerF = glui->add_spinner_to_panel(panel_PinceTypeF,
            "k2", GLUI_SPINNER_FLOAT, &(PinceArctanK2ValueF));
    pinceArctanK2SpinnerF -> set_float_limits(0.0, 20.0);
    GLUI_Spinner *pinceArctanK3SpinnerF = glui->add_spinner_to_panel(panel_PinceTypeF,
            "k3", GLUI_SPINNER_FLOAT, &(PinceArctanK3ValueF));
    pinceArctanK3SpinnerF -> set_float_limits(0.0, 1.0);

    /*ventilation*/
    GLUI_Rollout *panel_Vent =	glui->add_rollout("Ventilation holes",false);
    //GLUI_Panel *panel_Vent = glui->add_panel_to_panel(RolloutVent, "Vent holes");
    glui->add_checkbox_to_panel(panel_Vent, "Vent holes", &VentHoles);
    glui->add_checkbox_to_panel(panel_Vent, "double vents", &VentHolesDouble);
    glui->add_checkbox_to_panel(panel_Vent, "klapans", &LayoutKlapans);
    GLUI_Spinner *SpinVentDeb = glui->add_spinner_to_panel(panel_Vent, "Deb", GLUI_SPINNER_FLOAT, &(VentHolesDeb));
    SpinVentDeb -> set_float_limits(0.0, 100.0);
    GLUI_Spinner *SpinVentFin = glui->add_spinner_to_panel(panel_Vent, "Fin", GLUI_SPINNER_FLOAT, &(VentHolesFin));
    SpinVentFin -> set_float_limits(0.0, 100.0);

    //VentHolesNervsText = glui->add_edittext_to_panel(panel_Vent, "Nrvs:", GLUI_EDITTEXT_TEXT, "11111111111111111111111111111111111111111111111111111111111111111111111111111111111");
    //VentHolesNervsText -> set_w(350);
    //glui->add_checkbox_to_panel(panel_Vent, "central nerv", &VentCentralNerv);

    //printf ("\n before new VentHoles");
    FicVentHoles = glui->add_statictext_to_panel(panel_Vent, "???");
    FicVentHoles->set_text(fileNameVentHoles);
    //glui->add_column_to_panel(panel_Vent, false);
    
    //BlankST = glui->add_statictext_to_panel(panel_RepPoints, "");
    //BlankST->set_text("");
    
    GLUI_Button *boutonLoadVentHoles =
            glui->add_button_to_panel(panel_Vent, "Load", 0, &readVentHoles);
    boutonLoadVentHoles->set_w(10);
    
    /*3eme colonne*/
    glui->add_column(true);

    GLUI_Rollout *panel_marges = glui->add_rollout("Marges",false);

    /*Marge de couture Deb en cm*/
    GLUI_Spinner *SpinMargeDeb = glui->add_spinner_to_panel(panel_marges,
            "Marge Deb", GLUI_SPINNER_FLOAT, &(MargeDeb));
    SpinMargeDeb -> set_float_limits(0.0, 5.0);
    /*Marge de couture Fin en cm*/
    GLUI_Spinner *SpinMargeFin = glui->add_spinner_to_panel(panel_marges,
            "Marge Fin", GLUI_SPINNER_FLOAT, &(MargeFin));
    SpinMargeFin -> set_float_limits(0.0, 5.0);

    GLUI_Spinner *SpinExtMargeFin = glui->add_spinner_to_panel(panel_marges,
            "Marge Fin Ext", GLUI_SPINNER_FLOAT, &(margeFinExt));
    SpinExtMargeFin -> set_float_limits(0.0, 20.0);

    GLUI_Spinner *SpinIntMargeFin = glui->add_spinner_to_panel(panel_marges,
            "Marge Fin Int", GLUI_SPINNER_FLOAT, &(margeFinInt));
    SpinIntMargeFin -> set_float_limits(0.0, 20.0);

    GLUI_Spinner *SpinMargeFinNerv = glui->add_spinner_to_panel(panel_marges,
            "Marge Fin Nervs", GLUI_SPINNER_FLOAT, &(margeFinNerv));
    SpinMargeFinNerv -> set_float_limits(0.0, 20.0);

    GLUI_Spinner *SpinMargeFinDiagNerv = glui->add_spinner_to_panel(panel_marges,
            "Marge Fin Diag Nervs", GLUI_SPINNER_FLOAT, &(margeFinDiagNerv));
    SpinMargeFinDiagNerv -> set_float_limits(0.0, 20.0);

    /*reperes de couture*/
    /*
            glui->add_checkbox("reperes couture",&ReperesCoutures);
     */
    /*numerotation*/
    glui->add_checkbox("ventilation", &Ventilation, 0, &modifVentilation);

    GLUI_Rollout *RolloutNum=glui->add_rollout("Numerotation",false);
    //GLUI_Panel *panelNum = glui->add_panel_to_panel(RolloutNum, "Numerotation");

    glui->add_checkbox_to_panel(RolloutNum, "numerotation", &Numerotation);
    NumText = glui->add_edittext_to_panel(RolloutNum, "Text", GLUI_EDITTEXT_TEXT, "TEXT");
    NumText -> set_w(80);
    GLUI_Spinner *SpinTextX = glui->add_spinner_to_panel(RolloutNum, "Text X", GLUI_SPINNER_FLOAT, &(textX));
    SpinTextX -> set_float_limits(-2.0, 2.0);
    GLUI_Spinner *SpinTextY = glui->add_spinner_to_panel(RolloutNum, "Text Y", GLUI_SPINNER_FLOAT, &(textY));
    SpinTextY -> set_float_limits(-2.0, 2.0);

    /*mode de pincement des pinces !*/
/*    GLUI_Spinner *SpinModePinces = glui->add_spinner("mode pinces", GLUI_SPINNER_FLOAT, &(modePinces));
    SpinModePinces -> set_float_limits(0.0, 3.0);
    GLUI_Spinner *SpinModePinces2 = glui->add_spinner("mode pinces2", GLUI_SPINNER_FLOAT, &(modePinces2));
    SpinModePinces2 -> set_float_limits(0.01, 10.0);*/

    /*nom fichier forme*/
    GLUI_Panel *panel_forme = glui->add_panel("Form", GLUI_PANEL_EMBOSSED);
    FicForm = glui->add_statictext_to_panel(panel_forme, "???");
    FicForm->set_text(fileNameForm);
    glui->add_column_to_panel(panel_forme, false);
    GLUI_Button *bouton =
            glui->add_button_to_panel(panel_forme, "Load", 0, &readForm);
    bouton->set_w(10);
    glui->add_column_to_panel(panel_forme, false);
    GLUI_Button *bouton2 =
            glui->add_button_to_panel(panel_forme, "Load2", 0, &readForm2);
    bouton2->set_w(10);
    /*choix projection orthogonale ou perspective*/
    GLUI_Panel *panel_project = glui->add_panel("Project", GLUI_PANEL_EMBOSSED);
    FicProject = glui->add_statictext_to_panel(panel_project, "???");
    FicProject->set_text(fileNameProject);
    glui->add_column_to_panel(panel_project, false);
    GLUI_Button *boutonLoadProject =
            glui->add_button_to_panel(panel_project, "Load", 0, &readProject);
    boutonLoadProject -> set_w(10);
    glui->add_column_to_panel(panel_project, false);
    GLUI_Button *boutonsaveProject = 
			glui->add_button_to_panel(panel_project, "Save", 0, &saveProject);
	boutonsaveProject -> set_w(10);
    //
    GLUI_Panel *panel_proj = glui->add_panel("Projection");
    GLUI_RadioGroup *radio_proj =
            glui->add_radiogroup_to_panel(panel_proj, &ProjOrthoPers, 0, &modifProjection3d);
    glui->add_radiobutton_to_group(radio_proj, "Orthogonale");
    glui->add_radiobutton_to_group(radio_proj, "Perspective");


    GLUI_Rollout *RolloutPts=glui->add_rollout("Points", false);
	// GLUI_Panel *panelNum = glui->add_panel_to_panel(RolloutNum, "Numerotation");

    GLUI_Panel *panel_RepPoints = glui->add_panel_to_panel(RolloutPts, "Rep pts");
    glui->add_checkbox_to_panel(panel_RepPoints, "from file", &ReperPointsFromFile, 0, &modifRepPts);
    glui->add_checkbox_to_panel(RolloutPts, "pts suspentes", &ViewPtsSuspentes, 0, &modifViewSymetrique);
    glui->add_checkbox_to_panel(RolloutPts, "plotter format", &ReperPointsPlotterFormat, 0, &modifRepPts);
    GLUI_Spinner *SpinXMashtab = glui->add_spinner_to_panel(RolloutPts, "X mashtab", GLUI_SPINNER_FLOAT, &(XMashtab));
    SpinXMashtab -> set_float_limits(0.01, 100);


    FicRepPoints = glui->add_statictext_to_panel(panel_RepPoints, "???");
    FicRepPoints->set_text(fileNameRepPoints);
    glui->add_column_to_panel(panel_RepPoints, false);
    BlankST = glui->add_statictext_to_panel(panel_RepPoints, "");
    BlankST->set_text("");
    GLUI_Button *boutonLoadRepPoints =
            glui->add_button_to_panel(panel_RepPoints, "Load", 0, &readRepPoints);
    boutonLoadRepPoints->set_w(10);


    glui->add_checkbox("symetrique", &ViewSymetrique, 0, &modifViewSymetrique);

    /*choix visu des points de suspentage*/

    glui->add_checkbox("make pinces", &GoPince, 0, &modifViewSymetrique);
    /* bouton apply */

    GLUI_Spinner *SpinCoeff = glui->add_spinner("Coeff", GLUI_SPINNER_FLOAT, &(coeffMult));
    SpinCoeff -> set_float_limits(0.0, 2.0);

    glui->add_button("apply", 0, &apply);
    //glui->add_button("Layout1 to DXF", 0, &applyMagic);
    //GLUI_Panel *panel_Layout2 = glui->add_panel("Layout");
    
	GLUI_Rollout *panel_Layout2 = glui->add_rollout("Layout",false);
    glui->add_button_to_panel(panel_Layout2, "Layout to DXF", 0, &applyMagic2);
    //glui->add_checkbox_to_panel(panel_Layout2, "sym layout", &LayoutSymetrique, 0, &modifLayoutSymetrique);
    glui->add_checkbox_to_panel(panel_Layout2, "ventilation", &VentilationLayout, 0, &modifVentilationLayout);
	glui->add_checkbox_to_panel(panel_Layout2, "correct rep points", &CorrectRepPoints, 0, &modifVentilationLayout);

	//--------------- Design GUI -------------------
	GLUI_Rollout *panelDesign = glui->add_rollout("Design", false);
    glui->add_checkbox_to_panel(panelDesign, "Layout with design", &layoutWithDesign, 0, &modifLayoutWithDesign);

	GLUI_Button *btnLoadDesignExt = glui->add_button_to_panel(panelDesign, "Load design Ext", 0, &readFileDesignExt);
    btnLoadDesignExt->set_w(10);

	GLUI_Button *btnLoadDesignInt = glui->add_button_to_panel(panelDesign, "Load design Int", 0, &readFileDesignInt);
    btnLoadDesignInt->set_w(10);

    glui->add_column_to_panel(panelDesign, false);

	glui->add_statictext_to_panel(panelDesign, "");

	// load default design
	strcpy(fileNameDesignExt, "f17ext.wdn");
    FicDesignExt = glui->add_statictext_to_panel(panelDesign, "???");
    FicDesignExt->set_text(fileNameDesignExt);
	kiteDesignExt = readKiteDesignFromFile(fileNameDesignExt);

	strcpy(fileNameDesignInt, "f17int.wdn");
    FicDesignInt = glui->add_statictext_to_panel(panelDesign, "???");
    FicDesignInt->set_text(fileNameDesignInt);
	kiteDesignInt = readKiteDesignFromFile(fileNameDesignInt);

	//----------------------------------------------



    GLUI_Spinner *SpinTochnostLayout2 = glui->add_spinner_to_panel(panel_Layout2, "Accuracy", GLUI_SPINNER_FLOAT, &(tochnostLayout2));
    SpinTochnostLayout2 -> set_float_limits(0.0, 3.0);
    /* bouton save au format Texte */
    glui->add_button("save TXT", 0, &save);
    /* bouton save au format DXF */
    //glui->add_button("save DXF", 1, &save);
    /* bouton save au format Poly DXF */
    glui->add_button("save Poly DXF", 2, &save);
    glui->add_button( "3D->DXF", 0, &saveFichier3dDXF );

    //glui->add_button("Export to WPA", 0, &ExportWpa);
    /* bouton quit */
    glui->add_button("quit", 0, &quit);
    if (TEST_BUTTON) glui->add_button("Test", 0, &Test);
    /* Link windows to GLUI, and register idle callback */
    glui->set_main_gfx_window(windowPatron);
    /* We register the idle callback with GLUI, not with GLUT */
    GLUI_Master.set_glutIdleFunc(NULL);
    /*init Value boite de dialogue*/
    initValueDialogue();
    glui->sync_live();
    /*premiers calcs et trace des axes*/
    //calcVue3dEtPatron();
    apply(0);
    AxeSel = Axe3d;
    display();

    applyMagic2(0);
    /*boucle glut*/
    glutMainLoop();
    return 0;
}
