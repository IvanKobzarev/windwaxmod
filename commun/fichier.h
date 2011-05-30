#pragma once

#include "form.h"
#include "ballon.h"
#include "patternsproject.h"
#include "design.h"


typedef struct InfoForm TInfoForm;

struct InfoForm
{
	double surface, envergure;
	double surfaceProj, envergureProj;
	double allongement, allongementProj;
	double cordeMin, cordeMax;

	double largMin, largMax;

};

void readFichierProfil(const char* NomProf, Matrix** extrados, Matrix** intrados);

Ballonement* readBallonementFromFile(char* NomFic);

Form* readFichierForm(char* NomFic);

Ballonement* readBallonementFromFile(char* NomFic);

Form* readFichierForm2(char* NomFic);

void writeFichierForm(char *fileName, Form *f);

void writeFichierForm2(char *fileName, Form *f);

bool TrouveMotDansFichierTexte(FILE* fid, char* Mot);


void writeFichierFGen(char *fileName, Form *foil);

void calcInfoForm( Form* F, TInfoForm* info );

void AfficheInfoForm( TInfoForm info );

WindPatternsProject* readWindPatternsProject(char* NomFic);

int* readFichierVentHoles(char* NomFic, int* quant, int* central);

int* readFichierDiagNervs(char* NomFic, int* quant);

void writeWindPatternsProject(char *fileName, WindPatternsProject *wpp);

KiteDesign* readKiteDesignFromFile(const char* FilePath);