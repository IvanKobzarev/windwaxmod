#pragma once

#include "matrice.h"
#include "fichier.h"

//#include <afx.h>		//class CString
//#include <afxdlgs.h>	//class CFileDialog


class WindPatternsProject
{
    public:
        WindPatternsProject() {
			debug = true;
            for (int i = 0; i < 2; i++) {
                    PosPinceBA[i] = 0.0f;
                    PosPinceBF[i] = 0.0f;
                    AmpPinceBA[i] = 0.0f;
                    AmpPinceBF[i] = 0.0f;
            }
        }
        virtual ~WindPatternsProject();
        char name[255];
        char fromFile[255];
        int isPinces;
        int PincePowerA, PincePowerF, PinceArctanA, PinceArctanF;
        double PincePowerValueA, PincePowerValueF,PinceArctanK1ValueA, PinceArctanK2ValueA, PinceArctanK3ValueA, PinceArctanK1ValueF, PinceArctanK2ValueF, PinceArctanK3ValueF;
        float Deb[2];
        float Fin[2];

        float PosPinceBA[2];
        float PosPinceBF[2];
        float AmpPinceBA[2];
        float AmpPinceBF[2];
        int RaskladSymetrique;
        Matrice *IntProfCent, *ExtProfCent, *IntProfBout ,*ExtProfBout;
        int PinceNosRadio;
        int PinceHvostRadio;
        float PinceRadiusAlgKNos;
        float PinceRadiusAlgKHvost;
        //CString* diagNervText;
        Forme* Forme;
        int quantDiag;
        int *noNervD;
        int FenetrePatron;
        float PosDiagNerv2A, PosDiagNerv2F, PosDiagNerv1A, PosDiagNerv1F;
        float PosKlapanFin;
        float Marge[2];
        float MargeDeb, MargeFin;
        int ReperPointsFromFile;
        Matrice** ReperPoints;
        int ReperesSuspentes[2];
        int ReperesProfile[2];
        float XMashtab;
        int PinceNosEqualAmp;
        int PinceHvostEqualAmp;
        int PinceNos0Amp;
        int PinceHvost0Amp;
        int RaskladKlapans;
        int Ventilation;
        int VentilationLayout;
        int ReperPointsPlotterFormat;
		int CorrectRepPoints;

        float textX;
        float textY;
        float tochnostRasklad2;
        float modePinces;
        int  VentHoles;
        int RepPoints;
        int  DiagNervs;
        int  VentHolesDouble;
        int VentCentralNerv;
        int* noNervVH;
        float VentHolesDeb, VentHolesFin;
        float margeFinExt;
        float margeFinInt;
        float margeFinNerv;
        float margeFinDiagNerv;
        int quantVH;
        int* isPosNerv;
        float* PosNerv;
        char NomFichierRepPoints[255];
        char NomFichierVentHoles[255];
        char NomFichierDiagNerv[255];
        char NomFichierForme[255];
        char logFileName[255];
		bool debug;

        void print();

};
