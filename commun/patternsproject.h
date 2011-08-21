#pragma once

#include "matrice.h"
#include "fichier.h"
#include "design.h"

//#include <afx.h>		//class CString
//#include <afxdlgs.h>	//class CFileDialog


class WindPatternsProject
{
    public:
        WindPatternsProject() {
			debug = false;
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
        float PincePowerValueA, PincePowerValueF,PinceArctanK1ValueA, PinceArctanK2ValueA, PinceArctanK3ValueA, PinceArctanK1ValueF, PinceArctanK2ValueF, PinceArctanK3ValueF;
        float Deb[2];
        float Fin[2];

        float PosPinceBA[2];
        float PosPinceBF[2];
        float AmpPinceBA[2];
        float AmpPinceBF[2];
        int LayoutSymetrique;
        Matrix *IntProfCent, *ExtProfCent, *IntProfBout ,*ExtProfBout;
        int PinceNosRadio;
        int PinceHvostRadio;
        float PinceRadiusAlgKNos;
        float PinceRadiusAlgKHvost;
        //CString* diagNervText;
        Form* Form;
        int quantDiag;
        int *noNervD;
        int windowPatron;
        float PosDiagNerv2A, PosDiagNerv2F, PosDiagNerv1A, PosDiagNerv1F;
        float PosKlapanFin;
        float Marge[2];
        float MargeDeb, MargeFin;
        int ReperPointsFromFile;
        Matrix** ReperPoints;
        int ReperesSuspentes[2];
        int ReperesProfile[2];
        float XMashtab;
        int PinceNosEqualAmp;
        int PinceHvostEqualAmp;
        int PinceNos0Amp;
        int PinceHvost0Amp;
        int LayoutKlapans;
        int Ventilation;
        int VentilationLayout;
        int ReperPointsPlotterFormat;
		int CorrectRepPoints;
        int layoutWithDesign;

        float textX;
        float textY;
        float tochnostLayout2;
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
        char fileNameRepPoints[255];
        char fileNameVentHoles[255];
        char fileNameDiagNerv[255];
        char fileNameForm[255];
		char ballonementPath[255];
        char logFileName[255];
		bool debug;
		Ballonement* ballonement;
        KiteDesign* kiteDesignExt;
        KiteDesign* kiteDesignInt;
        void print();

};
