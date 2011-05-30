#pragma once

#include "matrice.h"
#include "fichier.h"
#include "patternsproject.h"

//#include <afx.h>		//class CString
//#include <afxdlgs.h>	//class CFileDialog


void calcForm3DBallonement
				(WindPatternsProject* gfd, Form *forme, int isPercent, double percent,
				   Matrix *ExtProfCent, Matrix *IntProfCent,
				   Matrix *ExtProfBout, Matrix *IntProfBout,
				   Matrix **XExt, Matrix **YExt, Matrix **ZExt,
				   Matrix **XInt, Matrix **YInt, Matrix **ZInt);
