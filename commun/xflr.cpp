

/*void Form::SerializeFoil (QDataStream &ar, Matrix* ExtProf, Matrix* IntProf) { 
	CFoil* foil = new CFoil();
	foil->InitFromWindWax(ExtProf, IntProf);
	foil->m_FoilName="DummyFoil";
	foil->Serialize(ar, true);
}*/

/*
void calcForm3D(Form *forme, int isPercent, double percent,
				   Matrix *ExtProfCent, Matrix *IntProfCent,
				   Matrix *ExtProfBout, Matrix *IntProfBout,
				   Matrix **XExt, Matrix **YExt, Matrix **ZExt,
				   Matrix **XInt, Matrix **YInt, Matrix **ZInt)

{
    Matrix *ExtProfCentN, *ExtProfBoutN;
	double LongNerv, EpaiRel, xp,yp, xo,yo,zo, a,v,m;
	double EpaiRelProfCent, EpaiRelProfBout;
	double coeffx, coeffyCent, coeffyBout;
	int i,j;
    bool isCenterPanel = (1 & forme->NbCaiss);

	EpaiRelProfCent = EpaisseurRelative(ExtProfCent, IntProfCent);
	EpaiRelProfBout = EpaisseurRelative(ExtProfBout, IntProfBout);
	
	*XExt = Zeros(forme->m_nbProfils, ExtProfCent->GetLignes());
	*YExt = Zeros(forme->m_nbProfils, ExtProfCent->GetLignes());
	*ZExt = Zeros(forme->m_nbProfils, ExtProfCent->GetLignes());

	*XInt = Zeros(forme->m_nbProfils, IntProfCent->GetLignes());
	*YInt = Zeros(forme->m_nbProfils, IntProfCent->GetLignes());
	*ZInt = Zeros(forme->m_nbProfils, IntProfCent->GetLignes());

	for (i=0; i<forme->m_nbProfils; i++)
	{
		LongNerv = forme->m_pProfils[i]->m_fLength;
		EpaiRel = forme->m_pProfils[i]->m_fWidth;
		xo = forme->m_pProfils[i]->m_fNezX;
		yo = forme->m_pProfils[i]->m_fNezY;
		zo = forme->m_pProfils[i]->m_fNezZ;
		a = forme->m_pProfils[i]->m_fInclin;
        if ((isCenterPanel) && (i == 0)) a = forme->m_pProfils[0]->m_fInclin - 0.33333333f * fabs(forme->m_pProfils[0]->m_fInclin - forme->m_pProfils[1]->m_fInclin);
                
		v = forme->m_pProfils[i]->m_fWash;
		m = forme->m_pProfils[i]->m_fMorph;

		coeffx = LongNerv/100.0f;
		coeffyCent = LongNerv*EpaiRel/(EpaiRelProfCent*100.0f);
		coeffyBout = LongNerv*EpaiRel/(EpaiRelProfBout*100.0f);

		for (j=0; j<IntProfCent->GetLignes(); j++)
		{
			xp = IntProfCent->Element(j,0)*coeffx;
			yp = IntProfCent->Element(j,1)*coeffyCent*m
				+ IntProfBout->Element(j,1)*coeffyBout*(1.0f-m);
			(*XInt)->SetElement(i,j, xo+(yp-xp*(double)sin(v))*(double)cos(a));
			(*YInt)->SetElement(i,j, yo+(yp-xp*(double)sin(v))*(double)sin(a));
			(*ZInt)->SetElement(i,j, zo-(xp*(double)cos(v)));

		}

		for (j=0; j<ExtProfCent->GetLignes(); j++)
		{
            xp = ExtProfCent->Element(j,0)*coeffx;
            yp = ExtProfCent->Element(j,1)*coeffyCent*m
		                        + ExtProfBout->Element(j,1)*coeffyBout*(1.0f-m);

			(*XExt)->SetElement(i,j, xo+(yp-xp*(double)sin(v))*(double)cos(a));
			(*YExt)->SetElement(i,j, yo+(yp-xp*(double)sin(v))*(double)sin(a));
			(*ZExt)->SetElement(i,j, zo-(xp*(double)cos(v)));
  		}
	}
}
*/

/*
void Form::WritePolars (QDataStream &ar) {
		printf ("\n WritePolars..");
	    ar << 100003;
        int m_oaFoilSize=m_nbProfils;
        ar << (int)m_oaFoilSize;

		//SerializeFoil(ar, ExtProfCent, IntProfCent);
        //for (int i=0; i<m_oaFoilSize; i++){
                //pFoil = (CFoil*)m_oaFoil.GetAt(i);
                //pFoil->Serialize(ar);
                //SerializeFoil(ar);
		//}
		double LongNerv, EpaiRel, xp,yp, xo,yo,zo, a,v,m, w, kw;
		double EpaiRelProfCent, EpaiRelProfBout;
		double coeffx, coeffyCent, coeffyBout;
		int i, j;
        //bool isCenterPanel = (1 & forme->NbCaiss);

		EpaiRelProfCent = EpaisseurRelative(ExtProfCent, IntProfCent);
		EpaiRelProfBout = EpaisseurRelative(ExtProfBout, IntProfBout);
		char foilName[50];
		printf ("\n Form->m_nbProfils=%d [0 -> %d]", m_nbProfils, (m_nbProfils-1));
		double w0 = m_pProfils[0]->m_fWidth;
		for (i=0; i<m_nbProfils; i++)
		{
			Matrix *ExtProf, *IntProf;
			ExtProf = Zeros(ExtProfCent->GetLignes(), 2);
			IntProf = Zeros(IntProfCent->GetLignes(), 2);

			//LongNerv = forme->m_pProfils[i]->m_fLength;
			w = m_pProfils[i]->m_fWidth;
			
			kw = w / w0;
			if (i == (m_nbProfils-1))  kw=0.1;
			//xo = forme->m_pProfils[i]->m_fNezX;
			//yo = forme->m_pProfils[i]->m_fNezY;
			//zo = forme->m_pProfils[i]->m_fNezZ;
			//a = forme->m_pProfils[i]->m_fInclin;
			//if ((isCenterPanel) && (i == 0)) a = forme->m_pProfils[0]->m_fInclin - 0.33333333f * fabs(forme->m_pProfils[0]->m_fInclin - forme->m_pProfils[1]->m_fInclin);
	                
			//v = forme->m_pProfils[i]->m_fWash;
			m = m_pProfils[i]->m_fMorph;
			//if (DEBUG) printf ("\n%3d -> a=%f v=%f m=%f", i, a * 180.0f/pi, v, m);
			coeffx = 1.0f;//LongNerv/100.0f;
			coeffyCent = 1.0f;// LongNerv*EpaiRel/(EpaiRelProfCent*100.0f);
			coeffyBout = 1.0f;// LongNerv*EpaiRel/(EpaiRelProfBout*100.0f);
			//printf ("\n i=%d kw=%f", i , kw);
			for (j=0; j<IntProfCent->GetLignes(); j++)
			{
				xp = IntProfCent->Element(j,0)*coeffx;
				yp = IntProfCent->Element(j,1)*coeffyCent*m
					+ IntProfBout->Element(j,1)*coeffyBout*(1.0f-m);
				if (kw != 1) yp = yp * kw;
				if (i == (m_nbProfils-1)) printf (" %f", yp);
				IntProf->SetElement(j,0,xp);
				IntProf->SetElement(j,1,yp);
			}

			for (j=0; j<ExtProfCent->GetLignes(); j++)
			{
				xp = ExtProfCent->Element(j,0)*coeffx;
				yp = ExtProfCent->Element(j,1)*coeffyCent*m
						+ ExtProfBout->Element(j,1)*coeffyBout*(1.0f-m);
				if (kw != 1) yp = yp * kw;
				ExtProf->SetElement(j,0,xp);
				ExtProf->SetElement(j,1,yp);
			}

			CFoil* foil = new CFoil();
			foil->InitFromWindWax(ExtProf, IntProf);
			sprintf(foilName, "P%03d", i);
			foil->m_FoilName = foilName;
			printf ("\n Serialize: %s", foilName);
			foil->Serialize(ar, true, 5);

			delete (foil);
			delete (ExtProf);
			delete (IntProf);
		}

        //then write polars
        ar << 0;
        ar << 0;
		printf ("\n ...WritePolars..");
}
*/

/*
bool writeFichierWpa(char *fileName, Form *forme)
{
	//CWaitCursor wait;
	CFileException fe;
    string nomFichierStr = fileName;
	QFile fp(fileName);

	if (!fp.open(QIODevice::WriteOnly))
	{
		printf ("\nCould not open the file for writing");
		return false;
	}

	QDataStream ar(&fp);//, CArchive::store);
	ar.setByteOrder(QDataStream::LittleEndian);
	printf ("\n writeFichierWpa");
	forme -> SerializeToWpa(ar);
	printf ("\n...writeFichierWpa");
	//ar.Close();
	fp.close();
	return true;
}*/



/* void Form::SerializeWing (QDataStream &ar) {
	printf ("\nin Form::SerializeWing()");

	CWing* wing = new CWing ();

	wing -> m_WingName = "Kite";
	char foilNameL[50];
	char foilNameR[50];

	double LongNerv, EpaiRel, xp,yp, xo,yo,zo, a,v,m;



        int isCenterPanel = (1 & NbCaiss);
		double positionSum = 0;
		double	xprev = m_pProfils[0]->m_fNezX;
		double	yprev = m_pProfils[0]->m_fNezY;

		wing->m_NPanel = isCenterPanel + m_nbProfils-1;	
		if (isCenterPanel) {
			LongNerv = m_pProfils[0]->m_fLength;
			EpaiRel = m_pProfils[0]->m_fWidth;

			xo = m_pProfils[0]->m_fNezX;
			yo = m_pProfils[0]->m_fNezY;
			zo = m_pProfils[0]->m_fNezZ;

			wing->m_RFoil[0]= "P000";
			wing->m_LFoil[0]= "P000";

			wing->m_TDihedral[0]=0;
			wing->m_TChord[0] = LongNerv;
			wing->m_TPos[0] = 0.0;

			wing->m_TOffset[0]=0.0;
			wing->m_TTwist[0]=0.0;

			wing->m_NXPanels[0]=2;
			wing->m_NYPanels[0]=2;
			wing->m_XPanelDist[0]=1;
			wing->m_YPanelDist[0]=0;

			xprev = 0;
		}

		double x1, y1, x2, y2;

		for (int i = 0; i < m_nbProfils; i++) {
			
			LongNerv = m_pProfils[i]->m_fLength;
			EpaiRel = m_pProfils[i]->m_fWidth;
			
			xo = m_pProfils[i]->m_fNezX;
			yo = m_pProfils[i]->m_fNezY;
			zo = m_pProfils[i]->m_fNezZ;
			positionSum += dist2d(xo, yo, xprev, yprev);
			xprev = xo;
			yprev = yo;

			wing->m_TDihedral[i+isCenterPanel]=0;
			v = m_pProfils[i]->m_fWash;
			printf ("\n vrillage=%f", (v*180/CHISLOPI));
			if (i < m_nbProfils-1) {
				x1 = m_pProfils[i]->m_fNezX;
				y1 = m_pProfils[i]->m_fNezY;

				x2 = m_pProfils[i+1]->m_fNezX;
				y2 = m_pProfils[i+1]->m_fNezY;
				wing->m_TDihedral[i+isCenterPanel] = -180/CHISLOPI*atan(abs(y1-y2)/abs(x1-x2));
			} 

			printf ("\ni=%d panel %d", i, (i+isCenterPanel));
			sprintf(foilNameL, "P%03d", i);
			sprintf(foilNameR, "P%03d", i);
			wing->m_RFoil[i+isCenterPanel]= foilNameL;
			wing->m_LFoil[i+isCenterPanel]= foilNameR;

	

			wing->m_TChord[i+isCenterPanel] = LongNerv;
			wing->m_TPos[i+isCenterPanel] = positionSum;
			//wing->m_Twist[i+isCenterPanel] = ;

			wing->m_TOffset[i+isCenterPanel]=m_pProfils[0]->m_fNezZ - zo;
			wing->m_TTwist[i+isCenterPanel]=v*180/CHISLOPI;

			wing->m_NXPanels[i+isCenterPanel]=2;
			wing->m_NYPanels[i+isCenterPanel]=2;
			wing->m_XPanelDist[i+isCenterPanel]=1;
			wing->m_YPanelDist[i+isCenterPanel]=0;
		}

		printf ("\n ...for");

		wing->m_bVLMAutoMesh = 1;
		wing->m_bSymetric = 1;
		printf ("\n in Form::SerializeWing() ComputeGeometry()");
		//wing->ComputeGeometry();
		
		wing -> SerializeWing(ar, true, 5);
 		printf ("\n...SerializeWing()");
}
*/


/* void Form::SerializeToWpa(QDataStream &ar) 
{
	printf ("\nSerializeToWpa()");
	int m_LengthUnit;
	int m_AreaUnit;
	int m_WeightUnit;
	int m_SpeedUnit;
	int m_ForceUnit;
	int m_MomentUnit;
	
	m_LengthUnit  = 3;
	m_AreaUnit    = 3;
	m_WeightUnit  = 1;
	m_SpeedUnit   = 0;
	m_ForceUnit   = 0;
	m_MomentUnit  = 0;

	int m_Type=0;
	int m_UnitType=2;//1= International, 2= English
	int m_AnalysisType=1; //0=LLT;1=VLM;2=Panel
	int m_RefAreaType=1; //1=classic or 2=projected on x-yplane

	bool m_bVLM1=true; //true if Classic, false if Quendez
	bool m_bThinSurfaces=true;//true if Plane Panel calcation on middle surface, false if on top & bottom
	bool m_bWakeRollUp=true;//true if wake roll up is to be taken into account in calcation
	bool m_bTiltedGeom=true;//true if calcation is performed on the tilted geometry, at alpha=0.0
	bool m_bViscous=true;
    bool m_bGround=true;

	double m_QInf = 9.0, m_Weight = 1.0, m_Alpha = 2, m_XCmRef = 0.25;

	double m_CoG_x=0.25;
	double m_CoG_y=0;
	double m_CoG_z=0;

	double m_Beta = 0;
	double m_Density = 1.225, m_Viscosity = 1.5e-5;
	double m_WingLoad = 1;
	double m_Height = 0;

    int wingSize = 1;
    int polarSize = 0;
    int bodySize = 0;
    int planeSize = 0;

	ar << 100013;
	ar << m_LengthUnit;
	ar << m_AreaUnit;
	ar << m_WeightUnit;
	ar << m_SpeedUnit;
	ar << m_ForceUnit;
	ar << m_MomentUnit;

	ar << m_Type;
	ar << (float)m_Weight;
	ar << (float)m_QInf;
	ar << (float)m_CoG_x;
	ar << (float)m_CoG_y;
	ar << (float)m_CoG_z;

	ar << (float)m_Density;
	ar << (float)m_Viscosity;
	ar << (float)m_Alpha;
	ar << (float)m_Beta;
	ar << m_AnalysisType;

	if (m_bVLM1) ar << 1; else ar << 0;
	ar << 1;
	if (m_bTiltedGeom) ar << 1; else ar << 0;
	if (m_bWakeRollUp) ar << 1; else ar << 0;

	ar << 1;
    SerializeWing(ar);

	// now store all the WPolars
	ar << (int)polarSize;

	// next store all the WOpps
    ar << 0;

	// then the foils,  polars and Opps
		
	WritePolars(ar);

	// next the bodies
	ar << (int)bodySize;

	// last write the planes...
	ar << (int)planeSize;

    ar << 0;

	CSF *m_pSF;
	CPF *m_pPF;

	m_pSF = new CSF();
	m_pSF->m_bmodified = false;
	m_pSF->InitSplineFoil();

	m_pPF = new CPF();
	m_pPF->m_bmodified = false;
	m_pPF->InitSplinedFoil();

	m_pSF->Serialize(ar, true);
	m_pPF->Serialize(ar, true);

	printf ("\n...SerializeToWpa()");
}
*/
