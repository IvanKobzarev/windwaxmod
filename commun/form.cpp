#include "form.h"
#include "matrice.h"

Form::Form() 
	: m_pProfils( NULL ), m_nbProfils(0)
{
	//Value par default
	CoeffProgGeom = 0.93f;
        CoeffExp = 0.0f;
	EpaiRelCent = 17.0f;
	
	EpaiRelProfCent=0;
	EpaiRelProfBout=0;

	NbCaiss = 21;

	//init ptr tableaux a NULL
	mCtrlNez = NULL;
	mCtrlFui = NULL;
	mCtrlA = NULL;
	mCtrlB = NULL;
	mCtrlC = NULL;
	mCtrlD = NULL;
	mCtrlE = NULL;
	mCtrlDiedre = NULL;
	mCtrlEpaiRel = NULL;
	mCtrlMorphing = NULL;
	mCtrlVrillage = NULL;
}


Matrix* Form::getExtProf(int nerv, bool realSize) {
	Matrix* res = new Matrix(ExtProfCent->GetLignes(), 2);
	double m = m_pProfils[nerv]->m_fMorph;
	double LongNerv = m_pProfils[nerv]->m_fLength;
	if (realSize == false) LongNerv = 100.0f;
	double EpaiRel = m_pProfils[nerv]->m_fWidth;
	double coeffx = LongNerv/100.0f;
	double coeffyCent = LongNerv*EpaiRel/(EpaiRelProfCent*100.0f);
	double coeffyBout = LongNerv*EpaiRel/(EpaiRelProfBout*100.0f);
	double xp=0, yp = 0;
	for (int j=0; j<ExtProfCent->GetLignes(); j++)
	{
	    xp = ExtProfCent->Element(j,0)*coeffx;
        yp = ExtProfCent->Element(j,1)*coeffyCent*m + ExtProfBout->Element(j,1)*coeffyBout*(1.0f-m);
		res->SetElement(j, 0, xp);
		res->SetElement(j, 1, yp);
    }
	return res;
}

Matrix* Form::getIntProf(int nerv, bool realSize) {
	Matrix* res = new Matrix(IntProfCent->GetLignes(), 2);
	double m = m_pProfils[nerv]->m_fMorph;
	double EpaiRel = m_pProfils[nerv]->m_fWidth;
	double LongNerv = m_pProfils[nerv]->m_fLength;
	if (realSize == false) LongNerv = 100.0f;
	double coeffx = LongNerv/100.0f;
	double coeffyCent = LongNerv*EpaiRel/(EpaiRelProfCent*100.0f);
	double coeffyBout = LongNerv*EpaiRel/(EpaiRelProfBout*100.0f);
	double xp=0, yp=0;
	for (int j=0; j<IntProfCent->GetLignes(); j++)
	{
		xp = IntProfCent->Element(j,0)*coeffx;
		yp = IntProfCent->Element(j,1)*coeffyCent*m	+ IntProfBout->Element(j,1)*coeffyBout*(1.0f-m);
		res->SetElement(j, 0, xp);
		res->SetElement(j, 1, yp);
    }
	return res;

}


Form3D::Form3D(){
}

Form3D::~Form3D(){
	delete this->XExt;
	delete this->YExt;
	delete this->ZExt;

	delete this->XInt;
	delete this->YInt;
	delete this->ZInt;
}

FormProjection::FormProjection(){
}

FormProjection::~FormProjection(){
	delete this->X;
	delete this->Y;
}

Ballonement::Ballonement(){
}

Ballonement::~Ballonement(){
	delete this->kChord;
	delete this->dyw;
	delete this->wN;
	delete this->kMf;
}

void Form::DeleteProfils()
{
	if ( m_pProfils != NULL )
	{
		for( int i = 0; i<m_nbProfils; i++ )
		{
			if ( m_pProfils[i] != NULL )
				delete m_pProfils[i];
			m_pProfils[i] = NULL;
		}
		delete m_pProfils; 
	}
	m_pProfils = NULL;
	m_nbProfils = 0;
}

void Form::Validate()
{
	for (int i = 0; i<m_nbProfils; i++) {
		if ((m_pProfils[i]->m_fMorph > 1.0) || (m_pProfils[i]->m_fMorph < 0.0)) {
		   char msg[200];
		   sprintf(msg, " bad forme, morph at %d profil == %f, should be in [0, 1]" ,
			   i, m_pProfils[i]->m_fMorph);
		   throw msg;
		}
	}
}


void Form::AllocateProfils( int nb )
{
	DeleteProfils();
	m_nbProfils = nb;
	m_pProfils = new Profil*[nb];
	for( int i = 0; i < m_nbProfils; i++ )
	{
		m_pProfils[i] = new Profil();
	}
}

Form::~Form() 
{ 
	DeleteProfils();

	if ( mCtrlNez != NULL ) 
		delete(mCtrlNez);
	mCtrlNez = NULL;

	if ( mCtrlFui != NULL ) 
		delete(mCtrlFui);
	mCtrlFui = NULL;

	if ( mCtrlA != NULL ) 
		delete(mCtrlA);
	mCtrlA = NULL;

	if ( mCtrlB != NULL ) 
		delete(mCtrlB);
	mCtrlB = NULL;

	if ( mCtrlC != NULL ) 
		delete(mCtrlC);
	mCtrlC = NULL;

	if ( mCtrlD != NULL ) 
		delete(mCtrlD);
	mCtrlD = NULL;

	if ( mCtrlE != NULL ) 
		delete(mCtrlE);
	mCtrlE = NULL;

	if ( mCtrlDiedre != NULL ) 
		delete(mCtrlDiedre);
	mCtrlDiedre = NULL;

	if ( mCtrlEpaiRel != NULL ) 
		delete(mCtrlEpaiRel);
	mCtrlEpaiRel = NULL;

	if ( mCtrlMorphing != NULL ) 
		delete(mCtrlMorphing);
	mCtrlMorphing = NULL;

	if ( mCtrlVrillage != NULL ) 
		delete(mCtrlVrillage);
	mCtrlVrillage = NULL;
}

void Ballonement::loadFromFile(const char* fileName) {
	FILE *fid;

	if( (fid = fopen( fileName, "rt" )) == NULL )
	{
		printf( "\nErreur ouverture fichier: '%s'", fileName );
		char ex[100];
		sprintf (ex, " could not open `%s`", fileName);
		throw ex;
	}
	else
	{
		//fscanf(fid,"%d %lf %lf %lf",&x);
	}
	 
}

/*****************/
/* calcForm3D */
/*****************/

void calcForm3D(Form *forme, int isPercent, double percent,
				   Matrix *ExtProfCent, Matrix *IntProfCent,
				   Matrix *ExtProfBout, Matrix *IntProfBout,
				   Matrix **XExt, Matrix **YExt, Matrix **ZExt,
				   Matrix **XInt, Matrix **YInt, Matrix **ZInt)

{
    //printf ("\n calcForm3D");
    Matrix *ExtProfCentN, *ExtProfBoutN;
	double LongNerv, EpaiRel, xp,yp, xo,yo,zo, a,v,m;
	double EpaiRelProfCent, EpaiRelProfBout;
	double coeffx, coeffyCent, coeffyBout;
	int i,j;
    bool isCenterPanel = (1 & forme->NbCaiss);

	/*calc epaisseur relative profil central et bout*/
	EpaiRelProfCent = EpaisseurRelative(ExtProfCent, IntProfCent);
	//printf ("\nEpaiRelProfCent=%f", EpaiRelProfCent);
	EpaiRelProfBout = EpaisseurRelative(ExtProfBout, IntProfBout);
	//printf ("\nEpaiRelProfBout=%f", EpaiRelProfBout);
	/*init matrices extrados*/
	//delete(*XExt); delete(*YExt); delete(*ZExt); 
	//*XExt = Zeros(5, 5);
	*XExt = Zeros(forme->m_nbProfils, ExtProfCent->GetLignes());
	*YExt = Zeros(forme->m_nbProfils, ExtProfCent->GetLignes());
	*ZExt = Zeros(forme->m_nbProfils, ExtProfCent->GetLignes());

	/*init matrices intrados*/
	//delete(*XInt); delete(*YInt); delete(*ZInt); 
	*XInt = Zeros(forme->m_nbProfils, IntProfCent->GetLignes());
	*YInt = Zeros(forme->m_nbProfils, IntProfCent->GetLignes());
	*ZInt = Zeros(forme->m_nbProfils, IntProfCent->GetLignes());
	/*boucle � partir de la 1ere nervure du centre vers l'extr�mit�*/


    if (isPercent) {
        //printf ("\n isPercent");
        makePosProfile(ExtProfCent, IntProfCent, percent, &ExtProfCentN);
        makePosProfile(ExtProfBout, IntProfBout, percent, &ExtProfBoutN);
    }
        
	for (i=0; i<forme->m_nbProfils; i++)
	{
		//longueur nervure courante
		LongNerv = forme->m_pProfils[i]->m_fLength;
		//epaisseur relative
		EpaiRel = forme->m_pProfils[i]->m_fWidth;
		//position xo,yo,zo du nez
		xo = forme->m_pProfils[i]->m_fNezX;
		yo = forme->m_pProfils[i]->m_fNezY;
		zo = forme->m_pProfils[i]->m_fNezZ;
        //printf ("\n %d long=%f width=%f (%f, %f, %f)", i, forme->m_pProfils[i]->m_fLength, forme->m_pProfils[i]->m_fWidth, xo, yo, zo);
		//inclinaison de la nervure par rapport a l'horizontale
		a = forme->m_pProfils[i]->m_fInclin;
        if ((isCenterPanel) && (i == 0)) a = forme->m_pProfils[0]->m_fInclin - 0.33333333f * fabs(forme->m_pProfils[0]->m_fInclin - forme->m_pProfils[1]->m_fInclin);
                
		//angle vrillage
		v = forme->m_pProfils[i]->m_fWash;
		//coeff morphing
		m = forme->m_pProfils[i]->m_fMorph;
               
		//printf ("\n%3d -> LNerv=%f EpRel=%f a=%f v=%f m=%f",
		//				i,  LongNerv, EpaiRel, a * 180.0f/pi, v, m);

		//calc coeffx et coeffy des points du profil en
		//fonction de l'�paisseur relative et de la longueur de nervure
		coeffx = LongNerv/100.0f;
		coeffyCent = LongNerv*EpaiRel/(EpaiRelProfCent*100.0f);
		coeffyBout = LongNerv*EpaiRel/(EpaiRelProfBout*100.0f);

        /*boucle sur les points du profil en intrados*/
		for (j=0; j<IntProfCent->GetLignes(); j++)
		{
			xp = IntProfCent->Element(j,0)*coeffx;
			yp = IntProfCent->Element(j,1)*coeffyCent*m
				+ IntProfBout->Element(j,1)*coeffyBout*(1.0f-m);
			(*XInt)->SetElement(i,j, xo+(yp-xp*(double)sin(v))*(double)cos(a));
			(*YInt)->SetElement(i,j, yo+(yp-xp*(double)sin(v))*(double)sin(a));
			(*ZInt)->SetElement(i,j, zo-(xp*(double)cos(v)));

		}

		/*boucle sur les points du profil en extrados*/
		for (j=0; j<ExtProfCent->GetLignes(); j++)
		{
            if (isPercent) {
                xp = ExtProfCentN->Element(j,0)*coeffx;
                yp = ExtProfCentN->Element(j,1)*coeffyCent*m
                        + ExtProfBoutN->Element(j,1)*coeffyBout*(1.0f-m);
            } else {
                xp = ExtProfCent->Element(j,0)*coeffx;
                yp = ExtProfCent->Element(j,1)*coeffyCent*m
                        + ExtProfBout->Element(j,1)*coeffyBout*(1.0f-m);
            }
			(*XExt)->SetElement(i,j, xo+(yp-xp*(double)sin(v))*(double)cos(a));
			(*YExt)->SetElement(i,j, yo+(yp-xp*(double)sin(v))*(double)sin(a));
			(*ZExt)->SetElement(i,j, zo-(xp*(double)cos(v)));
          /*  if (i == 0) {
                printf ("\n%d DELTA(%f, %f, %f)",j, (yp-xp*(double)sin(v))*(double)cos(a), (yp-xp*(double)sin(v))*(double)sin(a), -(xp*(double)cos(v)));
            } */
		}

	}
    //printf ("\n ...calcForm3D");

}

FormProjection* getFormProjection(Form3D* f3d) 
{
	//printf ("\n getFormProjection()");
	FormProjection* formePrj = new FormProjection();
	// we should X, 0, Z
	Matrix *X = Zeros(f3d->XExt->GetLignes(), 2); 
	Matrix *Y = Zeros(f3d->YExt->GetLignes(), 2);
		
	for (int i = 0; i < f3d->XExt->GetLignes(); i++)
	{
		X->SetElement(i, 0, f3d->XExt->Element(i, 0));
		Y->SetElement(i, 0, f3d->ZExt->Element(i, 0));
		//printf ("\n %d, 0: %f, %f",i, f3d->XExt->Element(i, 0), f3d->ZExt->Element(i, 0));
		X->SetElement(i, 1, f3d->XExt->Element(i, f3d->XExt->GetColonnes()-1));
		Y->SetElement(i, 1, f3d->ZExt->Element(i, f3d->ZExt->GetColonnes()-1));
		//printf ("\n %d, 1: %f, %f", i, f3d->XExt->Element(i, f3d->XExt->GetColonnes()-1), f3d->ZExt->Element(i, f3d->XExt->GetColonnes()-1));
	}

	formePrj->X = X;
	formePrj->Y = Y;
	return formePrj;
}

Form3D* getForm3D(Form *forme, int isPercent, double percent)
				   //Matrix *ExtProfCent, Matrix *IntProfCent,
				   //Matrix *ExtProfBout, Matrix *IntProfBout)
{
	Form3D* forme3D = new Form3D();
    Matrix *XExt, *YExt, *ZExt;
    Matrix *XInt, *YInt, *ZInt;

    calcForm3D(forme, isPercent, percent,
            forme->ExtProfCent, forme->IntProfCent, forme->ExtProfBout, forme->IntProfBout,
            &XExt, &YExt, &ZExt, &XInt, &YInt, &ZInt);

	forme3D->XExt=XExt;
	forme3D->YExt=YExt;
	forme3D->ZExt=ZExt;

	forme3D->XInt=XInt;
	forme3D->YInt=YInt;
	forme3D->ZInt=ZInt;

	//forme3D->ExtProfCent=ExtProfCent;
	//forme3D->IntProfCent=IntProfCent;

	//forme3D->ExtProfBout=ExtProfBout;
	//forme3D->IntProfBout=IntProfBout;

	forme3D->forme=forme;
	return forme3D;
}

double calcWidthNervs(WindPatternsProject* gfd, int noNerv1, int noNerv2, int face) {
    Matrix * Xd1[2], *Yd1[2]; //, *Xd1p[2], *Yd1p[2];
    Matrix * X1[2], *Y1[2], *Z1[2], *P1[2];
    calcPatron(gfd, noNerv1, false, face, face, 0.0f, 100.0f,
            noNerv2, false, face, face, 0.0f, 100.0f,
            &Xd1[0], &Yd1[0], &Xd1[1], &Yd1[1],
            &X1[0], &Y1[0], &Z1[0], &P1[0],
            &X1[1], &Y1[1], &Z1[1], &P1[1]);
    double width;
    calcWidth(Xd1[0], Yd1[0], Xd1[1], Yd1[1], &width);
    delete (Xd1[0]);
    delete (Yd1[0]);
    delete (Xd1[1]);
    delete (Yd1[1]);
    delete(X1[0]);
    delete(Y1[0]);
    delete(Z1[0]);
    delete(P1[0]);
    delete(X1[1]);
    delete(Y1[1]);
    delete(Z1[1]);
    delete(P1[1]);
    return width;
}
