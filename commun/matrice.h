#pragma once
/*************************************/
/* def type matrice et operateurs ...*/
/*************************************/

/****************/
/* def matrices */
/****************/

class Matrice
{
public:
	Matrice( int lig, int col);
	virtual ~Matrice();

	double Element( int l, int c );

	double* GetLigne(int l);
        void print( int c );
	void SetElement( int l, int c, double value );
	void AddElement( int l, int c, double value );
	void RemoveElement( int l, int c, double value );
	void MultiplyElement( int l, int c, double value );
	void DivideElement( int l, int c, double value );

	int GetLignes() { return  m_nLig; }
	int GetColonnes() { return m_nCol; }

protected:
	int m_nLig, m_nCol;	/*nombre de lignes et colonnes*/

	double **m_fTab;	/*tableau a 2 dimensions*/

};



/**************/

/* procedures */

/**************/

Matrice* AddMat(Matrice *m1,Matrice *m2);

Matrice* AddReelMat(Matrice *m, double reel);

Matrice* MultMat(Matrice *m1,Matrice *m2);

Matrice* MultReelMat(Matrice *m, double reel);

Matrice* LinSpace(double deb, double fin, int n);

Matrice* Zeros(int lig, int col);

Matrice* Ones(int lig, int col);

void VoirMat(Matrice *m);

void CopierMat(Matrice *dst, Matrice *src);

Matrice* CloneMat(Matrice *src);

void CopierMatInd(Matrice *dst, int lDebDst, int lFinDst,int cDebDst, int cFinDst,

				  Matrice *src, int lDebSrc, int lFinSrc,int cDebSrc, int cFinSrc);

/*void InterpLinMat(TMatrice *m, TMatrice *mi);*/

void InterpLinMat(Matrice *x, Matrice *y, Matrice *xi, Matrice *yi);

double InterpLinX(Matrice *m, double xi);

Matrice* ExtraitMat(Matrice *m, int lDeb, int lFin,int cDeb, int cFin);

void Ind(double val, Matrice *m, int *lig, int *col);

bool ValeurPresente(double val, Matrice *m);

void AjouteValeurCroissant(double val, Matrice **m);

void ResolutionGauss(Matrice *A, Matrice *B, Matrice **x);



