#pragma once

class Matrix
{
public:
	Matrix( int lig, int col);
	virtual ~Matrix();

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
	int m_nLig, m_nCol;

	double **m_fTab;

};

Matrix* AddMat(Matrix *m1,Matrix *m2);

Matrix* AddReelMat(Matrix *m, double reel);

Matrix* MultMat(Matrix *m1,Matrix *m2);

Matrix* MultReelMat(Matrix *m, double reel);

Matrix* LinSpace(double deb, double fin, int n);

Matrix* Zeros(int lig, int col);

Matrix* Ones(int lig, int col);

void VoirMat(Matrix *m);

void CopierMat(Matrix *dst, Matrix *src);

Matrix* CloneMat(Matrix *src);

void CopierMatInd(Matrix *dst, int lDebDst, int lFinDst,int cDebDst, int cFinDst,

				  Matrix *src, int lDebSrc, int lFinSrc,int cDebSrc, int cFinSrc);

/*void InterpLinMat(TMatrix *m, TMatrix *mi);*/

void InterpLinMat(Matrix *x, Matrix *y, Matrix *xi, Matrix *yi);

double InterpLinX(Matrix *m, double xi);

Matrix* ExtraitMat(Matrix *m, int lDeb, int lFin,int cDeb, int cFin);

void Ind(double val, Matrix *m, int *lig, int *col);

bool ValeurPresente(double val, Matrix *m);

void addeValeurCroissant(double val, Matrix **m);

void ResolutionGauss(Matrix *A, Matrix *B, Matrix **x);
