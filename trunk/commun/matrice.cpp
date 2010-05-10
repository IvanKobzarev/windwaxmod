/*************************************/

/* def type matrice et operateurs ...*/

/*************************************/

#pragma warning(disable:4514)
#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>
//#include "Matrice.h"
#include "matrice.h"
#ifndef DEBUG
    #define DEBUG false
#endif

/************************************/

/* allocation espace pour matrice   */

/************************************/

Matrice::Matrice( int lig, int col )
	: m_nLig(lig), m_nCol(col)
{
    double* tmp = new double[500000];
    delete tmp;
    //printf ("\n Matrice(%d, %d)", m_nLig, m_nCol);
	m_fTab = new double*[m_nLig];
	for( int l=0; l<m_nLig; l++) {
            if (DEBUG) printf (" l=%d",  l);
		m_fTab[l] = new double[m_nCol];
        }
    //printf ("\n... Matrice(%d, %d)", lig, col);
}

Matrice::~Matrice()
{ 
	if ( m_fTab != NULL )
	{
		for( int l=0; l<m_nLig; l++)
			delete [] m_fTab[l];
		delete [] m_fTab;
	}
	m_fTab = NULL;
	m_nLig = m_nCol = -1;
}

double Matrice::Element( int l, int c )
{
	if ( m_fTab == NULL || l < 0 || l > m_nLig || c < 0 || c > m_nCol )
		return 0.0f;

	return m_fTab[l][c];
}

double* Matrice::GetLigne( int l )
{
	if ( m_fTab == NULL || l < 0 || l > m_nLig )
		return NULL; ;

	return m_fTab[l];
}

void Matrice::SetElement( int l, int c, double value )
{
	if ( m_fTab == NULL || l < 0 || l > m_nLig || c < 0 || c > m_nCol )
		return;

	m_fTab[l][c] = value;
}

void Matrice::AddElement( int l, int c, double value )
{
	if ( m_fTab == NULL || l < 0 || l > m_nLig || c < 0 || c > m_nCol )
		return;

	m_fTab[l][c] += value;
}

void Matrice::RemoveElement( int l, int c, double value )
{
	if ( m_fTab == NULL || l < 0 || l > m_nLig || c < 0 || c > m_nCol )
		return;

	m_fTab[l][c] -= value;
}

void Matrice::MultiplyElement( int l, int c, double value )
{
	if ( m_fTab == NULL || l < 0 || l > m_nLig || c < 0 || c > m_nCol )
		return;

	m_fTab[l][c] *= value;
}

void Matrice::print( int c )
{
    for (int i = 0; i < m_nLig; i++) printf ("\n %d -> %f", i, m_fTab[i][c]);
}


void Matrice::DivideElement( int l, int c, double value )
{
	if ( m_fTab == NULL || l < 0 || l > m_nLig || c < 0 || c > m_nCol )
		return;

	m_fTab[l][c] /= value;
}

/********************************/

/* addition de deux matrices	*/

/********************************/

Matrice* AddMat(Matrice *m1,Matrice *m2)
{
	Matrice *res = NULL;

	if(m1==NULL)
		printf("erreur fonction AddMat: matrice m1 = NULL !!!");

	else if(m2==NULL)
		printf("erreur fonction AddMat: matrice m2 = NULL !!!");

	else
	{
		/*test coherence des tailles*/
		if ((m1->GetLignes() != m2->GetLignes())||(m1->GetColonnes() != m2->GetColonnes()))
		{
			res=NULL;
			printf("\n erreur fonction AddMat: tailles diff�rentes");
		}
		else
		{
			res = new Matrice(m1->GetLignes(),m1->GetColonnes());

			for( int l = 0; l < m1->GetLignes(); l++ )
			{
				for( int c = 0; c < m1->GetColonnes(); c++ )
				{ 
					res->SetElement( l, c, m1->Element(l,c) + m2->Element(l,c) );
				}
			}
		}
	}

	return res;
}



/************************************/

/* addition d'un r�el � une matrice */

/************************************/

Matrice* AddReelMat(Matrice *m, double value )
{
	Matrice *res = NULL;

	if (m==NULL)
		printf("erreur fonction AddReelMat: m = NULL !!!");

	else
	{
		res = new Matrice( m->GetLignes(), m->GetColonnes() );

		for( int l = 0; l < m->GetLignes(); l++ )
		{
			for( int c = 0; c < m->GetColonnes(); c++ )
			{
				res->SetElement( l, c, m->Element(l,c) + value );
			}
		}
	}

	return res;
}



/************************************/

/* multiplication de deux matrices	*/

/************************************/

Matrice* MultMat(Matrice *m1,Matrice *m2)
{
	Matrice *res = NULL;

	if(m1==NULL)
		printf("erreur fonction MultMat: m1 = NULL !!!");

	else if(m2==NULL)
		printf("erreur fonction MultMat: m2 = NULL !!!");

	else
	{
		if( m1->GetColonnes() != m2->GetLignes() )
			printf("erreur dimension procedure MultMat");

		else
		{
			res = new Matrice( m1->GetLignes(), m2->GetColonnes() );

			for( int l = 0; l < m1->GetLignes(); l++ )
			{
				for( int c = 0; c < m2->GetColonnes(); c++ )
				{
					double accu = 0.0;

					for( int c1 = 0; c1 < m1->GetColonnes(); c1++ )
					{
						accu += m1->Element(l,c1) * m2->Element(c1,c);
					}
					res->SetElement( l, c, accu );
				}
			}
		}
	}

	return res;
}





/********************************************/

/* multiplication d'une matrice par un reel */

/********************************************/

Matrice * MultReelMat(Matrice *m, double value)
{
	Matrice *res = NULL;

	if (m!=NULL)
	{
		res = new Matrice( m->GetLignes(), m->GetColonnes() );

		for( int l = 0; l < m->GetLignes(); l++ )
		{
			for( int c = 0; c < m->GetColonnes(); c++ )
			{
				res->SetElement( l, c, m->Element(l,c) * value );
			}
		}
	}

	return res;
}



/********************************************/

/* cre� une matrice ligne contenant des     */

/* reels lin�airement espac� entre 2 bornes */

/********************************************/

Matrice* LinSpace(double deb, double fin, int n)
{
	if ( n <= 1 )
		return NULL;

	Matrice *res = new Matrice( n, 1 );

	double pas = (fin-deb)/(n-1); 
	double accu = deb;

	for( int i = 0; i < n; i++ )
	{
		res->SetElement(i,0,accu);
		accu += pas;
	}

	res->SetElement(n-1,0,fin); /*si erreur d'arrondi*/

	return res;
}



/***********************************************/

/* cre� une matrice ne contenant que des zeros */

/***********************************************/

Matrice* Zeros(int lig, int col)
{
   if (DEBUG) printf ("\n Zeros(%d, %d)", lig, col);
	Matrice *res = new Matrice( lig, col );

	for( int l = 0; l < lig; l++ )
	{
            if (DEBUG) printf ("\n l=%d", l);
		for( int c = 0; c < col; c++ )
		{
                    if (DEBUG) printf (" c=%d", c);
			res->SetElement( l, c, 0.0 );
		}
	}
    if (DEBUG) printf ("\n ...Zeros()");
	return res;
}



/***********************************************/

/* cre� une matrice ne contenant que des uns   */

/***********************************************/

Matrice* Ones(int lig, int col)
{
	Matrice *res = new Matrice( lig, col );

	for( int l = 0; l < lig; l++ )
	{
		for( int c = 0; c < col; c++ )
		{
			res->SetElement( l, c, 1.0 );
		}
	}
	return res;
}


/*********************************************************/

/* copie tous les elements d'une matrice dans une autre  */

/*********************************************************/

void CopierMat(Matrice *dst, Matrice *src)
{
	if(src==NULL)
		printf("erreur fonction CopierMat: src = NULL !!!");

	else if(dst==NULL)
		printf("erreur fonction CopierMat: dst = NULL !!!");

	else
	{
		if ( (dst->GetLignes()!=src->GetLignes()) || (dst->GetColonnes()!=src->GetColonnes()) )
			printf("erreur dimension procedure CopierMat");

		else
		{
			for( int l = 0; l < src->GetLignes(); l++ )
			{
				for( int c = 0; c < src->GetColonnes(); c++ )
				{
					dst->SetElement( l, c, src->Element(l,c) );
				}
			}
		}
	}
}

/********************************/
/* cree une copie d'une matrice */
/********************************/

Matrice* CloneMat(Matrice *src)
{
	Matrice *dst = NULL;

	if(src==NULL)
		printf("erreur fonction CopierMat: src = NULL !!!");

	else
	{
		dst = new Matrice( src->GetLignes(), src->GetColonnes() );
		CopierMat( dst, src );
	}

	return dst;
}

/****************************************************/
/* copie les elements d'une matrice dans une autre  */
/* avec indices de debut et fin pour lig et col     */

/* pour la matrice destination et source

/****************************************************/

void CopierMatInd(Matrice *dst, int lDebDst, int lFinDst,int cDebDst, int cFinDst,

				  Matrice *src, int lDebSrc, int lFinSrc,int cDebSrc, int cFinSrc)

{
	if(src==NULL)
		printf("erreur fonction CopierMatInd: src = NULL !!!");

	else if(dst==NULL)
		printf("erreur fonction CopierMatInd: dst = NULL !!!");

	else
	{
		int dl = lFinDst - lDebDst; 
		int dc = cFinDst - cDebDst;

		/*test coherence dimensions*/

		if ((dl != lFinSrc - lDebSrc) || (dc != cFinSrc - cDebSrc)

			||(dst->GetLignes() < lFinDst)||(dst->GetColonnes() < cFinDst)

			||(src->GetLignes() < lFinSrc)||(src->GetColonnes() < lFinSrc)

			||(dl<0)||(dc<0))
		{
			printf("erreur dimension procedure CopierMatInd");
		}
		else
		{
			for( int l = 0; l < dl+1; l++ )
			{
				for( int c = 0; c < dc+1; c++ )
				{
					dst->SetElement( l+lDebDst, c+cDebDst, src->Element(l+lDebSrc, c+cDebSrc) );
				}
			}
		}
	}

}

/************************************/
/* visualise une matrice            */
/************************************/

void VoirMat(Matrice *m)
{
	if( m == NULL ) 
		printf("\n erreur fonction VoirMat: m = NULL !!!");

	else
	{
		int cpt = 0;

		for( int l = 0; l < m->GetLignes(); l++ )
		{
			for( int c = 0; c < m->GetColonnes(); c++ )
			{
				printf("\nm[%d][%d]=%f",l,c,m->Element(l,c));

				if(cpt++ == 20)
				{
					{printf("\n<pause - suite>"); fflush(stdin); getch(); cpt=0;}
				}
			}
		}
	}

	printf("\n<pause - fin>"); fflush(stdin); getch();
}

/*****************/
/* InterpLinMat  */
/*****************/
/*void InterpLinMat(Matrice *m, Matrice *mi)*/

void InterpLinMat(Matrice *x, Matrice *y, Matrice *xi, Matrice *yi)
{
	/*test matrices non nulle*/
	if(x==NULL) 
		printf("\n erreur fonction InterpLinMat: x = NULL !!!");
	else if(y==NULL) 
		printf("\n erreur fonction InterpLinMat: y = NULL !!!");
	else if(xi==NULL) 
		printf("\n erreur fonction InterpLinMat: xi = NULL !!!");
	else if(yi==NULL) 
		printf("\n erreur fonction InterpLinMat: yi = NULL !!!");
	else
	{
		int ipc = 0;
		for( int i = 0; i < xi->GetLignes(); i++ )
		{
			/*recherche indice pts avant/apres*/
			while( (x->Element(ipc,0) <= xi->Element(i,0))
					&& ( ipc < x->GetLignes()-1))
				ipc++;
			if(ipc>0) ipc--; 
			yi->SetElement(i,0,
					y->Element(ipc,0) + (y->Element(ipc+1,0) - y->Element(ipc,0))
					*
					( xi->Element(i,0) - x->Element(ipc,0))/(x->Element(ipc+1,0) - x->Element(ipc,0) ) 
				);
		}
	}
}

/***************/
/* InterpLinX  */
/***************/

double InterpLinX(Matrice *m, double xi)
{
	double yi = 0.0;
	/*test matrices non nulle*/
	if(m==NULL) 
		printf("\n erreur fonction InterpLinX: m = NULL !!!");
	else
	{
		/*recherche indice pts avant/apres*/
		int ipc=0; 
		while( (m->Element(ipc,0) <= xi) && (ipc<m->GetLignes()-1) ) 
			ipc++;

		if(ipc>0) ipc--; 

		/*interpole y*/
		yi = 
			m->Element(ipc,1) + 
				(m->Element(ipc+1,1) - m->Element(ipc,1))
				*
				(xi - m->Element(ipc,0))/(m->Element(ipc+1,0) - m->Element(ipc,0))
			;
	}

	return yi;
}

/******************************************************/
/* cree une sous-matrice extraite d'une autre matrice */
/******************************************************/

Matrice* ExtraitMat(Matrice *m, int lDeb, int lFin,int cDeb, int cFin)

{
	Matrice* res = NULL;

	if(m==NULL)
		printf("erreur fonction ExtraitMat: m = NULL !!!");

	else
	{
		int dl = lFin - lDeb; 
		int dc = cFin - cDeb;

		/*test coherence dimensions*/
		if ( (m->GetLignes() < lFin) || (m->GetColonnes() < cFin) || (dl<0) || (dc<0) )
			printf("erreur dimension procedure ExtraitMat");

		else
		{
			res = new Matrice(dl+1,dc+1);

			for( int l = 0; l < dl+1; l++ )
			{
				for( int c = 0; c < dc+1; c++ )
				{
					res->SetElement( l, c, m->Element(l+lDeb,c+cDeb) );
				}
			}
		}
	}

	return res;
}

/******************************************************/
/* donne l'indice d'une valeur d'un tableau ou de sa  */
/* valeur la plus proche                              */
/******************************************************/

void Ind(double val, Matrice *m, int *lig, int *col)
{
        //    printf ("\n Ind");
	double mini = 10000000000000000.0f;
	int lmini = 0;
	int cmini = 0;
	for( int l = 0; l < m->GetLignes(); l++ )
	{
		for( int c = 0; c < m->GetColonnes(); c++ )
		{
			double dist = (double)fabs(val - m->Element(l,c));
			if ( mini > dist )
			{ 
				mini=dist; 
				lmini=l; 
				cmini=c;
			}
		}
	}
	*lig = lmini; 
	*col = cmini;
        //printf ("\n Ind end");
}

/******************************************************/
/* indique si une valeur est pr�sente dans un tableau */
/******************************************************/

bool ValeurPresente(double val, Matrice *m)
{
	for( int l = 0; l < m->GetLignes(); l++ )
	{
		for( int c = 0; c < m->GetColonnes(); c++ )
		{
			if( val == m->Element(l,c) ) 
				return true;
		}
	}

	return false;
}

/*****************************************************************/
/* ajoute une valeur dans un vecteur colonne par ordre croissant */
/*****************************************************************/

void AjouteValeurCroissant(double val, Matrice **m)

{
	Matrice *res = new Matrice( (*m)->GetLignes()+1, (*m)->GetColonnes() );
	if(val < (*m)->Element(0,0)) //ajout devant ?
	{
		res->SetElement(0,0,val);
		for(int l=0; l<(*m)->GetLignes(); l++) 
			res->SetElement( l+1, 0, (*m)->Element(l,0) );
	}
	else if (val>(*m)->Element((*m)->GetLignes()-1,0)) //ajout derriere ?
	{
		for( int l=0; l<(*m)->GetLignes(); l++ ) res->SetElement( l, 0,(*m)->Element(l,0) );

		res->SetElement((*m)->GetLignes()-1,0,val);
	}

	else //ajout dans le tableau
	{
		int l = 0;
		while(val>(*m)->Element(l,0))
		{ 
			res->SetElement( l, 0, (*m)->Element(l,0) ); 
			l++;
		}

		res->SetElement(l,0,val);
		for( int i=l+1; i<res->GetLignes(); i++) 
			res->SetElement(i,0,(*m)->Element(i-1,0));
	}

	//mise a jour pointeurs ...
	delete *m;
	*m = res;
}

/*****************************************************************

* This is a gaussian elimination program with pivoting.          *

* The matrix is in an array called a; which is an n by n+1 array.*

* The choice of pivots is kept in an array called                *

* row[].                                                         *

* proc�dure r�cup�r�e sur internet et adapt�e � mes besoins pour *

* la resolution du systeme: A . x = B 

******************************************************************/

void ResolutionGauss(Matrice *A, Matrice *B, Matrice **x)
{
	//init var n et n1
	int n = A->GetLignes(); 
	int n1 = n-1;

	//init matrice a
	Matrice *a = Zeros(n, n+1);

	for( int i=0; i<n; i++)
	{
		a->SetElement(i,n,B->Element(i,0));

		for(int j=0; j<n; j++) 
			a->SetElement(i,j,A->Element(i,j));
	}

	/*initialize the array row*/
	Matrice *row = Zeros(n,1);

	for(int i=0;i<=n1;i++) 
		row->SetElement(i,0,(double)i);

	/*This outermost of 3 loops that performs the elimination step */
	for(int i=0; i<=n1;i++) 
	{
		/* First, find the choice of pivot */
		int p  = i;
		int i1 = i+1;

		double temp=(double)fabs(a->Element((int)row->Element(i,0),i));

		double temp1;
		for( int j=i1; j<n; j++ )
		{
			if( (temp1=(double)fabs(a->Element((int)row->Element(j,0),i))) > temp )
			{
				p=j; 
				temp=temp1;
			}
		}

		if(temp<1.E-9) 
		{
			printf("no solution is likely \n");
			delete a;
			delete row;
			return;
		}

		/* The largest pivot value is in temp and p is the
		corresponding index of the array row.
		If the largest pivot is small (< 10^(-9)), then
		a message is printed that no solution is likely, or can be trusted*/

		/*This step performs a simulated row switch by switching
		the values of row->t[i] with row->t[p] */

		int trow = (int)row->Element(i,0);
		row->SetElement(i,0,row->Element(p,0));
		row->SetElement(p,0,(double)trow);

		/* This is the elimination step */
		for( int j=i1; j<=n1; j++ )
		{
			int j1 = j+1;
			double mult = a->Element((int)row->Element(j,0),i) 
						/ a->Element((int)row->Element(i,0),i);
			for( int l=i1; l<=n; l++ )
			{
				a->SetElement(
					(int)row->Element(j,0),
					l,
					a->Element((int)row->Element(j,0),l) - mult*a->Element((int)row->Element(i,0),l)
					);
			}
		}
	}


	/* Now initialize the solution array, x->t[] */
	//for(k=0;k<=n1;k++) x->Element(k,0)=0.0;

	*x = Zeros(n,1);

	/* Now perform back substitution
	to calculate the solutions, x->Element(k,0)*/

	n1 = n-1;
	(*x)->SetElement(n1,0,a->Element((int)row->Element(n1,0),n) / a->Element((int)row->Element(n1,0),n1));

	for(int i=2;i<=n;i++)
	{
		double sum=0.0;
		for(int j=n-i+1;j<n;j++) 
			sum = sum + a->Element((int)row->Element(n-i,0),j) * (*x)->Element(j,0);

		(*x)->SetElement(n-i,0,
			(a->Element((int)row->Element(n-i,0),n)-sum) / a->Element((int)row->Element(n-i,0),n-i));
	}

	//liberation matrice de calcul
	delete row;
	delete a;
}

