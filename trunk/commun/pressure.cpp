
/*******************************************************************************
*
*     Linear vorticity surface panel method for airfoils.
*
*     Adapted from Kuethe and Chow 4th Edition
*     This version by I. Kroo  2-2-87
*
*     No attempt has been made to make this particularly fast or efficient.
*
*     Inputs:
*     -------
*     XB, YB     x and y coordinates of panel edges starting at trailing edge
*                proceeding forward on lower surface, wrapping around
*                leading edge and then running back to the trailingedge.
*     alphaD     Angle of attack in degrees
*
*     Outputs:
*     --------
*     X, Y       coordinates of panel centers
*     Cp         incompressible pressure coefficient at panel center
*     
********************************************************************************/

// version originale en fortran, traduit en C et Matlab par Thierry P�bayle, 2001
// function [X,Y,Cp]=LVFoil(XB,YB,alphaD);

void LVFoil(Matrix *XB, Matrix *YB, double alphaD, Matrix **X, Matrix **Y, Matrix **Cp )
{
	int n, np1, i, ip1, j;
	double alpha, SINA, COSA;
	Matrix *S=NULL, *THETA=NULL, *SINT=NULL, *COST=NULL;
	Matrix *CN1=NULL, *CN2=NULL, *CT1=NULL, *CT2=NULL, *AN=NULL, *AT=NULL, *RHS=NULL, *XX=NULL, *V=NULL;
	double A,B,C,D,E,F,G,P1,P2,P3,P4,P,Q;
	double ANmax;
	//
	//     calcation of geometric data:
	//     ------------------------------

	n=XB->GetLignes()-1; // n:Number of panels
	np1 = n+1;

	//init matrices
	*X=Zeros(n,1); *Y=Zeros(n,1); *Cp=Zeros(n,1);
	S=Zeros(n,1); THETA=Zeros(n,1); SINT=Zeros(n,1); COST=Zeros(n,1);

	//calc position centre et taille des panneaux 
	for(i=0;i<n; i++)
	{
		ip1=i+1;
		(*X)->SetElement(i,0, .5f * (XB->Element(i,0) + XB->Element(ip1,0)));
		(*Y)->SetElement(i,0, .5f * (YB->Element(i,0) + YB->Element(ip1,0)));
		S->SetElement(i,0, (double)sqrt( sqr(XB->Element(ip1,0)-XB->Element(i,0)) + sqr(YB->Element(ip1,0)-YB->Element(i,0)) ));
		THETA->SetElement(i,0, (double)atan2( (YB->Element(ip1,0)-YB->Element(i,0)), (XB->Element(ip1,0)-XB->Element(i,0)) ));
		SINT->SetElement(i,0, (double)sin(THETA->Element(i,0)));
		COST->SetElement(i,0, (double)cos(THETA->Element(i,0)));
	}	
	alpha = alphaD * DEG2RAD;

	SINA = (double)sin(alpha); 

	COSA = (double)cos(alpha);

	

	//

	//     calcation of influence coefficients

	//     -------------------------------------

	//

	

	// init tableaux

	CN1=Zeros(n,n); CN2=Zeros(n,n); CT1=Zeros(n,n); CT2=Zeros(n,n);

	AN=Zeros(np1,np1); AT=Zeros(np1,np1);

	

	for (i=0; i<n; i++)

	{

		for (j=0; j<n; j++)

		{

			if (i==j) 

			{

				CN1->SetElement(i,j, -1.0);

				CN2->SetElement(i,j, 1.0);

				CT1->SetElement(i,j, 0.5*pi);

				CT2->SetElement(i,j, CT1->Element(i,j));

			}

			else

			{

				A = -((*X)->Element(i,0)-XB->Element(j,0))*COST->Element(j,0)

					- ((*Y)->Element(i,0)-YB->Element(j,0))*SINT->Element(j,0);

				B = sqr((*X)->Element(i,0)-XB->Element(j,0)) + sqr((*Y)->Element(i,0)-YB->Element(j,0));

				C = SINT->Element(i,0)*COST->Element(j,0)-COST->Element(i,0)*SINT->Element(j,0);

				D = COST->Element(i,0)*COST->Element(j,0)+SINT->Element(i,0)*SINT->Element(j,0);

				E = ((*X)->Element(i,0)-XB->Element(j,0))*SINT->Element(j,0)

					- ((*Y)->Element(i,0)-YB->Element(j,0))*COST->Element(j,0);

				F = (double)log( 1. + S->Element(j,0)*(S->Element(j,0)+2.*A)/B );

				G = (double)atan2(E*S->Element(j,0), B+A*S->Element(j,0));

				P1 = 1.0f-2.0f*SINT->Element(j,0)*SINT->Element(j,0);

				P2 = 2.0f*SINT->Element(j,0)*COST->Element(j,0);

				P3 = SINT->Element(i,0)*P1 - COST->Element(i,0)*P2;

				P4 = COST->Element(i,0)*P1 + SINT->Element(i,0)*P2;

				P = ((*X)->Element(i,0)-XB->Element(j,0)) * P3 + ((*Y)->Element(i,0)-YB->Element(j,0)) * P4;

				Q = ((*X)->Element(i,0)-XB->Element(j,0)) * P4 - ((*Y)->Element(i,0)-YB->Element(j,0)) * P3;

				CN2->SetElement(i,j, D + ( .5f*Q*F - (A*C+D*E)*G )/S->Element(j,0));

				CN1->SetElement(i,j, .5f*D*F + C*G - CN2->Element(i,j));

				CT2->SetElement(i,j, C + ( .5f*P*F + (A*D-C*E)*G )/S->Element(j,0));

				CT1->SetElement(i,j, .5f*C*F - D*G - CT2->Element(i,j));

			}

		}

	}

	

	ANmax = 0.0;

	for (i=0; i<n; i++)

	{

		AN->SetElement(i,0, CN1->Element(i,0));

		AN->SetElement(i,n, CN2->Element(i,n-1));

		AT->SetElement(i,0, CT1->Element(i,0));

		AT->SetElement(i,n, CT2->Element(i,n-1));

		for (j=1; j<n; j++)

		{

			AN->SetElement(i,j, CN1->Element(i,j) + CN2->Element(i,j-1));

			AT->SetElement(i,j, CT1->Element(i,j) + CT2->Element(i,j-1));

			if(fabs(AN->Element(i,j))>ANmax)

			{

				ANmax = (double)fabs(AN->Element(i,j));

			}

		}

	}

				

	//     The Kutta Condition is imposed by the relation:

	//     A*gamma(1) + A*gamma(n+1) = 0.  A would be 1 but for matrix

	//     conditioning problems.

			

	AN->SetElement(n,0, 1.0);
	AN->SetElement(n,n, 1.0);

	for(j=1; j<n; j++)
	{
		AN->SetElement(n,j, 0.0);
	}
				

	//

	//     Decompose and solve the system:

	//     -------------------------------

	//

	//     AN.XX = RHS

	//     avec  RHS(i) = sin( THETA(i)-alpha )

	//			

	RHS=Zeros(np1,1);

	for(i=0; i<n; i++)

	{

		RHS->SetElement(i,0, SINT->Element(i,0)*COSA - COST->Element(i,0)*SINA);

	}

	//XX=inv(AN)*RHS;

	ResolutionGauss(AN, RHS, &XX);

				

	//

	//     Compute derived quantities:

	//     ---------------------------

	//

	//     V(i) = cos( THETA(i) - alpha )

	V=Zeros(n,1);

	*Cp=Zeros(n,1);

	for (i=0; i<n; i++)

	{

		V->SetElement(i,0,  COST->Element(i,0)*COSA + SINT->Element(i,0)*SINA);

		for (j=0; j<np1; j++)

		{

			V->SetElement(i,0, V->Element(i,0) + AT->Element(i,j) * XX->Element(j,0));

		}

		(*Cp)->SetElement(i,0, 1-V->Element(i,0)*V->Element(i,0));

	}



	//liberation matrices de calc

	delete(S); delete(THETA); delete(SINT); delete(COST);

	delete(CN1); delete(CN2); delete(CT1); delete(CT2);

	delete(AN); delete(AT); delete(RHS); delete(XX);

	delete(V);



}

