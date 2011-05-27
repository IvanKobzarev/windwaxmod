
void calcForm3DBallonement
				(WindPatternsProject* gfd, Form *forme, int isPercent, double percent,
				   Matrix *ExtProfCent, Matrix *IntProfCent,
				   Matrix *ExtProfBout, Matrix *IntProfBout,
				   Matrix **XExt, Matrix **YExt, Matrix **ZExt,
				   Matrix **XInt, Matrix **YInt, Matrix **ZInt)

{
	// !!! support of isPercent, percent not implemented yet!!!

    //printf ("\n calcForm3DBallonement");
    Matrix *ExtProfCentN, *ExtProfBoutN;
	double LongNerv, EpaiRel, xp, yp, xo, yo, zo, a, a0, v, m;
	double EpaiRelProfCent, EpaiRelProfBout;
	double coeffx, coeffy, coeffyCent, coeffyBout;
	int i,j;
    bool isCenterPanel = (1 & forme->NbCaiss);

	/*calc epaisseur relative profil central et bout*/
	EpaiRelProfCent = EpaisseurRelative(ExtProfCent, IntProfCent);
	EpaiRelProfBout = EpaisseurRelative(ExtProfBout, IntProfBout);

	*XExt = Zeros(forme->m_nbProfils*2 - 1, ExtProfCent->GetLignes());
	*YExt = Zeros(forme->m_nbProfils*2 - 1, ExtProfCent->GetLignes());
	*ZExt = Zeros(forme->m_nbProfils*2 - 1, ExtProfCent->GetLignes());

	*XInt = Zeros(forme->m_nbProfils*2 - 1, IntProfCent->GetLignes());
	*YInt = Zeros(forme->m_nbProfils*2 - 1, IntProfCent->GetLignes());
	*ZInt = Zeros(forme->m_nbProfils*2 - 1, IntProfCent->GetLignes());

    if (isPercent) {
        makePosProfile(ExtProfCent, IntProfCent, percent, &ExtProfCentN);
        makePosProfile(ExtProfBout, IntProfBout, percent, &ExtProfBoutN);
    }
	//printf ("\n calcForm3DBallonement go in FOR");        
	for (i=0; i<forme->m_nbProfils; i++)
	{
		// normal nervure profile
		//printf ("\n\n\n %d nervure normal", i);        
		LongNerv = forme->m_pProfils[i]->m_fLength;
		EpaiRel = forme->m_pProfils[i]->m_fWidth;
		xo = forme->m_pProfils[i]->m_fNezX;
		yo = forme->m_pProfils[i]->m_fNezY;
		zo = forme->m_pProfils[i]->m_fNezZ;
		a = forme->m_pProfils[i]->m_fInclin;
        if ((isCenterPanel) && (i == 0)) a = forme->m_pProfils[0]->m_fInclin - 0.33333333f * fabs(forme->m_pProfils[0]->m_fInclin - forme->m_pProfils[1]->m_fInclin);
              
		//angle vrillage
		v = forme->m_pProfils[i]->m_fWash;
		//coeff morphing
		m = forme->m_pProfils[i]->m_fMorph;
		coeffx = LongNerv/100.0f;
		coeffyCent = LongNerv*EpaiRel/(EpaiRelProfCent*100.0f);
		coeffyBout = LongNerv*EpaiRel/(EpaiRelProfBout*100.0f);
		//printf ("\n %d N longNerv=%f", i, LongNerv);
		//printf ("\n %d N coeffx=%f coeffy=(%f, %f)", i, coeffx, coeffyCent, coeffyBout);
		//printf ("\n %d N(init) EpaiRelProfCent=%f", i, EpaiRelProfCent);
		//printf ("\n %d N(init) EpaiRelProfBout=%f", i, EpaiRelProfBout);
		for (j=0; j<IntProfCent->GetLignes(); j++)
		{
			xp = IntProfCent->Element(j,0)*coeffx;
			yp = IntProfCent->Element(j,1)*coeffyCent*m
				+ IntProfBout->Element(j,1)*coeffyBout*(1.0f-m);
			(*XInt)->SetElement(2*i,j, xo+(yp-xp*(double)sin(v))*(double)cos(a));
			(*YInt)->SetElement(2*i,j, yo+(yp-xp*(double)sin(v))*(double)sin(a));
			(*ZInt)->SetElement(2*i,j, zo-(xp*(double)cos(v)));

		}
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
			(*XExt)->SetElement(2*i,j, xo+(yp-xp*(double)sin(v))*(double)cos(a));
			(*YExt)->SetElement(2*i,j, yo+(yp-xp*(double)sin(v))*(double)sin(a));
			(*ZExt)->SetElement(2*i,j, zo-(xp*(double)cos(v)));
		}
		//printf ("\n len (%d) ext(%d) int(%d)",i, ExtProfCent->GetLignes(), IntProfCent->GetLignes());
		// --  normal nervure profile
		if ( i < (forme->m_nbProfils - 1)) {
			// nervure ballonement
			LongNerv = (forme->m_pProfils[i]->m_fLength + forme->m_pProfils[i+1]->m_fLength) * 0.5f;
			EpaiRel = (forme->m_pProfils[i]->m_fWidth + forme->m_pProfils[i+1]->m_fWidth) * 0.5f;
			xo = (forme->m_pProfils[i]->m_fNezX + forme->m_pProfils[i+1]->m_fNezX) * 0.5f;
			yo = (forme->m_pProfils[i]->m_fNezY + forme->m_pProfils[i+1]->m_fNezY) * 0.5f;
			zo = (forme->m_pProfils[i]->m_fNezZ + forme->m_pProfils[i+1]->m_fNezZ) * 0.5f;
			a0 = (forme->m_pProfils[i]->m_fInclin + forme->m_pProfils[i+1]->m_fInclin) * 0.5f;
			if ((isCenterPanel) && (i == 0)) a0 = forme->m_pProfils[0]->m_fInclin - 0.33333333f * fabs(forme->m_pProfils[0]->m_fInclin - forme->m_pProfils[1]->m_fInclin);
			a = (a0 + forme->m_pProfils[i+1]->m_fInclin) * 0.5f;
			//angle vrillage
			v = (forme->m_pProfils[i]->m_fWash + forme->m_pProfils[i]->m_fWash) * 0.5f;
			//coeff morphing
			m = (forme->m_pProfils[i]->m_fMorph + forme->m_pProfils[i]->m_fMorph) * 0.5f;
			// time to calcate profile ballone
			ProfilGeom* pgCur = getProfile(gfd, forme, i);
			ProfilGeom* pgCurBal = getBalloneProfilGeom(pgCur, gfd->ballonement->kChord->Element(i, 0), gfd->ballonement->kMf->Element(i, 0), EpaiRel, gfd->ballonement->wN->Element(i, 0), gfd->ballonement->dyw->Element(i, 0));
			double l = abs (pgCur->ExtProf->Element(pgCur->ExtProf->GetLignes() - 1, 0) - pgCur->ExtProf->Element(0, 0));
			double xv = l * (100.0f - (gfd->PosPinceBF[0])) * 0.01f;
			ProfilGeom* pg = getProfilGeomTailDown(pgCurBal, pgCur, xv, gfd->ballonement->powerTail->Element(i, 0));
			//  -- time to calcate profile ballone
			coeffx = LongNerv/100.0f;
			double EpaiRelCur = EpaisseurRelative(pgCur->ExtProf, pgCur->IntProf);
			coeffy = LongNerv*EpaiRel/(EpaiRelCur*100.0f);
			//printf ("\n B longNerv=%f", LongNerv);
			//printf ("\n B coeffx=%f coeffy=(%f)", coeffx, coeffy);
			//printf ("\n B EpaiRelCur=%f", EpaiRelCur);
			double _x=0, _y=0, _z=0;
			for (j=0; j < pg->IntProf->GetLignes(); j++)
			{
				xp = pg->IntProf->Element(j, 0) * coeffx;
				yp = pg->IntProf->Element(j, 1) * coeffy;
				_x = xo+(yp-xp*(double)sin(v))*(double)cos(a);
				_y = yo+(yp-xp*(double)sin(v))*(double)sin(a);
				_z = zo-(xp*(double)cos(v));
				(*XInt)->SetElement(2*i+1, j, _x);
				(*YInt)->SetElement(2*i+1, j, _y);
				(*ZInt)->SetElement(2*i+1, j, _z);
			//	printf ("\n int %d (%f, %f, %f)", j, (*XInt)->Element(2*i+1, j), (*YInt)->Element(2*i+1, j), (*ZInt)->Element(2*i+1, j));
			}
			if (pg->IntProf->GetLignes() < IntProfCent->GetLignes()){
				for ( j = pg->IntProf->GetLignes(); j < IntProfCent->GetLignes(); j++) {
					(*XInt)->SetElement(2*i+1, j, _x);
					(*YInt)->SetElement(2*i+1, j, _y);
					(*ZInt)->SetElement(2*i+1, j, _z);
			//		printf ("\n +int %d (%f, %f, %f)", j, (*XInt)->Element(2*i+1, j), (*YInt)->Element(2*i+1, j), (*ZInt)->Element(2*i+1, j));
				}
			}

			//printf ("\n\n ext(%d) int(%d)", pg->ExtProf->GetLignes(), pg->IntProf->GetLignes());
			//printf ("\n\n BECOME ext(%d) int(%d)", ExtProfCent->GetLignes(), IntProfCent->GetLignes());
			for (j=0; j < pg->ExtProf->GetLignes(); j++)
			{
/*				if (isPercent) {
					xp = ExtProfCentN->Element(j,0)*coeffx;
					yp = ExtProfCentN->Element(j,1)*coeffyCent*m
							+ ExtProfBoutN->Element(j,1)*coeffyBout*(1.0f-m);
				} else { */
				xp = pg->ExtProf->Element(j, 0) * coeffx;
				yp = pg->ExtProf->Element(j, 1) * coeffy;
				//}
				_x = xo+(yp-xp*(double)sin(v))*(double)cos(a);
				_y = yo+(yp-xp*(double)sin(v))*(double)sin(a);
				_z = zo-(xp*(double)cos(v));
				(*XExt)->SetElement(2*i+1, j,_x);
				(*YExt)->SetElement(2*i+1, j, _y);
				(*ZExt)->SetElement(2*i+1, j, _z);
			//	printf ("\n ext %d (%f, %f, %f)",j, (*XExt)->Element(2*i+1, j), (*YExt)->Element(2*i+1, j), (*ZExt)->Element(2*i+1, j));
			}
			if (pg->ExtProf->GetLignes() < ExtProfCent->GetLignes()){
				for ( j = pg->ExtProf->GetLignes(); j < ExtProfCent->GetLignes(); j++) {
					(*XExt)->SetElement(2*i+1, j, _x);
					(*YExt)->SetElement(2*i+1, j, _y);
					(*ZExt)->SetElement(2*i+1, j, _z);
			//		printf ("\n +ext %d (%f, %f, %f)",j, (*XExt)->Element(2*i+1, j), (*YExt)->Element(2*i+1, j), (*ZExt)->Element(2*i+1, j));
				}
			}

			//printf ("\n --%d nervure ballone", i);        
		}
	}
}
