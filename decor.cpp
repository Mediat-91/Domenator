#include <string.h>
#include <stdio.h>
#include "stdafx.h"
#include <time.h>
//#include "globales.h"
#include "functions.h"
#include "Types.h"

#define _max(a, b) (a < b ? b : a)
#define _min(a, b) (a < b ? a : b)

extern DomeParameters myDome;
extern Tube* mySpiraleData;
extern Tube* myFriseData;
extern Tube* mySnakeData;
extern Tube* myCombData;
extern TestSlope TestSlopeFunctionListe[];
extern Parametric_4 FriseFunction;
#define sgn(a) (a < 0 ? -1 : (a > 0 ? 1 : 0))

float ScalarProduct(Vector3D* V, GLpoint *C, int i, int j, int nMaxI, int nMaxJ)
{
	double x, y, z;
	float u, v;
	u = (float)i / (float)(nMaxI - 1);
	u *= (float)(2 * _PI);
	v = (float)j / (float)(nMaxJ - 1);

	Calcul(u, v, &x, &y, &z, false);

	float p = (float)((x - C->x) * V->x + (y - C->y) * V->y + (z -C->z) * V->z);

//	printf("%15.12f\n", p / sqrtf((float)(V->x * V->x + V->y * V->y + V->z * V->z)));
	return fabsf(p) / sqrtf((float)(V->x * V->x + V->y * V->y + V->z * V->z));
}
double sPower(double x, double y)
{
	return sgn(x) * pow(sgn(x) * x, y);
}
double NextValue(double v, int j, int nMaxJ)
{
	static float fTotalLength;
	static float fPreviousLength;
	float x0, y0;
	float x1, y1;
	float fIncrement = 0.00001f;

	if (j == nMaxJ - 1)
	{
		fTotalLength = ComputeSlopeLength();
		fPreviousLength = 0;
		return v;
	}
	x0 = (float)SlopeFunction(0, myDome.fSlopeParameterXY, SLOPE_XY);
	y0 = (float)SlopeFunction(0, myDome.fSlopeParameterZ, SLOPE_Z);

	float fLimit = (float)j * fTotalLength / (float)(nMaxJ - 1);

	while (fPreviousLength < fLimit)
	{
		v -= fIncrement;
		x1 = (float)SlopeFunction(v, myDome.fSlopeParameterXY, SLOPE_XY);
		y1 = (float)SlopeFunction(v, myDome.fSlopeParameterZ, SLOPE_Z);
		fPreviousLength += sqrtf((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		x0 = x1;
		y0 = y1;

	}
	return v;
}


void SMA(Tube* Data, int nStart, int nMax)
{
	int nLargeur = 10;
	float x, y, z;
	Tube* Local;

	Local = (Tube*)malloc((nMax - nStart) * sizeof(Tube));
//	printf("%d --> %d\n", nStart, nMax);
	for (int i = 0; i < (nMax - nStart); i++)
	{
		Local[i].Centre.x = Data[nStart + i].Centre.x;
		Local[i].Centre.y = Data[nStart + i].Centre.y;
		Local[i].Centre.z = Data[nStart + i].Centre.z;
		Local[i].Radius = Data[nStart + i].Radius;
//		printf("(%f, %f, %f)\n", Local[i].Centre.x, Local[i].Centre.y, Local[i].Centre.z);
	}
//	printf("\ndeuxieme partie\n");

	for (int i = nStart; i < nMax; i++)
	{
		x = 0;
		y = 0;
		z = 0;
		int cc = _max(nStart, i - nLargeur);
		int dd = _min(nMax - 1, i + nLargeur);
		for (int j = _max(nStart, i - nLargeur); j <= _min(nMax - 1, i + nLargeur); j++)
		{
			x += Local[j - nStart].Centre.x;
			y += Local[j - nStart].Centre.y;
			z += Local[j - nStart].Centre.z;
		}
		int nb = (_min(nMax - 1, i + nLargeur) - _max(nStart, i - nLargeur) + 1);
		Data[i].Centre.x = x / nb;
		Data[i].Centre.y = y / nb;
		Data[i].Centre.z = z / nb;
//		printf("(%f, %f, %f)\n", Data[i].Centre.x, Data[i].Centre.y, Data[i].Centre.z);
	}
//	printf("\nNouvelle dent\n");
	Data[nStart].Centre.x = Local[0].Centre.x;
	Data[nStart].Centre.y = Local[0].Centre.y;
	Data[nStart].Centre.z = Local[0].Centre.z;

	Data[nMax - 1].Centre.x = Local[nMax - nStart - 1].Centre.x;
	Data[nMax - 1].Centre.y = Local[nMax - nStart - 1].Centre.y;
	Data[nMax - 1].Centre.z = Local[nMax - nStart - 1].Centre.z;

	free(Local);
}
////////////////////////////////////////////////
//               Variantes                    //
////////////////////////////////////////////////

int trueCalculComb()
{
	int nMaxI, nMaxJ;

	Vector3D AB;

	double u, v;
	double Xx0, Xy0, Xz0;

	int i, j;

	int n;

	double u_A, u_B, y;

	GLpoint c0, c1;
	float fShrink;

	int   nComb = myDome.nComb;
	int   nSteps = 5 * myDome.nCombStep;
	float fCombSpacing = myDome.fCombSpacing;
	float fCombRadius = myDome.fCombRadius;
	int   nCombShrink = myDome.nCombShrink;
	int   nCombToward = myDome.nCombToward;

	int nSideI = nSteps * myDome.nMeridianMain * (myDome.nMeridianMinor + 1);

	if (nComb == 0)
		return true;

	nMaxI = ((myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides);
	nMaxJ = (myDome.nParallelMinor + 1) * (myDome.nParallelMain - 1) + 1;
	nMaxI *= nSteps;
	nMaxJ *= nSteps;

	if (myCombData)
		free(myCombData);

	myCombData = (Tube*)malloc(10*nMaxI * nComb * sizeof(Tube));
	n = 0;
	for (int nS = 0; nS < myDome.nBaseSides; nS++)
	{
		i = nSideI * nS;
		u = (float)i/(float)nMaxI;
		u *= 2 * _PI;
		u_A = i;
		Calcul(u, 1, &Xx0, &Xy0, &Xz0, false);

		AB.x = -Xx0;
		AB.y = -Xy0;
		AB.z = -Xz0;

//		printf(" A = %f ->(%f, %f, %f)\n", u, AB.x, AB.y, AB.z);


		i = nSideI * (nS + 1);

		u = (float)(i) / (float)nMaxI;
		u *= 2 * _PI;
		u_B = (float)(i);
		Calcul(u, 1, &Xx0, &Xy0, &Xz0, false);

		AB.x += Xx0;
		AB.y += Xy0;
		AB.z += Xz0;
		// AB est calculé !
//		printf("AB = %f ->(%f, %f, %f)\n", u, AB.x, AB.y, AB.z);



		for (int k = 0; k < nComb; k++)
		{
			i = nSideI * nS + ((k + 1) * nSteps * myDome.nMeridianMain * (myDome.nMeridianMinor + 1)) / (nComb + 1);
			j = nMaxJ - 1;

			y = (2.0 * (float)i - (u_A + u_B)) / (u_B - u_A);
			i = (int)((sPower(y, fCombSpacing) * (u_B - u_A) + (u_A + u_B)) / 2);

			u = (float)i / (float)(nMaxI - 1);
			u *= 2 * _PI;

			v = 1;



			Calcul(u, v, &Xx0, &Xy0, &Xz0, false);
			c0.x = (float)Xx0;
			c0.y = (float)Xy0;
			c0.z = (float)Xz0;

//			printf("c0 ->%f(%f, %f, %f)\n", u, c0.x, c0.y, c0.z);

			/*
				ooo
				.xo
				...
			*/
			int I, J, Iprec, Jprec, Iprecprec, Jprecprec;
			float a, b, c, d, e, f;

			ChargeShrink(nCombShrink);

			e = 0;

			fShrink = (float)ShrinkComb(j, nMaxJ);

			myCombData[n].Centre = c0;
			myCombData[n].Radius = fCombRadius * fShrink;;
			n++;
			Iprec = -1;
			Jprec = -1;
			Iprecprec = -1;
			Jprecprec = -1;

//			printf("c0 ->(%f, %f, %f)\n", c0.x, c0.y, c0.z);

			while ((i != Iprecprec || j != Jprecprec) && e < 10.01 && i > nS * nSideI && i < (nS + 1) * nSideI && j > 0)
				{
				Iprecprec = Iprec;
				Jprecprec = Jprec;
				Iprec = i;
				Jprec = j;
				a = ScalarProduct(&AB, &c0, i - 1, j - 1, nMaxI, nMaxJ);
				b = ScalarProduct(&AB, &c0, i, j - 1, nMaxI, nMaxJ);
				c = ScalarProduct(&AB, &c0, i + 1, j - 1, nMaxI, nMaxJ);
				d = ScalarProduct(&AB, &c0, i + 1, j, nMaxI, nMaxJ);
				f = ScalarProduct(&AB, &c0, i - 1, j, nMaxI, nMaxJ);

				// Mais il faut empecher les boucles



				if (a < b)
				{
					e = a;
					I = i - 1;
					J = j - 1;
				}
				else
				{
					e = b;
					I = i;
					J = j - 1;
				}
				if (c < e)
				{
					e = c;
					I = i + 1;
					J = j - 1;
				}
				if (d < e)
				{
					e = d;
					I = i + 1;
					J = j;
				}
				if (f < e)
				{
					e = f;
					I = i - 1;
					J = j;
				}


				//			printf("(%3d, %3d), a = %15.12f, b = %15.12f, c = %15.12f, d = %15.12f, Min = %15.12f\n", i, j, a, b, c, d, e);
				i = I;
				j = J;

				u = (float)i / (float)(nMaxI - 1);
				u *= 2 * _PI;
				v = (float)j / (float)(nMaxJ - 1);

				Calcul(u, v, &Xx0, &Xy0, &Xz0, false);
				c1.x = (float)Xx0;
				c1.y = (float)Xy0;
				c1.z = (float)Xz0;
//				printf("c1 ->(%f, %f, %f)\n", c1.x, c1.y, c1.z);

				fShrink = (float)ShrinkComb(j, nMaxJ);

				myCombData[n].Centre = c1;
				myCombData[n].Radius = fCombRadius * fShrink;;
				n++;
			}
			
//			printf("(%15.12f, %15.12f, %15.12f)\n", myCombData[0].Centre.x, myCombData[0].Centre.y, myCombData[0].Centre.z);
//			printf("(%15.12f, %15.12f, %15.12f)\n", myCombData[1].Centre.x, myCombData[1].Centre.y, myCombData[1].Centre.z);

			myCombData[n].Centre = c0;
			myCombData[n].Radius = -1;
			n++;
		}
	}
	myCombData[n].Centre = c0;
	myCombData[n].Radius = -2;

//	printf("n = %d\n", n);
	return n;


}
int VarCalculSnake()
{
	int nMaxI;

	double u, v;
	double Xx0, Xy0, Xz0;
	GLpoint c0;
	float fStart;

	int   nSnake = myDome.nSnake;
	float fSnakeTour = (float)(myDome.fSnakeTour * 2 * _PI);
	int   nSnakeEvolve = myDome.nSnakeEvolve;
	float fSnakeShift = myDome.fSnakeShift;
	float fSnakeRadius = myDome.fSnakeRadius;
	int   nSnakeShrink = myDome.nSnakeShrink;
	int   nSnakeChirale = myDome.nSnakeChirale;
	int   nSnakeToward = myDome.nSnakeToward;
	float fSnakeOffset = myDome.fSnakeOffset;

	int   k;
	float fSens;
	float fShrink;

	if (nSnake == 0)
		return true;

	nMaxI = myDome.nSnakeStep * ((myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides + 1);


	if (mySnakeData)
		free(mySnakeData);

	int nTours = ((nSnake) / myDome.nBaseSides);

	mySnakeData = (Tube*)malloc((nSnakeChirale == CHIRAL_BOTH ? 2 : 1) * nMaxI * nSnake * sizeof(Tube));
	if (mySnakeData)
	{
		ChargeShrink(nSnakeShrink);

		for (int nbChirale = 0; nbChirale < (nSnakeChirale == CHIRAL_BOTH ? 2 : 1); nbChirale++)
		{
			for (int nK = 0; nK < nSnake; nK++)
			{
				nTours = ((nSnake) / myDome.nBaseSides);
				k = ((nK % myDome.nBaseSides) * myDome.nBaseSides) / (nK < nTours * myDome.nBaseSides ? myDome.nBaseSides : nSnake % myDome.nBaseSides);
				int nTourCourant = nK / myDome.nBaseSides;
				//				printf("%d    ", k);
				//				if (nK == floor((float)(k * nSnake) / (float)myDome.nBaseSides))
				{
					if (nSnakeChirale == CHIRAL_RIGHT || nbChirale == 1)
						fSens = -1;
					else
						fSens = 1;
					nTours = ((nSnake - 1) / myDome.nBaseSides);
					fSnakeOffset = myDome.fSnakeOffset + (float)(2 * nTourCourant) / (float)(nTours + 1);

					fStart = (float)(k * 2.0 * _PI) / (float)myDome.nBaseSides + (float)(fSnakeOffset * _PI) / (float)myDome.nBaseSides;

					for (int i = 0; i < nMaxI; i++)
					{
						u = ((sin((double)(fSens * fSnakeTour * i) / (double)(nMaxI - 1) + _PI + myDome.fSnakeShift * _PI / 2) + 1.0) / 2.0);
						u = 1-cos(_PI*(double)i / (double)(2*(nMaxI - 1)));
						u *= 2.0 * _PI / myDome.nBaseSides;
						u += fStart;


						v = 1 - TestSlopeFunctionListe[myDome.nSnakeSlope].Function((double)i / (double)(nMaxI - 1), myDome.fSnakeSlopeParameter);
						if (isnan(v))
						{
							v = 1;
							printf("isnan\n");
						}

						switch (myDome.nSnakeEvolve)//Upscale
						{
						case SHRINK_SQUARE_ROOT:
							v = sqrt(v);
							break;
						case SHRINK_SQUARE:
							v = v * v;
							break;
						}
						Calcul(u, v, &Xx0, &Xy0, &Xz0, false);
						c0.x = (float)Xx0;
						c0.y = (float)Xy0;
						c0.z = (float)Xz0;


						fShrink = (float)ShrinkSnake(i, nMaxI);


						mySnakeData[nbChirale * nSnake * nMaxI + nK * nMaxI + i].Centre = c0;
						mySnakeData[nbChirale * nSnake * nMaxI + nK * nMaxI + i].Radius = fSnakeRadius * fShrink;
					}
				}

			}
		}

		//		printf("\n");
		return nMaxI;
	}
	return 0;
}
int CalculSnakeold()
{
	int nMaxI;

	double u, v;
	double Xx0, Xy0, Xz0;
	GLpoint c0;
	float fStart;

	int   nSnake = myDome.nSnake;
	float fSnakeTour = (float)(myDome.fSnakeTour * 2 * _PI);
	int   nSnakeEvolve = myDome.nSnakeEvolve;
	float fSnakeShift = myDome.fSnakeShift;
	float fSnakeRadius = myDome.fSnakeRadius;
	int   nSnakeShrink = myDome.nSnakeShrink;
	int   nSnakeChirale = myDome.nSnakeChirale;
	int   nSnakeToward = myDome.nSnakeToward;
	float fSnakeOffset = myDome.fSnakeOffset;

	int   nK;
	float fSens;
	float fShrink;

	if (nSnake == 0)
		return true;
	if (nSnake > myDome.nBaseSides)
		nSnake = myDome.nBaseSides;

	nMaxI = myDome.nSnakeStep * ((myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides + 1);


	if (mySnakeData)
		free(mySnakeData);

	mySnakeData = (Tube*)malloc((nSnakeChirale == CHIRAL_BOTH ? 2 : 1) * nMaxI * nSnake * sizeof(Tube));
	if (mySnakeData)
	{
		ChargeShrink(nSnakeShrink);

		for (int nbChirale = 0; nbChirale < (nSnakeChirale == CHIRAL_BOTH ? 2 : 1); nbChirale++)
		{
			nK = 0;
			for (int k = 0; k < myDome.nBaseSides; k++)
			{
				if (nK == floor((float)(k * nSnake) / (float)myDome.nBaseSides))
				{
					if (nSnakeChirale == CHIRAL_RIGHT || nbChirale == 1)
						fSens = -1;
					else
						fSens = 1;
					fStart = (float)(k * 2.0 * _PI) / (float)myDome.nBaseSides + (float)(fSnakeOffset * _PI) / (float)myDome.nBaseSides;
					for (int i = 0; i < nMaxI; i++)
					{
						u = ((sin((double)(fSens * fSnakeTour * i) / (double)(nMaxI - 1) + _PI + myDome.fSnakeShift * _PI / 2) + 1.0) / 2.0);
						u *= 2.0 * _PI / myDome.nBaseSides;
						//					u = _PI / myDome.nBaseSides - u; //OK pour k = 0
						u += fStart;

						//						v = 1.0 - (double)i / (double)(nMaxI - 1);

						//						v = 1.0 - fabs(cos(2.0 * (double)i * _PI / ((double)nMaxI - 1.0)));

						//						v = 1 - SlopeCasablanca((double)i / (double)(nMaxI - 1), 2);

						v = 1 - TestSlopeFunctionListe[myDome.nSnakeSlope].Function((double)i / (double)(nMaxI - 1), myDome.fSnakeSlopeParameter);
						if (isnan(v))
							v = 1;

						switch (myDome.nSnakeEvolve)//Upscale
						{
						case SHRINK_SQUARE_ROOT:
							v = sqrt(v);
							break;
						case SHRINK_SQUARE:
							v = v * v;
							break;
						}
						Calcul(u, v, &Xx0, &Xy0, &Xz0, false);
						c0.x = (float)Xx0;
						c0.y = (float)Xy0;
						c0.z = (float)Xz0;


						fShrink = (float)ShrinkSnake(i, nMaxI);


						mySnakeData[nbChirale * nSnake * nMaxI + nK * nMaxI + i].Centre = c0;
						mySnakeData[nbChirale * nSnake * nMaxI + nK * nMaxI + i].Radius = fSnakeRadius * fShrink;
					}
					nK++;
				}

			}
		}

		return nMaxI;
	}
	return 0;
}

////////////////////////////////////////////////
//               Variantes                    //
////////////////////////////////////////////////


int CalculSpirale()
{
	float fSpiraleTour = myDome.fSpiraleTour;
	int   nSpirale = myDome.nSpirale;
	int   nSteer = (myDome.nSpiraleChirale == 2 ? 2 : 1);
	int   nSpiraleChirale = myDome.nSpiraleChirale;
	float fSpiraleDelay = myDome.fSpiraleDelay;
	float fSpiraleDecalage = myDome.fSpiraleDecalage;
	float fSpiraleSpeed = myDome.fSpiraleSpeed;
	int   nSpiraleShrink = myDome.nSpiraleShrink;	// none=0, sqrt=1, Id, carre
	float fSpiraleRadius = myDome.fSpiraleRadius;
	int   nSpiraleStep = myDome.nSpiraleStep;
	int   nSpiraleToward = myDome.nSpiraleToward;

	float fShrink;

	float u, v;
	double Xx0, Xy0, Xz0;
	GLpoint c0;

	float fStart;
	int   nStart;

	int nMaxI = nSpiraleStep * (myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides + 1;
	int nMaxSpirale = (int)((nMaxI - 1) * fSpiraleTour + 1);

	if (mySpiraleData)
		free(mySpiraleData);
	mySpiraleData = (Tube*)malloc(nSpirale * nSteer * nMaxSpirale * sizeof(Tube));
	if (mySpiraleData)
	{
		ChargeShrink(nSpiraleShrink);

		//		printf("%d - %d\n", (nSpiraleChirale % 2 == 0 ? 0 : nSpirale), (nSpiraleChirale == CHIRAL_LEFT ? nSpirale : 2 * nSpirale));

		for (int k = (nSpiraleChirale % 2 == 0 ? 0 : nSpirale); k < (nSpiraleChirale == CHIRAL_LEFT ? nSpirale : 2 * nSpirale); k++)
		{
			fStart = (float)((k % nSpirale) * 2 * _PI / nSpirale) + (k >= nSpirale ? (float)(2.0 * _PI) - fSpiraleDelay + (fSpiraleDecalage) * (float)(2 * _PI) / (float)nSpirale : fSpiraleDelay);
			nStart = nMaxSpirale * (k % nSpirale) + (nSpiraleChirale == CHIRAL_BOTH && k >= nSpirale ? nMaxSpirale * nSpirale : 0);
			u = fStart + (float)0.0;
			v = 0.0;

			Calcul(u, v, &Xx0, &Xy0, &Xz0, false);
			c0.x = (float)Xx0;
			c0.y = (float)Xy0;
			c0.z = (float)Xz0;

			mySpiraleData[nStart].Centre = c0;
			mySpiraleData[nStart].Radius = fSpiraleRadius * (float)ShrinkSpirale(0, nMaxSpirale);



			for (int i = 1; i < nMaxSpirale; i++)
			{
				//				u = fStart + (float)(2.0*_PI)*(float)(i % (nMaxI - 1)) / (float)(nMaxI - 1);

				u = fStart + (float)(2.0 * (float)_PI) * powf((float)(i % (nMaxI - 1)) / (float)(nMaxI - 1), fSpiraleSpeed);

				if (k >= nSpirale)
					u = myMod(-u, 2.0f * (float)_PI);
				//				v = (float)SlopeStraight((double)i / ((double)nMaxSpirale - 1), fSpiraleSpeed);

				v = (float)TestSlopeFunctionListe[myDome.nSpiraleSlope].Function((double)i / ((double)nMaxSpirale - 1), myDome.fSpiraleSlopeParameter);
				if (isnan(v))
					v = 0;

				/*				if (u < 0)
									printf("u = %lf\n", u);
								if (v < 0 || v > 1)
									printf("v = %lf\n", v);*/
				if (v > 1)
					v = 1;
				Calcul(u, v, &Xx0, &Xy0, &Xz0, false);

				c0.x = (float)Xx0;
				c0.y = (float)Xy0;
				c0.z = (float)Xz0;
				fShrink = (float)ShrinkSpirale(i, nMaxSpirale);
				mySpiraleData[nStart + i].Centre = c0;
				mySpiraleData[nStart + i].Radius = fSpiraleRadius * fShrink;

			}
		}
	}
	else
		return 0;
	return nMaxI;
}
int CalculFrise()
{
	int nMaxI = myDome.nFriseStep * ((myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides + 1);
	float u, v;
	float fShrink;
	double Xx0, Xy0, Xz0;
	GLpoint c0;

	int   nFrise = myDome.nFrise;
	float fFriseHeight = myDome.fFriseHeight;
	float fFrisePerSide = myDome.fFrisePerSide;
	float fFriseDepart = myDome.fFriseDepart;
	float fFriseShift = myDome.fFriseShift;
	float fFriseRadius = myDome.fFriseRadius;
	int   nFriseShrink = myDome.nFriseShrink;
	float fFrisePower = myDome.fFrisePower;
	float fFriseOffset = myDome.fFriseOffset * 2.0f * (float)_PI / myDome.nBaseSides;
	float fFriseGap = myDome.fFriseGap;
	int   nFriseVerse = 2 * myDome.bFriseInverse + myDome.bFriseReverse;
	float fFriseSpacing = myDome.fFriseSpacing;
	ChargeFrise(myDome.nFriseFunction);


	if (myFriseData)
		free(myFriseData);

	myFriseData = (Tube*)malloc(nMaxI * nFrise * sizeof(Tube));
	if (myFriseData)
	{
		ChargeShrink(nFriseShrink);



		for (int k = 0; k < nFrise; k++)
		{
			fFriseOffset = (myDome.fFriseOffset + (float)k * fFriseSpacing) * 2.0f * (float)_PI / myDome.nBaseSides;
			u = (float)0.0;
//			v = (float)1.0 - (fFriseHeight * (float)FriseFunction((u - fFriseOffset), myDome.nBaseSides, fFrisePerSide, fFrisePower, true) + (k + 1) * (fFriseDepart + fFriseShift * fabsf(sinf(myDome.nBaseSides * fFrisePerSide * (u - fFriseOffset) / 2))));

			switch (nFriseVerse)
			{
			case 0: v = (float)1.0 - (fFriseHeight * (float)FriseFunction((u - fFriseOffset), myDome.nBaseSides, fFrisePerSide, fFrisePower, true) + fFriseDepart + (k) * (fFriseGap + fFriseShift * /*fabsf */ (sinf(myDome.nBaseSides * fFrisePerSide * (u - fFriseOffset) / 2))));
				break;
			case 1:	v = (float)1.0 - (fFriseHeight * (1 - (float)FriseFunction((u - fFriseOffset), myDome.nBaseSides, fFrisePerSide, fFrisePower, true)) + fFriseDepart + (k) * (fFriseGap + fFriseShift * /*fabsf */(sinf(myDome.nBaseSides * fFrisePerSide * (u - fFriseOffset) / 2))));
				break;
			case 2:	v = (fFriseHeight * (float)FriseFunction((u - fFriseOffset), myDome.nBaseSides, fFrisePerSide, fFrisePower, true) + fFriseDepart + (k) * (fFriseGap + fFriseShift * /*fabsf */(sinf(myDome.nBaseSides * fFrisePerSide * (u - fFriseOffset) / 2))));
				break;
			case 3:	v = (fFriseHeight * (1 - (float)FriseFunction((u - fFriseOffset), myDome.nBaseSides, fFrisePerSide, fFrisePower, true)) + fFriseDepart + (k) * (fFriseGap + fFriseShift * /*fabsf */(sinf(myDome.nBaseSides * fFrisePerSide * (u - fFriseOffset) / 2))));
			}
			Calcul(u, v, &Xx0, &Xy0, &Xz0, false);
			c0.x = (float)Xx0;
			c0.y = (float)Xy0;
			c0.z = (float)Xz0;

			fShrink = (float)ShrinkFrise(0, nMaxI);

			myFriseData[k * nMaxI].Centre = c0;
			myFriseData[k * nMaxI].Radius = fFriseRadius * fShrink;

			for (int i = 1; i < nMaxI; i++)
			{
				u = (float)(2.0 * _PI) * (float)(i % (nMaxI - 1)) / (float)(nMaxI - 1);// +(float)k * fFriseShift / (float)(nMaxI - 1);
//				v = (float)1.0 - (fFriseHeight * (float)FriseFunction((u - fFriseOffset), myDome.nBaseSides, fFrisePerSide, fFrisePower, true) + (k + 1) * (fFriseDepart + fFriseShift * fabsf(sinf(myDome.nBaseSides * fFrisePerSide * (u - fFriseOffset) / 2))));
				switch (nFriseVerse)
				{
				case 0: 
					v = (float)1.0 - (fFriseHeight * (float)FriseFunction((u - fFriseOffset), myDome.nBaseSides, fFrisePerSide, fFrisePower, true) + fFriseDepart + (k) * (fFriseGap + fFriseShift * /*fabsf*/(sinf(myDome.nBaseSides * fFrisePerSide * (u - fFriseOffset) / 2))));
					break;
				case 1:
					v = (float)1.0 - (fFriseHeight * (1 - (float)FriseFunction((u - fFriseOffset), myDome.nBaseSides, fFrisePerSide, fFrisePower, true)) + fFriseDepart + (k) * (fFriseGap + fFriseShift * /*fabsf*/(sinf(myDome.nBaseSides * fFrisePerSide * (u - fFriseOffset) / 2))));
					break;
				case 2:
					v = (fFriseHeight * (float)FriseFunction((u - fFriseOffset), myDome.nBaseSides, fFrisePerSide, fFrisePower, true) + fFriseDepart + (k) * (fFriseGap + fFriseShift * /*fabsf*/(sinf(myDome.nBaseSides * fFrisePerSide * (u - fFriseOffset) / 2))));
					break;
				case 3:
					v = (fFriseHeight * (1 - (float)FriseFunction((u - fFriseOffset), myDome.nBaseSides, fFrisePerSide, fFrisePower, true)) + fFriseDepart + (k) * (fFriseGap + fFriseShift * /*fabsf*/(sinf(myDome.nBaseSides * fFrisePerSide * (u - fFriseOffset) / 2))));
				}
				Calcul(u, v, &Xx0, &Xy0, &Xz0, false);

				c0.x = (float)Xx0;
				c0.y = (float)Xy0;
				c0.z = (float)Xz0;
				fShrink = (float)ShrinkFrise(i, nMaxI);
				myFriseData[k * nMaxI + i].Centre = c0;
				myFriseData[k * nMaxI + i].Radius = fFriseRadius * fShrink;
			}
		}
		return nMaxI;
	}
	return 0;
}
int CalculSnake()
{
	int nMaxI;

	double u, v;
	double Xx0, Xy0, Xz0;
	GLpoint c0;
	float fStart;

	int   nSnake = myDome.nSnake;
	float fSnakeTour = (float)(myDome.fSnakeTour * 2 * _PI);
	int   nSnakeEvolve = myDome.nSnakeEvolve;
	float fSnakeShift = myDome.fSnakeShift;
	float fSnakeRadius = myDome.fSnakeRadius;
	int   nSnakeShrink = myDome.nSnakeShrink;
	int   nSnakeChirale = myDome.nSnakeChirale;
	int   nSnakeToward = myDome.nSnakeToward;
	float fSnakeOffset = myDome.fSnakeOffset;

	int   k;
	float fSens;
	float fShrink;

	if (nSnake == 0)
		return true;

	nMaxI = myDome.nSnakeStep * ((myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides + 1);


	if (mySnakeData)
		free(mySnakeData);

	int nTours = ((nSnake )/ myDome.nBaseSides);

	mySnakeData = (Tube*)malloc((nSnakeChirale == CHIRAL_BOTH ? 2 : 1) * nMaxI * nSnake * sizeof(Tube));
	if (mySnakeData)
	{
		ChargeShrink(nSnakeShrink);

		for (int nbChirale = 0; nbChirale < (nSnakeChirale == CHIRAL_BOTH ? 2 : 1); nbChirale++)
		{
			for (int nK = 0; nK < nSnake; nK++)
			{
				nTours = ((nSnake) / myDome.nBaseSides);
				k = ((nK % myDome.nBaseSides) * myDome.nBaseSides) / (nK < nTours * myDome.nBaseSides ? myDome.nBaseSides : nSnake % myDome.nBaseSides);
				int nTourCourant = nK / myDome.nBaseSides;
//				printf("%d    ", k);
//				if (nK == floor((float)(k * nSnake) / (float)myDome.nBaseSides))
				{
					if (nSnakeChirale == CHIRAL_RIGHT || nbChirale == 1)
						fSens = -1;
					else
						fSens = 1;
					nTours = ((nSnake-1) / myDome.nBaseSides);
					fSnakeOffset = myDome.fSnakeOffset + (float)(2*nTourCourant) / (float)(nTours+1);

					fStart = (float)(k * 2.0 * _PI) / (float)myDome.nBaseSides + (float)(fSnakeOffset * _PI) / (float)myDome.nBaseSides;

					for (int i = 0; i < nMaxI; i++)
					{
						u = ((sin((double)(fSens * fSnakeTour * i) / (double)(nMaxI - 1) + _PI + myDome.fSnakeShift * _PI / 2) + 1.0) / 2.0);
						u *= 2.0 * _PI / myDome.nBaseSides;
						u += fStart;


						v = 1 - TestSlopeFunctionListe[myDome.nSnakeSlope].Function((double)i / (double)(nMaxI - 1), myDome.fSnakeSlopeParameter);
						if (isnan(v))
						{
							v = 1;
							printf("isnan\n");
						}

						switch (myDome.nSnakeEvolve)//Upscale
						{
						case SHRINK_SQUARE_ROOT:
							v = sqrt(v);
							break;
						case SHRINK_SQUARE:
							v = v * v;
							break;
						}
						Calcul(u, v, &Xx0, &Xy0, &Xz0, false);
						c0.x = (float)Xx0;
						c0.y = (float)Xy0;
						c0.z = (float)Xz0;


						fShrink = (float)ShrinkSnake(i, nMaxI);


						mySnakeData[nbChirale * nSnake * nMaxI + nK * nMaxI + i].Centre = c0;
						mySnakeData[nbChirale * nSnake * nMaxI + nK * nMaxI + i].Radius = fSnakeRadius * fShrink;
					}
				}

			}
		}

//		printf("\n");
		return nMaxI;
	}
	return 0;
}
int CalculComb()
{
	int nMaxI, nMaxJ;

	Vector3D AB;

	double u, v;
	double Xx0, Xy0, Xz0;

	int i, j;

	int n;

	double u_A, u_B, y;

	GLpoint c0, c1;
	float fShrink;

	int   nGap = myDome.nCombThreadGap;
	int   nComb = myDome.nComb;
	int   nSteps = 5 * myDome.nCombStep;
	float fCombSpacing = myDome.fCombSpacing;
	float fCombRadius = myDome.fCombRadius;
	int   nCombShrink = myDome.nCombShrink;
	int   nCombToward = myDome.nCombToward;
	int   nMaxThread = (myDome.nCombType == COMB_UPRIGHT_BASE ? 1 : 2);
	int   nThread = myDome.nCombThread;
	int   nSideI = nSteps * myDome.nMeridianMain * (myDome.nMeridianMinor + 1);
	int nStart;
	float fDistanceMax = .01f;
	float fOffset = myDome.fCombOffset;

	if (nComb == 0)
		return true;

	nMaxI = ((myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides);
	nMaxJ = (myDome.nParallelMinor + 1) * (myDome.nParallelMain - 1) + 1;
	nMaxI *= nSteps;
	nMaxJ *= nSteps;
//	printf("\n---------------\n");
	if (myCombData)
		free(myCombData);

	myCombData = (Tube*)malloc(100 * nMaxI * nComb *nThread * sizeof(Tube));
	n = 0;
	for (int nS = 0; nS < myDome.nBaseSides; nS++)
	{
		for (int kT = 0; kT < nMaxThread; kT++)
		{
			i = nSideI * (nS + kT);
			u = (float)i / (float)nMaxI + _PI * (1 - 2 * kT) * fOffset / (float)myDome.nBaseSides;
			u *= 2 * _PI;
			u_A = (float)i;
			Calcul(u, 1, &Xx0, &Xy0, &Xz0, false);

			if (myDome.nCombType == COMB_PARALLEL_EDGE)
			{
				AB.x = Xy0;
				AB.y = -Xx0;
				AB.z = Xz0;
			}
			else
			{
				AB.x = -Xx0;
				AB.y = -Xy0;
				AB.z = -Xz0;
			}


			printf(" A = %f ->(%f, %f, %f)\n", u, AB.x, AB.y, AB.z);


			i = nSideI * (nS - kT + 1);

			u = (float)(i) / (float)nMaxI ;
			u *= 2 * _PI;
			u_B = (float)(i);
			if (myDome.nCombType == COMB_UPRIGHT_BASE)
			{
				Calcul(u, 1, &Xx0, &Xy0, &Xz0, false);

				AB.x += Xx0;
				AB.y += Xy0;
				AB.z += Xz0;
			}
			// AB est calculé !
	//		printf("AB = %f ->(%f, %f, %f)\n", u, AB.x, AB.y, AB.z);


			for (int nT = -nThread + 1; nT < nThread; nT += 2)
			{
				for (int k = 0; k < nComb; k++)
				{
					i = (nGap * nT) + nSideI * nS + ((k + 1) * nSteps * myDome.nMeridianMain * (myDome.nMeridianMinor + 1)) / (nComb + 1);
					j = nMaxJ - 1;

					y = (2.0 * (float)i - (u_A + u_B)) / (u_B - u_A);
					i = (int)((sPower(y, fCombSpacing) * (u_B - u_A) + (u_A + u_B)) / 2) ;
					u = (float)i / (float)(nMaxI - 1);
					u *= 2 * _PI;
					v = 1;
					v = NextValue(v, j, nMaxJ);



					Calcul(u, v, &Xx0, &Xy0, &Xz0, false);
					printf("n = %3d, (u, v) = (%f, %f)		(i, j) = (%d/%d, %d/%d)		X = (%f, %f, %f)\n", n, u, v, i, nMaxI, j, nMaxJ, Xx0, Xy0, Xz0);
					c0.x = (float)Xx0;
					c0.y = (float)Xy0;
					c0.z = (float)Xz0;

					//			printf("c0 ->%f(%f, %f, %f)\n", u, c0.x, c0.y, c0.z);

								/*
									ooo
									.xo
									...
								*/
					int I, J, Iprec, Jprec, Iprecprec, Jprecprec;
					float NordOuest, Nord, NordEst, Est, Ouest;
					float fDistance;
					ChargeShrink(nCombShrink);

					fShrink = (float)ShrinkComb(j, nMaxJ);
					nStart = n;

					myCombData[n].Centre = c0;
					myCombData[n].Radius = fCombRadius * fShrink;;
					n++;
					Iprec = -1;
					Jprec = -1;
					Iprecprec = -1;
					Jprecprec = -1;

					//			printf("c0 ->(%f, %f, %f)\n", c0.x, c0.y, c0.z);
					fDistance = fDistanceMax - 1;
					while (fDistance < fDistanceMax && i > nS * nSideI && i < (nS + 1) * nSideI && j > 0)
						{
						Iprecprec = Iprec;
						Jprecprec = Jprec;
						Iprec = i;
						Jprec = j;
						NordOuest = ScalarProduct(&AB, &c0, i - 1, j - 1, nMaxI, nMaxJ);
						Nord = ScalarProduct(&AB, &c0, i, j - 1, nMaxI, nMaxJ);
						NordEst = ScalarProduct(&AB, &c0, i + 1, j - 1, nMaxI, nMaxJ);
						Est = ScalarProduct(&AB, &c0, i + 1, j, nMaxI, nMaxJ);
						Ouest = ScalarProduct(&AB, &c0, i - 1, j, nMaxI, nMaxJ);

						// Mais il faut empecher les boucles

						fDistance = fDistanceMax;
						I = i;
						J = j;
						if (NordOuest < fDistance && (i - 1 != Iprecprec || j - 1 != Jprecprec))
						{
							fDistance = NordOuest;
							I = i - 1;
							J = j - 1;
						}
						if(Nord < fDistance && (i != Iprecprec || j - 1 != Jprecprec))
						{
							fDistance = Nord;
							I = i;
							J = j - 1;
						}
						if (NordEst < fDistance && (i + 1 != Iprecprec || j - 1 != Jprecprec))
						{
							fDistance = NordEst;
							I = i + 1;
							J = j - 1;
						}
						if (Est < fDistance && (i + 1!= Iprecprec || j != Jprecprec))
						{
							fDistance = Est;
							I = i + 1;
							J = j;
						}
						if (Ouest < fDistance && (i - 1 != Iprecprec || j != Jprecprec))
						{
							fDistance = Ouest;
							I = i - 1;
							J = j;
						}


//			printf("(%3d, %3d), a = %15.12f, b = %15.12f, c = %15.12f, d = %15.12f, Min = %15.12f\n", i, j, a, b, c, d, e);
//						printf("%f\n", fDistance);
						i = I;
						j = J;

						u = (float)i / (float)(nMaxI - 1);
						u *= 2 * _PI;
						v = (float)j / (float)(nMaxJ - 1);

						Calcul(u, v, &Xx0, &Xy0, &Xz0, false);
						c1.x = (float)Xx0;
						c1.y = (float)Xy0;
						c1.z = (float)Xz0;
						//				printf("c1 ->(%f, %f, %f)\n", c1.x, c1.y, c1.z);

						fShrink = (float)ShrinkComb(j, nMaxJ);
//						if (fDistance < fDistanceMax)
						{
							myCombData[n].Centre = c1;
							myCombData[n].Radius = fCombRadius * fShrink;
							n++;
						}
//						fDistance = 0;
					}
//					printf("fDistance = %f, fDistanceMax = %f , i = %d ,nS * nSideI = %d, (nS + 1) * nSideI = %d, j = %d\n", fDistance, fDistanceMax, i, nS * nSideI, (nS + 1) * nSideI, j);
					//			printf("(%15.12f, %15.12f, %15.12f)\n", myCombData[0].Centre.x, myCombData[0].Centre.y, myCombData[0].Centre.z);
					//			printf("(%15.12f, %15.12f, %15.12f)\n", myCombData[1].Centre.x, myCombData[1].Centre.y, myCombData[1].Centre.z);

//					if (nS == 1 && kT == 0)
					{
						SMA(myCombData, nStart, n);
//						printf("n = %5d\n", n);

					}
					myCombData[n].Centre = c0;
					myCombData[n].Radius = -1;
					n++;
//					printf("%d / %d\n", n, 100 * nMaxI * nComb * nThread);
				}
			}
		}
	}
	myCombData[n].Centre = c0;
	myCombData[n].Radius = -2;

	//	printf("n = %d\n", n);
	return n;


}

void DrawSpirale()
{
	int   nSpirale = myDome.nSpirale;
	float fSpiraleTour = myDome.fSpiraleTour;
	float fSpiraleSpeed = myDome.fSpiraleSpeed;
	float fSpiraleDelay = myDome.fSpiraleDelay;
	float fSpiraleDecalage = myDome.fSpiraleDecalage;
	int   nSpiraleChirale = myDome.nSpiraleChirale;
	int   nMaxI = myDome.nSpiraleStep * (myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides + 1;

	int   nMaxSpirale;
	int   nStart;


//	nMaxI = CalculSpirale();
	if (nMaxI > 0 && mySpiraleData)
	{

		glLineWidth(3);

		nMaxSpirale = (int)((nMaxI - 1) * fSpiraleTour + 1);
		for (int k = (nSpiraleChirale % 2 == 0 ? 0 : nSpirale); k < (nSpiraleChirale == CHIRAL_LEFT ? nSpirale : 2 * nSpirale); k++)
		{
			nStart = nMaxSpirale * (k % nSpirale) + (nSpiraleChirale == CHIRAL_BOTH && k >= nSpirale ? nMaxSpirale * nSpirale : 0);

			glBegin(GL_LINE_STRIP);

			glVertex3f(mySpiraleData[nStart].Centre.x, mySpiraleData[nStart].Centre.y, mySpiraleData[nStart].Centre.z);

			for (int i = 1; i < nMaxSpirale; i++)
				glVertex3f(mySpiraleData[nStart + i].Centre.x, mySpiraleData[nStart + i].Centre.y, mySpiraleData[nStart + i].Centre.z);

			glEnd();
		}
	}
}
void DrawFrise()
{
	int nMaxI = myDome.nFriseStep * ((myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides + 1);
	int nFrise = myDome.nFrise;

//	nMaxI = CalculFrise();
	if (nMaxI > 0 && nFrise > 0 && myFriseData != NULL)
	{
		glLineWidth(3);

		for (int k = 0; k < nFrise; k++)
		{

			glBegin(GL_LINE_STRIP);

			glVertex3f(myFriseData[nMaxI * k].Centre.x, myFriseData[nMaxI * k].Centre.y, myFriseData[nMaxI * k].Centre.z);

			for (int i = 1; i < nMaxI; i++)
				glVertex3f(myFriseData[nMaxI * k + i].Centre.x, myFriseData[nMaxI * k + i].Centre.y, myFriseData[nMaxI * k + i].Centre.z);

			glEnd();
		}
	}

}
void DrawSnake()
{
	int nSnake = myDome.nSnake;
	int nSnakeChirale = myDome.nSnakeChirale;
	int nMaxI = myDome.nSnakeStep * ((myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides + 1);


	if (nSnake > 0)
	{
//		nMaxI = CalculSnake();
		glLineWidth(3);

		if (nMaxI > 0 && mySnakeData != NULL)
		{
			for (int nbChirale = 0; nbChirale < (nSnakeChirale == CHIRAL_BOTH ? 2 : 1); nbChirale++)
			{
				for (int k = 0; k < nSnake; k++)
				{
					glBegin(GL_LINE_STRIP);
					for (int i = 0; i < nMaxI; i++)
						glVertex3f(mySnakeData[nbChirale * nSnake * nMaxI + k * nMaxI + i].Centre.x, mySnakeData[nbChirale * nSnake * nMaxI + k * nMaxI + i].Centre.y, mySnakeData[nbChirale * nSnake * nMaxI + k * nMaxI + i].Centre.z);
					glEnd();
				}
			}
		}
	}
}
void DrawComb()
{
	int nComb = myDome.nComb;
	int i = 0;

	if (nComb > 0)
	{
//		nMaxI = CalculComb();
		glLineWidth(3);

		if ( myCombData != NULL)
		{
			while (myCombData[i].Radius > -2)
			{
				glBegin(GL_LINE_STRIP);
				while (myCombData[i].Radius > 0)
				{
					glVertex3f(myCombData[i].Centre.x, myCombData[i].Centre.y, myCombData[i].Centre.z);
					i++;
				}
				glEnd();
				i++;
			}
		}
	}
}

int PrintSpirale(FILE* Dome, int nbObjects, PRINT_TYPE Type)
{
	int   nSpirale = myDome.nSpirale;
	float fSpiraleTour = myDome.fSpiraleTour;
	int   nSpiraleChirale = myDome.nSpiraleChirale;
	char  sSpiraleMaterial[400];
	strcpy_s(sSpiraleMaterial, myDome.sSpiraleMaterial);
	GLpoint c0, c1;
	int nMaxI, nMaxSpirale;
	int nStart;

	if (Type == PRINT_STL_BINARY)
		nbObjects = 0;
	nMaxI = CalculSpirale();
	if (nMaxI > 0)
	{
		nMaxSpirale = (int)((nMaxI - 1) * fSpiraleTour + 1);

		for (int k = (nSpiraleChirale % 2 == 0 ? 0 : nSpirale); k < (nSpiraleChirale == CHIRAL_LEFT ? nSpirale : 2 * nSpirale); k++)
		{
			nStart = nMaxSpirale * (k % nSpirale) + (nSpiraleChirale == CHIRAL_BOTH && k >= nSpirale ? nMaxSpirale * nSpirale : 0);

			c0.x = 2 * mySpiraleData[nStart].Centre.x;
			c0.y = 2 * mySpiraleData[nStart].Centre.y;
			c0.z = 2 * mySpiraleData[nStart].Centre.z;

			if (Type == PRINT_STL_BINARY)
				nbObjects += PrintSTLSphereBinary(c0, mySpiraleData[nStart].Radius, GetColor(sSpiraleMaterial), Dome);
			else if (Type == PRINT_OBJ)
				nbObjects += PrintOBJSphere(c0, mySpiraleData[nStart].Radius, Dome, nbObjects, "Spirale Sphere", sSpiraleMaterial);
			else if (Type == PRINT_POV)
				PrintPOVSphere(c0, mySpiraleData[nStart].Radius, Dome, sSpiraleMaterial);

			for (int i = 1; i < nMaxSpirale; i++)
			{
				c1.x = 2 * mySpiraleData[nStart + i].Centre.x;
				c1.y = 2 * mySpiraleData[nStart + i].Centre.y;
				c1.z = 2 * mySpiraleData[nStart + i].Centre.z;

				if (Type == PRINT_STL_BINARY)
				{
					nbObjects += PrintSTLSphereBinary(c1, mySpiraleData[nStart + i].Radius, GetColor(sSpiraleMaterial), Dome);
					nbObjects += PrintSTLConeBinary(c0, c1, mySpiraleData[nStart + i - 1].Radius, mySpiraleData[nStart + i].Radius, GetColor(sSpiraleMaterial), Dome);
				}
				else if (Type == PRINT_OBJ)
				{
					nbObjects += PrintOBJSphere(c1, mySpiraleData[nStart + i].Radius, Dome, nbObjects, "Spirale Sphere", sSpiraleMaterial);
					nbObjects += PrintOBJCone(c0, c1, mySpiraleData[nStart + i - 1].Radius, mySpiraleData[nStart + i].Radius, Dome, nbObjects, "Spirale Cone", sSpiraleMaterial);
				}
				else if (Type == PRINT_POV)
					PrintPOVRoundedCone(c0, c1, mySpiraleData[nStart + i - 1].Radius, mySpiraleData[nStart + i].Radius, Dome, sSpiraleMaterial);

				c0.x = c1.x;
				c0.y = c1.y;
				c0.z = c1.z;
			}
		}
	}
	return nbObjects;
}
int PrintFrise(FILE* Dome, int nbObjects, PRINT_TYPE Type)
{
	int   nFrise = myDome.nFrise;
	char  sFriseMaterial[400];

	strcpy_s(sFriseMaterial, myDome.sFriseMaterial);
	GLpoint c0, c1;
	int nMaxI;

	if (Type == PRINT_STL_BINARY)
		nbObjects = 0;

	nMaxI = CalculFrise();
	if (nMaxI > 0)
	{
		for (int k = 0; k < nFrise; k++)
		{

			c0.x = 2 * myFriseData[k * nMaxI].Centre.x;
			c0.y = 2 * myFriseData[k * nMaxI].Centre.y;
			c0.z = 2 * myFriseData[k * nMaxI].Centre.z;

			if (Type == PRINT_STL_BINARY)
				nbObjects += PrintSTLSphereBinary(c0, myFriseData[k * nMaxI].Radius, GetColor(sFriseMaterial), Dome);
			else if (Type == PRINT_OBJ)
				nbObjects += PrintOBJSphere(c0, myFriseData[k * nMaxI].Radius, Dome, nbObjects, "Frise Sphere", sFriseMaterial);
			else if (Type == PRINT_POV)
				PrintPOVSphere(c0, myFriseData[k * nMaxI].Radius, Dome, sFriseMaterial);


			for (int i = 1; i < nMaxI; i++)
			{

				c1.x = 2 * myFriseData[k * nMaxI + i].Centre.x;
				c1.y = 2 * myFriseData[k * nMaxI + i].Centre.y;
				c1.z = 2 * myFriseData[k * nMaxI + i].Centre.z;

				if (Type == PRINT_STL_BINARY)
				{
					nbObjects += PrintSTLSphereBinary(c1, myFriseData[k * nMaxI + i].Radius, GetColor(sFriseMaterial), Dome);
					nbObjects += PrintSTLConeBinary(c0, c1, myFriseData[k * nMaxI + i - 1].Radius, myFriseData[k * nMaxI + i].Radius, GetColor(sFriseMaterial), Dome);
				}
				else if (Type == PRINT_OBJ)
				{
					nbObjects += PrintOBJSphere(c1, myFriseData[k * nMaxI + i].Radius, Dome, nbObjects, "Frise Sphere", sFriseMaterial);
					nbObjects += PrintOBJCone(c0, c1, myFriseData[k * nMaxI + i - 1].Radius, myFriseData[k * nMaxI + i].Radius, Dome, nbObjects, "Frise Cone", sFriseMaterial);
				}
				else if (Type == PRINT_POV)
				{
					PrintPOVSphere(c1, myFriseData[k * nMaxI + i].Radius, Dome, sFriseMaterial);
					PrintPOVCone(c0, c1, myFriseData[k * nMaxI + i - 1].Radius, myFriseData[k * nMaxI + i].Radius, Dome, sFriseMaterial);
				}


				c0.x = c1.x;
				c0.y = c1.y;
				c0.z = c1.z;

			}
		}
	}
	return nbObjects;
}
int PrintSnake(FILE* Dome, int nbObjects, PRINT_TYPE Type)
{
	int nSnake = myDome.nSnake;
	int nSnakeChirale = myDome.nSnakeChirale;
	int nMaxI;
	char sSnakeMaterial[400];

	GLpoint c0, c1;
	int nStart;

	if (Type == PRINT_STL_BINARY)
		nbObjects = 0;
	strcpy_s(sSnakeMaterial, myDome.sSnakeMaterial);

//	if (nSnake > myDome.nBaseSides)
//		nSnake = myDome.nBaseSides;


	if (nSnake > 0)
	{
		nMaxI = CalculSnake();
		if (nMaxI > 0)
		{
			for (int nbChirale = 0; nbChirale < (nSnakeChirale == CHIRAL_BOTH ? 2 : 1); nbChirale++)
			{
				for (int k = 0; k < nSnake; k++)
				{
					nStart = nbChirale * nSnake * nMaxI + k * nMaxI;

					c0.x = 2 * mySnakeData[nStart].Centre.x;
					c0.y = 2 * mySnakeData[nStart].Centre.y;
					c0.z = 2 * mySnakeData[nStart].Centre.z;

					if (Type == PRINT_STL_BINARY)
						nbObjects += PrintSTLSphereBinary(c0, mySnakeData[nStart].Radius, GetColor(sSnakeMaterial), Dome);
					else if (Type == PRINT_OBJ)
						nbObjects += PrintOBJSphere(c0, mySnakeData[nStart].Radius, Dome, nbObjects, "Snake Sphere", sSnakeMaterial);
					else if (Type == PRINT_POV)
						PrintPOVSphere(c0, mySnakeData[nStart].Radius, Dome, sSnakeMaterial);
					for (int i = 1; i < nMaxI; i++)
					{

						c1.x = 2 * mySnakeData[nStart + i].Centre.x;
						c1.y = 2 * mySnakeData[nStart + i].Centre.y;
						c1.z = 2 * mySnakeData[nStart + i].Centre.z;
						if (Type == PRINT_STL_BINARY)
						{
							nbObjects += PrintSTLSphereBinary(c1, mySnakeData[nStart + i].Radius, GetColor(sSnakeMaterial), Dome);
							nbObjects += PrintSTLConeBinary(c0, c1, mySnakeData[nStart + i - 1].Radius, mySnakeData[nStart + i].Radius, GetColor(sSnakeMaterial), Dome);
						}
						else if (Type == PRINT_OBJ)
						{
							nbObjects += PrintOBJSphere(c1, mySnakeData[nStart + i].Radius, Dome, nbObjects, "Snake Sphere", sSnakeMaterial);
							nbObjects += PrintOBJCone(c0, c1, mySnakeData[nStart + i - 1].Radius, mySnakeData[nStart + i].Radius, Dome, nbObjects, "Spirale Cone", sSnakeMaterial);
						}
						else if (Type == PRINT_POV)
							PrintPOVRoundedCone(c0, c1, mySnakeData[nStart + i - 1].Radius, mySnakeData[nStart + i].Radius, Dome, sSnakeMaterial);

						c0.x = c1.x;
						c0.y = c1.y;
						c0.z = c1.z;
					}
				}
			}
		}
		return nbObjects;
	}
	return 0;
}
int PrintComb(FILE* Dome, int nbObjects, PRINT_TYPE Type)
{
	int nComb = myDome.nComb;
	int nMaxI;
	char sCombMaterial[400];

	GLpoint c0, c1;

	if (Type == PRINT_STL_BINARY)
		nbObjects = 0;
	strcpy_s(sCombMaterial, myDome.sCombMaterial);

	int i = 0;

	if (nComb > 0)
	{
		nMaxI = CalculComb();
		if (nMaxI > 0)
		{


			while (myCombData[i].Radius > -2)
			{
				c0.x = 2 * myCombData[i].Centre.x;
				c0.y = 2 * myCombData[i].Centre.y;
				c0.z = 2 * myCombData[i].Centre.z;

				if (Type == PRINT_STL_BINARY)
					nbObjects += PrintSTLSphereBinary(c0, myCombData[i].Radius, GetColor(sCombMaterial), Dome);
				else if (Type == PRINT_OBJ)
					nbObjects += PrintOBJSphere(c0, myCombData[i].Radius, Dome, nbObjects, "Comb Sphere", sCombMaterial);
				else if (Type == PRINT_POV)
					PrintPOVSphere(c0, myCombData[i].Radius, Dome, sCombMaterial);

				while (myCombData[i].Radius > 0)
				{
					c1.x = 2 * myCombData[i].Centre.x;
					c1.y = 2 * myCombData[i].Centre.y;
					c1.z = 2 * myCombData[i].Centre.z;
					if (Type == PRINT_STL_BINARY)
					{
						nbObjects += PrintSTLSphereBinary(c1, myCombData[i].Radius, GetColor(sCombMaterial), Dome);
						nbObjects += PrintSTLConeBinary(c0, c1, myCombData[i - 1].Radius, myCombData[i].Radius, GetColor(sCombMaterial), Dome);
					}
					else if (Type == PRINT_OBJ)
					{
						nbObjects += PrintOBJSphere(c1, myCombData[i].Radius, Dome, nbObjects, "Comb Sphere", sCombMaterial);
						nbObjects += PrintOBJCone(c0, c1, myCombData[i - 1].Radius, myCombData[i].Radius, Dome, nbObjects, "Comb Cone", sCombMaterial);
					}
					else if (Type == PRINT_POV)
						PrintPOVRoundedCone(c0, c1, myCombData[i - 1].Radius, myCombData[i].Radius, Dome, sCombMaterial);

					c0.x = c1.x;
					c0.y = c1.y;
					c0.z = c1.z;
					i++;
				}
				i++;
			}

			
		}
		return nbObjects;
	}
	return 0;
}

