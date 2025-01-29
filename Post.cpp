#include <string.h>
#include <stdio.h>
#include "stdafx.h"
#include <time.h>
//#include "globales.h"
#include "functions.h"
#include "Types.h"

extern DomeParameters myDome;
extern DomeDisplay myDisplay;
extern GLUI_StaticText* stTest;
extern Filtre FiltreListe[];
extern GLUI_Listbox* lbFiltreType;
extern GLUI_Spinner* spLissageSigma1;
extern GLUI_Spinner* spLissageSigma2;
extern GLUI_Listbox* lbLissageDirection;
extern GLUI_Spinner* spPoidsSelf;
extern GLUI_Spinner* spPoidsVertical;
extern GLUI_Spinner* spPoidsHorizontal;
extern GLUI_Spinner* spPoidsDiagonales;
extern TestSlope TestSlopeFunctionListe[];

float MarqueeFunction(float x, float y)
{
	float a = (float)myDome.nMeridianMinor / 2;
	float b = (float)myDome.nParallelMinor / 2;
	float alpha  = myDome.fMarqueeAlpha;
	float beta  = myDome.fMarqueeBeta;
	Parametric_2 f1;
	float Param1 = myDome.fMarqueeParameter;
	int nFunction = myDome.nMarqueeFunction;


	f1 = TestSlopeFunctionListe[nFunction].Function;
//	printf("%s\n", TestSlopeFunctionListe[nFunction].sNom);

	double z =  ((1.0f - alpha - beta) * f1((x / a)* (x / a), Param1) * f1((y / b)* (y / b), Param1) + alpha * f1((x / a)*(x / a), Param1) + beta * f1((y / b)* (y / b), Param1));
//	z = ((1.0f - alpha - beta) * x * x * y * y / (a * a * b * b) + alpha * x * x / (a * a) + beta * y * y / (b * b));

	return (float)z;
}

float MarqueeFunction(int i, int j)
{
	float a = (float)(myDome.nMeridianMinor + 1) / 2;
	float b = (float)(myDome.nParallelMinor + 1) / 2;
	float alpha = myDome.fMarqueeAlpha;
	float beta = myDome.fMarqueeBeta;
	Parametric_2 f1;
	float Param1 = myDome.fMarqueeParameter;
	int nFunction = myDome.nMarqueeFunction;
	float x, y;

	f1 = TestSlopeFunctionListe[nFunction].Function;
	x = ((float)i - a) / a;
	y = ((float)j - b) / b;

	double z = ((1.0f - alpha - beta) * f1(x * x, Param1) * f1(y * y, Param1) + alpha * f1(x * x, Param1) + beta * f1(y * y, Param1));

	return (float)z;
}
void FilterInput(TYPE_FILTER nFiltre, bool bDefault = true)
{
	switch (nFiltre)
	{
	case FILTER_DoG:
		spLissageSigma1->enable();
		spLissageSigma2->enable();
		lbLissageDirection->disable();
		spPoidsSelf->disable();
		spPoidsVertical->disable();
		spPoidsHorizontal->disable();
		spPoidsDiagonales->disable();
		break;
	case FILTER_MDIF:
		spLissageSigma1->enable();
		lbLissageDirection->enable();

		spLissageSigma2->disable();
		spPoidsSelf->disable();
		spPoidsVertical->disable();
		spPoidsHorizontal->disable();
		spPoidsDiagonales->disable();
		break;
	case FILTER_Gauss:
	case FILTER_Circular:
	case FILTER_LoG:
	case FILTER_Exponential:
	case FILTER_HighPass:
	case FILTER_LowPass:
		spLissageSigma1->enable();
		spLissageSigma2->disable();
		lbLissageDirection->disable();
		spPoidsSelf->disable();
		spPoidsVertical->disable();
		spPoidsHorizontal->disable();
		spPoidsDiagonales->disable();
		break;
	case FILTER_Prewitt:
	case FILTER_Sobel:
		lbLissageDirection->enable();
		spLissageSigma1->disable();
		spLissageSigma2->disable();
		spPoidsSelf->disable();
		spPoidsVertical->disable();
		spPoidsHorizontal->disable();
		spPoidsDiagonales->disable();
		break;
	case FILTER_Mean:
	case FILTER_Binomial:
	case FILTER_Pyramidal:
	case FILTER_Laplacian:
	case FILTER_Conic:
	case FILTER_Enhancer:
		spLissageSigma1->disable();
		spLissageSigma2->disable();
		lbLissageDirection->disable();
		spPoidsSelf->disable();
		spPoidsVertical->disable();
		spPoidsHorizontal->disable();
		spPoidsDiagonales->disable();
		break;
	case FILTER_Perso:
		spPoidsSelf->enable();
		spPoidsVertical->enable();
		spPoidsHorizontal->enable();
		spPoidsDiagonales->enable();
		spLissageSigma1->disable();
		spLissageSigma2->disable();
		lbLissageDirection->disable();

	}
	if (bDefault)
	{
		myDome.fLissagePuissance = FiltreListe[nFiltre].fPower;
		myDome.fLissageSigma1 = FiltreListe[nFiltre].fSigma1;
		myDome.fLissageSigma2 = FiltreListe[nFiltre].fSigma2;
		myDome.nLissageDirection = FiltreListe[nFiltre].nDirection;
		myDome.nLissageDistanceType = FiltreListe[nFiltre].nDistance;
		myDome.fLissagePoidsSelf = FiltreListe[nFiltre].fSelf;
		myDome.fLissagePoidsHorizontal = FiltreListe[nFiltre].fHorizontal;
		myDome.fLissagePoidsVertical = FiltreListe[nFiltre].fVertical;
		myDome.fLissagePoidsDiagonales = FiltreListe[nFiltre].fDiagonal;
	}
}

unsigned long CoefBinomial(int n, int k)
{
	if (k == 0 || k == n)
		return 1;
	return CoefBinomial(n - 1, k - 1) + CoefBinomial(n - 1, k);
}
void MatricesMultiplication(float* A, float* B, float* C, int nSize)
{
	for (int i = 0; i < nSize; i++)
		for (int j = 0; j < nSize; j++)
			C[i * nSize + j] = 0;

	for (int i = 0; i < nSize; i++)
		for (int j = 0; j < nSize; j++)
			for (int k = 0; k < nSize; k++)
				C[i * nSize + j] += A[i * nSize + k] * B[k * nSize + j];


}

float Distance(int i, int j, DISTANCE_Type nDistance)
{
	float fDistance;
	float x = (float)i;
	float y = (float)j;

	switch (nDistance)
	{
	case DISTANCE_EUCLIDE:
	default:
		fDistance = sqrtf(x * x + y * y);
		break;
	case DISTANCE_MANHATTAN:
		fDistance = fabsf(x) + fabsf(y);
		break;
	case DISTANCE_ULTRA:
		fDistance = __max(fabsf(x), fabsf(y));
		break;
	case DISTANCE_3:
		fDistance = powf(fabsf(x * x * x) + fabsf(y * y * y), 1.0f / 3.0f);
		break;
	case DISTANCE_DISCRETE:
		fDistance = (x == 0 ? (y == 0.0f ? 0 : 1.0f) : 1);
	}
	return fDistance;
}
float Max(float* Matrix, int nCote)
{
	float fMax = -1000000;
	int nSize = 2 * nCote + 1;

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			if (fMax < Matrix[(i + nCote) * nSize + (j + nCote)])
				fMax = Matrix[(i + nCote) * nSize + (j + nCote)];
		}
	return fMax;
}

void PrintMatrix(float* Matrix, int nCote, float fTotal, int nType)
{
	float fMin = 1000000;
	int nSize = 2 * nCote + 1;
	float fEcart = 100;
	float fMax = -1;
	float fDiv = -1;

	printf("\n");

	for (int j = -nCote; j <= nCote; j++)
	{
		for (int i = -nCote; i <= nCote; i++)
		{
			if (fabsf(Matrix[(j + nCote) * nSize + (i + nCote)]) < fMin && Matrix[(j + nCote) * nSize + (i + nCote)] != 0)
				fMin = fabsf(Matrix[(j + nCote) * nSize + (i + nCote)]);
			if (fabsf(Matrix[(j + nCote) * nSize + (i + nCote)]) > fMax)
				fMax = fabsf(Matrix[(j + nCote) * nSize + (i + nCote)]);
			if (fabsf(Matrix[(j + nCote) * nSize + (i + nCote)]) != 0 && fMax / fabsf(Matrix[(j + nCote) * nSize + (i + nCote)]) < fEcart && fMax / fabsf(Matrix[(j + nCote) * nSize + (i + nCote)]) > fDiv)
				fDiv = fMax / fabsf(Matrix[(j + nCote) * nSize + (i + nCote)]);

		}
	}
	fMin = (fMin == 0 ? 1 : fMin);
	fTotal = (fTotal == 0 ? 1 : fTotal);

	for (int i = -nCote; i <= nCote; i++)
	{
		for (int j = -nCote; j <= nCote; j++)
		{
			switch (nType)
			{
			case TYPE_Brut:
				printf("%8.5f;", Matrix[(j + nCote) * nSize + (i + nCote)]);
				break;
			case TYPE_Norme:
				printf("%8.5f;", Matrix[(j + nCote) * nSize + (i + nCote)] / fTotal);
				break;
			case TYPE_Integer:
				printf("%8d;", (int)round(Matrix[(j + nCote) * nSize + (i + nCote)] / fMin));
				break;
			case TYPE_Round:
				printf("%8d;", (int)round(Matrix[(j + nCote) * nSize + (i + nCote)]));
				break;
			case TYPE_Standard:
				printf("%8d;  ", (int)round(Matrix[(j + nCote) * nSize + (i + nCote)] * fDiv / fMax));
				break;
			default:
				break;
			}
		}
		printf("\n");
	}
	printf("\n");
}

float Perso(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	int nSize = 2 * nCote + 1;
	float fDistance;
	// Ajouter un coef dist dans le calcul de fCoef, et le choix de la distance
	float fTotal = 0;
	float fSelf = Parameters->fSelf;
	float fDiag = Parameters->fDiagonal;
	float fVertic = Parameters->fVertical;
	float fHoriz = Parameters->fHorizontal;

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			//			nDistance = __max(abs(i), abs(j));
			//			fDistance = sqrtf((float)(i * i + j * j));
			fDistance = Distance(i, j, nDistance);
			float fCoef = 1.0f / powf(fDistance, fPower);

			if (i == 0 && j == 0)
				Matrix[(i + nCote) * nSize + (j + nCote)] = fSelf;
			else if (abs(i) == abs(j))
				Matrix[(i + nCote) * nSize + (j + nCote)] = fDiag * abs(i) * fCoef;
			else if (abs(i) > abs(j))
				Matrix[(i + nCote) * nSize + (j + nCote)] = (((abs(j) * fDiag) + ((abs(i) - abs(j)) * fHoriz)) * fCoef);
			else
				Matrix[(i + nCote) * nSize + (j + nCote)] = (((abs(i) * fDiag) + ((abs(j) - abs(i)) * fVertic)) * fCoef);

			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];
		}

	return fTotal;

}
float Gauss(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	float fSigma = Parameters->fSigma1;
	float fCoef = 1.0f / (2.0f * (float)_PI * fSigma * fSigma);
	int nSize = 2 * nCote + 1;
	float fTotal = 0;
	float fDistance;

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			fDistance = powf(Distance(i, j, nDistance), fPower);
			Matrix[(i + nCote) * nSize + (j + nCote)] = fCoef*expf(-(fDistance) / (2 * fSigma * fSigma));
			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];
		}

	return fTotal;
}
float Mean(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	int nSize = 2 * nCote + 1;
	float fDis;

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			fDis = 1 - powf(Distance(i, j, nDistance) / Distance(nCote + 1, nCote + 1, nDistance), fPower);
			Matrix[(i + nCote) * nSize + (j + nCote)] = fDis;
		}
	return (float)(nSize * nSize);
}
float Binomial(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	int nSize = 2 * nCote + 1;
	float fTotal = 0;

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			float fDis = powf(Distance(i, j, nDistance), fPower);
			fDis = (fDis == 0 ? 1 : fDis);
			if (abs(i) == nCote)
				Matrix[(i + nCote) * nSize + (j + nCote)] = (float)CoefBinomial(2 * nCote, nCote - abs(j)) / fDis;
			else if (abs(j) == nCote)
				Matrix[(i + nCote) * nSize + (j + nCote)] = (float)CoefBinomial(2 * nCote, nCote - abs(i)) / fDis;
			else
				Matrix[(i + nCote) * nSize + (j + nCote)] = (float)CoefBinomial(2 * nCote, nCote - abs(i)) * (float)CoefBinomial(2 * nCote, nCote - abs(j)) / fDis;

			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];

		}
	return fTotal;
}
float Pyramidal(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	int nSize = 2 * nCote + 1;
	float fTotal = 0;

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			float fDis = powf(Distance(i, j, nDistance), fPower);
			fDis = (fDis == 0 ? 1 : fDis);

			if (abs(i) == nCote)
				Matrix[(i + nCote) * nSize + (j + nCote)] = (float)(nCote - abs(j) + 1) / fDis;
			else if (abs(j) == nCote)
				Matrix[(i + nCote) * nSize + (j + nCote)] = (float)(nCote - abs(i) + 1) / fDis;
			else
				Matrix[(i + nCote) * nSize + (j + nCote)] = (float)((nCote - abs(i) + 1) * (nCote - abs(j) + 1)) / fDis;
			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];

		}
	return fTotal;
}
float Circular(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	int nSize = 2 * nCote + 1;
	float fTotal = 0;
	float fDistance;
	float fRayon = Parameters->fSigma1;

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			fDistance = powf(Distance(i, j, nDistance), fPower);
			if (fDistance < fRayon)
				Matrix[(i + nCote) * nSize + (j + nCote)] = 1.0f;
			else
				Matrix[(i + nCote) * nSize + (j + nCote)] = 0.0f;
			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];

		}
	return fTotal;
}
float LoG(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	float fSigma = Parameters->fSigma1;
	float fCoef = 1.0f / ((float)_PI * fSigma * fSigma * fSigma * fSigma);
	int nSize = 2 * nCote + 1;
	float fTotal = 0;
	float fDis, fDistance;
	float fNorme;

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			fDistance = powf(Distance(i, j, nDistance), fPower);
			fDis = (fDistance) / (2 * fSigma * fSigma);

			Matrix[(i + nCote) * nSize + (j + nCote)] = -fCoef * (1 - fDis) * expf(-fDis);
			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];
		}
	fTotal /= Matrix[(0 + nCote) * nSize + (0 + nCote)];
	fNorme = Matrix[(0 + nCote) * nSize + (0 + nCote)];

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
			Matrix[(i + nCote) * nSize + (j + nCote)] /= fNorme;
	return fTotal;
}
float Exponential(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	float fGamma = Parameters->fSigma1;
	float fCoef = fGamma * fGamma;
	int nSize = 2 * nCote + 1;
	float fTotal = 0;
	float fDis;

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			fDis = powf(Distance(i, j, nDistance), fPower) * fGamma;
			//			printf("%d ", (i + nCote) * nSize + (j + nCote));
			Matrix[(i + nCote) * nSize + (j + nCote)] = fCoef * expf(-fDis);
			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];
		}

	return fTotal;

}
float Sobel(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	int nSize = 2 * nCote + 1;
	float fTotal = 0;
	float fDis;
	float fMax = -1000000;
	int nSens = Parameters->nDirection;
	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			fDis = Distance(i, j, nDistance);

			if ((nSens == DIRECTION_NORTH ? i : j) == 0)
				Matrix[(i + nCote) * nSize + (j + nCote)] = 0;
			else
				//	Matrix[(i + nCote) * nSize + (j + nCote)] = (float)(nSens == DIRECTION_NORTH ? i : j) / (fDis * abs(i));
				Matrix[(i + nCote) * nSize + (j + nCote)] = (float)(nSens == DIRECTION_NORTH ? i : j) / (powf((fDis), fPower) * (nSens == DIRECTION_NORTH ? abs(i) : abs(j)));
			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];
			if (fMax < abs(Matrix[(i + nCote) * nSize + (j + nCote)]))
				fMax = abs(Matrix[(i + nCote) * nSize + (j + nCote)]);
		}
	if (fMax > 0)
		for (int i = -nCote; i <= nCote; i++)
			for (int j = -nCote; j <= nCote; j++)
				Matrix[(i + nCote) * nSize + (j + nCote)] /= fMax;
	Matrix[(0 + nCote) * nSize + (0 + nCote)] = 1;
	return fMax;

}
float Laplacian(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	int nSize = 2 * nCote + 1;
	float fTotal = 0;
	float* A, * B;

	A = (float*)malloc(nSize * nSize * sizeof(float));
	B = (float*)malloc(nSize * nSize * sizeof(float));

	Parameters->nDirection = DIRECTION_NORTH;
	Sobel(B, nCote, nDistance, fPower, Parameters);
	Parameters->nDirection = DIRECTION_WEST;
	Sobel(A, nCote, nDistance, fPower, Parameters);

	MatricesMultiplication(A, B, Matrix, nSize);
	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];

	free(A);
	free(B);
	return fTotal;
}
float MDIF(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	int nSize = 2 * nCote + 1;
	float fTotal = 0;
	float fSigma = Parameters->fSigma1;
	float fCoe = fSigma * fSigma;
	float fMax = -1000000;
	int nDirection = Parameters->nDirection;
	int k;

	Gauss(Matrix, nCote, nDistance, fPower, Parameters);

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			k = (nDirection == SENS_Vertical ? i : j);
			Matrix[(i + nCote) * nSize + (j + nCote)] *= -k / fCoe;
			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];
			if (fMax < abs(Matrix[(i + nCote) * nSize + (j + nCote)]))
				fMax = abs(Matrix[(i + nCote) * nSize + (j + nCote)]);
		}
	Matrix[(0 + nCote) * nSize + (0 + nCote)] = 1;
	return fMax;
}
float Conic(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters) 
{
	int nSize = 2 * nCote + 1;
	float fTotal = 0;
	float fDis;

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			fDis = powf(Distance(i, j, nDistance), fPower);

			Matrix[(i + nCote) * nSize + (j + nCote)] = roundf(__max(nCote + 1 - fDis, 0) * __max(nCote + 1 - fDis, 0) / 2);

			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];
		}

	return fTotal;

}
float HighPass(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	int nSize = 2 * nCote + 1;
	float fTotal = 0;
	float fDis;
	float fGamma = Parameters->fSigma1;

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			//			fDis = (float)(abs(i) + abs(j));
			fDis = Distance(i, j, nDistance);
			if (i == 0 && j == 0)
				Matrix[(i + nCote) * nSize + (j + nCote)] = 0;
			else
				Matrix[(i + nCote) * nSize + (j + nCote)] = -fGamma / (powf(fDis, fPower));

			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];
		}
	Matrix[(nCote)*nSize + (nCote)] = 1 - fTotal;
	return fTotal;
}
float Prewitt(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	int nSize = 2 * nCote + 1;
	float fTotal = 0;
	float fVal;
	float fMax = -1000000;
	int nDirection = Parameters->nDirection;

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			float fDis = powf(Distance(i, j, nDistance), fPower);
			fDis = (fDis == 0 ? 1 : fDis);

			switch (nDirection)
			{
			case DIRECTION_WEST:
				if (i < 0)
					fVal = -1;
				else if (i > 0)
					fVal = 1;
				else
					fVal = 0;
				break;
			case DIRECTION_EAST:
				if (i < 0)
					fVal = 1;
				else if (i > 0)
					fVal = -1;
				else
					fVal = 0;
				break;
			case DIRECTION_NORTH:
				if (j < 0)
					fVal = -1;
				else if (j > 0)
					fVal = 1;
				else
					fVal = 0;
				break;
			case DIRECTION_SOUTH:
				if (j < 0)
					fVal = 1;
				else if (j > 0)
					fVal = -1;
				else
					fVal = 0;
				break;
			case DIRECTION_NORTH_WEST:
				if (i + j < 0)
					fVal = -1;
				else if (i + j > 0)
					fVal = 1;
				else
					fVal = 0;
				break;
			case DIRECTION_SOUTH_EAST:
				if (i + j > 0)
					fVal = -1;
				else if (i + j < 0)
					fVal = 1;
				else
					fVal = 0;
				break;
			case DIRECTION_SOUTH_WEST:
				if (i - j < 0)
					fVal = -1;
				else if (i - j > 0)
					fVal = 1;
				else
					fVal = 0;
				break;
			case DIRECTION_NORTH_EAST:
				if (i - j > 0)
					fVal = -1;
				else if (i - j < 0)
					fVal = 1;
				else
					fVal = 0;
				break;

			}
			fVal /= fDis;
			Matrix[(i + nCote) * nSize + (j + nCote)] = fVal;
//			Matrix[(0 + nCote) * nSize + (0 + nCote)] = 1;
			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];
			if (fMax < abs(Matrix[(i + nCote) * nSize + (j + nCote)]))
				fMax = abs(Matrix[(i + nCote) * nSize + (j + nCote)]);

		}

	return fMax;

}
float LowPass(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	int nSize = 2 * nCote + 1;
	float fTotal = 0;
	float fDis;
	float fAlpha = Parameters->fSigma1;

	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			//			fDis = (float)(abs(i) + abs(j));
			fDis = Distance(i, j, nDistance);
			Matrix[(i + nCote) * nSize + (j + nCote)] = fAlpha * fAlpha / (1 + powf(fDis, fPower));
			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];
		}
	return fTotal;
}
float DoG(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	int nSize = 2 * nCote + 1;
	float fTotal = 0;
	float* A, * B;
	float fTotal1, fTotal2;
	float fSigma1 = Parameters->fSigma1;
	float fSigma2 = Parameters->fSigma2;
	float fMax = -1000000;

	A = (float*)malloc(nSize * nSize * sizeof(float));
	B = (float*)malloc(nSize * nSize * sizeof(float));

	fTotal1 = Gauss(A, nCote, nDistance, fPower, Parameters);
	fTotal1 = Max(A, nCote);
	Parameters->fSigma1 = Parameters->fSigma2;
	fTotal2 = Gauss(B, nCote, nDistance, fPower, Parameters);
	fTotal2 = Max(B, nCote);
	for (int i = -nCote; i <= nCote; i++)
		for (int j = -nCote; j <= nCote; j++)
		{
			Matrix[(i + nCote) * nSize + (j + nCote)] = A[(i + nCote) * nSize + (j + nCote)] / fTotal1 - B[(i + nCote) * nSize + (j + nCote)] / fTotal2;
			fTotal += Matrix[(i + nCote) * nSize + (j + nCote)];
			if (fMax < abs(Matrix[(i + nCote) * nSize + (j + nCote)]))
				fMax = abs(Matrix[(i + nCote) * nSize + (j + nCote)]);
		}

	free(A);
	free(B);
	if (fTotal < 0.000001)
		return fMax;
	return fTotal;
}
float Enhancer(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters)
{
	float fTotal;
	int nSize = 2 * nCote + 1;

	fTotal = Laplacian(Matrix, nCote, nDistance, fPower, Parameters);
	Matrix[(0 + nCote) * nSize + (0 + nCote)]++;

	return fTotal + 1;
}


void Normalize(Vector3D* V)
{
	double fNorme = sqrt(V->x * V->x + V->y * V->y + V->z * V->z);

	if (fNorme != 0.0)
	{
		V->x /= fNorme;
		V->y /= fNorme;
		V->z /= fNorme;
	}
}

float Poids(int a, int b, float p0, float p1, float p2, float p3, int puissance)
{
	return (__max(a, b) == 0 ? p0 : (__min(a, b) * p2 + (__max(a, b) - __min(a, b)) * p1) / (__max(a, b) * __max(a, b)) + p3 / powf((float)(a * a + b * b), (float)puissance));
}
void Coord(int i, int j, int a, int b, int nMaxI, int nMaxJ, int* n_I, int* n_J)
{
	*n_I = (myDome.nDomeEnd == 360 ? (i + a + nMaxI - 1) % (nMaxI - 1) : (i + a >= nMaxI ? -1 : (i + a < 0 ? -1 : i + a)));
	*n_J = (j + b >= nMaxJ ? -1 : (j + b < 0 ? -1 : j + b));
}

void EvoluteOldOld(GLpoint* Vertices, int nMaxI, int nMaxJ)
{
	bool bEvolute = (myDome.bDiffGeomMeridien || myDome.bDiffGeomParallel);
	float Lambda = myDome.fDiffGeomLambda;
	GLpoint* Evolute;
	GLpoint Temp;
	int n_I1, n_I2;

	if (myDome.nDomeEnd != 360 && false)
		nMaxI = (int)((float)(nMaxI * myDome.nDomeEnd) / 360.0) + 1;

	if (bEvolute && myDisplay.bDisplayMinor)
	{
		Evolute = (GLpoint*)malloc(nMaxI * nMaxJ * sizeof(GLpoint));
		if (Evolute == 0)
		{
			Message("     Memory Allocation Evolute Problem   ");
			return;
		}

		memcpy_s(Evolute, nMaxI * nMaxJ * sizeof(GLpoint), Vertices, nMaxI * nMaxJ * sizeof(GLpoint));
		for (int i = 0; i < nMaxI; i++)
			for (int j = 1; j < nMaxJ - 1; j++)
			{
				Temp.x = 0.0;
				Temp.y = 0.0;
				Temp.z = 0.0;

				if (myDome.bDiffGeomMeridien)
				{
					Temp.x = Evolute[i * nMaxJ + j - 1].x + Evolute[i * nMaxJ + j + 1].x - 2 * Evolute[i * nMaxJ + j].x;
					Temp.y = Evolute[i * nMaxJ + j - 1].y + Evolute[i * nMaxJ + j + 1].y - 2 * Evolute[i * nMaxJ + j].y;
					Temp.z = Evolute[i * nMaxJ + j - 1].z + Evolute[i * nMaxJ + j + 1].z - 2 * Evolute[i * nMaxJ + j].z;
				}

				if (myDome.bDiffGeomParallel)
				{
					n_I1 = (myDome.nDomeEnd == 360 ? (i - 2 + nMaxI) % (nMaxI - 1) : (i - 1 >= nMaxI ? -1 : (i - 1 < 0 ? -1 : i - 1)));
					n_I2 = (myDome.nDomeEnd == 360 ? (i + 1) % (nMaxI - 1) : (i + 1 >= nMaxI ? -1 : (i + 1 < 0 ? -1 : i + 1)));

					if (n_I1 != -1 && n_I2 != -1)
					{
						Temp.x += Evolute[n_I1 * nMaxJ + j].x + Evolute[n_I2 * nMaxJ + j].x - 2 * Evolute[i * nMaxJ + j].x;
						Temp.y += Evolute[n_I1 * nMaxJ + j].y + Evolute[n_I2 * nMaxJ + j].y - 2 * Evolute[i * nMaxJ + j].y;
						Temp.z += Evolute[n_I1 * nMaxJ + j].z + Evolute[n_I2 * nMaxJ + j].z - 2 * Evolute[i * nMaxJ + j].z;
					}
				}


				Vertices[i * nMaxJ + j].x = Evolute[i * nMaxJ + j].x + (1 - Lambda) * Temp.x;
				Vertices[i * nMaxJ + j].y = Evolute[i * nMaxJ + j].y + (1 - Lambda) * Temp.y;
				Vertices[i * nMaxJ + j].z = Evolute[i * nMaxJ + j].z + (1 - Lambda) * Temp.z;

			}


		free(Evolute);
	}

}
float Norme(GLpoint* Evolute, int nMaxJ, int i, int j)
{
	return (float)sqrt(Evolute[i * nMaxJ + j].x * Evolute[i * nMaxJ + j].x + Evolute[i * nMaxJ + j].y * Evolute[i * nMaxJ + j].y);
}
void PedalOld(GLpoint* Vertices, int nMaxI, int nMaxJ)
{
	bool bPedal = (myDome.bDiffGeomMeridien || myDome.bDiffGeomParallel);
	float Lambda = myDome.fDiffGeomLambda;
	GLpoint* Pedal;
	GLpoint Temp;


	if (bPedal && myDisplay.bDisplayMinor)
	{
		Pedal = (GLpoint*)malloc(nMaxI * nMaxJ * sizeof(GLpoint));
		if (Pedal == 0)
		{
			Message("     Memory Allocation Pedal Problem   ");
			return;
		}

		memcpy_s(Pedal, nMaxI * nMaxJ * sizeof(GLpoint), Vertices, nMaxI * nMaxJ * sizeof(GLpoint));

		for (int i = 0; i < nMaxI; i++)
			for (int j = 1; j < nMaxJ; j++)
			{
				Temp.x = 0.0;
				Temp.y = 0.0;
				Temp.z = 0.0;

				if (myDome.bDiffGeomMeridien)
				{
					float x = Norme(Pedal, nMaxJ, i, j);
					float y = Pedal[i * nMaxJ + j].z;

					float dx, dy;
					float ddx, ddy;

					if (j == 1 || j == nMaxJ - 2)
					{
						//central  1
						dx = -Norme(Pedal, nMaxJ, i, j - 1) / 2 + Norme(Pedal, nMaxJ, i, j + 1) / 2;
						dy = -Pedal[i * nMaxJ + j - 1].z / 2 + Pedal[i * nMaxJ + j + 1].z / 2;

						ddx = Norme(Pedal, nMaxJ, i, j - 1) - 2 * Norme(Pedal, nMaxJ, i, j) + Norme(Pedal, nMaxJ, i, j + 1);
						ddy = Pedal[i * nMaxJ + j - 1].z - 2 * Pedal[i * nMaxJ + j].z + Pedal[i * nMaxJ + j + 1].z;
					}
					else if (j == nMaxJ - 1)
					{
						//backward 2
						dx = Norme(Pedal, nMaxJ, i, j - 2) / 2 - 2 * Norme(Pedal, nMaxJ, i, j - 1) + 3 * Norme(Pedal, nMaxJ, i, j) / 2;
						dy = Pedal[i * nMaxJ + j - 2].z / 2 - 2 * Pedal[i * nMaxJ + j - 1].z + 3 * Pedal[i * nMaxJ + j].z / 2;

						ddx = -Norme(Pedal, nMaxJ, i, j - 3) + 4 * Norme(Pedal, nMaxJ, i, j - 2) - 5 * Norme(Pedal, nMaxJ, i, j - 1) + 2 * Norme(Pedal, nMaxJ, i, j);
						ddy = -Pedal[i * nMaxJ + j - 3].z + 4 * Pedal[i * nMaxJ + j - 2].z - 5 * Pedal[i * nMaxJ + j - 1].z + 2 * Pedal[i * nMaxJ + j].z;


					}
					else
					{
						//central 2
						dx = Norme(Pedal, nMaxJ, i, j - 2) / 12 - 2 * Norme(Pedal, nMaxJ, i, j - 1) / 3 + 2 * Norme(Pedal, nMaxJ, i, j + 1) / 3 - Norme(Pedal, nMaxJ, i, j + 2) / 12;
						dy = Pedal[i * nMaxJ + j - 2].z / 12 - 2 * Pedal[i * nMaxJ + j - 1].z / 3 + 2 * Pedal[i * nMaxJ + j + 1].z / 3 - Pedal[i * nMaxJ + j + 2].z / 12;

						ddx = -Norme(Pedal, nMaxJ, i, j - 2) / 12 + 4 * Norme(Pedal, nMaxJ, i, j - 1) / 3 - 5 * Norme(Pedal, nMaxJ, i, j) / 2 + 4 * Norme(Pedal, nMaxJ, i, j + 1) / 3 - Norme(Pedal, nMaxJ, i, j + 2) / 12;
						ddy = -Pedal[i * nMaxJ + j - 2].z / 12 + 4 * Pedal[i * nMaxJ + j - 1].z / 3 - 5 * Pedal[i * nMaxJ + j].z / 2 + 4 * Pedal[i * nMaxJ + j + 1].z / 3 - Pedal[i * nMaxJ + j + 2].z / 12;

					}

					float X, Y;

					X = dy * (x * dy - dx * y) / (dx * dx + dy * dy);
					Y = -dx * (x * dy - dx * y) / (dx * dx + dy * dy);

					int n = 1;
					//				printf("(%3d, %3d) --> %15.5f   ", i, j, X);
					X = (X > 2 ? 2 : X < -2 ? -2 : X);
					Y = (Y > 2 ? 2 : Y < -2 ? -2 : Y);
					switch (n)
					{
					case 1:
						Temp.x = Pedal[i * nMaxJ + j].x * X / x;
						Temp.y = Pedal[i * nMaxJ + j].y * X / x;
						Temp.z = Y;
						break;

					case 2:
						Temp.x = Pedal[i * nMaxJ + j - 1].x + Pedal[i * nMaxJ + j + 1].x - 2 * Pedal[i * nMaxJ + j].x;
						Temp.y = Pedal[i * nMaxJ + j - 1].y + Pedal[i * nMaxJ + j + 1].y - 2 * Pedal[i * nMaxJ + j].y;
						Temp.z = Pedal[i * nMaxJ + j - 1].z + Pedal[i * nMaxJ + j + 1].z - 2 * Pedal[i * nMaxJ + j].z;
						break;

					}
				}

				if (myDome.bDiffGeomParallel)
				{
					float dx = Pedal[myMod(i - 2, nMaxI - 1) * nMaxJ + j].x / 12 - 2 * Pedal[myMod(i - 1, nMaxI - 1) * nMaxJ + j].x / 3 +
						   2 * Pedal[myMod(i + 1, nMaxI - 1) * nMaxJ + j].x / 3      - Pedal[myMod(i + 2, nMaxI - 1) * nMaxJ + j].x / 12;
					float dy = Pedal[myMod(i - 2, nMaxI - 1) * nMaxJ + j].y / 12 - 2 * Pedal[myMod(i - 1, nMaxI - 1) * nMaxJ + j].y / 3 +
					       2 * Pedal[myMod(i + 1, nMaxI - 1) * nMaxJ + j].y / 3      - Pedal[myMod(i + 2, nMaxI - 1) * nMaxJ + j].y / 12;

					float ddx = -Pedal[myMod(i - 2, nMaxI - 1) * nMaxJ + j].x / 12 + 4 * Pedal[myMod(i - 1, nMaxI - 1) * nMaxJ + j].x / 3 -
						     5 * Pedal[i * nMaxJ + j].x / 2 +
					     	 4 * Pedal[myMod(i + 1, nMaxI - 1) * nMaxJ + j].x / 3      - Pedal[myMod(i + 2, nMaxI - 1) * nMaxJ + j].x / 12;
					float ddy = -Pedal[myMod(i - 2, nMaxI - 1) * nMaxJ + j].y / 12 + 4 * Pedal[myMod(i - 1, nMaxI - 1) * nMaxJ + j].y / 3 -
						     5 * Pedal[i * nMaxJ + j].y / 2 +
						     4 * Pedal[myMod(i + 1, nMaxI - 1) * nMaxJ + j].y / 3      - Pedal[myMod(i + 2, nMaxI - 1) * nMaxJ + j].y / 12;

					float x = Pedal[i * nMaxJ + j].x;
					float y = Pedal[i * nMaxJ + j].y;

					float X, Y;

					X =  dy * (x * dy - dx * y) / (dx * dx + dy * dy);
					Y = -dx * (x * dy - dx * y) / (dx * dx + dy * dy);

//					X = (X > 1 ? 1 : X < -1 ? -1 : X);
//					Y = (Y > 1 ? 1 : Y < -1 ? -1 : Y);


					Temp.x += X;
					Temp.y += Y;
					Temp.z += Pedal[i * nMaxJ + j].z;


				}

				int n = 3;
				switch (n)
				{
				case 1:
				{
					Vertices[i * nMaxJ + j].x = Lambda * Temp.x;
					Vertices[i * nMaxJ + j].y = Lambda * Temp.y;
					Vertices[i * nMaxJ + j].z = Lambda * Temp.z;
					break;
				}
				case 2:
				{
					Vertices[i * nMaxJ + j].x = Pedal[i * nMaxJ + j].x + Lambda * Temp.x;
					Vertices[i * nMaxJ + j].y = Pedal[i * nMaxJ + j].y + Lambda * Temp.y;
					Vertices[i * nMaxJ + j].z = Pedal[i * nMaxJ + j].z + Lambda * Temp.z;
					break;
				}
				case 3:
				{
					Vertices[i * nMaxJ + j].x = Pedal[i * nMaxJ + j].x + Lambda * (Temp.x - Pedal[i * nMaxJ + j].x);
					Vertices[i * nMaxJ + j].y = Pedal[i * nMaxJ + j].y + Lambda * (Temp.y - Pedal[i * nMaxJ + j].y);
					Vertices[i * nMaxJ + j].z = Pedal[i * nMaxJ + j].z + Lambda * (Temp.z - Pedal[i * nMaxJ + j].z);
					break;
				}
				}


			}


		free(Pedal);
	}

}
void RadialOld(GLpoint* Vertices, int nMaxI, int nMaxJ)
{
	bool bRadial = (myDome.bDiffGeomMeridien || myDome.bDiffGeomParallel);
	float Lambda = myDome.fDiffGeomLambda;
	GLpoint* Radial;
	GLpoint Temp;


	if (bRadial && myDisplay.bDisplayMinor)
	{
		Radial = (GLpoint*)malloc(nMaxI * nMaxJ * sizeof(GLpoint));
		if (Radial == 0)
		{
			Message("     Memory Allocation Radial Problem   ");
			return;
		}

		memcpy_s(Radial, nMaxI * nMaxJ * sizeof(GLpoint), Vertices, nMaxI * nMaxJ * sizeof(GLpoint));

		for (int i = 0; i < nMaxI; i++)
			for (int j = 1; j < nMaxJ; j++)
			{
				Temp.x = 0.0;
				Temp.y = 0.0;
				Temp.z = 0.0;

				if (myDome.bDiffGeomMeridien)
				{
					float x = Norme(Radial, nMaxJ, i, j);
					float y = Radial[i * nMaxJ + j].z;

					float dx, dy;
					float ddx, ddy;

					if (j == 1 || j == nMaxJ - 2)
					{
						//central  1
						dx = -Norme(Radial, nMaxJ, i, j - 1) / 2 + Norme(Radial, nMaxJ, i, j + 1) / 2;
						dy = -Radial[i * nMaxJ + j - 1].z / 2 + Radial[i * nMaxJ + j + 1].z / 2;

						ddx = Norme(Radial, nMaxJ, i, j - 1) - 2 * Norme(Radial, nMaxJ, i, j) + Norme(Radial, nMaxJ, i, j + 1);
						ddy = Radial[i * nMaxJ + j - 1].z - 2 * Radial[i * nMaxJ + j].z + Radial[i * nMaxJ + j + 1].z;
					}
					else if (j == nMaxJ - 1)
					{
						//backward 2
						dx = Norme(Radial, nMaxJ, i, j - 2) / 2 - 2 * Norme(Radial, nMaxJ, i, j - 1) + 3 * Norme(Radial, nMaxJ, i, j) / 2;
						dy = Radial[i * nMaxJ + j - 2].z / 2 - 2 * Radial[i * nMaxJ + j - 1].z + 3 * Radial[i * nMaxJ + j].z / 2;

						ddx = -Norme(Radial, nMaxJ, i, j - 3) + 4 * Norme(Radial, nMaxJ, i, j - 2) - 5 * Norme(Radial, nMaxJ, i, j - 1) + 2 * Norme(Radial, nMaxJ, i, j);
						ddy = -Radial[i * nMaxJ + j - 3].z + 4 * Radial[i * nMaxJ + j - 2].z - 5 * Radial[i * nMaxJ + j - 1].z + 2 * Radial[i * nMaxJ + j].z;


					}
					else
					{
						//central 2
						dx = Norme(Radial, nMaxJ, i, j - 2) / 12 - 2 * Norme(Radial, nMaxJ, i, j - 1) / 3 + 2 * Norme(Radial, nMaxJ, i, j + 1) / 3 - Norme(Radial, nMaxJ, i, j + 2) / 12;
						dy = Radial[i * nMaxJ + j - 2].z / 12 - 2 * Radial[i * nMaxJ + j - 1].z / 3 + 2 * Radial[i * nMaxJ + j + 1].z / 3 - Radial[i * nMaxJ + j + 2].z / 12;

						ddx = -Norme(Radial, nMaxJ, i, j - 2) / 12 + 4 * Norme(Radial, nMaxJ, i, j - 1) / 3 - 5 * Norme(Radial, nMaxJ, i, j) / 2 + 4 * Norme(Radial, nMaxJ, i, j + 1) / 3 - Norme(Radial, nMaxJ, i, j + 2) / 12;
						ddy = -Radial[i * nMaxJ + j - 2].z / 12 + 4 * Radial[i * nMaxJ + j - 1].z / 3 - 5 * Radial[i * nMaxJ + j].z / 2 + 4 * Radial[i * nMaxJ + j + 1].z / 3 - Radial[i * nMaxJ + j + 2].z / 12;

					}

					float X, Y;

					X = -dy * (dx * dx + dy * dy) / (dx * ddy - ddx * dy);
					Y =  dx * (dx * dx + dy * dy) / (dx * ddy - ddx * dy);

					int n = 1;
					//				printf("(%3d, %3d) --> %15.5f   ", i, j, X);
					X = (X > 2 ? 2 : X < -2 ? -2 : X);
					Y = (Y > 2 ? 2 : Y < -2 ? -2 : Y);
					switch (n)
					{
					case 1:
						Temp.x = Radial[i * nMaxJ + j].x * X / x;
						Temp.y = Radial[i * nMaxJ + j].y * X / x;
						Temp.z = Y;
						break;

					case 2:
						Temp.x = Radial[i * nMaxJ + j - 1].x + Radial[i * nMaxJ + j + 1].x - 2 * Radial[i * nMaxJ + j].x;
						Temp.y = Radial[i * nMaxJ + j - 1].y + Radial[i * nMaxJ + j + 1].y - 2 * Radial[i * nMaxJ + j].y;
						Temp.z = Radial[i * nMaxJ + j - 1].z + Radial[i * nMaxJ + j + 1].z - 2 * Radial[i * nMaxJ + j].z;
						break;

					}
				}

				if (myDome.bDiffGeomParallel)
				{
					float dx = Radial[myMod(i - 2, nMaxI - 1) * nMaxJ + j].x / 12 - 2 * Radial[myMod(i - 1, nMaxI - 1) * nMaxJ + j].x / 3 +
						2 * Radial[myMod(i + 1, nMaxI - 1) * nMaxJ + j].x / 3 - Radial[myMod(i + 2, nMaxI - 1) * nMaxJ + j].x / 12;
					float dy = Radial[myMod(i - 2, nMaxI - 1) * nMaxJ + j].y / 12 - 2 * Radial[myMod(i - 1, nMaxI - 1) * nMaxJ + j].y / 3 +
						2 * Radial[myMod(i + 1, nMaxI - 1) * nMaxJ + j].y / 3 - Radial[myMod(i + 2, nMaxI - 1) * nMaxJ + j].y / 12;

					float ddx = -Radial[myMod(i - 2, nMaxI - 1) * nMaxJ + j].x / 12 + 4 * Radial[myMod(i - 1, nMaxI - 1) * nMaxJ + j].x / 3 -
						5 * Radial[i * nMaxJ + j].x / 2 +
						4 * Radial[myMod(i + 1, nMaxI - 1) * nMaxJ + j].x / 3 - Radial[myMod(i + 2, nMaxI - 1) * nMaxJ + j].x / 12;
					float ddy = -Radial[myMod(i - 2, nMaxI - 1) * nMaxJ + j].y / 12 + 4 * Radial[myMod(i - 1, nMaxI - 1) * nMaxJ + j].y / 3 -
						5 * Radial[i * nMaxJ + j].y / 2 +
						4 * Radial[myMod(i + 1, nMaxI - 1) * nMaxJ + j].y / 3 - Radial[myMod(i + 2, nMaxI - 1) * nMaxJ + j].y / 12;

					float x = Radial[i * nMaxJ + j].x;
					float y = Radial[i * nMaxJ + j].y;

					float X, Y;

					X = -dy * (dx * dx + dy * dy) / (dx * ddy - ddx * dy);
					Y =  dx * (dx * dx + dy * dy) / (dx * ddy - ddx * dy);

					//					X = (X > 1 ? 1 : X < -1 ? -1 : X);
					//					Y = (Y > 1 ? 1 : Y < -1 ? -1 : Y);


					Temp.x += X;
					Temp.y += Y;
					Temp.z += Radial[i * nMaxJ + j].z;


				}

				int n = 3;
				switch (n)
				{
				case 1:
				{
					Vertices[i * nMaxJ + j].x = Lambda * Temp.x;
					Vertices[i * nMaxJ + j].y = Lambda * Temp.y;
					Vertices[i * nMaxJ + j].z = Lambda * Temp.z;
					break;
				}
				case 2:
				{
					Vertices[i * nMaxJ + j].x = Radial[i * nMaxJ + j].x + Lambda * Temp.x;
					Vertices[i * nMaxJ + j].y = Radial[i * nMaxJ + j].y + Lambda * Temp.y;
					Vertices[i * nMaxJ + j].z = Radial[i * nMaxJ + j].z + Lambda * Temp.z;
					break;
				}
				case 3:
				{
					Vertices[i * nMaxJ + j].x = Radial[i * nMaxJ + j].x + Lambda * (Temp.x - Radial[i * nMaxJ + j].x);
					Vertices[i * nMaxJ + j].y = Radial[i * nMaxJ + j].y + Lambda * (Temp.y - Radial[i * nMaxJ + j].y);
					Vertices[i * nMaxJ + j].z = Radial[i * nMaxJ + j].z + Lambda * (Temp.z - Radial[i * nMaxJ + j].z);
					break;
				}
				}


			}


		free(Radial);
	}

}
void EvoluteOld(GLpoint* Vertices, int nMaxI, int nMaxJ)
{
	bool bEvolute = (myDome.bDiffGeomMeridien || myDome.bDiffGeomParallel);
	float Lambda = myDome.fDiffGeomLambda;
	GLpoint* Evolute;
	GLpoint Temp;

	switch(myDome.nDifferentialGeometry)
	{
	case DIFF_GEOM_EVOLUTE:
		break;
	case DIFF_GEOM_PEDAL:
		PedalOld(Vertices, nMaxI, nMaxJ);
		return;
	case DIFF_GEOM_RADIAL:
		RadialOld(Vertices, nMaxI, nMaxJ);
		return;
	}

	if (bEvolute && myDisplay.bDisplayMinor)
	{
		Evolute = (GLpoint*)malloc(nMaxI * nMaxJ * sizeof(GLpoint));
		if (Evolute == 0)
		{
			Message("     Memory Allocation Evolute Problem   ");
			return;
		}

		memcpy_s(Evolute, nMaxI * nMaxJ * sizeof(GLpoint), Vertices, nMaxI * nMaxJ * sizeof(GLpoint));

		for (int i = 0; i < nMaxI; i++)
			for (int j = 1; j < nMaxJ; j++)
			{
				Temp.x = 0.0;
				Temp.y = 0.0;
				Temp.z = 0.0;

				if (myDome.bDiffGeomMeridien)
				{
					float x = Norme(Evolute, nMaxJ, i, j);
					float y = Evolute[i * nMaxJ + j].z;

					float dx, dy;
					float ddx, ddy;

					if (j == 1 || j == nMaxJ - 2)
					{
						//central  1
						dx = -Norme(Evolute, nMaxJ, i, j - 1) / 2 + Norme(Evolute, nMaxJ, i, j + 1) / 2;
						dy = -Evolute[i * nMaxJ + j - 1].z / 2 + Evolute[i * nMaxJ + j + 1].z / 2;

						ddx = Norme(Evolute, nMaxJ, i, j - 1) - 2 * Norme(Evolute, nMaxJ, i, j) + Norme(Evolute, nMaxJ, i, j + 1);
						ddy = Evolute[i * nMaxJ + j - 1].z - 2 * Evolute[i * nMaxJ + j].z + Evolute[i * nMaxJ + j + 1].z;
					}
					else if (j == nMaxJ - 1)
					{
						//backward 2
						dx = Norme(Evolute, nMaxJ, i, j - 2) / 2 - 2 * Norme(Evolute, nMaxJ, i, j - 1) + 3 * Norme(Evolute, nMaxJ, i, j) / 2;
						dy = Evolute[i * nMaxJ + j - 2].z / 2 - 2 * Evolute[i * nMaxJ + j - 1].z + 3 * Evolute[i * nMaxJ + j].z / 2;

						ddx = -Norme(Evolute, nMaxJ, i, j - 3) + 4 * Norme(Evolute, nMaxJ, i, j - 2) - 5 * Norme(Evolute, nMaxJ, i, j - 1) + 2 * Norme(Evolute, nMaxJ, i, j);
						ddy = -Evolute[i * nMaxJ + j - 3].z + 4 * Evolute[i * nMaxJ + j - 2].z - 5 * Evolute[i * nMaxJ + j - 1].z + 2 * Evolute[i * nMaxJ + j].z;


					}
					else
					{
						//central 2
						dx = Norme(Evolute, nMaxJ, i, j - 2) / 12 - 2 * Norme(Evolute, nMaxJ, i, j - 1) / 3 + 2 * Norme(Evolute, nMaxJ, i, j + 1) / 3 - Norme(Evolute, nMaxJ, i, j + 2) / 12;
						dy = Evolute[i * nMaxJ + j - 2].z / 12 - 2 * Evolute[i * nMaxJ + j - 1].z / 3 + 2 * Evolute[i * nMaxJ + j + 1].z / 3 - Evolute[i * nMaxJ + j + 2].z / 12;

						ddx = -Norme(Evolute, nMaxJ, i, j - 2) / 12 + 4 * Norme(Evolute, nMaxJ, i, j - 1) / 3 - 5 * Norme(Evolute, nMaxJ, i, j) / 2 + 4 * Norme(Evolute, nMaxJ, i, j + 1) / 3 - Norme(Evolute, nMaxJ, i, j + 2) / 12;
						ddy = -Evolute[i * nMaxJ + j - 2].z / 12 + 4 * Evolute[i * nMaxJ + j - 1].z / 3 - 5 * Evolute[i * nMaxJ + j].z / 2 + 4 * Evolute[i * nMaxJ + j + 1].z / 3 - Evolute[i * nMaxJ + j + 2].z / 12;

					}

					float X, Y;

					X = x - (dy * (dx * dx + dy * dy)) / (dx * ddy - ddx * dy);
					Y = y + (dx * (dx * dx + dy * dy)) / (dx * ddy - ddx * dy);
					int n = 1;
					//				printf("(%3d, %3d) --> %15.5f   ", i, j, X);
					X = (X > 1 ? 1 : X < -1 ? -1 : X);
					Y = (Y > 1 ? 1 : Y < -1 ? -1 : Y);
					switch (n)
					{
					case 1:
						Temp.x = Evolute[i * nMaxJ + j].x * X / x;
						Temp.y = Evolute[i * nMaxJ + j].y * X / x;
						Temp.z = Y;
						break;

					case 2:
						Temp.x = Evolute[i * nMaxJ + j - 1].x + Evolute[i * nMaxJ + j + 1].x - 2 * Evolute[i * nMaxJ + j].x;
						Temp.y = Evolute[i * nMaxJ + j - 1].y + Evolute[i * nMaxJ + j + 1].y - 2 * Evolute[i * nMaxJ + j].y;
						Temp.z = Evolute[i * nMaxJ + j - 1].z + Evolute[i * nMaxJ + j + 1].z - 2 * Evolute[i * nMaxJ + j].z;
						break;

					}
				}

				if (myDome.bDiffGeomParallel)
				{
					float dx = Evolute[myMod(i - 2, nMaxI - 1) * nMaxJ + j].x / 12 - 2 * Evolute[myMod(i - 1, nMaxI - 1) * nMaxJ + j].x / 3 +
						2 * Evolute[myMod(i + 1, nMaxI - 1) * nMaxJ + j].x / 3 - Evolute[myMod(i + 2, nMaxI - 1) * nMaxJ + j].x / 12;
					float dy = Evolute[myMod(i - 2, nMaxI - 1) * nMaxJ + j].y / 12 - 2 * Evolute[myMod(i - 1, nMaxI - 1) * nMaxJ + j].y / 3 +
						2 * Evolute[myMod(i + 1, nMaxI - 1) * nMaxJ + j].y / 3 - Evolute[myMod(i + 2, nMaxI - 1) * nMaxJ + j].y / 12;

					float ddx = -Evolute[myMod(i - 2, nMaxI - 1) * nMaxJ + j].x / 12 + 4 * Evolute[myMod(i - 1, nMaxI - 1) * nMaxJ + j].x / 3 -
						5 * Evolute[i * nMaxJ + j].x / 2 +
						4 * Evolute[myMod(i + 1, nMaxI - 1) * nMaxJ + j].x / 3 - Evolute[myMod(i + 2, nMaxI - 1) * nMaxJ + j].x / 12;
					float ddy = -Evolute[myMod(i - 2, nMaxI - 1) * nMaxJ + j].y / 12 + 4 * Evolute[myMod(i - 1, nMaxI - 1) * nMaxJ + j].y / 3 -
						5 * Evolute[i * nMaxJ + j].y / 2 +
						4 * Evolute[myMod(i + 1, nMaxI - 1) * nMaxJ + j].y / 3 - Evolute[myMod(i + 2, nMaxI - 1) * nMaxJ + j].y / 12;

					float x = Evolute[i * nMaxJ + j].x;
					float y = Evolute[i * nMaxJ + j].y;

					float X, Y;

					X = x - (dy * (dx * dx + dy * dy)) / (dx * ddy - ddx * dy);
					Y = y + (dx * (dx * dx + dy * dy)) / (dx * ddy - ddx * dy);

					X = (X > 2 ? 2 : X < -2 ? -2 : X);
					Y = (Y > 2 ? 2 : Y < -2 ? -2 : Y);


					Temp.x += X;
					Temp.y += Y;
					Temp.z += Evolute[i * nMaxJ + j].z;


				}

				int n = 3;
				switch (n)
				{
				case 1:
				{
					Vertices[i * nMaxJ + j].x = Lambda * Temp.x;
					Vertices[i * nMaxJ + j].y = Lambda * Temp.y;
					Vertices[i * nMaxJ + j].z = Lambda * Temp.z;
					break;
				}
				case 2:
				{
					Vertices[i * nMaxJ + j].x = Evolute[i * nMaxJ + j].x + Lambda * Temp.x;
					Vertices[i * nMaxJ + j].y = Evolute[i * nMaxJ + j].y + Lambda * Temp.y;
					Vertices[i * nMaxJ + j].z = Evolute[i * nMaxJ + j].z + Lambda * Temp.z;
					break;
				}
				case 3:
				{
					Vertices[i * nMaxJ + j].x = Evolute[i * nMaxJ + j].x + Lambda * (Temp.x - Evolute[i * nMaxJ + j].x);
					Vertices[i * nMaxJ + j].y = Evolute[i * nMaxJ + j].y + Lambda * (Temp.y - Evolute[i * nMaxJ + j].y);
					Vertices[i * nMaxJ + j].z = Evolute[i * nMaxJ + j].z + Lambda * (Temp.z - Evolute[i * nMaxJ + j].z);
					break;
				}
				}


			}


		free(Evolute);
	}

}

void DifferentialGeometry(GLpoint* Vertices, int nMaxI, int nMaxJ)
{
	bool bDifferentialGeometry = (myDome.bDiffGeomMeridien || myDome.bDiffGeomParallel);
	float Lambda = myDome.fDiffGeomLambda;
	GLpoint* DifferentialGeometry;
	GLpoint Temp;
	int n;

	if (bDifferentialGeometry && myDisplay.bDisplayMinor)
	{
		DifferentialGeometry = (GLpoint*)malloc(nMaxI * nMaxJ * sizeof(GLpoint));
		if (DifferentialGeometry == 0)
		{
			Message("     Memory Allocation DifferentialGeometry Problem   ");
			return;
		}

		memcpy_s(DifferentialGeometry, nMaxI * nMaxJ * sizeof(GLpoint), Vertices, nMaxI * nMaxJ * sizeof(GLpoint));

		for (int i = 0; i < nMaxI; i++)
			for (int j = 1; j < nMaxJ; j++)
			{
				Temp.x = 0.0;
				Temp.y = 0.0;
				Temp.z = 0.0;

				if (myDome.bDiffGeomMeridien)
				{
					float x = Norme(DifferentialGeometry, nMaxJ, i, j);
					float y = DifferentialGeometry[i * nMaxJ + j].z;

					float dx, dy;
					float ddx, ddy;

					if (j == 1 || j == nMaxJ - 2)
					{
						//central  1
						dx = -Norme(DifferentialGeometry, nMaxJ, i, j - 1) / 2 + Norme(DifferentialGeometry, nMaxJ, i, j + 1) / 2;
						dy = -DifferentialGeometry[i * nMaxJ + j - 1].z / 2 + DifferentialGeometry[i * nMaxJ + j + 1].z / 2;

						ddx = Norme(DifferentialGeometry, nMaxJ, i, j - 1) - 2 * Norme(DifferentialGeometry, nMaxJ, i, j) + Norme(DifferentialGeometry, nMaxJ, i, j + 1);
						ddy = DifferentialGeometry[i * nMaxJ + j - 1].z - 2 * DifferentialGeometry[i * nMaxJ + j].z + DifferentialGeometry[i * nMaxJ + j + 1].z;
					}
					else if (j == nMaxJ - 1)
					{
						//backward 2
						dx = Norme(DifferentialGeometry, nMaxJ, i, j - 2) / 2 - 2 * Norme(DifferentialGeometry, nMaxJ, i, j - 1) + 3 * Norme(DifferentialGeometry, nMaxJ, i, j) / 2;
						dy = DifferentialGeometry[i * nMaxJ + j - 2].z / 2 - 2 * DifferentialGeometry[i * nMaxJ + j - 1].z + 3 * DifferentialGeometry[i * nMaxJ + j].z / 2;

						ddx = -Norme(DifferentialGeometry, nMaxJ, i, j - 3) + 4 * Norme(DifferentialGeometry, nMaxJ, i, j - 2) - 5 * Norme(DifferentialGeometry, nMaxJ, i, j - 1) + 2 * Norme(DifferentialGeometry, nMaxJ, i, j);
						ddy = -DifferentialGeometry[i * nMaxJ + j - 3].z + 4 * DifferentialGeometry[i * nMaxJ + j - 2].z - 5 * DifferentialGeometry[i * nMaxJ + j - 1].z + 2 * DifferentialGeometry[i * nMaxJ + j].z;


					}
					else
					{
						//central 2
						dx = Norme(DifferentialGeometry, nMaxJ, i, j - 2) / 12 - 2 * Norme(DifferentialGeometry, nMaxJ, i, j - 1) / 3 + 2 * Norme(DifferentialGeometry, nMaxJ, i, j + 1) / 3 - Norme(DifferentialGeometry, nMaxJ, i, j + 2) / 12;
						dy = DifferentialGeometry[i * nMaxJ + j - 2].z / 12 - 2 * DifferentialGeometry[i * nMaxJ + j - 1].z / 3 + 2 * DifferentialGeometry[i * nMaxJ + j + 1].z / 3 - DifferentialGeometry[i * nMaxJ + j + 2].z / 12;

						ddx = -Norme(DifferentialGeometry, nMaxJ, i, j - 2) / 12 + 4 * Norme(DifferentialGeometry, nMaxJ, i, j - 1) / 3 - 5 * Norme(DifferentialGeometry, nMaxJ, i, j) / 2 + 4 * Norme(DifferentialGeometry, nMaxJ, i, j + 1) / 3 - Norme(DifferentialGeometry, nMaxJ, i, j + 2) / 12;
						ddy = -DifferentialGeometry[i * nMaxJ + j - 2].z / 12 + 4 * DifferentialGeometry[i * nMaxJ + j - 1].z / 3 - 5 * DifferentialGeometry[i * nMaxJ + j].z / 2 + 4 * DifferentialGeometry[i * nMaxJ + j + 1].z / 3 - DifferentialGeometry[i * nMaxJ + j + 2].z / 12;

					}

					float X, Y;

					switch (myDome.nDifferentialGeometry)
					{
					case DIFF_GEOM_EVOLUTE:
						X = x - (dy * (dx * dx + dy * dy)) / (dx * ddy - ddx * dy);
						Y = y + (dx * (dx * dx + dy * dy)) / (dx * ddy - ddx * dy);
						break;
					case DIFF_GEOM_PEDAL:
						X =  dy * (x * dy - dx * y) / (dx * dx + dy * dy);
						Y = -dx * (x * dy - dx * y) / (dx * dx + dy * dy);
						break;
					case DIFF_GEOM_RADIAL:
						X = -dy * (dx * dx + dy * dy) / (dx * ddy - ddx * dy);
						Y =  dx * (dx * dx + dy * dy) / (dx * ddy - ddx * dy);
						break;
					}

					if (X > 2 || X < -2 || Y > 2 || Y < -2)
						printf("Dépassement Méridien\n");
					n = 1;
					//				printf("(%3d, %3d) --> %15.5f   ", i, j, X);
					X = (X > 2 ? 2 : X < -2 ? -2 : X);
					Y = (Y > 2 ? 2 : Y < -2 ? -2 : Y);
					switch (n)
					{
					case 1:
						Temp.x = DifferentialGeometry[i * nMaxJ + j].x * Lambda * X / x;
						Temp.y = DifferentialGeometry[i * nMaxJ + j].y * Lambda * X / x;
						Temp.z = Y * Lambda;
						break;

					case 2:
						Temp.x = DifferentialGeometry[i * nMaxJ + j - 1].x + DifferentialGeometry[i * nMaxJ + j + 1].x - 2 * DifferentialGeometry[i * nMaxJ + j].x;
						Temp.y = DifferentialGeometry[i * nMaxJ + j - 1].y + DifferentialGeometry[i * nMaxJ + j + 1].y - 2 * DifferentialGeometry[i * nMaxJ + j].y;
						Temp.z = DifferentialGeometry[i * nMaxJ + j - 1].z + DifferentialGeometry[i * nMaxJ + j + 1].z - 2 * DifferentialGeometry[i * nMaxJ + j].z;
						break;

					}
				}

				if (myDome.bDiffGeomParallel)
				{
					float dx = DifferentialGeometry[myMod(i - 2, nMaxI - 1) * nMaxJ + j].x / 12 - 2 * DifferentialGeometry[myMod(i - 1, nMaxI - 1) * nMaxJ + j].x / 3 +
						2 * DifferentialGeometry[myMod(i + 1, nMaxI - 1) * nMaxJ + j].x / 3 - DifferentialGeometry[myMod(i + 2, nMaxI - 1) * nMaxJ + j].x / 12;
					float dy = DifferentialGeometry[myMod(i - 2, nMaxI - 1) * nMaxJ + j].y / 12 - 2 * DifferentialGeometry[myMod(i - 1, nMaxI - 1) * nMaxJ + j].y / 3 +
						2 * DifferentialGeometry[myMod(i + 1, nMaxI - 1) * nMaxJ + j].y / 3 - DifferentialGeometry[myMod(i + 2, nMaxI - 1) * nMaxJ + j].y / 12;

					float ddx = -DifferentialGeometry[myMod(i - 2, nMaxI - 1) * nMaxJ + j].x / 12 + 4 * DifferentialGeometry[myMod(i - 1, nMaxI - 1) * nMaxJ + j].x / 3 -
						5 * DifferentialGeometry[i * nMaxJ + j].x / 2 +
						4 * DifferentialGeometry[myMod(i + 1, nMaxI - 1) * nMaxJ + j].x / 3 - DifferentialGeometry[myMod(i + 2, nMaxI - 1) * nMaxJ + j].x / 12;
					float ddy = -DifferentialGeometry[myMod(i - 2, nMaxI - 1) * nMaxJ + j].y / 12 + 4 * DifferentialGeometry[myMod(i - 1, nMaxI - 1) * nMaxJ + j].y / 3 -
						5 * DifferentialGeometry[i * nMaxJ + j].y / 2 +
						4 * DifferentialGeometry[myMod(i + 1, nMaxI - 1) * nMaxJ + j].y / 3 - DifferentialGeometry[myMod(i + 2, nMaxI - 1) * nMaxJ + j].y / 12;

					float x = DifferentialGeometry[i * nMaxJ + j].x;
					float y = DifferentialGeometry[i * nMaxJ + j].y;

					float X, Y;

					switch (myDome.nDifferentialGeometry)
					{
					case DIFF_GEOM_EVOLUTE:
						X = x - (dy * (dx * dx + dy * dy)) / (dx * ddy - ddx * dy);
						Y = y + (dx * (dx * dx + dy * dy)) / (dx * ddy - ddx * dy);
						break;
					case DIFF_GEOM_PEDAL:
						X =  dy * (x * dy - dx * y) / (dx * dx + dy * dy);
						Y = -dx * (x * dy - dx * y) / (dx * dx + dy * dy);
						break;
					case DIFF_GEOM_RADIAL:
						X = -dy * (dx * dx + dy * dy) / (dx * ddy - ddx * dy);
						Y =  dx * (dx * dx + dy * dy) / (dx * ddy - ddx * dy);
						break;
					}

					//					X = (X > 1 ? 1 : X < -1 ? -1 : X);
					//					Y = (Y > 1 ? 1 : Y < -1 ? -1 : Y);

					if (X > 2 || X < -2 || Y > 2 || Y < -2)
						printf("Dépassement Parallele\n");

					X = (X > 2 ? 2 : X < -2 ? -2 : X);
					Y = (Y > 2 ? 2 : Y < -2 ? -2 : Y);

					Temp.x += X;
					Temp.y += Y;
					Temp.z += DifferentialGeometry[i * nMaxJ + j].z;


				}


				n = 3;
				switch (n)
				{
				case 1:
					Vertices[i * nMaxJ + j].x = Lambda * Temp.x;
					Vertices[i * nMaxJ + j].y = Lambda * Temp.y;
					Vertices[i * nMaxJ + j].z = Lambda * Temp.z;
					break;
				case 2:
					Vertices[i * nMaxJ + j].x = DifferentialGeometry[i * nMaxJ + j].x + Lambda * Temp.x;
					Vertices[i * nMaxJ + j].y = DifferentialGeometry[i * nMaxJ + j].y + Lambda * Temp.y;
					Vertices[i * nMaxJ + j].z = DifferentialGeometry[i * nMaxJ + j].z + Lambda * Temp.z;
					break;
				case 3:
					Vertices[i * nMaxJ + j].x = DifferentialGeometry[i * nMaxJ + j].x + Lambda * (Temp.x - DifferentialGeometry[i * nMaxJ + j].x);
					Vertices[i * nMaxJ + j].y = DifferentialGeometry[i * nMaxJ + j].y + Lambda * (Temp.y - DifferentialGeometry[i * nMaxJ + j].y);
					Vertices[i * nMaxJ + j].z = DifferentialGeometry[i * nMaxJ + j].z + Lambda * (Temp.z - DifferentialGeometry[i * nMaxJ + j].z);
					break;
				}

			}


		free(DifferentialGeometry);
	}

}

void Lissage(GLpoint* Vertices, int nMaxI, int nMaxJ)
{
	int bLissage = myDome.bLissage;
	int nPuissance = myDome.nLissageDistanceType;
	int nDistance = myDome.nLissageDistanceType;

	int nTaille = myDome.nLissageTaille;
	float fPower = myDome.fLissagePuissance;
	FilterParameter Parameters;
	int nSize = 2 * nTaille + 1;
	int n_I, n_J;
	float PoidsTotal;
	GLpoint* Lissage, Temp;
	float* Matrix;
	float fTotal;

	if (myDome.nDomeEnd != 360)
		nMaxI = (int)((float)(nMaxI * myDome.nDomeEnd) / 360.0) + 1;


	if (bLissage)
	{
		Lissage = (GLpoint*)malloc(nMaxI * nMaxJ * sizeof(GLpoint));
		Matrix = (float*)malloc(nSize * nSize * sizeof(float));

		if (Lissage == NULL || Matrix == NULL)
		{
			Message("     Memory Allocation Lissage Problem   ");
			return;
		}

		Parameters.fDiagonal = myDome.fLissagePoidsDiagonales;
		Parameters.fHorizontal = myDome.fLissagePoidsHorizontal;
		Parameters.fVertical = myDome.fLissagePoidsVertical;
		Parameters.fSelf = myDome.fLissagePoidsSelf;
		Parameters.fSigma1 = myDome.fLissageSigma1;
		Parameters.fSigma2 = myDome.fLissageSigma2;
		Parameters.nDirection = (DIRECTION_Type)myDome.nLissageDirection;

		fTotal = FiltreListe[myDome.nLissageTypeFiltre].Function(Matrix, nTaille, (DISTANCE_Type)nDistance, fPower, &Parameters);
//		fTotal = Gauss(Matrix, nTaille, nDistance, fPower, &Parameters);
		PrintMatrix(Matrix, nTaille, fTotal, TYPE_Standard);

		memcpy_s(Lissage, nMaxI * nMaxJ * sizeof(GLpoint), Vertices, nMaxI * nMaxJ * sizeof(GLpoint));
		for (int i = 0; i < nMaxI; i++)
			for (int j = 1; j < nMaxJ; j++)
			{
				Temp.x = 0.0;
				Temp.y = 0.0;
				Temp.z = 0.0;
				PoidsTotal = 0;;

				for (int a = -nTaille; a <= nTaille; a++)
//					for (int b = (j == nMaxJ - 1 ? 0 : -nTaille); b <= (j == nMaxJ - 1 ? 0 : nTaille); b++)
					for (int b = -nTaille; b <= nTaille; b++)
					{
						Coord(i, j, a, b, nMaxI, nMaxJ, &n_I, &n_J);
						if (n_J != -1 && n_I != -1)
						{

							Temp.x += Lissage[n_I * nMaxJ + n_J].x * Matrix[(a + nTaille) * nSize + (b + nTaille)];
							Temp.y += Lissage[n_I * nMaxJ + n_J].y * Matrix[(a + nTaille) * nSize + (b + nTaille)];
							Temp.z += Lissage[n_I * nMaxJ + n_J].z * Matrix[(a + nTaille) * nSize + (b + nTaille)];
							PoidsTotal += Matrix[(a + nTaille) * nSize + (b + nTaille)];
							if (i == 14 && j == 38 && i==j)
							{
//								printf("L=%f  ;   ", Lissage[n_I * nMaxJ + n_J].x);
								if (b == nTaille)
									printf("\n");
							}

						}
					}
				if (i == 14 && j == 38 && i ==j)
					printf("Temp=%f ;Total =%f ; Poids = %f", Temp.x, fTotal, PoidsTotal);

				if (PoidsTotal == 0)
					PoidsTotal = 1;
				Vertices[i * nMaxJ + j].x = Temp.x / PoidsTotal;
				Vertices[i * nMaxJ + j].y = Temp.y / PoidsTotal;
				Vertices[i * nMaxJ + j].z = Temp.z / PoidsTotal;

//				if (Vertices[i * nMaxJ + j].y > 1 || Vertices[i * nMaxJ + j].y < -1)
//					printf("i = %d, j = %d, .x = %f\n", i, j, Vertices[i * nMaxJ + j].y);
			}

		free(Lissage);
	}

}
void ComputeVector(GLpoint* Vertices, int nMaxI, int nMaxJ, int i, int j, GLpoint *V, int nMode, int nNormale)
{
	float N = sqrtf(Vertices[i * nMaxJ + j].x * Vertices[i * nMaxJ + j].x + Vertices[i * nMaxJ + j].y * Vertices[i * nMaxJ + j].y + Vertices[i * nMaxJ + j].z * Vertices[i * nMaxJ + j].z);
	float fFacteur;
	float fLambda = myDome.fEnhanceLambda;
	float fDamping = myDome.fEnhanceDamping;
	int n = myDome.nMeridianMain * myDome.nBaseSides * (myDome.nMeridianMinor + 1) ;

	float x = (float)j / (float)(nMaxJ - 1);

	if (nMode == CALCUL_Somme)
		fFacteur = (1 + fLambda / (10 * N));
	else
		fFacteur = fLambda;


//	fFacteur *= 1 - (1 - x) * fDamping;
	fFacteur = fFacteur + (1 - fFacteur) * (1 - x) * fDamping;

	if ("cylindical" && nNormale == ENHANCE_CYLINDRICAL)
	{
		V->x = Vertices[i * nMaxJ + j].x * fFacteur - Vertices[i * nMaxJ + j].x;
		V->y = Vertices[i * nMaxJ + j].y * fFacteur - Vertices[i * nMaxJ + j].y;
		V->z = Vertices[i * nMaxJ + j].z - Vertices[i * nMaxJ + j].z;
	}

	if ("vertical" && nNormale == ENHANCE_VERTICAL)
	{
		V->x = Vertices[i * nMaxJ + j].x - Vertices[i * nMaxJ + j].x;
		V->y = Vertices[i * nMaxJ + j].y - Vertices[i * nMaxJ + j].y;
		V->z = Vertices[i * nMaxJ + j].z * fFacteur - Vertices[i * nMaxJ + j].z;
	}

	if ("radial" && nNormale == ENHANCE_RADIAL)
	{
		V->x = Vertices[i * nMaxJ + j].x * fFacteur - Vertices[i * nMaxJ + j].x;
		V->y = Vertices[i * nMaxJ + j].y * fFacteur - Vertices[i * nMaxJ + j].y;
		V->z = Vertices[i * nMaxJ + j].z * fFacteur - Vertices[i * nMaxJ + j].z;
	}

	/*
			A
			|
		  	|
	D------ O ------B
			|
			|
			C
	
	*/

	if ("OB OA" && nNormale == ENHANCE_OB_OA)
	{
		GLpoint v0, v1, v2;
		fFacteur *= 300;

		v0.x = Vertices[i * nMaxJ + j - 1].x - Vertices[i * nMaxJ + j].x;
		v0.y = Vertices[i * nMaxJ + j - 1].y - Vertices[i * nMaxJ + j].y;
		v0.z = Vertices[i * nMaxJ + j - 1].z - Vertices[i * nMaxJ + j].z;

		v1.x = Vertices[myMod(i + 1, n) * nMaxJ + j].x - Vertices[i * nMaxJ + j].x;
		v1.y = Vertices[myMod(i + 1, n) * nMaxJ + j].y - Vertices[i * nMaxJ + j].y;
		v1.z = Vertices[myMod(i + 1, n) * nMaxJ + j].z - Vertices[i * nMaxJ + j].z;

		CrossProduct(&v0, &v1, &v2);

		V->x = Vertices[i * nMaxJ + j].x - v2.x * fFacteur - Vertices[i * nMaxJ + j].x;
		V->y = Vertices[i * nMaxJ + j].y - v2.y * fFacteur - Vertices[i * nMaxJ + j].y;
		V->z = Vertices[i * nMaxJ + j].z - v2.z * fFacteur - Vertices[i * nMaxJ + j].z;


	}
	if ("DB AC" && nNormale == ENHANCE_DB_AC)
	{
		GLpoint v0, v1, v2;
		fFacteur *= 100;
		if (j < nMaxJ - 1)
		{
			v0.x = Vertices[i * nMaxJ + j - 1].x - Vertices[i * nMaxJ + j + 1].x;
			v0.y = Vertices[i * nMaxJ + j - 1].y - Vertices[i * nMaxJ + j + 1].y;
			v0.z = Vertices[i * nMaxJ + j - 1].z - Vertices[i * nMaxJ + j + 1].z;

			v1.x = Vertices[myMod(i + 1, n) * nMaxJ + j].x - Vertices[myMod(i - 1, n) * nMaxJ + j].x;
			v1.y = Vertices[myMod(i + 1, n) * nMaxJ + j].y - Vertices[myMod(i - 1, n) * nMaxJ + j].y;
			v1.z = Vertices[myMod(i + 1, n) * nMaxJ + j].z - Vertices[myMod(i - 1, n) * nMaxJ + j].z;

			CrossProduct(&v0, &v1, &v2);

			V->x = - v2.x * fFacteur;
			V->y = - v2.y * fFacteur;
			V->z = - v2.z * fFacteur;
		}
		else
		{
			V->x = 0;
			V->y = 0;
			V->z = 0;
		}

	}
	
	if ("(D + B)/2 O" && nNormale == ENHANCE_D_B_O)
	{
		GLpoint v0, v1;
		fFacteur *= 50;
		v0.x = (Vertices[myMod(i + 1, n) * nMaxJ + j].x + Vertices[myMod(i - 1, n) * nMaxJ + j].x) / 2;
		v0.y = (Vertices[myMod(i + 1, n) * nMaxJ + j].y + Vertices[myMod(i - 1, n) * nMaxJ + j].y) / 2;
		v0.z = (Vertices[myMod(i + 1, n) * nMaxJ + j].z + Vertices[myMod(i - 1, n) * nMaxJ + j].z) / 2;

//		printf("(%3d, %3d) --> (%3d, %3d) --> (%3d, %3d)\n", i, j, myMod(i + 1, n), nMaxI, myMod(i - 1, n), n);
		v1.x = Vertices[i * nMaxJ + j].x - v0.x;
		v1.y = Vertices[i * nMaxJ + j].y - v0.y;
		v1.z = Vertices[i * nMaxJ + j].z - v0.z;

//		printf("(%d, %d) --> (%f, %f, %f)\n", i, j, v0.x, v0.y, v0.z);
		V->x = Vertices[i * nMaxJ + j].x + v1.x * fFacteur - Vertices[i * nMaxJ + j].x;
		V->y = Vertices[i * nMaxJ + j].y + v1.y * fFacteur - Vertices[i * nMaxJ + j].y;
		V->z = Vertices[i * nMaxJ + j].z + v1.z * fFacteur - Vertices[i * nMaxJ + j].z;


	}
	if ("(A + C)/2 O" && nNormale == ENHANCE_A_C_O)
	{
		GLpoint v0, v1;
		fFacteur *= 300;
		if (j < nMaxJ - 1 && j > 0)

		{
			v0.x = (Vertices[i * nMaxJ + j - 1].x + Vertices[i * nMaxJ + j + 1].x) / 2;
			v0.y = (Vertices[i * nMaxJ + j - 1].y + Vertices[i * nMaxJ + j + 1].y) / 2;
			v0.z = (Vertices[i * nMaxJ + j - 1].z + Vertices[i * nMaxJ + j + 1].z) / 2;

			//		printf("(%3d, %3d) --> (%3d, %3d) --> (%3d, %3d)\n", i, j, myMod(i + 1, n), nMaxI, myMod(i - 1, n), n);
			v1.x = Vertices[i * nMaxJ + j].x - v0.x;
			v1.y = Vertices[i * nMaxJ + j].y - v0.y;
			v1.z = Vertices[i * nMaxJ + j].z - v0.z;

			//		printf("(%d, %d) --> (%f, %f, %f)\n", i, j, v0.x, v0.y, v0.z);
			V->x = v1.x * fFacteur;
			V->y = v1.y * fFacteur;
			V->z = v1.z * fFacteur;
		}
		else
		{
			V->x = 0;
			V->y = 0;
			V->z = 0;
		}

	}	
		
	if ("ABCD_O" && nNormale == ENHANCE_ABCD_O)
	{
		GLpoint v0, v1;
		fFacteur *= 200;

		if (j < nMaxJ - 1 && j > 0)
		{
			v0.x = (Vertices[i * nMaxJ + j - 1].x + Vertices[i * nMaxJ + j + 1].x + Vertices[myMod(i + 1, n) * nMaxJ + j].x + Vertices[myMod(i - 1, n) * nMaxJ + j].x) / 4;
			v0.y = (Vertices[i * nMaxJ + j - 1].y + Vertices[i * nMaxJ + j + 1].y + Vertices[myMod(i + 1, n) * nMaxJ + j].y + Vertices[myMod(i - 1, n) * nMaxJ + j].y) / 4;
			v0.z = (Vertices[i * nMaxJ + j - 1].z + Vertices[i * nMaxJ + j + 1].z + Vertices[myMod(i + 1, n) * nMaxJ + j].z + Vertices[myMod(i - 1, n) * nMaxJ + j].z) / 4;

			//		printf("(%3d, %3d) --> (%3d, %3d) --> (%3d, %3d)\n", i, j, myMod(i + 1, n), nMaxI, myMod(i - 1, n), n);
			v1.x = Vertices[i * nMaxJ + j].x - v0.x;
			v1.y = Vertices[i * nMaxJ + j].y - v0.y;
			v1.z = Vertices[i * nMaxJ + j].z - v0.z;

			//		printf("(%d, %d) --> (%f, %f, %f)\n", i, j, v0.x, v0.y, v0.z);
			V->x = v1.x * fFacteur;
			V->y = v1.y * fFacteur;
			V->z = v1.z * fFacteur;

		}
		else
		{
			V->x = 0;
			V->y = 0;
			V->z = 0;
		}
	}
	if("(OA + OB) * (OA + OD)" && nNormale == ENHANCE_OA_OB_OA_OD)
	{
		GLpoint vA, vB, vD, v0, v1, v2;
		fFacteur *= 500;

		vA.x = Vertices[i * nMaxJ + j - 1].x - Vertices[i * nMaxJ + j].x;
		vA.y = Vertices[i * nMaxJ + j - 1].y - Vertices[i * nMaxJ + j].y;
		vA.z = Vertices[i * nMaxJ + j - 1].z - Vertices[i * nMaxJ + j].z;

		vB.x = Vertices[myMod(i + 1, n) * nMaxJ + j].x - Vertices[i * nMaxJ + j].x;
		vB.y = Vertices[myMod(i + 1, n) * nMaxJ + j].y - Vertices[i * nMaxJ + j].y;
		vB.z = Vertices[myMod(i + 1, n) * nMaxJ + j].z - Vertices[i * nMaxJ + j].z;

		vD.x = Vertices[myMod(i - 1, n) * nMaxJ + j].x - Vertices[i * nMaxJ + j].x;
		vD.y = Vertices[myMod(i - 1, n) * nMaxJ + j].y - Vertices[i * nMaxJ + j].y;
		vD.z = Vertices[myMod(i + 1, n) * nMaxJ + j].z - Vertices[i * nMaxJ + j].z;

		v0.x = (vA.x + vB.x) / 2;
		v0.y = (vA.y + vB.y) / 2;
		v0.z = (vA.z + vB.z) / 2;

		v1.x = (vA.x + vD.x) / 2;
		v1.y = (vA.y + vD.y) / 2;
		v1.z = (vA.z + vD.z) / 2;

		CrossProduct(&v0, &v1, &v2);

		V->x = Vertices[i * nMaxJ + j].x + v2.x * fFacteur - Vertices[i * nMaxJ + j].x;
		V->y = Vertices[i * nMaxJ + j].y + v2.y * fFacteur - Vertices[i * nMaxJ + j].y;
		V->z = Vertices[i * nMaxJ + j].z + v2.z * fFacteur - Vertices[i * nMaxJ + j].z;


	}

	if("None" && nNormale == ENHANCE_NONE)
	{
		V->x = 0;
		V->y = 0;
		V->z = 0;
	}

}
void Enhance(GLpoint* Vertices, int nMaxI, int nMaxJ)
{
	GLpoint V;
	int nMerge = myDome.nEnhanceMerge; //0 = aucun, 1 = 1, 2 = 2 (+)
	int nMode = myDome.nEnhanceCalcul; // 0 = +, 1 = *
	int nNormale = myDome.nEnhanceType;
	bool bMeridian = myDome.bEnhanceMeridian;
	bool bParallel = myDome.bEnhanceParallel;

	GLpoint* LocalVertices;
	if (!myDisplay.bDisplayMinor)
		return;
	LocalVertices  = (GLpoint*)malloc(nMaxI * nMaxJ * sizeof(GLpoint));
	if (LocalVertices == 0)
	{
		Message("     Memory Allocation LocalVertices Problem   ");
		return;
	}

	memcpy_s(LocalVertices, nMaxI * nMaxJ * sizeof(GLpoint), Vertices, nMaxI * nMaxJ * sizeof(GLpoint));

	if (bParallel)
	{
		for (int j = 1; j < myDome.nParallelMain; j++)
			for (int i = 0; i < nMaxI; i++)
			{
				if ((myMod(i, myDome.nMeridianMinor + 1) == 0) && bMeridian && (nMerge != MERGE_SOMME) )
					;
				else
				{
					ComputeVector(LocalVertices, nMaxI, nMaxJ, i, j * (myDome.nParallelMinor + 1), &V, nMode, nNormale);
					Vertices[i * nMaxJ + j * (myDome.nParallelMinor + 1)].x += V.x;
					Vertices[i * nMaxJ + j * (myDome.nParallelMinor + 1)].y += V.y;
					Vertices[i * nMaxJ + j * (myDome.nParallelMinor + 1)].z += V.z;
				}
			}
	}

	if (bMeridian)
	{
		for (int i = 0 ; i < myDome.nMeridianMain * myDome.nBaseSides + 1; i++)
			for (int j = 0; j < nMaxJ; j++)
			{
				if ((myMod(j, myDome.nParallelMinor + 1) == 0) && bParallel && (nMerge == MERGE_ZERO))
				;
				else
				{
					ComputeVector(LocalVertices, nMaxI, nMaxJ, i * (myDome.nMeridianMinor + 1), j, &V, nMode, nNormale);
					Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j].x += V.x;
					Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j].y += V.y;
					Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j].z += V.z;
				}

			}
	}
	free(LocalVertices);
}
void NormaleSharpening(GLpoint* Vertices, int nMaxI, int nMaxJ)
{
	bool bNormale = myDome.bNormale;
	float Lambda = myDome.fNormaleLambda;
	float Facteur;
	GLpoint* Normale;
	GLpoint Temp, U, V, W;
	int n_I1, n_I2;
	int nNombre;


	if (myDome.nDomeEnd != 360)
		nMaxI = (int)((float)(nMaxI * myDome.nDomeEnd) / 360.0) + 1;

	Facteur = sqrtf((float)(nMaxI * nMaxI + nMaxJ * nMaxJ)) / 4;
	if (bNormale && myDisplay.bDisplayMinor)
	{
		Normale = (GLpoint*)malloc(nMaxI * nMaxJ * sizeof(GLpoint));
		if (Normale == 0)
		{
			Message("     Memory Allocation Evolute Problem   ");
			return;
		}

		memcpy_s(Normale, nMaxI * nMaxJ * sizeof(GLpoint), Vertices, nMaxI * nMaxJ * sizeof(GLpoint));
		for (int i = 0; i < nMaxI; i++)
			for (int j = 1; j < nMaxJ - 1; j++)
			{
				Temp.x = 0.0;
				Temp.y = 0.0;
				Temp.z = 0.0;
				nNombre = 0;

				for (int k = 0; k < 8; k++)
				{
					int ai = (int)round(sin(_PI * (k + 2) / 4));
					int aj = (int)round(sin(_PI * (k + 4) / 4));

					n_I1 = (myDome.nDomeEnd == 360 ? (i + ai + nMaxI - 1) % (nMaxI - 1) : (i + ai < 0 || i + ai > nMaxI - 1 ? -1 : i + ai));
					//					if (i == 0 && j == 1)
					//						printf("i + %d, j +%d\n", ai, aj);
					if (n_I1 != -1)
					{
						U.x = Normale[n_I1 * nMaxJ + j + aj].x - Normale[i * nMaxJ + j].x;
						U.y = Normale[n_I1 * nMaxJ + j + aj].y - Normale[i * nMaxJ + j].y;
						U.z = Normale[n_I1 * nMaxJ + j + aj].z - Normale[i * nMaxJ + j].z;
					}

					ai = (int)round(sin(_PI * (k + 3) / 4));
					aj = (int)round(sin(_PI * (k + 5) / 4));
					n_I2 = (myDome.nDomeEnd == 360 ? (i + ai + nMaxI - 1) % (nMaxI - 1) : (i + ai < 0 || i + ai>nMaxI - 1 ? -1 : i + ai));
					//					if (i == 0 && j == 1)
					//						printf("i + %d, j +%d\n\n", ai, aj);
					if (n_I2 != -1)
					{
						V.x = Normale[n_I2 * nMaxJ + j + aj].x - Normale[i * nMaxJ + j].x;
						V.y = Normale[n_I2 * nMaxJ + j + aj].y - Normale[i * nMaxJ + j].y;
						V.z = Normale[n_I2 * nMaxJ + j + aj].z - Normale[i * nMaxJ + j].z;
					}

					if (n_I1 != -1 && n_I2 != -1)
					{
						CrossProduct(&U, &V, &W);

						Temp.x += W.x;
						Temp.y += W.y;
						Temp.z += W.z;
						nNombre++;
					}
				}


				Vertices[i * nMaxJ + j].x = Normale[i * nMaxJ + j].x + Facteur * Lambda * Temp.x / (nNombre == 0 ? 1 : nNombre);
				Vertices[i * nMaxJ + j].y = Normale[i * nMaxJ + j].y + Facteur * Lambda * Temp.y / (nNombre == 0 ? 1 : nNombre);
				Vertices[i * nMaxJ + j].z = Normale[i * nMaxJ + j].z + Facteur * Lambda * Temp.z / (nNombre == 0 ? 1 : nNombre);

			}
		free(Normale);
	}

}

bool PostValide(int iMain, int jMain)
{
	int nMeridianEvery = myDome.nPostMeridianEvery;
	int nMeridianShift = myDome.nPostMeridianShift;
	int nMeridianEvolve = myDome.nPostParallelEvolve;
	int nParallelFirst = myDome.nPostParallelFirst - 1;
	int nParallelEvery = myDome.nPostParallelEvery;
	int nMeridianPerSide = myDome.bPostMeridianPerSide;
	bool bPost = true;
	int ii, jj;

	//				printf("(%2d, %2d)", iMain, jMain);

	jj = (myDome.nParallelMain - 1) - jMain;
	if (jj < nParallelFirst)
		bPost = false;
	if (((jj - nParallelFirst) % nParallelEvery) != 0)
		bPost = false;

	ii = (nMeridianPerSide == 0 ? iMain : iMain % myDome.nMeridianMain) + nMeridianShift;
	nMeridianShift = 0;
//	if (ii < nMeridianDelay)
//		bPost = false;
	if ((ii - nMeridianShift - (nMeridianEvolve * (jj - nParallelFirst)) / nParallelEvery) % nMeridianEvery != 0)
		bPost = false;

	return bPost;
}
void BumpsOld(int nMaxI, int nMaxJ, GLpoint* Vertices, bool bDraw)
{
	float fLambda = 10 * myDome.fBumpLambda;
	float fDistanceMin = myDome.fBumpDistance;
	bool bBumps = (myDome.bBump && myDome.nMeridianMinor > 0 && myDome.nParallelMinor > 0);
	float fPower = myDome.fBumpPower;
	float fDamping = myDome.fBumpDamping;
	int nType = myDome.nBumpNormal;
	int nDistanceType = myDome.nBumpDistanceType;
	int bInverse = myDome.nBumpDistanceInverse;
	Vector3D v0, v1, v2, v3, w0, w1, w2, w3, V;
	float fDistance, fDistanceEdge;
	bool bTest;
	float dMap;

	if (bBumps && (myDisplay.bDisplayMinor || bDraw))
	{
		for (int i = 0; i < myDome.nMeridianMain * myDome.nBaseSides; i++)
			for (int j = 0; j < myDome.nParallelMain - 0; j++)
			{
				//calcul Vector
				v0.x = Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].x - Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].x;
				v0.y = Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].y - Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].y;
				v0.z = Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].z - Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].z;

				v1.x = Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].x - Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].x;
				v1.y = Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].y - Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].y;
				v1.z = Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].z - Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].z;

				CrossProduct(&v0, &v1, &w0);

				v2.x = Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].x - Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].x;
				v2.y = Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].y - Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].y;
				v2.z = Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nMeridianMinor + 1)].z - Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].z;

				v3.x = Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].x - Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].x;
				v3.y = Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].y - Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].y;
				v3.z = Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].z - Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].z;

//				printf("i = %d, j = %d, v0 = (%d) - (%d)     v1 = (%d) - (%d)\n", i, j, i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1), i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1), (i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1), i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1));

				CrossProduct(&v2, &v3, &w1);


				CrossProduct(&v1, &v2, &w2);
				CrossProduct(&v3, &v0, &w3);


				/*		printf("(%3d, %3d) - %3d, %3d, %3d, %3d\n",i, j, i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1),
													 i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1),
											   (i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1),
						(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1));*/

				switch (nType)
				{
				case BUMP_MEAN:
					V.x = (Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].x +
						Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].x +
						Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].x +
						Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].x) / (4 * myDome.nParallelMain);
					V.y = (Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].y +
						Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].y +
						Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].y +
						Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].y) / (4 * myDome.nParallelMain);
					V.z = (Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].z +
						Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].z +
						Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].z +
						Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].z) / (4 * myDome.nParallelMain);
					break;
				case BUMP_1VECTEUR:
					V.x = w0.x;
					V.y = w0.y;
					V.z = w0.z;
					break;
				case BUMP_2VECTEUR:
					V.x = (w0.x + w1.x) / 2;
					V.y = (w0.y + w1.y) / 2;
					V.z = (w0.z + w1.z) / 2;
					break;
					//				case 3:
					//					V.x = (w0.x + w1.x + w2.x + w3.x) / 4;
					//					V.y = (w0.y + w1.y + w2.y + w3.y) / 4;
					//					V.z = (w0.z + w1.z + w2.z + w3.z) / 4;
					//					break;
				default:
					V.x = (w0.x + w1.x) / 2;
					V.y = (w0.y + w1.y) / 2;
					V.z = (w0.z + w1.z) / 2;
					break;

				}
				if (myDome.nBumpNormalize)
					Normalize(&V);
				V.x /= myDome.nParallelMain;
				V.y /= myDome.nParallelMain;
				V.z /= myDome.nParallelMain;

				for (int kM = 1; kM <= myDome.nMeridianMinor; kM++)  // peut-être =0; <= xxx +1
					for (int kP = 1; kP <= myDome.nParallelMinor; kP++)
					{
						float x, y, fNormeMax;
						float x1, y1;

						if (bInverse)
						{
							x = (float)__min(kM - 1, myDome.nMeridianMinor - kM);
							y = (float)__min(kP - 1, myDome.nParallelMinor - kP);
							x1 = ceilf(myDome.nMeridianMinor / 2.0f) - 1.0f;
							y1 = ceilf(myDome.nParallelMinor / 2.0f) - 1.0f;
						}
						else
						{
							x = floorf(abs((float)(myDome.nMeridianMinor + 1) / 2.0f - (float)kM));
							y = floorf(abs((float)(myDome.nParallelMinor + 1) / 2.0f - (float)kP));
							x1 = floorf(abs((float)(myDome.nMeridianMinor + 1) / 2.0f - 1.0f));
							y1 = floorf(abs((float)(myDome.nParallelMinor + 1) / 2.0f - 1.0f));
						}

						switch (nDistanceType)
						{
						case DISTANCE_EUCLIDE:
						default:
							fDistance = sqrtf(x * x + y * y);
							fNormeMax = sqrtf(x1 * x1 + y1 * y1);
							break;
						case DISTANCE_MANHATTAN:
							fDistance = x + y;
							fNormeMax = x1 + y1;
							break;
						case DISTANCE_ULTRA:
							fDistance = __max(x, y);
							fNormeMax = __max(x1, y1);
							break;
						case DISTANCE_3:
							fDistance = powf(x * x * x + y * y * y, 1.0f / 3.0f);
							fNormeMax = powf(x1 * x1 * x1 + y1 * y1 * y1, 1.0f / 3.0f);
							break;
						case DISTANCE_DISCRETE:
							fDistance = (x == 0 ? (y == 0.0f ? 0 : 1.0f) : 1);
							fNormeMax = 1;
						}
						if (fNormeMax > 0)
							fDistance /= fNormeMax;
						if (fDistance < 0)
							printf("%f\n", fDistance);

						//						fNormeMax = (float)sqrt(ceil((float)myDome.nMeridianMinor / 2.0f)*ceil((float)myDome.nMeridianMinor / 2.0f) + ceil((float)myDome.nParallelMinor / 2.0)*ceil((float)myDome.nParallelMinor / 2.0));

						/*						fDistanceEdge = sqrt(x*x + y*y) / fNormeMax;
												fDistanceCenter = sqrt((ceil((float)myDome.nMeridianMinor / 2.0f) - x) *(ceil((float)myDome.nMeridianMinor / 2.0f) - x) +  (ceil((float)myDome.nParallelMinor / 2.0f) - y) *(ceil((float)myDome.nParallelMinor / 2.0f) - y));

												if (sqrt((ceil((float)myDome.nParallelMinor / 2.0f) - 1) * (ceil((float)myDome.nParallelMinor / 2.0f) - 1) + (ceil((float)myDome.nMeridianMinor / 2.0f) - 1) * (ceil((float)myDome.nMeridianMinor / 2.0f) - 1)))
													fDistanceCenter /= sqrt((ceil((float)myDome.nParallelMinor / 2.0f) - 1) * (ceil((float)myDome.nParallelMinor / 2.0f) - 1) + (ceil((float)myDome.nMeridianMinor / 2.0f) - 1) * (ceil((float)myDome.nMeridianMinor / 2.0f) - 1));
						*/


						if (bInverse)
							bTest = (fDistance > 1 - fDistanceMin);
						else
							bTest = (fDistance <= fDistanceMin);

						//						printf("[%2d, %2d] -- (%5.3f, %5.3f) -- %5.3f\n", kM, kP, x, y, fDistance);

						if (bTest && fDistanceMin > 0)
						{
							if (myDome.nBumpType == BUMP_FUNNEL)
							{
								x = (float)__min(kM - 1, myDome.nMeridianMinor - kM);
								y = (float)__min(kP - 1, myDome.nParallelMinor - kP);
								x1 = ceilf(myDome.nMeridianMinor / 2.0f) - 1.0f;
								y1 = ceilf(myDome.nParallelMinor / 2.0f) - 1.0f;

								fDistanceEdge = sqrtf(x * x + y * y);
								if (sqrtf(x1 * x1 + y1 * y1 > 0))
									fDistanceEdge /= sqrtf(x1 * x1 + y1 * y1);
							}
							else
							{
								x = floorf(abs((float)(myDome.nMeridianMinor + 1) / 2.0f - (float)kM));
								y = floorf(abs((float)(myDome.nParallelMinor + 1) / 2.0f - (float)kP));
								x1 = floorf(abs((float)(myDome.nMeridianMinor + 1) / 2.0f - 1.0f));
								y1 = floorf(abs((float)(myDome.nParallelMinor + 1) / 2.0f - 1.0f));
								switch (myDome.nBumpType)
								{
								case BUMP_WAVE:

									dMap = 1.0f + sinf(5.0f * (float)_PI * fDistance);


									if (dMap < 0)
										dMap = 0;
									else
										dMap = sqrtf(dMap);
									break;
								case BUMP_FOURPEAKS:
									x = 2 * x - x1; 
									y = 2 * y - y1;
									dMap = 1 - (x / x1) * (x / x1) - (y / y1) * (y / y1);

									if (dMap < 0)
										dMap = 0;
									else
										dMap = sqrtf(dMap);
									break;
								case BUMP_ASTROID:
									dMap = 1 - powf(x / x1, 2.0f / 3.0f) - powf(y / y1, 2.0f / 3.0f);
									if (dMap < 0)
										dMap = 0;
									else
										dMap = powf(dMap, 3 / 2);
									break;

								case BUMP_ELLIPSOID:
									dMap = 1 - (x / x1) * (x / x1) - (y / y1) * (y / y1);
									if (dMap < 0)
										dMap = 0;
									else
										dMap = sqrtf(dMap);
									break;
								case BUMP_TRACTOID:
									dMap = (acoshf((float)(1.0f / (fDistance < .01 ? .01 : fDistance))) - sqrtf(1.0f - fDistance * fDistance / (fNormeMax * fNormeMax)));
									if (dMap < 0)
										dMap = 0;

									break;
								case BUMP_CONE:
								default:
									dMap = 0;
									break;
								}
								//								dMap = 1 - (x / x1) * (x / x1) + (y / y1) * (y / y1);
								fDistanceEdge = dMap;

							}


							fDistanceEdge = powf(fDistanceEdge, fPower) * (1 - fDamping + j * fDamping / (float)(myDome.nParallelMain - 1));

							Vertices[(i * (myDome.nMeridianMinor + 1) + kM) * nMaxJ + j * (myDome.nParallelMinor + 1) + kP].x += (float)(fLambda * fDistanceEdge * V.x);
							Vertices[(i * (myDome.nMeridianMinor + 1) + kM) * nMaxJ + j * (myDome.nParallelMinor + 1) + kP].y += (float)(fLambda * fDistanceEdge * V.y);
							Vertices[(i * (myDome.nMeridianMinor + 1) + kM) * nMaxJ + j * (myDome.nParallelMinor + 1) + kP].z += (float)(fLambda * fDistanceEdge * V.z);

						}
					}
			}

	}
}
void Marquees(int nMaxI, int nMaxJ, GLpoint* Vertices, bool bDraw)
{
	bool bMarquee = myDome.bMarquee;
	float fMarqueeLambda = myDome.fMarqueeLambda;
	float Lambda = myDome.fMarqueeLambda;
	float fTop = myDome.fMarqueeTop;
	float fDamping = myDome.fMarqueeDamping;
	GLpoint* TempVertices;
	Vector3D v0, v1, V;
	int i, j;
	int nMain = myDome.nMeridianMain * (myDome.nMeridianMinor + 1) * myDome.nBaseSides;
	nMain *= (myDome.nParallelMain - 1) * (myDome.nParallelMinor + 1) + 1;
	TempVertices = (GLpoint*)malloc((nMaxI * nMaxJ) * sizeof(GLpoint));
	float fCoeffDamping;

//	printf("Taille totale = %d\n", nMaxI * nMaxJ);

	float fCoeff;
	nMaxJ = myDome.nParallelMain + (myDome.nParallelMain - 1) * myDome.nParallelMinor;
	if (bMarquee && (myDisplay.bDisplayMinor || bDraw))
	{

		// Minor struts

		for (int iMain = 0; iMain < myDome.nMeridianMain * myDome.nBaseSides; iMain++)
		{
			for (int jMain = 1; jMain < myDome.nParallelMain; jMain++)
			{
				fCoeffDamping = fDamping * jMain / (myDome.nParallelMain - 1) + 1 - fDamping;

				for (int iMinor = 0; iMinor < myDome.nMeridianMinor + 1; iMinor++)
				{
					for (int jMinor = 0; jMinor < myDome.nParallelMinor + 1; jMinor++)
					{
						i = iMain * (myDome.nMeridianMinor + 1) + iMinor;
						j = jMain * (myDome.nParallelMinor + 1) - jMinor;
						int ii, jj;

						if (i == 0)
							ii = (myDome.nMeridianMain * myDome.nBaseSides) * (myDome.nMeridianMinor + 1) - 1;
						else
							ii = i - 1;

						v0.x = Vertices[(i + 1) * nMaxJ + j].x - Vertices[ii * nMaxJ + j].x;
						v0.y = Vertices[(i + 1) * nMaxJ + j].y - Vertices[ii * nMaxJ + j].y;
						v0.z = Vertices[(i + 1) * nMaxJ + j].z - Vertices[ii * nMaxJ + j].z;

						//						printf("(%4d, %4d) = %4d/%4d pour (%2d, %2d) = %3d\n", i, j, i * nMaxJ + j, nMaxI * nMaxJ, iMinor, jMinor, iMinor * (myDome.nParallelMinor + 1) + jMinor);
						if (j == nMaxJ - 1)
							jj = j;
						else
							jj = j + 1;
						v1.x = Vertices[i * nMaxJ + (j - 1)].x - Vertices[i * nMaxJ + jj].x;
						v1.y = Vertices[i * nMaxJ + (j - 1)].y - Vertices[i * nMaxJ + jj].y;
						v1.z = Vertices[i * nMaxJ + (j - 1)].z - Vertices[i * nMaxJ + jj].z;



						CrossProduct(&v0, &v1, &V);


						if (j == nMaxJ - 1)
						{
							V.x *= 2;
							V.y *= 2;
							V.z *= 2;
						}
						//							printf("(%d/%d, %d/%d)  : (%f, %f, %f) ", i, nMaxI, j, nMaxJ, Vertices[i * nMaxJ + j].x, Vertices[i * nMaxJ + j].y, Vertices[i * nMaxJ + j].z);
						//							printf("Marquee = %f\n", sqrt(V.x * V.x + V.y * V.y + V.z * V.z));

						if (myDome.bMarqueeNormalize)
						{
							float fNorme = 10.0f*(float)sqrt(V.x * V.x + V.y * V.y + V.z * V.z);
							if (fNorme != 0)
							{
								V.x /= fNorme;
								V.y /= fNorme;
								V.z /= fNorme;
							}
						}

						if (jMain == 1 && iMain == 0 && iMinor == 0 && jMinor == 0)
							fCoeff = (float)(V.x * V.x + V.y * V.y + V.z * V.z);


						if (myDome.nSpiraleStep == -2)
						{
							printf("float %d\n", myDome.nSpiraleStep);
							TempVertices[i * nMaxJ + j].x = Vertices[i * nMaxJ + j].x + (float)(Lambda * V.x * (myDome.nMeridianMinor + myDome.nParallelMinor)) * MarqueeFunction((float)(iMinor - myDome.nMeridianMinor / 2), (float)(jMinor - myDome.nParallelMinor / 2));
							TempVertices[i * nMaxJ + j].y = Vertices[i * nMaxJ + j].y + (float)(Lambda * V.y * (myDome.nMeridianMinor + myDome.nParallelMinor)) * MarqueeFunction((float)(iMinor - myDome.nMeridianMinor / 2), (float)(jMinor - myDome.nParallelMinor / 2));
							TempVertices[i * nMaxJ + j].z = Vertices[i * nMaxJ + j].z + (float)(Lambda * V.z * (myDome.nMeridianMinor + myDome.nParallelMinor)) * MarqueeFunction((float)(iMinor - myDome.nMeridianMinor / 2), (float)(jMinor - myDome.nParallelMinor / 2));
						}
//							printf("int %d\n", myDome.nSpiraleStep);
						TempVertices[i * nMaxJ + j].x = Vertices[i * nMaxJ + j].x + (float)(Lambda * fCoeffDamping * V.x * (myDome.nMeridianMinor + myDome.nParallelMinor)) * MarqueeFunction(iMinor, jMinor);
						TempVertices[i * nMaxJ + j].y = Vertices[i * nMaxJ + j].y + (float)(Lambda * fCoeffDamping * V.y * (myDome.nMeridianMinor + myDome.nParallelMinor)) * MarqueeFunction(iMinor, jMinor);
						TempVertices[i * nMaxJ + j].z = Vertices[i * nMaxJ + j].z + (float)(Lambda * fCoeffDamping * V.z * (myDome.nMeridianMinor + myDome.nParallelMinor)) * MarqueeFunction(iMinor, jMinor);
						
					}

				}

			}
		}

		// cas j = 0

		fCoeff = fTop * sqrtf(fCoeff);
//		printf("%f\n",fCoeff);
		for (int i = 0; i < nMaxI; i++)
		{
			TempVertices[i * nMaxJ + 0].x = Vertices[i * nMaxJ + 0].x + 0.0f * (float)(Lambda * (1 - fDamping) * fCoeff * (myDome.nMeridianMinor + myDome.nParallelMinor));
			TempVertices[i * nMaxJ + 0].y = Vertices[i * nMaxJ + 0].y + 0.0f * (float)(Lambda * (1 - fDamping) * fCoeff * (myDome.nMeridianMinor + myDome.nParallelMinor));
			TempVertices[i * nMaxJ + 0].z = Vertices[i * nMaxJ + 0].z + (float)(Lambda * (1 - fDamping) * fCoeff * (myDome.nMeridianMinor + myDome.nParallelMinor));

		}
		 
		// dernier meridien
		for (int j = 0; j <= (myDome.nParallelMain - 1) * (myDome.nParallelMinor + 1); j++)
		{
			TempVertices[(myDome.nMeridianMain * myDome.nBaseSides) * (myDome.nMeridianMinor + 1) * nMaxJ + j].x = TempVertices[0 + j].x;
			TempVertices[(myDome.nMeridianMain * myDome.nBaseSides) * (myDome.nMeridianMinor + 1) * nMaxJ + j].y = TempVertices[0 + j].y;
			TempVertices[(myDome.nMeridianMain * myDome.nBaseSides) * (myDome.nMeridianMinor + 1) * nMaxJ + j].z = TempVertices[0 + j].z;
		}

		for (int i = 1000000; i < nMaxI * nMaxJ; i++)
		{
			Vertices[i ].x = TempVertices[i ].x;
			Vertices[i ].y = TempVertices[i ].y;
			Vertices[i ].z = TempVertices[i ].z;
		}

		bool bMarquee = true;
//		int ii, jj;
//		printf("nMeridianEvery = %d	nMeridianDelay = %d	nMeridianEvolve = %d	nParallelFirst = %d	nParallelEvery = %d\n", 
//			nMeridianEvery, nMeridianDelay, nMeridianEvolve, nParallelFirst + 1, nParallelEvery);
		for (int jMain = 1; jMain < myDome.nParallelMain; jMain++)
		{
			for (int iMain = 0; iMain < myDome.nMeridianMain * myDome.nBaseSides; iMain++)
			{

				bMarquee = PostValide(iMain, jMain);

				if (bMarquee)
				{
					for (int iMinor = 0; iMinor < myDome.nMeridianMinor + 2; iMinor++)
					{
						for (int jMinor = 0; jMinor < myDome.nParallelMinor + 1; jMinor++)
						{
							i = iMain * (myDome.nMeridianMinor + 1) + iMinor;
							j = jMain * (myDome.nParallelMinor + 1) - jMinor;
							Vertices[i * nMaxJ + j].x = TempVertices[i * nMaxJ + j].x;
							Vertices[i * nMaxJ + j].y = TempVertices[i * nMaxJ + j].y;
							Vertices[i * nMaxJ + j].z = TempVertices[i * nMaxJ + j].z;

						}

					}
				}
//				printf(" %s", (bMarquee ? "X " : "  "));
			}
//			printf("\n");
		}

	}

	free(TempVertices);

}
void Bumps(int nMaxI, int nMaxJ, GLpoint* Vertices, bool bDraw)
{
	float fLambda = 10*myDome.fBumpLambda;
	float fDistanceMin = myDome.fBumpDistance;
	bool bBumps = (myDome.bBump && myDome.nMeridianMinor > 0 && myDome.nParallelMinor > 0);
	float fPower = myDome.fBumpPower;
	float fDamping = myDome.fBumpDamping;
	int nType = myDome.nBumpNormal;
	int nDistanceType = myDome.nBumpDistanceType;
	int bInverse = myDome.nBumpDistanceInverse;
	Vector3D v0, v1, v2, v3, w0, w1, w2, w3, V;
	float fDistance, fDistanceEdge;
	bool bTest;
	float dMap;
	int nMeridianEvery = myDome.nPostMeridianEvery;
	int nMeridianShift = myDome.nPostMeridianShift % myDome.nMeridianMain;
	int nParallelEvolve = myDome.nPostParallelEvolve;
	int nParallelFirst = myDome.nPostParallelFirst - 1;
	int nParallelEvery = myDome.nPostParallelEvery;
	bool bBump;
	int nInverse; 
	if (bBumps && (myDisplay.bDisplayMinor || bDraw))
	{
		for (int i = 0; i < myDome.nMeridianMain * myDome.nBaseSides; i++)
		{
			for (int j = 0; j < myDome.nParallelMain - 1; j++)
			{
/*				nInverse = myDome.nParallelMain - j - 2;
				bBump = true;
				if (nInverse < nParallelFirst)
					bBump = false;
				bBump = bBump && ((nInverse - nParallelFirst) % nParallelEvery == 0);

//				bBump = bBump && (((i % myDome.nMeridianMain) % nMeridianEvery) == ((nMeridianDelay + nMeridianEvolve) % myDome.nMeridianMain));
				bBump = bBump && (((i % myDome.nMeridianMain) % nMeridianEvery) - nMeridianDelay+100* nMeridianEvery) % nMeridianEvery == ((nInverse * nParallelEvolve) % myDome.nMeridianMain) % nMeridianEvery;
*/
				bBump = PostValide(i, j + 1);

				if (i < myDome.nMeridianMain && false)
					printf("(%d, %d)-->%s, (%d, %d)\n", i, nInverse, (bBump ? "true" : "false"), (((i % myDome.nMeridianMain) % nMeridianEvery) - nMeridianShift) % nMeridianEvery, ((nInverse * nParallelEvolve) % myDome.nMeridianMain) % nMeridianEvery);

				if (bBump)
				{	//calcul Vector
					v0.x = Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].x - Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].x;
					v0.y = Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].y - Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].y;
					v0.z = Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].z - Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].z;

					v1.x = Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].x - Vertices[(i) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].x;
					v1.y = Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].y - Vertices[(i) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].y;
					v1.z = Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].z - Vertices[(i) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].z;

					CrossProduct(&v0, &v1, &w0);

					v2.x = Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].x - Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].x;
					v2.y = Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].y - Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].y;
					v2.z = Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nMeridianMinor + 1)].z - Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].z;

					v3.x = Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 0) * (myDome.nParallelMinor + 1)].x - Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 0) * (myDome.nParallelMinor + 1)].x;
					v3.y = Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 0) * (myDome.nParallelMinor + 1)].y - Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 0) * (myDome.nParallelMinor + 1)].y;
					v3.z = Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 0) * (myDome.nParallelMinor + 1)].z - Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 0) * (myDome.nParallelMinor + 1)].z;

					//				printf("i = %d, j = %d, v0 = (%d) - (%d)     v1 = (%d) - (%d)\n", i, j, i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1), i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1), (i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1), i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1));
					//				printf("i = %d, j = %d, v0 =( %f, %f, %f)   v1 = (%f, %f, %f)\n", i, j, v0.x, v0.y, v0.z, v1.x, v1.y, v1.z);
					/*				printf("i = %d, j = %d, v2 = (%d) - (%d)     v3 = (%d) - (%d)\n", i, j, (i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1),
																											(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1),
																											i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 0) * (myDome.nParallelMinor + 1),
																											(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 0) * (myDome.nParallelMinor + 1));
					*/
					CrossProduct(&v3, &v2, &w1);


					CrossProduct(&v1, &v2, &w2);
					CrossProduct(&v3, &v0, &w3);


					/*		printf("(%3d, %3d) - %3d, %3d, %3d, %3d\n",i, j, i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1),
														 i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1),
												   (i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1),
							(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1));*/

					switch (nType)
					{
					case BUMP_MEAN:
						V.x = (Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].x +
							Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].x +
							Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].x +
							Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].x) / (4 * myDome.nParallelMain);
						V.y = (Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].y +
							Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].y +
							Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].y +
							Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].y) / (4 * myDome.nParallelMain);
						V.z = (Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].z +
							Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + j * (myDome.nParallelMinor + 1)].z +
							Vertices[i * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].z +
							Vertices[(i + 1) * (myDome.nMeridianMinor + 1) * nMaxJ + (j + 1) * (myDome.nParallelMinor + 1)].z) / (4 * myDome.nParallelMain);
						break;
					case BUMP_1VECTEUR:
						V.x = w0.x;
						V.y = w0.y;
						V.z = w0.z;
						break;
					case BUMP_2VECTEUR:
						V.x = (w0.x + w1.x) / 2;
						V.y = (w0.y + w1.y) / 2;
						V.z = (w0.z + w1.z) / 2;
						break;
						//				case 3:
						//					V.x = (w0.x + w1.x + w2.x + w3.x) / 4;
						//					V.y = (w0.y + w1.y + w2.y + w3.y) / 4;
						//					V.z = (w0.z + w1.z + w2.z + w3.z) / 4;
						//					break;
					default:
						V.x = (w0.x + w1.x) / 2;
						V.y = (w0.y + w1.y) / 2;
						V.z = (w0.z + w1.z) / 2;
						break;

					}
					if (myDome.nBumpNormalize)
						Normalize(&V);
					V.x /= myDome.nParallelMain;
					V.y /= myDome.nParallelMain;
					V.z /= myDome.nParallelMain;

					for (int kM = 1; kM <= myDome.nMeridianMinor; kM++)  // peut-être =0; <= xxx +1
						for (int kP = 1; kP <= myDome.nParallelMinor; kP++)
						{
							float x, y, fNormeMax;
							float x1, y1;

							if (bInverse)
							{
								x = (float)__min(kM - 1, myDome.nMeridianMinor - kM);
								y = (float)__min(kP - 1, myDome.nParallelMinor - kP);
								x1 = ceilf(myDome.nMeridianMinor / 2.0f) - 1.0f;
								y1 = ceilf(myDome.nParallelMinor / 2.0f) - 1.0f;
							}
							else
							{
								x = floorf(abs((float)(myDome.nMeridianMinor + 1) / 2.0f - (float)kM));
								y = floorf(abs((float)(myDome.nParallelMinor + 1) / 2.0f - (float)kP));
								x1 = floorf(abs((float)(myDome.nMeridianMinor + 1) / 2.0f - 1.0f));
								y1 = floorf(abs((float)(myDome.nParallelMinor + 1) / 2.0f - 1.0f));
							}

							switch (nDistanceType)
							{
							case DISTANCE_EUCLIDE:
							default:
								fDistance = sqrtf(x * x + y * y);
								fNormeMax = sqrtf(x1 * x1 + y1 * y1);
								break;
							case DISTANCE_MANHATTAN:
								fDistance = x + y;
								fNormeMax = x1 + y1;
								break;
							case DISTANCE_ULTRA:
								fDistance = __max(x, y);
								fNormeMax = __max(x1, y1);
								break;
							case DISTANCE_3:
								fDistance = powf(x * x * x + y * y * y, 1.0f / 3.0f);
								fNormeMax = powf(x1 * x1 * x1 + y1 * y1 * y1, 1.0f / 3.0f);
								break;
							case DISTANCE_DISCRETE:
								fDistance = (x == 0 ? (y == 0.0f ? 0 : 1.0f) : 1);
								fNormeMax = 1;
							}
							if (fNormeMax > 0)
								fDistance /= fNormeMax;
							if (fDistance < 0)
								printf("%f\n", fDistance);

							//						fNormeMax = (float)sqrt(ceil((float)myDome.nMeridianMinor / 2.0f)*ceil((float)myDome.nMeridianMinor / 2.0f) + ceil((float)myDome.nParallelMinor / 2.0)*ceil((float)myDome.nParallelMinor / 2.0));

							/*						fDistanceEdge = sqrt(x*x + y*y) / fNormeMax;
													fDistanceCenter = sqrt((ceil((float)myDome.nMeridianMinor / 2.0f) - x) *(ceil((float)myDome.nMeridianMinor / 2.0f) - x) +  (ceil((float)myDome.nParallelMinor / 2.0f) - y) *(ceil((float)myDome.nParallelMinor / 2.0f) - y));

													if (sqrt((ceil((float)myDome.nParallelMinor / 2.0f) - 1) * (ceil((float)myDome.nParallelMinor / 2.0f) - 1) + (ceil((float)myDome.nMeridianMinor / 2.0f) - 1) * (ceil((float)myDome.nMeridianMinor / 2.0f) - 1)))
														fDistanceCenter /= sqrt((ceil((float)myDome.nParallelMinor / 2.0f) - 1) * (ceil((float)myDome.nParallelMinor / 2.0f) - 1) + (ceil((float)myDome.nMeridianMinor / 2.0f) - 1) * (ceil((float)myDome.nMeridianMinor / 2.0f) - 1));
							*/


							if (bInverse)
								bTest = (fDistance > 1 - fDistanceMin);
							else
								bTest = (fDistance <= fDistanceMin);

							//						printf("[%2d, %2d] -- (%5.3f, %5.3f) -- %5.3f\n", kM, kP, x, y, fDistance);

							if (bTest && fDistanceMin > 0)
							{
								if (myDome.nBumpType == BUMP_FUNNEL)
								{
									x = (float)__min(kM - 1, myDome.nMeridianMinor - kM);
									y = (float)__min(kP - 1, myDome.nParallelMinor - kP);
									x1 = ceilf(myDome.nMeridianMinor / 2.0f) - 1.0f;
									y1 = ceilf(myDome.nParallelMinor / 2.0f) - 1.0f;

									fDistanceEdge = fDistance;
								}
								else
								{
									x = floorf(abs((float)(myDome.nMeridianMinor + 1) / 2.0f - (float)kM));
									y = floorf(abs((float)(myDome.nParallelMinor + 1) / 2.0f - (float)kP));
									x1 = floorf(abs((float)(myDome.nMeridianMinor + 1) / 2.0f - 1.0f));
									y1 = floorf(abs((float)(myDome.nParallelMinor + 1) / 2.0f - 1.0f));
									switch (myDome.nBumpType)
									{
									case BUMP_WAVE:

										dMap = 1.0f + sinf(5.0f * (float)_PI * fDistance);


										if (dMap < 0)
											dMap = 0;
										else
											dMap = sqrtf(dMap);
										break;
									case BUMP_FOURPEAKS:
										x = 2 * x - x1; //cross 2
										y = 2 * y - y1;
										dMap = 1 - (x / x1) * (x / x1) - (y / y1) * (y / y1);

										dMap = 1 - Distance((int)x, (int)y, (DISTANCE_Type)nDistanceType) / fNormeMax;
										if (dMap < 0)
											dMap = 0;
										else
											dMap = sqrtf(dMap);
										break;
									case BUMP_ASTROID:
										dMap = 1 - powf(x / x1, 2.0f / 3.0f) - powf(y / y1, 2.0f / 3.0f);
										dMap = 1 - powf(Distance((int)x, (int)y, (DISTANCE_Type)nDistanceType) / fNormeMax, 2.0f / 3.0f);

										if (dMap < 0)
											dMap = 0;
										else
											dMap = powf(dMap, 3.0f); // ou 3/2


										break;

									case BUMP_ELLIPSOID:
										dMap = 1 - (x / x1) * (x / x1) - (y / y1) * (y / y1);
										//									dMap = sqrtf((x / x1) * (x / x1) - (y / y1) * (y / y1)) - (Distance((int)x, (int)y, (DISTANCE_Type)nDistanceType) / fNormeMax) * (Distance((int)x, (int)y, (DISTANCE_Type)nDistanceType) / fNormeMax);
										dMap = 1 - (Distance((int)x, (int)y, (DISTANCE_Type)nDistanceType) / fNormeMax) * (Distance((int)x, (int)y, (DISTANCE_Type)nDistanceType) / fNormeMax);

										//									printf("%f\n", Distance((int)x, (int)y, (DISTANCE_Type)nDistanceType) / fNormeMax);
										if (dMap < 0)
											dMap = 0;

										//									else
										//										dMap = sqrtf(dMap);
										break;
									case BUMP_TRACTOID:
										dMap = (acoshf((float)(1.0f / (fDistance < .01 ? .01 : fDistance))) - sqrtf(1.0f - fDistance * fDistance / (fNormeMax * fNormeMax)));
										if (dMap < 0)
											dMap = 0;

										break;
									case BUMP_CONE:
										dMap = 1 - fDistance;
										break;
									default:
										dMap = 0;
										break;
									}
									fDistanceEdge = dMap;

								}


								fDistanceEdge = powf(fDistanceEdge, fPower) * (1 - fDamping + j * fDamping / (float)(myDome.nParallelMain - 1));

								Vertices[(i * (myDome.nMeridianMinor + 1) + kM) * nMaxJ + j * (myDome.nParallelMinor + 1) + kP].x += (float)(fLambda * fDistanceEdge * V.x);
								Vertices[(i * (myDome.nMeridianMinor + 1) + kM) * nMaxJ + j * (myDome.nParallelMinor + 1) + kP].y += (float)(fLambda * fDistanceEdge * V.y);
								Vertices[(i * (myDome.nMeridianMinor + 1) + kM) * nMaxJ + j * (myDome.nParallelMinor + 1) + kP].z += (float)(fLambda * fDistanceEdge * V.z);

							}
						}
				}
			}
		}
	}
}

