#include <string.h>
#include <stdio.h>
#include "stdafx.h"
#include <math.h>
#include <time.h>
#include <direct.h>

#include <GL\glui.h>
//#include "vc_01.h"
//#include "Functions.h"
#include "Globales.h"
#include "filebrowser.h"
#include "Ear.h"
//#define GLUI_FileBrowser FileBrowser

double MinU = 100, MaxU = -100, MinV = 100, MaxV = -100;
GLUI_Spinner* spCrenelPerSide;
GLUI_Spinner* spCrenelWide;
GLUI_Spinner* spCrenelDelay;
GLUI_Spinner* spWindowLower;
GLUI_Checkbox* ckWindowBlind;

float MaxRadius(char **sMaterial, SPHERE_Type Type)
{
	float fRadius;
	
	switch(Type)
	{
	case SPHERE_TOP:
		if (myDome.fSphereMPRadius >= myDome.fSphereMpRadius)
		{
			fRadius = myDome.fSphereMPRadius;
			*sMaterial = myDome.sSphereMPMaterial;
		}
		else
		{
			fRadius = myDome.fSphereMpRadius;
			*sMaterial = myDome.sSphereMpMaterial;
		}
		if (myDome.fSpheremPRadius > fRadius)
		{
			fRadius = myDome.fSpheremPRadius;
			*sMaterial = myDome.sSpheremPMaterial;
		}
		if (myDome.fSpherempRadius > fRadius)
		{
			fRadius = myDome.fSpherempRadius;
			*sMaterial = myDome.sSpherempMaterial;
		}
		if (myDome.fDiagonalDRRadius > fRadius)
		{
			fRadius = myDome.fDiagonalDRRadius;
			*sMaterial = myDome.sDiagonalDRMaterial;
		}
		if (myDome.fDiagonalURRadius > fRadius)
		{
			fRadius = myDome.fDiagonalURRadius;
			*sMaterial = myDome.sDiagonalURMaterial;
		}
		break;
	case SPHERE_MP:
	case SPHERE_Mp:
	case SPHERE_mP:
	case SPHERE_mp:
		if (Type == SPHERE_MP && myDome.fSphereMPRadius > 0)
		{
			fRadius = myDome.fSphereMPRadius;
			*sMaterial = myDome.sSphereMPMaterial;
		}
		else if (Type == SPHERE_Mp && myDome.fSphereMpRadius > 0)
		{
			fRadius = myDome.fSphereMpRadius;
			*sMaterial = myDome.sSphereMpMaterial;
		}
		else if (Type == SPHERE_mP && myDome.fSpheremPRadius > 0)
		{
			fRadius = myDome.fSpheremPRadius;
			*sMaterial = myDome.sSpheremPMaterial;
		}
		else if (Type == SPHERE_mp && myDome.fSpherempRadius > 0)
		{
			fRadius = myDome.fSpherempRadius;
			*sMaterial = myDome.sSpherempMaterial;
		}
		else if (myDome.fDiagonalDRRadius > myDome.fDiagonalURRadius)
		{
			fRadius = myDome.fDiagonalDRRadius;
			*sMaterial = myDome.sDiagonalDRMaterial;
		}
		else
		{
			fRadius = myDome.fDiagonalURRadius;
			*sMaterial = myDome.sDiagonalURMaterial;
		}
		break;
	}
	return fRadius;
}


void ProcessMessage(int Control)
{

	glMessage->hide();
}

void Message(const char *sMessage, int Type)
{
	stMessage->set_text(sMessage);

	glMessage->show();
}

void CrossProduct(GLpoint *V, GLpoint *W, GLpoint *Resultat)
{
	Resultat->x = (V->y * W->z - V->z * W->y);
	Resultat->y = (V->z * W->x - V->x * W->z);
	Resultat->z = (V->x * W->y - V->y * W->x);
}

void CrossProduct(Vector3D *V, Vector3D *W, Vector3D *Resultat)
{
	Resultat->x = (V->y * W->z - V->z * W->y);
	Resultat->y = (V->z * W->x - V->x * W->z);
	Resultat->z = (V->x * W->y - V->y * W->x);
}

float Colineaire(GLpoint O, GLpoint M, GLpoint N)
{
	//	return (b*z-c*y)*(b*z-c*y) + (c*x-a*z)*(c*x-a*z) + (a*y-b*x)*(a*y-b*x);
	return ((M.y - O.y)*(N.z - O.z) - (M.z - O.z)*(N.y - O.y)) * ((M.y - O.y)*(N.z - O.z) - (M.z - O.z)*(N.y - O.y)) +
		((M.z - O.z)*(N.x - O.x) - (M.x - O.x)*(N.z - O.z)) * ((M.z - O.z)*(N.x - O.x) - (M.x - O.x)*(N.z - O.z)) +
		((M.x - O.x)*(N.y - O.y) - (M.y - O.y)*(N.x - O.x)) * ((M.x - O.x)*(N.y - O.y) - (M.y - O.y)*(N.x - O.x));
}

float GetSegmentLength(float v0, float v1, float v2)
{
	//  Compute the length of a small section of a parametric curve
	//  from t0 to t2, with t1 being the mid point.
	//  The 3 points are precomputed.
	float epsilon  = 0.00001f;
	float epsilon2 = 0.000001f;
	float kMaxArc = 1.05f;
	float kLenRatio = 1.2f;

	GLpoint Point[3];
	Point[0].x = (float)SlopeFunction(v0, myDome.fSlopeParameterXY, SLOPE_XY);
	Point[0].y = (float)SlopeFunction(v0, myDome.fSlopeParameterZ, SLOPE_Z);
	Point[1].x = (float)SlopeFunction(v1, myDome.fSlopeParameterXY, SLOPE_XY);
	Point[1].y = (float)SlopeFunction(v1, myDome.fSlopeParameterZ, SLOPE_Z);
	Point[2].x = (float)SlopeFunction(v2, myDome.fSlopeParameterXY, SLOPE_XY);
	Point[2].y = (float)SlopeFunction(v2, myDome.fSlopeParameterZ, SLOPE_Z);

	float d1 = sqrtf((Point[0].x - Point[2].x) * (Point[0].x - Point[2].x) + (Point[0].y - Point[2].y) * (Point[0].y - Point[2].y));
	float da = sqrtf((Point[0].x - Point[1].x) * (Point[0].x - Point[1].x) + (Point[0].y - Point[1].y) * (Point[0].y - Point[1].y));
	float db = sqrtf((Point[1].x - Point[2].x) * (Point[1].x - Point[2].x) + (Point[1].y - Point[2].y) * (Point[1].y - Point[2].y));
	float d2 = da + db;
		//  if we’re in a region of high curvature, recurse, otherwise return the
		//  length as ( d2 + ( d2 – d1 ) / 3 )
	if (d2 < epsilon)
		return (d2 + (d2 - d1) / 3);
	else
		if (d1 < epsilon  || d2 / d1 > kMaxArc ||
			da < epsilon2 || db / da > kLenRatio ||
			db < epsilon2 || da / db > kLenRatio)
		{
			//  Recurse
			/*
			epsilon and epsilon2 in the above expression are there to guard against division
			by zero and underflow. The value of epsilon2 should be less than half that of
			epsilon otherwise unnecessary recursion occurs : values of 1e-5 and 1e-6 work
			satisfactorily.
			kMaxArc implicitly refers to the maximum allowable angle that can be
			subtended by a circular arc : a value of 1.05 gives good results in practice.
			kLenRatio refers to the maximum allowable ratio in the distance between
			successive data points : a value of 1.2 gives reasonable results.
			*/
			//GetPoint((t0 + t1) / 2, mid_pt)
			float d3 = GetSegmentLength(v0, (v0 + v1) / 2, v1);
			//GetPoint((t1 + t2) / 2, mid_pt)
			float d4 = GetSegmentLength(v1, (v1 + v2) / 2, v2);
			return d3 + d4;
		}
		else
			return (d2 + (d2 - d1) / 3);
}

float CurveLength()
{
//P: array of n_pts
//T : corresponding array of t values
	float Length = 0;

	int nbPoints = 1001;

	for (int i = 0; i < nbPoints - 2; i += 2)
		Length += GetSegmentLength((float)i / (float)(nbPoints - 1), (float)(i + 1) / (float)(nbPoints - 1), (float)(i + 2) / (float)(nbPoints - 1));

	return Length;
}

double CalculSegment(double v0, double v1)
{
	double x0, x1, y0, y1, xm, ym, z;
	double da, db, d;
	double epsilon = 0.0000001;
	Calcul(v0, 1, &x0, &y0, &z, false);
	Calcul(v1, 1, &x1, &y1, &z, false);
	Calcul((v0 + v1) / 2, 1, &xm, &ym, &z, false);

	da = sqrt((x0 - xm) * (x0 - xm) + (y0 - ym) * (y0 - ym));
	db = sqrt((xm - x1) * (xm - x1) + (ym - y1) * (ym - y1));
	d = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));

	if (da + db - d < epsilon)
		return d;
	else
		return CalculSegment(v0, (v0 + v1) / 2) + CalculSegment((v0 + v1) / 2, v1);
}

float ComputeBaseLength()
{
	double dLength = 0;
	double x0, y0, x1, y1, z;
	int nMax;

//	(myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor)* myDome.nBaseSides + 1
	nMax = myDome.nMeridianMain * (30*myDome.nMeridianMinor + 1) * myDome.nBaseSides ;
	;
	Calcul(0, 1, &x0, &y0, &z, false);
	for (int i = 0; i < nMax; i++)
	{
		Calcul(2.0 * _PI * (double)(i + 1) / (double)nMax, 1, &x1, &y1, &z, false);
		dLength += sqrt((x0-x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
		x0 = x1;
		y0 = y1;
	}
/*	double dLength2 = 0;
	nMax = myDome.nMeridianMain * (myDome.nMeridianMinor + 1) * myDome.nBaseSides -1;
	for (int i = 0; i < nMax; i++)
	{
		dLength2 += CalculSegment(2.0 * _PI * (double)(i) / (double)nMax, 2.0 * _PI * (double)(i + 1) / (double)nMax);
	}
*/
	return (float)dLength;
}

float ComputeSlopeLength()
{
	float x0, y0;
	float x1, y1;

//	ChargeSlope(myDome.nSlope);
	x0 = (float)SlopeFunction(0, myDome.fSlopeParameterXY, SLOPE_XY);
	y0 = (float)SlopeFunction(0, myDome.fSlopeParameterZ, SLOPE_Z);
	fSlopeLength = 0.0;



	for (float v = (float)0.0001; v < (float)1.00001; v += (float)0.0001)
	{
		if (v > 1)
			v = 1;
		x1 = (float)SlopeFunction(v, myDome.fSlopeParameterXY, SLOPE_XY);
		y1 = (float)SlopeFunction(v, myDome.fSlopeParameterZ, SLOPE_Z);
		fSlopeLength += sqrtf((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		x0 = x1;
		y0 = y1;
	}
/*	float Length = CurveLength();
	printf("Slope = %s, Curve  = %lf\n", SlopeFunctionListe[myDome.nSlope].sNom, Length);*/
	x0 = (float)SlopeFunction(0, myDome.fSlopeParameterXY, SLOPE_XY);
	y0 = (float)SlopeFunction(0, myDome.fSlopeParameterZ, SLOPE_Z);
	fSlopeLength = 0.0;
	for (float v = (float)0.00001; v < (float)1.000001; v += (float)0.00001)
	{
		if (v > 1)
			v = 1;
		x1 = (float)SlopeFunction(v, myDome.fSlopeParameterXY, SLOPE_XY);
		y1 = (float)SlopeFunction(v, myDome.fSlopeParameterZ, SLOPE_Z);
		fSlopeLength += sqrtf((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		x0 = x1;
		y0 = y1;
	}

	return fSlopeLength;
}

void BaseDefault()
{
	myDome.nBase = BASE_CIRCLE;
	myDome.nBaseSecondary = BASE_CONVEX_POLYGON;
	myDome.nBaseSides = 5;

	myDome.fBaseParameter = 1.0;
	myDome.fBaseSmoothing = 0.0;
	myDome.fBaseExpand = 1.0;
	myDome.fBaseDamping = 1.0;
	myDome.fBaseRound = 0.0;
	myDome.fBaseBaryCentre = 0.0;
	myDome.fBaseMorphing = 0.0;
	myDome.fBaseSpiral = 0.0;
	myDome.nBaseEvenly = RB_ALONG_PARAMETER;
	myDome.fBaseSecondaryRotation = 0.0;
	myDome.fBaseSecondaryExpansion = 1.0;
	myDome.bBaseMorphingDivision = 0;

}
void SlopeDefault()
{
	myDome.nSlope = SLOPE_HALF_CIRCLE;
	myDome.nSlopeTwist = 0;
	myDome.fSlopeDelay = 0;
	myDome.nSlopeEvenly = 0;
	myDome.fSlopeParameterZ = 2.0;
	myDome.fSlopeParameterXY = 1.0;
	myDome.bSlopeInverse = false;
	myDome.nCrenelPerSide = 0;
	myDome.fCrenelFront = 0.1f;
	myDome.fCrenelThick = 0.2f;
	myDome.fCrenelWide = 0.4f;
	myDome.fCrenelHeight = 0.5f;
	myDome.fCrenelPower = 0.30f;
	myDome.fCrenelAngle = 0.0f;
	myDome.fCrenelDelay = 0.0f;

}
void WindowDefault()
{
	myDome.nWindow = WINDOW_NONE;
	myDome.fWindowPower = 1.0;
	myDome.fWindowHeight = 0.5;
	myDome.fWindowDamping = 0.0;
	myDome.fWindowPerSide = 1.0;
	myDome.fWindowShift = 0;
	myDome.bWindowInverse = 0;
	myDome.fWindowDelay = 0.0;
	myDome.fWindowAttack = 0.0;
	myDome.fWindowAwning = 0.0;
	myDome.fWindowParameter = 0.0;
	myDome.bWindowBlind = 0;
	myDome.fWindowLower = 0.0;
}
void SpikeDefault()
{
	myDome.fSpikeMP = 0.0f;
	myDome.fSpikeMp = 0.0f;
	myDome.fSpikemP = 0.0f;
	myDome.fSpikemp = 0.0f;
	myDome.fSpikeAngle = 1.0f;
}
void LissageDefault()
{
	myDome.bLissage = false;
	myDome.fLissagePoidsSelf = 1.0;
	myDome.fLissagePoidsHorizontal = 0.0;
	myDome.fLissagePoidsVertical = 0.0;
	myDome.fLissagePoidsDiagonales = 0.0;
	myDome.fLissagePuissance = 1.0;
	myDome.nLissageDistanceType = DISTANCE_EUCLIDE;
	myDome.nLissageNombre = 1;
	myDome.nLissageTaille = 1;
	myDome.nLissageTypeFiltre = FILTER_Gauss;
	myDome.nLissageDirection = DIRECTION_NORTH;
	myDome.fLissageSigma1 = 1.0;
	myDome.fLissageSigma2 = 1.0;
}
void NoiseDefault()
{
	myDome.nNoiseOctave = 4;
	myDome.fNoiseDensity = 4.0;
	myDome.nNoiseSeed = 0;
	myDome.fNoiseFrequencyU = 1.0;
	myDome.fNoiseFrequencyV = 1.0;
	myDome.fNoiseAmplitudeR = 0.0;
	myDome.fNoiseAmplitudeZ = 0.0;
	myDome.fNoiseProgressiveR = 0.0;
	myDome.fNoiseProgressiveZ = 0.0;

}

float ComputeWindowCoefficient(float u, float a, float d)
{
	float x;
	float fCoeff = 1;

	if (myDome.bWindowInverse == 0)
		x = fmodf((float)(200 * (float)_PI) + u * a * d, (float)2 * (float)_PI);
	else
		x = fmodf((float)(200 * (float)_PI) + u * a * d + (float)_PI, (float)2 * (float)_PI);

	if (x <= myDome.fWindowDelay)
		fCoeff = 0;
	if ((myDome.fWindowDelay < x) && (x < myDome.fWindowDelay + myDome.fWindowAttack))
		fCoeff = (float)((1 - cos(_PI*(x - myDome.fWindowDelay) / myDome.fWindowAttack)) / 2);
	if ((myDome.fWindowDelay + myDome.fWindowAttack <= x) && (x <= 2 * _PI - myDome.fWindowDelay - myDome.fWindowAttack))
		fCoeff = 1;
	if ((2 * _PI - myDome.fWindowDelay - myDome.fWindowAttack < x) && (x < 2 * _PI - myDome.fWindowDelay))
		fCoeff = (float)((1 - cos(_PI*(2 * _PI - myDome.fWindowDelay - x) / myDome.fWindowAttack)) / 2);
	if (2 * _PI - myDome.fWindowDelay <= x)
		fCoeff = 0;
	if (myDisplay.nDisplay == DISPLAY_Slope)
		fCoeff = 0;
		
	return fCoeff;
}

/* Noise Perlin */

void PerlinNoiseNormalize2(float v[2])
{
	float s;

	s = (float)sqrt(v[0] * v[0] + v[1] * v[1]);
	s = 1.0f / s;
	v[0] = v[0] * s;
	v[1] = v[1] * s;
}
void PerlinNoiseInit(void)
{
	int i, j, k;

	for (i = 0; i < SAMPLE_SIZE; i++)
	{
		PerlinNoise_p[i] = i;
		PerlinNoise_g1[i] = (float)((rand() % (SAMPLE_SIZE + SAMPLE_SIZE)) - SAMPLE_SIZE) / SAMPLE_SIZE;
		for (j = 0; j < 2; j++)
			PerlinNoise_g2[i][j] = (float)((rand() % (SAMPLE_SIZE + SAMPLE_SIZE)) - SAMPLE_SIZE) / SAMPLE_SIZE;
		PerlinNoiseNormalize2(PerlinNoise_g2[i]);

		/*		for (j = 0 ; j < 3 ; j++)
		PerlinNoise_g3[i][j] = (float)((rand() % (SAMPLE_SIZE + SAMPLE_SIZE)) - SAMPLE_SIZE) / SAMPLE_SIZE;
		normalize3(PerlinNoise_g3[i]);*/
	}

	while (--i)
	{
		k = PerlinNoise_p[i];
		PerlinNoise_p[i] = PerlinNoise_p[j = rand() % SAMPLE_SIZE];
		PerlinNoise_p[j] = k;
	}

	for (i = 0; i < SAMPLE_SIZE + 2; i++)
	{
		PerlinNoise_p[SAMPLE_SIZE + i] = PerlinNoise_p[i];
		PerlinNoise_g1[SAMPLE_SIZE + i] = PerlinNoise_g1[i];
		for (j = 0; j < 2; j++)
			PerlinNoise_g2[SAMPLE_SIZE + i][j] = PerlinNoise_g2[i][j];
		for (j = 0; j < 3; j++)
			PerlinNoise_g3[SAMPLE_SIZE + i][j] = PerlinNoise_g3[i][j];
	}

}
float PerlinNoise2(float vec[2])
{
	int bx0, bx1, by0, by1, b00, b10, b01, b11;
	float rx0, rx1, ry0, ry1, *q, sx, sy, a, b, t, u, v;
	int i, j;

	if (myNoise.bStart)
	{
		srand(myNoise.nSeed);
		myNoise.bStart = false;
		PerlinNoiseInit();
	}

	PerlinNoiseSetup(0, bx0, bx1, rx0, rx1);
	PerlinNoiseSetup(1, by0, by1, ry0, ry1);

	i = PerlinNoise_p[bx0];
	j = PerlinNoise_p[bx1];

	b00 = PerlinNoise_p[i + by0];
	b10 = PerlinNoise_p[j + by0];
	b01 = PerlinNoise_p[i + by1];
	b11 = PerlinNoise_p[j + by1];

	sx = PerlinNoiseS_curve(rx0);
	sy = PerlinNoiseS_curve(ry0);

#define at2(rx,ry) ( rx * q[0] + ry * q[1] )

	q = PerlinNoise_g2[b00];
	u = at2(rx0, ry0);
	q = PerlinNoise_g2[b10];
	v = at2(rx1, ry0);
	a = PerlinNoiseLerp(sx, u, v);

	q = PerlinNoise_g2[b01];
	u = at2(rx0, ry1);
	q = PerlinNoise_g2[b11];
	v = at2(rx1, ry1);
	b = PerlinNoiseLerp(sx, u, v);

	return PerlinNoiseLerp(sy, a, b);
}
float PerlinNoise2D(float vec[2])
{
	int terms = myDome.nNoiseOctave; //myNoise.nOctaves;
	float freq = myDome.fNoiseDensity; //myNoise.dFrequency;
	float result = 0.0f;
	float amp = myNoise.fAmplitude;

	vec[0] *= myDome.fNoiseDensity;
	vec[1] *= myDome.fNoiseDensity;

	for (int i = 0; i < terms; i++)
	{
		result += PerlinNoise2(vec)*amp;
		vec[0] *= 2.0f;
		vec[1] *= 2.0f;
		amp *= 0.5f;
	}

	return result;
}
float PerlinNoiseGet2D(float x, float y)
{
	float vec[2];
	vec[0] = x;
	vec[1] = y;
	return PerlinNoise2D(vec);
};





void InitParameters()
{


	memset(&myDome, '0', sizeof(myDome));



	myDome.nDomeVersion = 13;
	strcpy_s(myDome.sDomeName, "my Name");
	myDome.fDomeScaleZ = 1.0;
	strcpy_s(myDome.sDomeMaterial, S_NONE);
	myDome.nDomeEnd = 360;
	myDome.fDomeSocle = 0;
	myDome.fDomeInFolding = 0;
	myDome.fTopPart = 0;
	myDome.fTopParameter = 1;
	myDome.nTopType = TOP_NONE;

	if (lbDomeFile)
		lbDomeFile->set_int_val(0);

	BaseDefault();

	if (myDome.nBaseSides % 2 && ckDivision)
		ckDivision->disable();

//	if (spBaseRound)
//		spBaseRound->disable();

	SlopeDefault();

	WindowDefault();
	myDome.nTopType = 0;
	myDome.fTopParameter = 1;
	myDome.fTopPart = .5;

	myDome.bFaceGenerate = 1;
	myDome.bFaceSmooth = 0;
	strcpy_s(myDome.sFaceCenterMaterial, S_NONE);
	strcpy_s(myDome.sFaceEdgeMaterial, S_NONE);
	myDome.nMeridianEdgeSize = 0;
	myDome.nParallelEdgeSize = 0;

	myDome.bWireGenerate = 1;
	myDome.fWireShrinking = 0.0;
	myDome.nWireShrinkingType = SHRINK_NONE;
	myDome.fWireShrinkingSpeed = 0.0;

	myDome.nMeridianMain = 4;
	myDome.fMeridianMainRadius = (float)0.03;
	strcpy_s(myDome.sMeridianMainMaterial, S_NONE);

	myDome.nMeridianMinor = 5;
	myDome.fMeridianMinorRadius = (float)0.01;
	strcpy_s(myDome.sMeridianMinorMaterial, S_NONE);

	myDome.nParallelMain = 20;
	myDome.fParallelMainRadius = (float)0.02;
	strcpy_s(myDome.sParallelMainMaterial, S_NONE);

	myDome.nParallelMinor = 1;
	myDome.fParallelMinorRadius = (float)0.005;
	strcpy_s(myDome.sParallelMinorMaterial, S_NONE);

	myDome.fDiagonalURRadius = (float)0.005;
	strcpy_s(myDome.sDiagonalURMaterial, S_NONE);

	myDome.fDiagonalDRRadius = (float)0.005;
	strcpy_s(myDome.sDiagonalDRMaterial, S_NONE);

	myDome.fSphereMPRadius = (float)0.03;
	strcpy_s(myDome.sSphereMPMaterial, S_NONE);

	myDome.fSphereMpRadius = (float)0.03;
	strcpy_s(myDome.sSphereMpMaterial, S_NONE);

	myDome.fSpheremPRadius = (float)0.02;
	strcpy_s(myDome.sSpheremPMaterial, S_NONE);

	myDome.fSpherempRadius = (float)0.01;
	strcpy_s(myDome.sSpherempMaterial, S_NONE);

	NoiseDefault();

	LissageDefault();

	myDome.fDiffGeomLambda = 0.5f;
	myDome.nDiffGeomNombre = 1;
	myDome.bDiffGeomMeridien = false;
	myDome.bDiffGeomParallel = false;

	myDome.bNormale = false;
	myDome.fNormaleLambda = 1.0;
	myDome.nNormaleNombre = 1;

	myDome.nSpirale = 0;
	myDome.nSpiraleShrink = SHRINK_NONE;
	myDome.fSpiraleRadius = (float).03;
	strcpy_s(myDome.sSpiraleMaterial, S_NONE);

	myDome.fSpiraleTour = 1;
	myDome.fSpiraleSpeed = 1.0;
	myDome.fSpiraleDelay = 0.0;
	myDome.fSpiraleDecalage= 0.0;
	myDome.nSpiraleChirale = CHIRAL_LEFT;
	myDome.nSpiraleStep = 1;
	myDome.nSpiraleToward = TOWARD_TOP;
	myDome.nSpiraleSlope = SLOPE_STRAIGHT;
	myDome.fSpiraleSlopeParameter = 1;

	myDome.nFrise = 0;
	myDome.fFriseHeight = (float)0.3;
	myDome.fFrisePerSide = (float)1;
	myDome.fFriseDepart = (float)0;
	myDome.fFriseShift = (float).05;
	myDome.fFriseRadius = (float).03;
	myDome.nFriseShrink = SHRINK_NONE;
	myDome.nFriseFunction = WINDOW_NONE;
	strcpy_s(myDome.sFriseMaterial, S_NONE);
	myDome.nFriseToward = TOWARD_TOP;
	myDome.nFriseStep = 1;
	myDome.fFrisePower = 1;
	myDome.fFriseOffset = 0;
	myDome.fFriseParameter = 1;
	myDome.fFriseGap = 0.0;
	myDome.bFriseInverse = 0;
	myDome.bFriseReverse = 0;
	myDome.fFriseSpacing = 0.0;

	myDome.nSnake = 0;;
	myDome.fSnakeTour = 1;
	myDome.nSnakeEvolve = SHRINK_LINEAR;
	myDome.fSnakeShift = 0;
	myDome.fSnakeRadius = (float).03;
	myDome.nSnakeShrink = SHRINK_NONE;
	strcpy_s(myDome.sSnakeMaterial, S_NONE);
	myDome.nSnakeChirale = CHIRAL_LEFT;
	myDome.nSnakeStep = 1;
	myDome.nSnakeToward = TOWARD_TOP;
	myDome.fSnakeOffset = 0;
	myDome.nSnakeSlope = SLOPE_STRAIGHT;
	myDome.fSnakeSlopeParameter = 1;

	myDome.fDerivative = 0;
	myDome.bDerivativeNorme = 0;

	myDome.fBumpDamping = 0;
	myDome.fBumpPower = 1;
	myDome.nBumpNormalize = 0;
	myDome.nBumpType = BUMP_FUNNEL;
	myDome.nBumpNormal = BUMP_MEAN;
	myDome.fBumpLambda = 1;
	myDome.bBump = 0;
	myDome.fBumpDistance = 1;
	myDome.nBumpDistanceType= DISTANCE_EUCLIDE;
	myDome.nBumpDistanceInverse = 0;
	myDome.nPostMeridianEvery = 1;
	myDome.nPostMeridianShift = 0;
	myDome.nPostParallelEvolve =0;
	myDome.nPostParallelFirst = 0;
	myDome.nPostParallelEvery = 1;

	myDome.bMarquee = 0;
	myDome.nMarqueeFunction = SLOPE_STRAIGHT;
	myDome.fMarqueeParameter = 1;
	myDome.fMarqueeLambda = 1;
	myDome.fMarqueeDamping = 0;
	myDome.fMarqueeAlpha = 0.5f;
	myDome.fMarqueeBeta = 0.5f;
	myDome.bMarqueeNormalize = 0;
	myDome.fMarqueeTop = 1.0f;
	myDome.bPostMeridianPerSide = 0;
	myDome.nDifferentialGeometry = DIFF_GEOM_EVOLUTE;
	myDome.nEnhanceType = ENHANCE_NONE;
	myDome.bEnhanceParallel = 0;
	myDome.bEnhanceMeridian = 0;
	myDome.nEnhanceMerge = MERGE_SOMME;
	myDome.nEnhanceCalcul = CALCUL_Produit;
	myDome.fEnhanceLambda = 1.1f;
	myDome.fEnhanceDamping = 0.0f;

	myDome.nComb = 0;
	myDome.fCombRadius = 0.03f;
	myDome.nCombShrink = SHRINK_NONE;
	strcpy_s(myDome.sCombMaterial, S_NONE);	
	myDome.fCombSpacing = 1.0;
	myDome.nCombStep = 1;
	myDome.nCombToward = TOWARD_TOP;
	myDome.nCombThread = 1;
	myDome.nCombType = COMB_PARALLEL_EDGE;
	myDome.nCombThreadGap = 2;
	myDome.fCombOffset = 0;
	SpikeDefault();

	myDome.bLid = false;
	strcpy_s(myDome.sLidMaterial, S_NONE);


	myNoise.nOctaves = myDome.nNoiseOctave;
	myNoise.fFrequency = myDome.fNoiseDensity;
	myNoise.fAmplitude = 1.0;
	myNoise.nSeed = myDome.nNoiseSeed;
	myNoise.bStart = true;



	myDisplay.nDisplayType = RB_WIRE;
	myDisplay.bDisplayMinor = 0;
	myDisplay.fScaleDisplay = 1.0;
	myDisplay.nLight_1 = 1;
	myDisplay.nLight_2 = 1;
	myDisplay.fSpeedCoefficient = 1.0;
	myDisplay.nMaterial = MAT_GOLD;
	myDisplay.bWireFrame = true;
	myDisplay.bSmooth = false;
	myDisplay.bFlatFace = false;
	if (ckWireFrame)
		ckWireFrame->disable();

	//fSlopeParameterXY = myDome.fSlopeParameterXY;
	ChargeSlope(myDome.nSlope);
	ComputeSlopeLength();
	bAllocation = 1;
	bDomeModifie = false;

	if (tbBase)
	{
		tbBase->set_text(BaseFunctionListe[myDome.nBase].sNom);
		tbWindow->set_text(WindowFunctionListe[myDome.nWindow].sNom);
		tbSlope->set_text(SlopeFunctionListe[myDome.nSlope].sNom);
		tbTop->set_text(TopListe[myDome.nTopType].sNom);
	}

	ChargeMesh();
	if (lbDomeFile)
		lbDomeFile->set_int_val(0);

}

void LoadDome()
{
	FILE *DomeFile;
	errno_t Erreur;
	char FileName[400];
	int  Retour = 0;


	strcpy_s(FileName, myDome.sDomeName);
	strcat_s(FileName, ".DOM");

	Erreur = fopen_s(&DomeFile, FileName, "rb");

	//	Retour = fread(&myDome, sizeof(myDome), 1, DomeFile);

	if (Erreur == 0)
	{
		Retour = fread(&myDome, 1, sizeof(myDome), DomeFile);
		fclose(DomeFile);
	}
	printf("Read %d \n", Retour);
	SpikeDefault();
	if (myDome.nDomeVersion <= 13)
	{
		if (myDome.nDomeVersion < 2)
		{
			myDome.fWindowDamping = 0.0;
			myDome.fBaseExpand = 1.0;
		}
		if (myDome.nDomeVersion < 3)
		{
			DefaultSlopeParameters();
		}
		if (myDome.nDomeVersion < 4)
		{
			myDome.fBaseDamping = 1.0;
		}
		if (myDome.nDomeVersion < 5)
		{
//			float  fSlopeParameterXY; //Attention load and save
//			float  fWindowPower;
//			int    nDomeEnd;
			myDome.fDomeSocle = 0.0;
			LissageDefault();

			myDome.fDiffGeomLambda = 0.5f;
			myDome.nDiffGeomNombre = 1;
			myDome.bDiffGeomMeridien = false;
			myDome.bDiffGeomParallel = false;

			myDome.bNormale = false;
			myDome.fNormaleLambda = 1.0;
			myDome.nNormaleNombre = 1;

		}
		if (myDome.nDomeVersion < 6)
		{
			myDome.nSpirale = 0;
			myDome.nSpiraleShrink = SHRINK_NONE;
			myDome.fSpiraleRadius = (float).03;
			strcpy_s(myDome.sSpiraleMaterial, S_NONE);
			myDome.fSpiraleTour = 1;
			myDome.fSpiraleSpeed = 1.0;
			myDome.fSpiraleDelay = 0.0;
			myDome.fSpiraleDecalage = 0.0;
			myDome.nSpiraleChirale = CHIRAL_LEFT;
			myDome.nSpiraleStep = 1;
			myDome.nSpiraleToward = TOWARD_TOP;

			myDome.nFrise = 0;
			myDome.fFriseHeight = (float).3;
			myDome.fFrisePerSide = (float)1;
			myDome.fFriseDepart = (float)0;
			myDome.fFriseShift = (float).05;
			myDome.fFriseRadius = (float).03;
			myDome.nFriseShrink = SHRINK_NONE;
			myDome.nFriseFunction = WINDOW_NONE;
			strcpy_s(myDome.sFriseMaterial, S_NONE);
			myDome.nFriseToward = TOWARD_TOP;
			myDome.nFriseStep = 1;
			myDome.fFrisePower = 1;
			myDome.fFriseOffset = 0;

			myDome.fDerivative = 0;
			myDome.bDerivativeNorme = 0;

			myDome.nSnake = 0;;
			myDome.fSnakeTour = 1;
			myDome.nSnakeEvolve = SHRINK_LINEAR;
			myDome.fSnakeShift = 0;
			myDome.fSnakeRadius = (float).03;
			myDome.nSnakeShrink = SHRINK_NONE;
			strcpy_s(myDome.sSnakeMaterial, S_NONE);
			myDome.nSnakeChirale = CHIRAL_LEFT;
			myDome.nSnakeStep = 1;
			myDome.nSnakeToward = 0;
			myDome.fSnakeOffset = 0;
			strcpy_s(myDome.sFaceEdgeMaterial, myDome.sFaceCenterMaterial);
			myDome.nMeridianEdgeSize = 0;
			myDome.nParallelEdgeSize = 0;
			myDome.fBumpDamping = 0;
			myDome.fBumpPower = 1;
			myDome.nBumpNormalize = 0;
			myDome.nBumpNormal = BUMP_MEAN;
			myDome.nBumpType = BUMP_FUNNEL;
			myDome.fBumpLambda = 1;
			myDome.bBump = 0;
			myDome.fBumpDistance = 1;
			myDome.nBumpDistanceType = DISTANCE_EUCLIDE;
			myDome.nBumpDistanceInverse = 0;
			myDome.nPostMeridianEvery = 1;
			myDome.nPostMeridianShift = 0;
			myDome.nPostParallelEvolve = 0;
			myDome.nPostParallelFirst = 0;
			myDome.nPostParallelEvery = 1;

			SpikeDefault();
			myDome.bLid = false;
			strcpy_s(myDome.sLidMaterial, S_NONE);

		}
		if (myDome.nDomeVersion < 7)
		{
			myDome.fTopPart = 0.0;
			myDome.bSlopeInverse = false;
			myDome.fBaseRound = 0.0;
			myDome.fBaseBaryCentre = 0.0;
		}
		if (myDome.nDomeVersion < 8)
		{
			myDome.nCrenelPerSide = 0;
			myDome.fCrenelFront = 0.1f;
			myDome.fCrenelThick = 0.2f;
			myDome.fCrenelWide = 0.4f;
			myDome.fCrenelHeight = 0.5f;
			myDome.fCrenelPower = 0.3f;
			myDome.fWindowAwning = 0.0f;
			myDome.fWindowParameter = 0.0;
			myDome.fBaseMorphing = 0.0;


			myDome.fBaseSecondaryRotation = 0.0;
			myDome.fBaseSecondaryExpansion = 1.0;
			myDome.bBaseMorphingDivision = 0;
			myDome.bWindowBlind = 0;
			myDome.fFriseParameter = 1;
			myDome.fFriseGap = 0.0;

			GetFunctions(myDome.nSlope);

		}
		if (myDome.nDomeVersion < 9)
		{
			GetFunctions(myDome.nSlope);
			myDome.fBaseSpiral = 0.0;
			myDome.nSnakeSlope = SLOPE_STRAIGHT;
			myDome.fSnakeSlopeParameter = 1;
			myDome.nSpiraleSlope = SLOPE_STRAIGHT;
			myDome.fSpiraleSlopeParameter = 1;
			myDome.fLissagePoidsVertical=0;
			myDome.nLissageTypeFiltre = FILTER_Gauss;
			myDome.nLissageDirection = DIRECTION_NORTH;
			myDome.fLissageSigma1 = 1.0;
			myDome.fLissageSigma2 = 1.0;

		}
		if (myDome.nDomeVersion < 10)
		{
			myDome.nTopType = TOP_NONE;
			myDome.fTopParameter = 1;
			myDome.fDomeInFolding = 0;
			myDome.fCrenelAngle = 0.0f;
			myDome.fCrenelDelay = 0.0f;
			myDome.nBaseEvenly = RB_ALONG_PARAMETER;
		}
		if (myDome.nDomeVersion < 11)
		{
			myDome.fSpikeAngle = 1.0f;
			myDome.bFriseInverse = 0;
			myDome.bFriseReverse = 0;
		}
		if (myDome.nDomeVersion < 12)
		{

			myDome.bMarquee = 0;
			myDome.nMarqueeFunction = SLOPE_STRAIGHT;
			myDome.fMarqueeParameter = 1;
			myDome.fMarqueeLambda = 1;
			myDome.fMarqueeDamping = 0;
			myDome.fMarqueeAlpha = 0.5f;
			myDome.fMarqueeBeta = 0.5f;
			myDome.bMarqueeNormalize = 0;
			myDome.fMarqueeTop = 1.0f;
			myDome.bPostMeridianPerSide = 0;
			myDome.nDifferentialGeometry = DIFF_GEOM_EVOLUTE;

			myDome.nEnhanceType = ENHANCE_NONE;
			myDome.bEnhanceParallel = 0;
			myDome.bEnhanceMeridian = 0;
			myDome.nEnhanceMerge = MERGE_SOMME;
			myDome.nEnhanceCalcul = CALCUL_Produit;
			myDome.fEnhanceLambda = 1.1f;
			myDome.fEnhanceDamping = 0.0f;
			myDome.fWindowLower = 0.0f;

			myDome.nComb = 0;
			myDome.fCombRadius = 0.03f;
			myDome.nCombShrink = SHRINK_NONE;
			strcpy_s(myDome.sCombMaterial, S_NONE);
			myDome.fCombSpacing = 1.0;
			myDome.nCombStep = 1;
			myDome.nCombToward = TOWARD_TOP;
			myDome.nCombThread = 1;
			myDome.nCombType = COMB_PARALLEL_EDGE;
			myDome.nCombThreadGap = 2;
			myDome.fCombOffset = 0;

			myDome.fFriseSpacing = 0.0;

		}
		myDome.nDomeVersion = 13;
		ChargeShrink(myDome.nWireShrinkingType);
		ChargeTestSlope(myDome.nTestSlopeXY, myDome.nTestSlopeZ);

		ComputeSlopeLength();
		bAllocation = 1;
	}
	else
	{
		Message("     Not a Valid DOM File   ");
		InitParameters();
		glParameter->sync_live();
		glDisplay->sync_live();
	}

	FilterInput((TYPE_FILTER)myDome.nLissageTypeFiltre, false);
	if (myDome.nBase == BASE_FOURIER)
		BaseFourierRho(0, 0, true);

	tbBase->set_text(BaseFunctionListe[myDome.nBase].sNom);
	tbWindow->set_text(WindowFunctionListe[myDome.nWindow].sNom);
	tbSlope->set_text(SlopeFunctionListe[myDome.nSlope].sNom);
	tbTop->set_text(TopListe[myDome.nTopType].sNom);

	glParameter->sync_live();
	glDecor->sync_live();
	glPostProcess->sync_live();
	glPreProcess->sync_live();
	glBase->sync_live();
	glSlope->sync_live();
	glTop->sync_live();
	glWindow->sync_live();
}

void SaveDome()
{
	FILE *DomeFile;
	errno_t Erreur;
	char FileName[400];
	int  Retour;


	strcpy_s(FileName, myDome.sDomeName);
	strcat_s(FileName, ".DOM");


	Erreur = fopen_s(&DomeFile, FileName, "wb");
	if (Erreur != 0)
		Message("     Problem Saving DOM File   ");
	else
	{
		Retour = fwrite(&myDome, 1, sizeof(myDome), DomeFile);
		fflush(DomeFile);
		fclose(DomeFile);

		if (Retour != sizeof(myDome))
			Message("     Problem saving the file   ");

		printf("Written %d \n", Retour);

		fill_lbfile(lbDomeFile, (char *)".DOM", FileName);

	}
}


double CalculSegmentX(double u0, double u1)
{
	double x0, x1, y0, y1, xm, ym;
	double da, db, d;
	double epsilon = 0.000000001;

	x0 = BaseFunction(u0, myDome.nBaseSides, myDome.fBaseParameter, BASE_X, 0);
	x1 = BaseFunction(u1, myDome.nBaseSides, myDome.fBaseParameter, BASE_X, 0);
	xm = BaseFunction((u0 + u1) / 2.0, myDome.nBaseSides, myDome.fBaseParameter, BASE_X, 0);
	y0 = BaseFunction(u0, myDome.nBaseSides, myDome.fBaseParameter, BASE_Y, 0);
	y1 = BaseFunction(u1, myDome.nBaseSides, myDome.fBaseParameter, BASE_Y, 0);
	ym = BaseFunction((u0 + u1) / 2.0, myDome.nBaseSides, myDome.fBaseParameter, BASE_Y, 0);

	da = sqrt((x0 - xm) * (x0 - xm) + (y0 - ym) * (y0 - ym));
	db = sqrt((xm - x1) * (xm - x1) + (ym - y1) * (ym - y1));
	d = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));

	if (da + db - d < epsilon)
		return da + db;
	else
		return CalculSegmentX(u0, (u0 + u1) / 2.0) + CalculSegmentX((u0 + u1) / 2.0, u1);
}



double CalculBaseLength(int nMax)
{
	double dLength = 0;

	for (int i = 0; i < nMax; i++)
	{
		dLength += CalculSegmentX(2.0 * _PI * (double)(i) / (double)nMax, 2.0 * _PI * (double)(i + 1) / (double)nMax);
	}
	return dLength;
}

double GetNextAngle(double dTest, double u, int nMax)
{
	double dAngle = 0;
	double x1, y1;
	double dPlus = 2.0 * _PI / (20 * nMax);

	for (; dAngle <= dTest && u < 2*_PI; u += dPlus)
	{
		if (u > 2 * _PI)
			u = 2 * _PI;
		x1 = BaseFunction(u , myDome.nBaseSides, myDome.fBaseParameter, BASE_X, 0);
		y1 = BaseFunction(u , myDome.nBaseSides, myDome.fBaseParameter, BASE_Y, 0);

//FuncX(u * 2.0 * _PI, nBaseSides, BaseFunctionListe[i].fBaseParameter);
//FuncY(u * 2.0 * _PI, nBaseSides, BaseFunctionListe[i].fBaseParameter);

		dAngle = myMod((float)(2.0 * _PI + atan2(y1, x1)), (float)(2.0 * _PI));

	}
	return u;
}

double GetNextLength1(double dTest, double u, double* dLength, int nMax)
{
	double x0, y0, x1, y1;

	x0 = BaseFunction(u, myDome.nBaseSides, myDome.fBaseParameter, BASE_X, 0);
	y0 = BaseFunction(u, myDome.nBaseSides, myDome.fBaseParameter, BASE_Y, 0);
	for (; *dLength <= dTest && u < 2.0 * _PI; )
	{
		u += (2.0 * _PI) / (50.0 * nMax);
		if (u > 2.0 * _PI)
			u = 2.0 * _PI;
		x1 = BaseFunction(u, myDome.nBaseSides, myDome.fBaseParameter, BASE_X, 0);
		y1 = BaseFunction(u, myDome.nBaseSides, myDome.fBaseParameter, BASE_Y, 0);
		*dLength += sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
		x0 = x1;
		y0 = y1;

	}

	return u;
}

double GetNextLength(double dTest, double u, double* dLength, int nMax)
{
	double dInc, epsilon;

	double dEcart;

	dInc = (2.0 * _PI) / (double)(50 * nMax);
	epsilon = 0.000001;
	while (*dLength < dTest && u < 2.0 * _PI)
	{
		if (u + dInc > 2.0 * _PI)
			u = 2.0 * _PI - dInc;

		dEcart = CalculSegmentX(u, u + dInc);
		*dLength += dEcart;
		u += dInc;
	}

	return u;
}

void CalculMesh(int nMaxI, int nMaxJ, GLpoint *Vertices)
{
	int i, j;

	double dX, dY, dZ;
	double dTest;

	int nmi;

/*	//Length of Base
	double X1, Y1, Len;
	double X = BaseFunctionX(0, myDome.nBaseSides, myDome.fBaseParameter);
	double Y = BaseFunctionY(0, myDome.nBaseSides, myDome.fBaseParameter);
	Len = 0;
	for (double u = 0; u < 2 * _PI;u += 0.0001)
	{
		X1 = BaseFunctionX(u, myDome.nBaseSides, myDome.fBaseParameter);
		Y1 = BaseFunctionY(u, myDome.nBaseSides, myDome.fBaseParameter);
		Len += sqrt((X1 - X)*(X1 - X) + (Y1 - Y)*(Y1 - Y));
		X = X1;
		Y = Y1;
	}
	printf("Longueur = %lf\n", Len);
	*/
//	if (myDome.nDomeEnd == 360)
		nmi = nMaxI;
//	else
//		nmi = (int)((float)(nMaxI * myDome.nDomeEnd) / 360.0) + 1;

//	printf("\nNombre de sides = ;%d;%d;%d;%d\n", myDome.nBaseSides, (myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor * myDisplay.bDisplayMinor), myDome.nMeridianMain, myDome.nMeridianMinor);
//	printf("Side;Variable;Fraction;Corrigé\n");

	double dBaseLength = CalculBaseLength(nmi);
	double dCurrentLength = 0;
//	printf("%16.14f\n", dBaseLength);

	float u = 0;
	double dStep;
//	printf("\n\nnmi = %d\n", nmi);
	for (i = 0; i < nmi; i++)
	{

		switch (myDome.nBaseEvenly)
		{
		case RB_ALONG_LENGTH:
			dStep = dBaseLength / (double)(nmi - 1);
			dTest = dStep * (double)i;
			
			u = (float)GetNextLength(dTest, u, &dCurrentLength, nmi);
			break;
		case RB_ALONG_ANGLE:
			dStep = 2.0 * _PI / (double)(nmi - 1);
		    dTest = dStep * (double)i;

			u = (float)GetNextAngle(dTest, u, nmi);
//			printf("%f\n", u);
			break;
		case RB_ALONG_PARAMETER:
		default:
			u = (float)(2 * _PI * (float)i / (float)(nMaxI - 1));
			break;
		}
		float fLength = 0;
		float w1 = 0;
		float w;
		float x0 = (float)SlopeFunction(0, myDome.fSlopeParameterZ, SLOPE_Z);
		float y0 = (float)SlopeFunction(0, myDome.fSlopeParameterXY, SLOPE_XY);

		for (j = 0; j < nMaxJ; j++)
		{
			float x1, y1;
			float v = ((float)j / (float)(nMaxJ - 1));

			if (myDome.nSlopeEvenly == RB_ALONG_LENGTH)
			{
				w = w1;
				for (; fLength < (float)(j + 1) * fSlopeLength / (float)(nMaxJ - 1); w1 += (float)0.0001)
				{
					x1 = (float)SlopeFunction(w1, myDome.fSlopeParameterZ, SLOPE_Z);
					y1 = (float)SlopeFunction(w1, myDome.fSlopeParameterXY, SLOPE_XY);
					fLength += sqrtf((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
					x0 = x1;
					y0 = y1;
				}
				if ((w1 > 1.0) || (j == nMaxJ - 2))  //1
					w1 = 1.0;
				//				w = 1.0;
			}
			else if (myDome.nSlopeEvenly == RB_ALONG_HEIGHT)
			{
				w = w1;
				x1 = x0;
				if (j < nMaxJ - 1)
				{
					for (float x1 = x0; x1 < (float)(j + 1) / (float)(nMaxJ - 1); w1 += (float)0.0001)
					{
						x1 = (float)SlopeFunction(w1, myDome.fSlopeParameterZ, SLOPE_Z);
						x0 = x1;
					}
					if ((w1 > 1.0) || (j == nMaxJ - 1))
						w1 = 1.0;
				}
				else
					w = 1.0;
			}

			else
				w = v;

			Calcul(u, w, &dX, &dY, &dZ, true);

			//		if (myDome.nBase == BASE_ROOF)
			//			fprintf_s(DomFile, "u = %lf, v = %lf, X = %lf, Y= %lf %s", u, w, dX, dY, CRLF);

			if ((i == nMaxI - 1) && (myDome.nBase != BASE_ROOF))
			{
				Vertices[i * nMaxJ + j].x = Vertices[0 * nMaxJ + j].x;
				Vertices[i * nMaxJ + j].y = Vertices[0 * nMaxJ + j].y;
				Vertices[i * nMaxJ + j].z = Vertices[0 * nMaxJ + j].z;
			}
			else if ((j == 0) && (i > 0) && (myDome.nBase != BASE_ROOF))//£
			{
				Vertices[i * nMaxJ + j].x = Vertices[0 * nMaxJ + 0].x;
				Vertices[i * nMaxJ + j].y = Vertices[0 * nMaxJ + 0].y;
				Vertices[i * nMaxJ + j].z = Vertices[0 * nMaxJ + 0].z;
			}
			else
			{
				Vertices[i * nMaxJ + j].x = (float)dX;
				Vertices[i * nMaxJ + j].y = (float)dY;
				Vertices[i * nMaxJ + j].z = (float)dZ;
			}

		}
	}
//	printf("\nFin\n");
//	fflush(DomFile);
//	fclose(DomFile);

//	printf("\n%lf %lf \n%lf %lf\n", MinU, MaxU, MinV, MaxV);
}

void Normales(int i, int j, int nMaxI, int nMaxJ, GLpoint *Vertices, GLpoint *Normals)
{
	GLpoint Vecteur[4];
	GLpoint Temp;
	float   fNorme;  

	if (i == 0)
		i = (nMaxI - 1);

	if (j == 0)
	{
		Vecteur[0].x = Vertices[(i - 1) * nMaxJ + 1].x - Vertices[i * nMaxJ + j].x;
		Vecteur[0].y = Vertices[(i - 1) * nMaxJ + 1].y - Vertices[i * nMaxJ + j].y;
		Vecteur[0].z = Vertices[(i - 1) * nMaxJ + 1].z - Vertices[i * nMaxJ + j].z;

		Vecteur[1].x = Vertices[i * nMaxJ + 1].x - Vertices[i * nMaxJ + j].x;
		Vecteur[1].y = Vertices[i * nMaxJ + 1].y - Vertices[i * nMaxJ + j].y;
		Vecteur[1].z = Vertices[i * nMaxJ + 1].z - Vertices[i * nMaxJ + j].z;

		CrossProduct(&Vecteur[0], &Vecteur[1], &Temp);
		Normals->x = -Temp.x;
		Normals->y = -Temp.y;
		Normals->z = -Temp.z;

		//		Normals->x = 0;
		//		Normals->y = 0;
		//		Normals->z = -1;
		return;
	}


	if (i == 0)
		i = (nMaxI - 1);

	Vecteur[0].x = Vertices[(i - 1) * nMaxJ + j].x - Vertices[i * nMaxJ + j].x;
	Vecteur[0].y = Vertices[(i - 1) * nMaxJ + j].y - Vertices[i * nMaxJ + j].y;
	Vecteur[0].z = Vertices[(i - 1) * nMaxJ + j].z - Vertices[i * nMaxJ + j].z;

	Vecteur[1].x = Vertices[i * nMaxJ + (j - 1)].x - Vertices[i * nMaxJ + j].x;
	Vecteur[1].y = Vertices[i * nMaxJ + (j - 1)].y - Vertices[i * nMaxJ + j].y;
	Vecteur[1].z = Vertices[i * nMaxJ + (j - 1)].z - Vertices[i * nMaxJ + j].z;

/*		Vecteur[2].x = Vertices[((i + 1)% nMaxI) * nMaxJ + j].x  -	Vertices[i * nMaxJ + j].x;
		Vecteur[2].y = Vertices[((i + 1)% nMaxI) * nMaxJ + j].y  -	Vertices[i * nMaxJ + j].y;
		Vecteur[2].z = Vertices[((i + 1)% nMaxI) * nMaxJ + j].z  -	Vertices[i * nMaxJ + j].z;
*/
	Vecteur[2].x = Vertices[(i % (nMaxI - 1) + 1) * nMaxJ + j].x - Vertices[i * nMaxJ + j].x;
	Vecteur[2].y = Vertices[(i % (nMaxI - 1) + 1) * nMaxJ + j].y - Vertices[i * nMaxJ + j].y;
	Vecteur[2].z = Vertices[(i % (nMaxI - 1) + 1) * nMaxJ + j].z - Vertices[i * nMaxJ + j].z;


	CrossProduct(&Vecteur[0], &Vecteur[1], &Temp);
	Normals->x = Temp.x;
	Normals->y = Temp.y;
	Normals->z = Temp.z;
	CrossProduct(&Vecteur[1], &Vecteur[2], &Temp);
	Normals->x += Temp.x;
	Normals->y += Temp.y;
	Normals->z += Temp.z;
	if (j < nMaxJ - 1)
	{
		Vecteur[3].x = Vertices[i * nMaxJ + (j + 1)].x - Vertices[i * nMaxJ + j].x;
		Vecteur[3].y = Vertices[i * nMaxJ + (j + 1)].y - Vertices[i * nMaxJ + j].y;
		Vecteur[3].z = Vertices[i * nMaxJ + (j + 1)].z - Vertices[i * nMaxJ + j].z;

		CrossProduct(&Vecteur[2], &Vecteur[3], &Temp);
		Normals->x += Temp.x;
		Normals->y += Temp.y;
		Normals->z += Temp.z;
		CrossProduct(&Vecteur[3], &Vecteur[1], &Temp);
		Normals->x += Temp.x;
		Normals->y += Temp.y;
		Normals->z += Temp.z;
	}

	fNorme = (sqrtf((Normals->x)*(Normals->x) + (Normals->y)*(Normals->y) + (Normals->z)*(Normals->z)));

	if (fNorme > 0.0)
	{
		Normals->x = Normals->x / fNorme;
		Normals->y = Normals->y / fNorme;
		Normals->z = Normals->z / fNorme;
	}

}

void CalculNormales(int nMaxI, int nMaxJ, GLpoint *Vertices, GLpoint *Normals)
{
	GLpoint Vecteur[4];
	GLpoint Temp;

	for (int i = 1; i < nMaxI; i++)
	{
		for (int j = 1; j < nMaxJ; j++)
		{
			Vecteur[0].x = Vertices[(i - 1) * nMaxJ + j].x - Vertices[i * nMaxJ + j].x;
			Vecteur[0].y = Vertices[(i - 1) * nMaxJ + j].y - Vertices[i * nMaxJ + j].y;
			Vecteur[0].z = Vertices[(i - 1) * nMaxJ + j].z - Vertices[i * nMaxJ + j].z;

			Vecteur[1].x = Vertices[i * nMaxJ + (j - 1)].x - Vertices[i * nMaxJ + j].x;
			Vecteur[1].y = Vertices[i * nMaxJ + (j - 1)].y - Vertices[i * nMaxJ + j].y;
			Vecteur[1].z = Vertices[i * nMaxJ + (j - 1)].z - Vertices[i * nMaxJ + j].z;

			Vecteur[2].x = Vertices[((i + 1) % nMaxI) * nMaxJ + j].x - Vertices[i * nMaxJ + j].x;
			Vecteur[2].y = Vertices[((i + 1) % nMaxI) * nMaxJ + j].y - Vertices[i * nMaxJ + j].y;
			Vecteur[2].z = Vertices[((i + 1) % nMaxI) * nMaxJ + j].z - Vertices[i * nMaxJ + j].z;

			CrossProduct(&Vecteur[0], &Vecteur[1], &Temp);
			Normals[i * nMaxJ + j].x = Temp.x;
			Normals[i * nMaxJ + j].y = Temp.y;
			Normals[i * nMaxJ + j].z = Temp.z;
			CrossProduct(&Vecteur[1], &Vecteur[2], &Temp);
			Normals[i * nMaxJ + j].x += Temp.x;
			Normals[i * nMaxJ + j].y += Temp.y;
			Normals[i * nMaxJ + j].z += Temp.z;

			if (j < nMaxJ - 1)
			{
				Vecteur[3].x = Vertices[i * nMaxJ + (j + 1)].x - Vertices[i * nMaxJ + j].x;
				Vecteur[3].y = Vertices[i * nMaxJ + (j + 1)].y - Vertices[i * nMaxJ + j].y;
				Vecteur[3].z = Vertices[i * nMaxJ + (j + 1)].z - Vertices[i * nMaxJ + j].z;

				CrossProduct(&Vecteur[2], &Vecteur[3], &Temp);
				Normals[i * nMaxJ + j].x += Temp.x;
				Normals[i * nMaxJ + j].y += Temp.y;
				Normals[i * nMaxJ + j].z += Temp.z;
				CrossProduct(&Vecteur[3], &Vecteur[1], &Temp);
				Normals[i * nMaxJ + j].x += Temp.x;
				Normals[i * nMaxJ + j].y += Temp.y;
				Normals[i * nMaxJ + j].z += Temp.z;
			}
		}
		Normals[i * nMaxJ + 0].x = 0;
		Normals[i * nMaxJ + 0].y = 0; // à vérifier
		Normals[i * nMaxJ + 0].z = -1; // à vérifier
	}
	for (int j = 0; j < nMaxJ; j++)
	{
		Normals[0 * nMaxJ + j].x = Normals[(nMaxI - 1) * nMaxJ + j].x;
		Normals[0 * nMaxJ + j].y = Normals[(nMaxI - 1) * nMaxJ + j].y;
		Normals[0 * nMaxJ + j].z = Normals[(nMaxI - 1) * nMaxJ + j].z;

	}

}

int myMod(int x, int y)
{
	int z;
	z = x % y;
	while (z < 0)
	{
		z += y;
	}
	return z;
}



float myMod(float x, float y)
{
	float z;
	z = (float)fmod(x, y);
	while (z < 0)
	{
		z = (float)fmod(z + y, y);
	}
	return z;
}
double myMod(double x, double y)
{
	double z;
	z = fmod(x, y);
	while (z < 0)
	{
		z = fmod(z + y, y);
	}
	return z;
}


int SpikeLength(GLpoint *Vertice, int i, int j, int nMaxI, int nMaxJ, Vector3D *V )
{
	float fSpikeMP = myDome.fSpikeMP;
	float fSpikeMp = myDome.fSpikeMp;
	float fSpikemP = myDome.fSpikemP;
	float fSpikemp = myDome.fSpikemp;
	float fSpike;
	float a = myDome.fSpikeAngle;
	float fModule;

	Vector3D W, T;
	
	if (i == 0 && j == 0)
	{
		V->x = 2 * Vertice[0].x * (fSpikeMP * ShrinkStruts(0, nMaxJ));
		V->y = 2 * Vertice[0].y * (fSpikeMP * ShrinkStruts(0, nMaxJ));
		V->z = 2 * Vertice[0].z * (fSpikeMP * ShrinkStruts(0, nMaxJ));
	}
	else
	{
		if ((j % (myDome.nParallelMinor + 1) == 0) && (i % (myDome.nMeridianMinor + 1) == 0))
			fSpike = fSpikeMP;
		if ((j % (myDome.nParallelMinor + 1) != 0) && (i % (myDome.nMeridianMinor + 1) == 0))
			fSpike = fSpikeMp;
		if ((j % (myDome.nParallelMinor + 1) == 0) && (i % (myDome.nMeridianMinor + 1) != 0))
			fSpike = fSpikemP;
		if ((j % (myDome.nParallelMinor + 1) != 0) && (i % (myDome.nMeridianMinor + 1) != 0))
			fSpike = fSpikemp;
		if(myDisplay.bDisplayMinor == false)
			fSpike = fSpikeMP;

		CalculDerivee(Vertice, i, j, V, nMaxI, nMaxJ);
		CalculDerivee(Vertice, i, j, &W, nMaxI, nMaxJ);

		fModule = (float)sqrt(W.x * W.x + W.y * W.y + W.z * W.z);

		if (W.z == 0)
		{
			V->x = 0;
			V->y = 0;
			V->z = 0;
			return 1;
		}

		T.x = 0;
		T.y = 0;
		T.z = W.z * fModule / fabs(W.z);

		V->x = fModule * ((1 - a * a) * T.x + a * W.x) / sqrt((1 - a * a) * (1 - a * a) * fModule * fModule + a * a * fModule * fModule + 2 * (1 - a*a) * a * T.z * W.z);
		V->y = fModule * ((1 - a * a) * T.y + a * W.y) / sqrt((1 - a * a) * (1 - a * a) * fModule * fModule + a * a * fModule * fModule + 2 * (1 - a*a) * a * T.z * W.z);
		V->z = fModule * ((1 - a * a) * T.z + a * W.z) / sqrt((1 - a * a) * (1 - a * a) * fModule * fModule + a * a * fModule * fModule + 2 * (1 - a*a) * a * T.z * W.z);

/*		float t = (float)sqrt(V->x * V->x + V->y * V->y + V->z * V->z);
		V->x = V->x * fModule / t;
		V->y = V->y * fModule / t;
		V->z = V->z * fModule / t;


		if (i == 2 && j == 3)
			printf(" a = %f, W : %f, T : %f, V : %f, t : %f\n", a, W.x * W.x + W.y * W.y + W.z * W.z, T.x * T.x + T.y * T.y + T.z * T.z, V->x * V->x + V->y * V->y + V->z * V->z, t);
*/

		V->x *= (GLfloat)(fSpike * ShrinkStruts(j, nMaxJ));
		V->y *= (GLfloat)(fSpike * ShrinkStruts(j, nMaxJ));
		V->z *= (GLfloat)(fSpike * ShrinkStruts(j, nMaxJ));
	}
	return 1;
}

int PrintSpike(GLpoint *Vertice, FILE *Dome, int nMaxI, int nMaxJ, int nbObjects, PRINT_TYPE Type)
{
	GLpoint Centre0, Centre1;
	Vector3D V;
	float fSphereRadius;
	char *sSphereMaterial;
	
	float fSpikeMP = myDome.fSpikeMP;
	float fSpikeMp = myDome.fSpikeMp;
	float fSpikemP = myDome.fSpikemP;
	float fSpikemp = myDome.fSpikemp;
	float fSpike;

	if (Type == PRINT_STL_BINARY)
		nbObjects = 0;

	
	// Top Sphere							
	fSphereRadius = MaxRadius(&sSphereMaterial, SPHERE_TOP);
	if (fSphereRadius > 0.0)
	{
		if (!myDome.bWireGenerate)
		{
			Centre0.x = 2 * Vertice[0].x;
			Centre0.y = 2 * Vertice[0].y;
			Centre0.z = 2 * Vertice[0].z;

			if (Type == PRINT_STL_BINARY)
				nbObjects += PrintSTLSphereBinary(Centre0, fSphereRadius*ShrinkStruts(0, nMaxJ), GetColor(sSphereMaterial), Dome);
			else if (Type == PRINT_OBJ)
				nbObjects += PrintOBJSphere(Centre0, fSphereRadius*ShrinkStruts(0, nMaxJ), Dome, nbObjects,  "Top Sphere Spike", sSphereMaterial);
			else if (Type == PRINT_POV)
				PrintPOVSphere(Centre0, fSphereRadius*ShrinkStruts(0, nMaxJ), Dome, sSphereMaterial);

		}
		Centre0.x = 2 * Vertice[0].x;
		Centre0.y = 2 * Vertice[0].y;
		Centre0.z = 2 * Vertice[0].z;
		SpikeLength(Vertice, 0, 0, nMaxI, nMaxJ, &V);

		Centre1.x = Centre0.x + (float)V.x;
		Centre1.y = Centre0.y + (float)V.y;
		Centre1.z = Centre0.z + (float)V.z;

		if (fSpikeMP * ShrinkStruts(0, nMaxJ) > 0)

			if (Type == PRINT_STL_BINARY)
				nbObjects += PrintSTLConeBinary(Centre0, Centre1, fSphereRadius*ShrinkStruts(0, nMaxJ), 0, GetColor(sSphereMaterial), Dome, false);
			else if (Type == PRINT_OBJ)
				nbObjects += PrintOBJCone(Centre0, Centre1, fSphereRadius*ShrinkStruts(0, nMaxJ), 0, Dome, nbObjects, "Top Spikes Cone", sSphereMaterial, false);
			else if (Type == PRINT_POV)
				PrintPOVCone(Centre0, Centre1, fSphereRadius*ShrinkStruts(0, nMaxJ), 0.0f, Dome, sSphereMaterial);
	}

	for (int i = 0; i < nMaxI - 1; i++)
	{
		for (int j = 1; j < nMaxJ; j++)
		{
			//Spheres					
			sSphereMaterial = myDome.sSphereMPMaterial;
			fSphereRadius = 0.0f;
			if ((j % (myDome.nParallelMinor + 1) == 0) && (i % (myDome.nMeridianMinor + 1) == 0))
			{
				sSphereMaterial = myDome.sSphereMPMaterial;
				fSpike = fSpikeMP;
				fSphereRadius = (fSpike > 0.0f ? myDome.fSphereMPRadius : 0);
			}
			if ((j % (myDome.nParallelMinor + 1) != 0) && (i % (myDome.nMeridianMinor + 1) == 0))
			{
				sSphereMaterial = myDome.sSphereMpMaterial;
				fSpike = fSpikeMp;
				fSphereRadius = (fSpike > 0.0f ? myDome.fSphereMpRadius : 0);
			}
			if ((j % (myDome.nParallelMinor + 1) == 0) && (i % (myDome.nMeridianMinor + 1) != 0))
			{
				sSphereMaterial = myDome.sSpheremPMaterial;
				fSpike = fSpikemP;
				fSphereRadius = (fSpike > 0.0f ? myDome.fSpheremPRadius : 0);
			}
			if ((j % (myDome.nParallelMinor + 1) != 0) && (i % (myDome.nMeridianMinor + 1) != 0))
			{
				sSphereMaterial = myDome.sSpherempMaterial;
				fSpike = fSpikemp;
				fSphereRadius = (fSpike > 0.0f ? myDome.fSpherempRadius : 0);
			}
			if (fSphereRadius > 0.0)
			{
				if (!myDome.bWireGenerate)
				{
					Centre0.x = 2 * Vertice[i * nMaxJ + j].x;
					Centre0.y = 2 * Vertice[i * nMaxJ + j].y;
					Centre0.z = 2 * Vertice[i * nMaxJ + j].z;

					if (Type == PRINT_STL_BINARY)
						nbObjects += PrintSTLSphereBinary(Centre0, fSphereRadius*ShrinkStruts(j, nMaxJ), GetColor(sSphereMaterial), Dome);
					else if (Type == PRINT_OBJ)
						nbObjects += PrintOBJSphere(Centre0, fSphereRadius*ShrinkStruts(j, nMaxJ), Dome, nbObjects, "Sphere Spike", sSphereMaterial);
					else if (Type == PRINT_POV)
						PrintPOVSphere(Centre0, fSphereRadius*ShrinkStruts(j, nMaxJ), Dome, sSphereMaterial);
				}
				Centre0.x = 2 * Vertice[i * nMaxJ + j].x;
				Centre0.y = 2 * Vertice[i * nMaxJ + j].y;
				Centre0.z = 2 * Vertice[i * nMaxJ + j].z;

				SpikeLength(Vertice, i, j, nMaxI, nMaxJ, &V);

				Centre1.x = Centre0.x + (float)V.x;
				Centre1.y = Centre0.y + (float)V.y;
				Centre1.z = Centre0.z + (float)V.z;


				if (Type == PRINT_STL_BINARY)
					nbObjects += PrintSTLConeBinary(Centre0, Centre1, fSphereRadius*ShrinkStruts(j, nMaxJ), 0, GetColor(sSphereMaterial), Dome, false);
				else if (Type == PRINT_OBJ)
					nbObjects += PrintOBJCone(Centre0, Centre1, fSphereRadius*ShrinkStruts(j, nMaxJ), 0, Dome, nbObjects, "Spikes Cone", sSphereMaterial, false);
				else if (Type == PRINT_POV)
					PrintPOVCone(Centre0, Centre1, fSphereRadius*ShrinkStruts(j, nMaxJ), 0.0f, Dome, sSphereMaterial);


			}
		}
	}
	return nbObjects;
}

int PrintSpike_deprecated(GLpoint *Vertice, FILE *Dome, int nMaxI, int nMaxJ, int nbObjects, PRINT_TYPE Type)
{
	GLpoint Centre0, Centre1;
	Vector3D V;
	float fSphereRadius;
	char *sSphereMaterial;

	float fSpikeMP = myDome.fSpikeMP;
	float fSpikeMp = myDome.fSpikeMp;
	float fSpikemP = myDome.fSpikemP;
	float fSpikemp = myDome.fSpikemp;
	float fSpike;

	if (Type == PRINT_STL_BINARY)
		nbObjects = 0;


	// Top Sphere							
	fSphereRadius = MaxRadius(&sSphereMaterial, SPHERE_TOP);
	if (fSphereRadius > 0.0)
	{
		if (!myDome.bWireGenerate)
		{
			Centre0.x = 2 * Vertice[0].x;
			Centre0.y = 2 * Vertice[0].y;
			Centre0.z = 2 * Vertice[0].z;

			if (Type == PRINT_STL_BINARY)
				nbObjects += PrintSTLSphereBinary(Centre0, fSphereRadius*ShrinkStruts(0, nMaxJ), GetColor(sSphereMaterial), Dome);
			else if (Type == PRINT_OBJ)
				nbObjects += PrintOBJSphere(Centre0, fSphereRadius*ShrinkStruts(0, nMaxJ), Dome, nbObjects, "Top Sphere Spike", sSphereMaterial);
			else if (Type == PRINT_POV)
				PrintPOVSphere(Centre0, fSphereRadius*ShrinkStruts(0, nMaxJ), Dome, sSphereMaterial);

		}
		Centre0.x = 2 * Vertice[0].x;
		Centre0.y = 2 * Vertice[0].y;
		Centre0.z = 2 * Vertice[0].z;

		Centre1.x = 2 * Vertice[0].x * (1.0f + fSpikeMP * ShrinkStruts(0, nMaxJ));
		Centre1.y = 2 * Vertice[0].y * (1.0f + fSpikeMP * ShrinkStruts(0, nMaxJ));
		Centre1.z = 2 * Vertice[0].z * (1.0f + fSpikeMP * ShrinkStruts(0, nMaxJ));
		if (fSpikeMP * ShrinkStruts(0, nMaxJ) > 0)

			if (Type == PRINT_STL_BINARY)
				nbObjects += PrintSTLConeBinary(Centre0, Centre1, fSphereRadius*ShrinkStruts(0, nMaxJ), 0, GetColor(sSphereMaterial), Dome);
			else if (Type == PRINT_OBJ)
				nbObjects += PrintOBJCone(Centre0, Centre1, fSphereRadius*ShrinkStruts(0, nMaxJ), 0, Dome, nbObjects, "Top Spikes Cone", sSphereMaterial);
			else if (Type == PRINT_POV)
				PrintPOVCone(Centre0, Centre1, fSphereRadius*ShrinkStruts(0, nMaxJ), 0.0f, Dome, sSphereMaterial);
	}

	for (int i = 0; i < nMaxI - 1; i++)
	{
		for (int j = 1; j < nMaxJ; j++)
		{
			//Spheres					
			sSphereMaterial = myDome.sSphereMPMaterial;
			fSphereRadius = 0.0f;
			if ((j % (myDome.nParallelMinor + 1) == 0) && (i % (myDome.nMeridianMinor + 1) == 0))
			{
				sSphereMaterial = myDome.sSphereMPMaterial;
				fSpike = fSpikeMP;
				fSphereRadius = (fSpike > 0.0f ? myDome.fSphereMPRadius : 0);
			}
			if ((j % (myDome.nParallelMinor + 1) != 0) && (i % (myDome.nMeridianMinor + 1) == 0))
			{
				sSphereMaterial = myDome.sSphereMpMaterial;
				fSpike = fSpikeMp;
				fSphereRadius = (fSpike > 0.0f ? myDome.fSphereMpRadius : 0);
			}
			if ((j % (myDome.nParallelMinor + 1) == 0) && (i % (myDome.nMeridianMinor + 1) != 0))
			{
				sSphereMaterial = myDome.sSpheremPMaterial;
				fSpike = fSpikemP;
				fSphereRadius = (fSpike > 0.0f ? myDome.fSpheremPRadius : 0);
			}
			if ((j % (myDome.nParallelMinor + 1) != 0) && (i % (myDome.nMeridianMinor + 1) != 0))
			{
				sSphereMaterial = myDome.sSpherempMaterial;
				fSpike = fSpikemp;
				fSphereRadius = (fSpike > 0.0f ? myDome.fSpherempRadius : 0);
			}
			if (fSphereRadius > 0.0)
			{
				if (!myDome.bWireGenerate)
				{
					Centre0.x = 2 * Vertice[i * nMaxJ + j].x;
					Centre0.y = 2 * Vertice[i * nMaxJ + j].y;
					Centre0.z = 2 * Vertice[i * nMaxJ + j].z;

					if (Type == PRINT_STL_BINARY)
						nbObjects += PrintSTLSphereBinary(Centre0, fSphereRadius*ShrinkStruts(j, nMaxJ), GetColor(sSphereMaterial), Dome);
					else if (Type == PRINT_OBJ)
						nbObjects += PrintOBJSphere(Centre0, fSphereRadius*ShrinkStruts(j, nMaxJ), Dome, nbObjects, "Sphere Spike", sSphereMaterial);
					else if (Type == PRINT_POV)
						PrintPOVSphere(Centre0, fSphereRadius*ShrinkStruts(j, nMaxJ), Dome, sSphereMaterial);
				}
				Centre0.x = 2 * Vertice[i * nMaxJ + j].x;
				Centre0.y = 2 * Vertice[i * nMaxJ + j].y;
				Centre0.z = 2 * Vertice[i * nMaxJ + j].z;

				if (CalculDerivee(Vertice, i, j, &V, nMaxI, nMaxJ) > 0 && fSpike > 0)
				{
					Centre1.x = 2 * Vertice[i * nMaxJ + j].x + (GLfloat)(fSpike*ShrinkStruts(j, nMaxJ) * V.x);
					Centre1.y = 2 * Vertice[i * nMaxJ + j].y + (GLfloat)(fSpike*ShrinkStruts(j, nMaxJ) * V.y);
					Centre1.z = 2 * Vertice[i * nMaxJ + j].z + (GLfloat)(fSpike*ShrinkStruts(j, nMaxJ) * V.z);


					if (Type == PRINT_STL_BINARY)
						nbObjects += PrintSTLConeBinary(Centre0, Centre1, fSphereRadius*ShrinkStruts(j, nMaxJ), 0, GetColor(sSphereMaterial), Dome);
					else if (Type == PRINT_OBJ)
						nbObjects += PrintOBJCone(Centre0, Centre1, fSphereRadius*ShrinkStruts(j, nMaxJ), 0, Dome, nbObjects, "Spikes Cone", sSphereMaterial);
					else if (Type == PRINT_POV)
						PrintPOVCone(Centre0, Centre1, fSphereRadius*ShrinkStruts(j, nMaxJ), 0.0f, Dome, sSphereMaterial);


				}
			}
		}
	}
	return nbObjects;
}

int PrintStrutsCone(GLpoint *Vertices, FILE *Dome, int nMaxI, int nMaxJ, int nbObjects, PRINT_TYPE Type)
{
	GLpoint Centre1, Centre2;
	char *sMeridianMaterial, *sParallelMaterial;
	float fMeridianRadius, fParallelRadius;

	if (Type == PRINT_STL_BINARY)
		nbObjects = 0;


	for (int i = 0; i < nMaxI - 1; i++)
	{
		fMeridianRadius = (i % (myDome.nMeridianMinor + 1) == 0) ? myDome.fMeridianMainRadius : myDome.fMeridianMinorRadius;
		sMeridianMaterial = ((i % (myDome.nMeridianMinor + 1) == 0) ? myDome.sMeridianMainMaterial : myDome.sMeridianMinorMaterial);

		for (int j = 1; j < nMaxJ; j++)
		{
			fParallelRadius = (j % (myDome.nParallelMinor + 1) == 0) ? myDome.fParallelMainRadius : myDome.fParallelMinorRadius;
			sParallelMaterial = (j % (myDome.nParallelMinor + 1) == 0) ? myDome.sParallelMainMaterial : myDome.sParallelMinorMaterial;
			//Meridians		
			Centre1.x = 2 * Vertices[i * nMaxJ + j].x;
			Centre1.y = 2 * Vertices[i * nMaxJ + j].y;
			Centre1.z = 2 * Vertices[i * nMaxJ + j].z;

			Centre2.x = 2 * Vertices[i * nMaxJ + j - 1].x;
			Centre2.y = 2 * Vertices[i * nMaxJ + j - 1].y;
			Centre2.z = 2 * Vertices[i * nMaxJ + j - 1].z;
			
			if (Type == PRINT_STL_BINARY)
				nbObjects += PrintSTLConeBinary(Centre1, Centre2, fMeridianRadius*ShrinkStruts(j, nMaxJ), fMeridianRadius*ShrinkStruts(j - 1, nMaxJ), GetColor(sMeridianMaterial), Dome);
			else if (Type == PRINT_OBJ)
				nbObjects += PrintOBJCone(Centre1, Centre2, fMeridianRadius*ShrinkStruts(j, nMaxJ), fMeridianRadius*ShrinkStruts(j - 1, nMaxJ), Dome, nbObjects, ((i % (myDome.nMeridianMinor + 1) == 0) ? "Main Meridian Cone" : "Minor Meridian Cone"), sMeridianMaterial);
			else if (Type == PRINT_POV)
				PrintPOVCone(Centre1, Centre2, fMeridianRadius*ShrinkStruts(j, nMaxJ), fMeridianRadius*ShrinkStruts(j - 1, nMaxJ), Dome, sMeridianMaterial);


			//Parallels		
			Centre2.x = 2 * Vertices[(i + 1)* nMaxJ + j].x;
			Centre2.y = 2 * Vertices[(i + 1)* nMaxJ + j].y;
			Centre2.z = 2 * Vertices[(i + 1)* nMaxJ + j].z;

			if (Type == PRINT_STL_BINARY)
				nbObjects += PrintSTLConeBinary(Centre2, Centre1, fParallelRadius*ShrinkStruts(j, nMaxJ), fParallelRadius*ShrinkStruts(j, nMaxJ), GetColor(sParallelMaterial), Dome);
			else if (Type == PRINT_OBJ)
				nbObjects += PrintOBJCone(Centre1, Centre2, fParallelRadius*ShrinkStruts(j, nMaxJ), fParallelRadius*ShrinkStruts(j, nMaxJ), Dome, nbObjects, (j % (myDome.nParallelMinor + 1) == 0) ? "Main Parallel Cone" : "Minor Parallel Cone", sParallelMaterial);
			else if (Type == PRINT_POV)
				PrintPOVCone(Centre1, Centre2, fParallelRadius*ShrinkStruts(j, nMaxJ), fParallelRadius*ShrinkStruts(j - 1, nMaxJ), Dome, sParallelMaterial);
			

			//Down Left to Up Right Diagonal		
			Centre2.x = 2 * Vertices[(i + 1)* nMaxJ + j - 1].x;
			Centre2.y = 2 * Vertices[(i + 1)* nMaxJ + j - 1].y;
			Centre2.z = 2 * Vertices[(i + 1)* nMaxJ + j - 1].z;

			if (Type == PRINT_STL_BINARY)
				nbObjects += PrintSTLConeBinary(Centre1, Centre2, myDome.fDiagonalURRadius*ShrinkStruts(j, nMaxJ), myDome.fDiagonalURRadius*ShrinkStruts(j - 1, nMaxJ), GetColor(myDome.sDiagonalURMaterial), Dome);
			else if (Type == PRINT_OBJ)
				nbObjects += PrintOBJCone(Centre1, Centre2, myDome.fDiagonalURRadius*ShrinkStruts(j, nMaxJ), myDome.fDiagonalURRadius*ShrinkStruts(j-1, nMaxJ), Dome, nbObjects, "Diagonal UR Cone", myDome.sDiagonalURMaterial);
			else if (Type == PRINT_POV)
				PrintPOVCone(Centre1, Centre2, myDome.fDiagonalURRadius*ShrinkStruts(j, nMaxJ), myDome.fDiagonalURRadius*ShrinkStruts(j - 1, nMaxJ), Dome, myDome.sDiagonalURMaterial);


			//Up Left to Down Right Diagonal	
			float fDistance = DistancePlan(Vertices[i * nMaxJ + j], Vertices[(i + 1)* nMaxJ + j - 1], Vertices[(i + 1) * nMaxJ + j], Vertices[i * nMaxJ + j - 1]);
			
				
			Centre1.x = 2 * Vertices[(i + 1) * nMaxJ + j].x;
			Centre1.y = 2 * Vertices[(i + 1) * nMaxJ + j].y;
			Centre1.z = 2 * Vertices[(i + 1) * nMaxJ + j].z;

			Centre2.x = 2 * Vertices[i * nMaxJ + j - 1].x;
			Centre2.y = 2 * Vertices[i * nMaxJ + j - 1].y;
			Centre2.z = 2 * Vertices[i * nMaxJ + j - 1].z;
			if (fDistance > (myDome.fDiagonalDRRadius*ShrinkStruts(j, nMaxJ) + myDome.fDiagonalDRRadius*ShrinkStruts(j - 1, nMaxJ)) / 2)
			{
				float Radius;
				GLpoint Centre3;

				Centre3.x = Vertices[(i + 1)* nMaxJ + j - 1].x + Vertices[i * nMaxJ + j].x;
				Centre3.y = Vertices[(i + 1)* nMaxJ + j - 1].y + Vertices[i * nMaxJ + j].y;
				Centre3.z = Vertices[(i + 1)* nMaxJ + j - 1].z + Vertices[i * nMaxJ + j].z;
				Radius = (myDome.fDiagonalDRRadius*ShrinkStruts(j, nMaxJ) + myDome.fDiagonalDRRadius*ShrinkStruts(j - 1, nMaxJ)) / 2.0f;
				if (Type == PRINT_STL_BINARY)
				{
					nbObjects += PrintSTLConeBinary(Centre1, Centre3, myDome.fDiagonalDRRadius*ShrinkStruts(j, nMaxJ), Radius, GetColor(myDome.sDiagonalDRMaterial), Dome);
					nbObjects += PrintSTLSphereBinary(Centre3, Radius, GetColor(myDome.sDiagonalDRMaterial), Dome);
					nbObjects += PrintSTLConeBinary(Centre3, Centre2, Radius, myDome.fDiagonalDRRadius*ShrinkStruts(j - 1, nMaxJ), GetColor(myDome.sDiagonalDRMaterial), Dome);
				}
				else if (Type == PRINT_OBJ)
				{
					nbObjects += PrintOBJCone(Centre1, Centre3, myDome.fDiagonalDRRadius*ShrinkStruts(j, nMaxJ), Radius, Dome, nbObjects, "Diagonal DR Cone", myDome.sDiagonalDRMaterial);
					nbObjects += PrintOBJSphere(Centre3, Radius, Dome, nbObjects, "Top Sphere", myDome.sDiagonalDRMaterial);
					nbObjects += PrintOBJCone(Centre3, Centre2, Radius, myDome.fDiagonalDRRadius*ShrinkStruts(j - 1, nMaxJ), Dome, nbObjects, "Diagonal DR Cone", myDome.sDiagonalDRMaterial);
				}
				else if (Type == PRINT_POV)
				{
					PrintPOVCone(Centre1, Centre3, myDome.fDiagonalDRRadius*ShrinkStruts(j, nMaxJ), Radius, Dome, myDome.sDiagonalDRMaterial);
					PrintPOVSphere(Centre3, Radius, Dome, myDome.sDiagonalDRMaterial);
					PrintPOVCone(Centre3, Centre2, Radius, myDome.fDiagonalDRRadius*ShrinkStruts(j - 1, nMaxJ), Dome, myDome.sDiagonalDRMaterial);
				}
			}
			else
			{
				if (Type == PRINT_STL_BINARY)
					nbObjects += PrintSTLConeBinary(Centre1, Centre2, myDome.fDiagonalDRRadius*ShrinkStruts(j, nMaxJ), myDome.fDiagonalDRRadius*ShrinkStruts(j - 1, nMaxJ), GetColor(myDome.sDiagonalDRMaterial), Dome);
				else if (Type == PRINT_OBJ)
					nbObjects += PrintOBJCone(Centre1, Centre2, myDome.fDiagonalDRRadius*ShrinkStruts(j, nMaxJ), myDome.fDiagonalDRRadius*ShrinkStruts(j - 1, nMaxJ), Dome, nbObjects, "Diagonal DR Cone", myDome.sDiagonalDRMaterial);
				else if (Type == PRINT_POV)
					PrintPOVCone(Centre1, Centre2, myDome.fDiagonalDRRadius*ShrinkStruts(j, nMaxJ), myDome.fDiagonalDRRadius*ShrinkStruts(j - 1, nMaxJ), Dome, myDome.sDiagonalDRMaterial);
			}
		}
	}
	return nbObjects;
}

int PrintStrutsSphere(GLpoint *Vertices, FILE *Dome, int nMaxI, int nMaxJ, int nbObjects, PRINT_TYPE Type)
{
	float fSphereRadius;
	char *sSphereMaterial;
	GLpoint Centre;

	char sTitre[100];

	if (Type == PRINT_STL_BINARY)
		nbObjects = 0;
	fSphereRadius = MaxRadius(&sSphereMaterial, SPHERE_TOP);
	if (fSphereRadius > 0.0)
	{
		Centre.x = 2 * Vertices[0].x;
		Centre.y = 2 * Vertices[0].y;
		Centre.z = 2 * Vertices[0].z;

		if (Type == PRINT_STL_BINARY)
			nbObjects += PrintSTLSphereBinary(Centre, fSphereRadius*ShrinkStruts(0, nMaxJ), GetColor(sSphereMaterial), Dome);
		else if (Type == PRINT_OBJ)
			nbObjects += PrintOBJSphere(Centre, fSphereRadius*ShrinkStruts(0, nMaxJ), Dome, nbObjects, "Top Sphere", sSphereMaterial);
		else if (Type == PRINT_POV)
			PrintPOVSphere(Centre, fSphereRadius*ShrinkStruts(0, nMaxJ), Dome, sSphereMaterial);
	}

	for (int i = 0; i < nMaxI - 1; i++)
	{
		for (int j = 1; j < nMaxJ; j++)
		{

			//Spheres			
			if ((j % (myDome.nParallelMinor + 1) == 0) && (i % (myDome.nMeridianMinor + 1) == 0))
			{
				fSphereRadius = MaxRadius(&sSphereMaterial, SPHERE_MP);
				strcpy_s(sTitre, 100, "Sphere MP");
			}
			if ((j % (myDome.nParallelMinor + 1) != 0) && (i % (myDome.nMeridianMinor + 1) == 0))
			{
				fSphereRadius = MaxRadius(&sSphereMaterial, SPHERE_Mp);
				strcpy_s(sTitre, 100, "Sphere Mp");
			}
			if ((j % (myDome.nParallelMinor + 1) == 0) && (i % (myDome.nMeridianMinor + 1) != 0))
			{
				fSphereRadius = MaxRadius(&sSphereMaterial, SPHERE_mP);
				strcpy_s(sTitre, 100, "Sphere mP");
			}
			if ((j % (myDome.nParallelMinor + 1) != 0) && (i % (myDome.nMeridianMinor + 1) != 0))
			{
				fSphereRadius = MaxRadius(&sSphereMaterial, SPHERE_mp);
				strcpy_s(sTitre, 100, "Sphere mp");
			}
			if (fSphereRadius > 0.0)
			{
				Centre.x = 2 * Vertices[i * nMaxJ + j].x;
				Centre.y = 2 * Vertices[i * nMaxJ + j].y;
				Centre.z = 2 * Vertices[i * nMaxJ + j].z;

				if (Type == PRINT_STL_BINARY)
					nbObjects += PrintSTLSphereBinary(Centre, fSphereRadius*ShrinkStruts(j, nMaxJ), GetColor(sSphereMaterial), Dome);
				else if (Type == PRINT_OBJ)
					nbObjects += PrintOBJSphere(Centre, fSphereRadius*ShrinkStruts(j, nMaxJ), Dome, nbObjects, "Top Sphere", sSphereMaterial);
				else if (Type == PRINT_POV)
					PrintPOVSphere(Centre, fSphereRadius*ShrinkStruts(j, nMaxJ), Dome, sSphereMaterial);
			}
		}
	}
	return nbObjects;
}


int PrintCloture(FILE *Dome, int nMaxI, int nMaxJ, GLpoint *Poly, int nbObjects, PRINT_TYPE Type)
{
	int i;
	GLpoint C0;
	int nCouleur = GetColor(myDome.sLidMaterial);
	char *sCoverMaterial = myDome.sLidMaterial;

	if (Type == PRINT_STL_BINARY)
		nbObjects = 0;

	C0.x = 0;
	C0.y = 0;
	C0.z = myDome.nWindow == WINDOW_NONE ? 0 : myDome.fWindowHeight / 2;

	if (!myDome.bLid)
		return nbObjects;
	if (Type == PRINT_OBJ)
	{
		fprintf(Dome, "o Mesh %s", CRLF);
		fprintf(Dome, "usemtl color_%s%s", myDome.sLidMaterial, CRLF);
		fprintf(Dome, "v %lf %lf %lf %s", 4*C0.x, 4*C0.y, 4*C0.z, CRLF);
	}

	if(Type == PRINT_POV)
		fprintf(Dome, "mesh {  // Cloture %s", CRLF);


	for (i = 0;i < nMaxI-1;i++)
	{
		if (Type == PRINT_STL_BINARY)
		{
			PrintSTLTriangle(Dome, C0, Poly[i * nMaxJ + nMaxJ - 1], Poly[((i + 1)) * nMaxJ + nMaxJ - 1], nCouleur);
			nbObjects++;
		}
		if (Type == PRINT_OBJ)
			fprintf(Dome, "v %lf %lf %lf %s", 4*Poly[i * nMaxJ + nMaxJ - 1].x, 4*Poly[i * nMaxJ + nMaxJ - 1].y, 4*Poly[i * nMaxJ + nMaxJ - 1].z, CRLF);
		if (Type == PRINT_POV)
			fprintf(Dome, "triangle { <%lf, %lf, %lf>, <%lf, %lf, %lf>, <%lf, %lf, %lf> }%s",
			C0.x, C0.z, C0.y,
			2 * Poly[i * nMaxJ + nMaxJ - 1].x, 2 * Poly[i * nMaxJ + nMaxJ - 1].z, 2 * Poly[i * nMaxJ + nMaxJ - 1].y,
			2 * Poly[((i + 1)) * nMaxJ + nMaxJ - 1].x, 2 * Poly[((i + 1)) * nMaxJ + nMaxJ - 1].z, 2 * Poly[((i + 1)) * nMaxJ + nMaxJ - 1].y,
			CRLF);

	}
	
	if (Type == PRINT_OBJ)
	{
		for (i = 0;i < nMaxI-1;i++)
			fprintf(Dome, "f %d %d %d %s", nbObjects + 1, nbObjects + i + 2, nbObjects + (i + 1)%(nMaxI-1) + 2, CRLF);

		nbObjects += nMaxI;
	}

	if (Type == PRINT_POV)
	{
		if (strcmp(myDome.sLidMaterial, S_NONE))
			fprintf(Dome, " material{%s}} %s", sCoverMaterial, CRLF);
		else
			fprintf(Dome, "}%s", CRLF);
	}

	return nbObjects;

}


unsigned long PrintColor(FILE *DomeSTL)
{
	GLpoint V1;
	unsigned long nbTriangles = 0;
	unsigned short Couleur;

	Couleur = atoi(myDome.sDomeName + 6);
	V1.x = -1;
	V1.y = 3;
	V1.z = 0;
	nbTriangles += PrintSTLSphereBinary(V1, .5, Couleur, DomeSTL);
	for(int i = 0;i<8;i++)
		for (int j = 0; j < 16;j++)
		{
			V1.x = (float)i;
			V1.y = (float)j;
			V1.z = 0;
			Couleur = i * 4 + 32 * j * 2 + Couleur;
			nbTriangles += PrintSTLSphereBinary(V1, .5, Couleur, DomeSTL);
		}

	return nbTriangles;

}


int PrintDome(GLpoint *Vertice, FILE *Dome, int nMaxI, int nMaxJ, PRINT_TYPE Type)
{
	clock_t StartProcess, EndProcess;
	double duration;

	bool bSpike = myDome.fSpikeMP > 0 || myDome.fSpikeMp > 0 || myDome.fSpikemP > 0 || myDome.fSpikemp > 0;

	int Start = 0;

	switch (Type)
	{
	case PRINT_OBJ:
		printf("OBJ ");
		break;
	case PRINT_POV:
		printf("POV ");
		break;
	case PRINT_STL_BINARY:
		printf("STL ");
		break;
	}
	printf("File\nDur%ce : Sphere      Cones       Faces       Spikes      Spirales      Frises      Snakes    Combs\n", 130);

	StartProcess = clock();
	if (myDome.bWireGenerate)
		Start = PrintStrutsSphere(Vertice, Dome, nMaxI, nMaxJ, Start, Type);
	EndProcess = clock();
	duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
	printf("      %8.3lf", duration);

	StartProcess = clock();
	if (myDome.bWireGenerate)
		Start = PrintStrutsCone(Vertice, Dome, nMaxI, nMaxJ, Start, Type);
	EndProcess = clock();
	duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
	printf("   %8.3lf", duration);

	StartProcess = clock();
	if (myDome.bFaceGenerate)
	{
		switch (Type)
		{
		case PRINT_OBJ:
			Start = PrintOBJMesh(Dome, Vertice, Start, nMaxI, nMaxJ);
			break;
		case PRINT_POV:
			printf("POV ");
			break;
		case PRINT_STL_BINARY:
			printf("STL ");
			break;
		}

	}
		
	EndProcess = clock();
	duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
	printf("       %5.3lf", duration);

	StartProcess = clock();
	if (bSpike)
		Start = PrintSpike(Vertice, Dome, nMaxI, nMaxJ, Start, Type);
	EndProcess = clock();
	duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
	printf("     %8.3lf", duration);

	StartProcess = clock();
	if (myDome.nSpirale > 0)
		Start = PrintSpirale(Dome, Start, Type);
	EndProcess = clock();
	duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
	printf("      %8.3lf", duration);

	StartProcess = clock();
	if (myDome.nFrise > 0)
		Start = PrintFrise(Dome, Start, Type);
	EndProcess = clock();
	duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
	printf("    %8.3lf", duration);

	StartProcess = clock();
	if (myDome.nSnake > 0)
		PrintSnake(Dome, Start, Type);
	EndProcess = clock();
	duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
	printf("    %8.3lf\n", duration);

	StartProcess = clock();
	if (myDome.nComb > 0)
		PrintComb(Dome, Start, Type);
	EndProcess = clock();
	duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
	printf("    %8.3lf\n", duration);

	return 1;
}


GLpoint SphereFormula(GLpoint Centre, double Rayon, double u, double v)
{
	GLpoint Retour;

	Retour.x = (GLfloat)(Centre.x + Rayon * cos(u) * sin(v));
	Retour.y = (GLfloat)(Centre.y + Rayon * sin(u) * sin(v));
	Retour.z = (GLfloat)(Centre.z + Rayon * cos(v));

	return Retour;
}




GLpoint ConeFormula(GLpoint Centre1, GLpoint Centre2, double Rayon1, double Rayon2, double u, double v, bool bStraight)
{
	GLpoint Retour;
	GLpoint AB, v0, v1, v2;
	float R;
	AB.x = Centre2.x - Centre1.x;
	AB.y = Centre2.y - Centre1.y;
	AB.z = Centre2.z - Centre1.z;
	R = sqrtf(AB.x*AB.x + AB.y*AB.y + AB.z*AB.z);

	v0.x = AB.x / R;
	v0.y = AB.y / R;
	v0.z = AB.z / R;

	if (v0.x == 0 && v0.y == 0)
	{
		v1.x = 1;
		v1.y = 0;
		v1.z = 0;
	}
	else
	{
		v1.x = v0.y;
		v1.y = -v0.x;
		v1.z = 0;
		R = sqrtf(v1.x*v1.x + v1.y*v1.y);
		v1.x /= R;
		v1.y /= R;
	}
	CrossProduct(&v0, &v1, &v2);

	if(bStraight)
		R = (float)(Rayon1 + v * (Rayon2 - Rayon1));
	else
		R = (float)(Rayon2 + (Rayon2 - Rayon1)*(6 * v*v*v*v*v - 15 * v*v*v*v + 10 * v*v*v - 1));

//		R = (float)(Rayon2 + (Rayon1 - Rayon2)*(2.0f*v*v*v - 3.0f*v*v + 1));

	Retour.x = (GLfloat)(Centre1.x + R * (v1.x*cos(u) + v2.x*sin(u)) + AB.x*v);
	Retour.y = (GLfloat)(Centre1.y + R * (v1.y*cos(u) + v2.y*sin(u)) + AB.y*v);
	Retour.z = (GLfloat)(Centre1.z + R * (v1.z*cos(u) + v2.z*sin(u)) + AB.z*v);

	return Retour;

}



int ChargeCouleur()

{
	FILE *ColorFile;
	errno_t Erreur;
	char  FileName[400];
	char sBuffer[1200];
	char *sRetour;
	char sCouleur[120];
	int  nbCouleur;
	strcpy_s(FileName, "..\\");
	strcat_s(FileName, "Colors.mtl");

	Erreur = fopen_s(&ColorFile, FileName, "r");
	if (Erreur != 0)
	{
		Message("     Error Creating Color File   ");
		printf("Erreur = %d\n", Erreur);
		return 0;
	}
	nbCouleur = 0;
	do
	{
		strcpy_s(sBuffer, 10, "#Fin");
		sRetour = fgets(sBuffer, 1200, ColorFile);
		strcpy_s(sCouleur, "NULL");
		if (sBuffer[0] != '#')
			if (!strncmp(sBuffer, "newmtl color_", strlen("newmtl color_")))
			{
				strcpy_s(sCouleur, sBuffer + strlen("newmtl color_"));
				strcpy_s(ListeCouleur[nbCouleur], sCouleur);
//				printf("%s\n", sCouleur);
				nbCouleur++;
			}
	} while (sRetour);

	fclose(ColorFile);
	return nbCouleur;
}
	


double SlopeFunction_old(double v, double fParam, double dAxe)
{
	double Socle = myDome.fDomeSocle;

	if (v > 1 - Socle)
		return (dAxe == SLOPE_XY ? 1 : v);

	v = (v) / (1 - Socle);
	if (dAxe == SLOPE_XY)
		return SlopeFunctionXY(v, fParam);
	else
		return (1 - Socle) * SlopeFunctionZ(v, fParam);

}


double SlopeFunction(double v, double fParam, double dAxe)
{
	double Socle = myDome.fDomeSocle;
	double dPart = myDome.fTopPart;
	double dParam = myDome.fTopParameter;

	double Retour = 0;
	int nTop = myDome.nTopType;

	double alpha = 3;
	double beta = myDome.fTopParameter; //2
	double g, h;

//	double dPart = myDome.fWindowShift; //entre 0 et 1;
	double dFootXY = SlopeFunctionXY(dPart, myDome.fSlopeParameterXY);
	double dFootZ = 1 - SlopeFunctionZ(dPart, myDome.fSlopeParameterZ);

	double dMiddle;
	double dLength;

	double a, b, c, d;
	double dSplit;
	double Rayon = dFootXY * dParam / 3;

	double Oxy = dFootXY - Rayon;
	double Oz = dFootZ;

	double w;

	switch (nTop)
	{
	case TOP_NONE:
		dPart = 0;
		dMiddle = 0;
		break;
	case TOP_ARROW:
		dPart = myDome.fTopPart;
		dMiddle = 0;
		break;
	case TOP_LANTERN:
		dPart = myDome.fTopPart;
		dMiddle = dPart / 2;
		break;
	case TOP_MERINGUE:
		dPart = myDome.fTopPart;
		dMiddle = dPart / 2;
		break;
	case TOP_FOOT:
		dPart = myDome.fTopPart;
		dLength = dParam * dFootXY + sqrt((1 - dFootZ) * (1 - dFootZ) + ((1 - dParam) * (1 - dParam) * dFootXY * dFootXY));
		dMiddle = dPart * dParam * dFootXY / dLength;
		break;
	}
	w = (_PI / 2) * (v - dPart) / (dMiddle - dPart);

	/* Arrow */
//	dHight = myDome.fWindowDelay; // entre 0 et 3

	if (v > 1 - Socle)
		Retour = (dAxe == SLOPE_XY ? 1 : v);
	else if (v < dMiddle)
	{
		switch (nTop)
		{
		case TOP_LANTERN:
			h = (dPart - 2 * v) / dPart;
			g = pow(h, alpha);
			Retour = (dAxe == SLOPE_XY ? (1 - h) * SlopeFunctionXY(dPart, fParam) : -g * beta * dPart + SlopeFunctionZ(dPart, fParam) - beta * dPart);
			break;
		case TOP_FOOT:
			Retour = (dAxe == SLOPE_XY ? v * dParam * dFootXY / (dMiddle) : 0);
			break;
		case TOP_MERINGUE:
			Retour = (dAxe == SLOPE_XY ? v * (Oxy) / dMiddle : 1 - (-v * (1 - Oz) / dMiddle + 1 + Rayon));
			break;
		}
	}
	else if (v < dPart)
	{
		switch (nTop)
		{
		case TOP_ARROW:
			Retour = (dAxe == SLOPE_XY ? v * dFootXY / dPart : v * (dParam - dFootZ) / dPart + 1 - dParam);
			break;
		case TOP_LANTERN:
			Retour = (dAxe == SLOPE_XY ? SlopeFunctionXY(dPart, fParam) : SlopeFunctionZ(dPart, fParam) - 2 * (dPart - v) * beta);
			break;
		case TOP_FOOT:
			dSplit = dMiddle / dPart;
			a = dFootXY * (dParam - 1) / (dPart * (dSplit - 1));
			b = dFootXY * (dSplit - dParam) / (dSplit - 1);
			c = (1 - dFootZ) / (dPart * (1 - dSplit));
			d = (dFootZ - 1) * dSplit / (1 - dSplit);
			//		printf("a = %f   b = %f   c = %f   d = %f\n", a, b, c, d);
			Retour = (dAxe == SLOPE_XY ? a * v + b : c * v + d);
			break;
		case TOP_MERINGUE:
			Retour = (dAxe == SLOPE_XY ? Rayon * cos(w) + Oxy : 1 - (Rayon * sin(w) + Oz));
			break;
		}
	}
	else
	{
		v = (v) / (1 - Socle);

		if (dAxe == SLOPE_XY)
			Retour = SlopeFunctionXY(v, fParam);
		else
			Retour = (1 - Socle) * SlopeFunctionZ(v, fParam);
	}

	double fTopMax;
	switch(nTop)
	{
	case TOP_ARROW:
		fTopMax = dParam;
		break;
	case TOP_LANTERN:
		fTopMax = 1 + beta * dPart - SlopeFunctionZ(dPart, fParam) + beta * dPart;
		break;
	case TOP_NONE:
	case TOP_FOOT:
		fTopMax = 1;
		break;
	case TOP_MERINGUE:
		fTopMax =  (1 + Rayon);
		break;
	}


	float fSeuil = 1 - myDome.fDomeInFolding;
	if (isnan(Retour))
		return 0;
	else if (1 - Retour >  fSeuil * fTopMax && dAxe == SLOPE_Z && true)
		return 1-(2 * fSeuil * fTopMax - 1 + Retour); // tenir compte du top

	else
		return Retour;

}

double BaseFunctionOld(double u, double a, double b, double dAxe, double w)
{

	double Round = myDome.fBaseRound;
	double Inter = myDome.fBaseBaryCentre;
	double v;
	double Alpha;
	double Coefficient;
	double X, Y, X1, Y1;
	double dMorphing;
	double dMorphingMax = myDome.fBaseMorphing;

	int nSideNumber;
	double Uc;

	static int nn = 0;
	double XR, YR;
	if (dMorphingMax == 0)
		dMorphing = 1;
	else if (w < dMorphingMax)
		dMorphing = 0;
	else
		dMorphing = (w - dMorphingMax) / (1 - dMorphingMax);




	
	XR = BaseFunctionX(u, a, b);
	YR = BaseFunctionY(u, a, b);




	BaseFunctionSecondaryX = BaseFunctionListe[myDome.nBaseSecondary].FunctionX;
	BaseFunctionSecondaryY = BaseFunctionListe[myDome.nBaseSecondary].FunctionY;

//	BaseFunctionSecondaryX = BaseCircle_X;
//	BaseFunctionSecondaryY = BaseCircle_Y;

//	X = dMorphing * BaseFunctionX(u, a, b) + (1 - dMorphing) * BaseFunctionSecondaryX(u, a, b);
//	Y = dMorphing * BaseFunctionY(u, a, b) + (1 - dMorphing) * BaseFunctionSecondaryY(u, a, b);

	double fraction = fmod((double)myDome.nBaseSides * u, 2 * _PI);

	double XS, YS;
	double nStruts = (myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor * myDisplay.bDisplayMinor);
	nSideNumber = (int)round((((double)myDome.nBaseSides * u - fraction) / (2 * _PI)));
	nSideNumber = (int)floor(((double)myDome.nBaseSides * u + (_PI / a)) / (2 * _PI));
	if (nSideNumber % 2)
		Uc = (nSideNumber + 1) * 2 * _PI / a;
	else
		Uc = 2 * _PI * (nSideNumber * nStruts + 2 * (u * nStruts * (double)myDome.nBaseSides / (2 * _PI) - nSideNumber * nStruts)) / ((double)myDome.nBaseSides * nStruts);
	if (w == 0 && dAxe == 17)
	{
		nn++;
		if (fabs(fabs(fraction - _PI) - _PI) < _PI / (2 * a))
			printf("\n");
		printf("%3d; %8.7f;%8.7f; %8.7f\n", nSideNumber, u, fraction, Uc);
	}

	if ((int)a % 2 == 0)
	{
		XS = BaseFunctionSecondaryX(fmod(Uc + 2 * _PI * myDome.fCrenelFront, 2 * _PI), a / 2, b);
		YS = BaseFunctionSecondaryY(fmod(Uc + 2 * _PI * myDome.fCrenelFront, 2 * _PI), a / 2, b);
	}
	else
	{
		XS = BaseFunctionSecondaryX(fmod(u + 2 * _PI * myDome.fCrenelFront, 2 * _PI), a, b);
		YS = BaseFunctionSecondaryY(fmod(u + 2 * _PI * myDome.fCrenelFront, 2 * _PI), a, b);
	}


	X = (dMorphing)*XR + (1 - dMorphing) * XS;
	Y = (dMorphing)*YR + (1 - dMorphing) * YS;


	if(Round == 0 && Inter == 0)
		return (dAxe == BASE_X ? X : Y);

	if (Round == 0 && Inter != 0)
		return (dAxe == BASE_X ? ((1 - Inter) * X + Inter * BaseFunctionSecondaryX(u, a, b)) :
			                     ((1 - Inter) * Y + Inter * BaseFunctionSecondaryY(u, a, b)));

	v = floor(a * u / (2 * _PI)) * 2 * _PI + _PI * (fmod(a * u, 2 * _PI) - Round) / (_PI - Round);
	v = v / a;
	Alpha = Round * (fmod(a * u, 2 * _PI) - _PI) / (a*(Round - _PI));
	Coefficient = 1 / sqrt(BaseFunctionX(0, a, b) * BaseFunctionX(0, a, b) + BaseFunctionY(0, a, b) * BaseFunctionY(0, a, b));



	if ((fmod(a * u, 2 * _PI) < Round) || (fmod(a * u, 2 * _PI) > 2 * _PI - Round))
		return (dAxe == BASE_X ? BaseFunctionSecondaryX(u, a, b) : BaseFunctionSecondaryY(u, a, b));
	
	
	X = BaseFunctionX(v, a, b);
	Y = BaseFunctionY(v, a, b);


	X1 = X * cos(Alpha) - Y * sin(Alpha);;
	Y1 = X * sin(Alpha) + Y * cos(Alpha);

	X1 *= Coefficient;
	Y1 *= Coefficient;
	if (Inter == 0)
		return (dAxe == BASE_X ? X1 : Y1);
	else
		return (dAxe == BASE_X ? ((1 - Inter) * X1 + Inter * BaseFunctionSecondaryX(u, a, b)) :
			                     ((1 - Inter) * Y1 + Inter * BaseFunctionSecondaryY(u, a, b)));

}

double BaseFunction(double u, double a, double b, double dAxe, double w)
{

	double Round = myDome.fBaseRound;
	double Inter = myDome.fBaseBaryCentre;
	double v;
	double Alpha;
	double Coefficient;
	double X, Y, X1, Y1;
	double dMorphing;
	double dMorphingMax = myDome.fBaseMorphing;
	int nSideNumber;
	double Uc;
	double fRotation = myDome.fBaseSecondaryRotation;
	double fExpansion = myDome.fBaseSecondaryExpansion;
	bool bDivision = myDome.bBaseMorphingDivision;

	static int nn = 0;
	double XR, YR;
	if (dMorphingMax == 0)
		dMorphing = 1;
	else if (w <= dMorphingMax)
		dMorphing = 0;
	else
		dMorphing = (w - dMorphingMax) / (1 - dMorphingMax);


//	printf("(%8.7lf, %8.7lf)    (%8.7lf, %8.7lf)\n", BaseLotus_X(u, a, b), BaseLotus_Y(u, a, b), BaseLotus(u, a, b, BASE_X), BaseLotus(u, a, b, BASE_Y));
	BaseFunctionSecondaryX = BaseFunctionListe[myDome.nBaseSecondary].FunctionX;
	BaseFunctionSecondaryY = BaseFunctionListe[myDome.nBaseSecondary].FunctionY;


	double XS, YS, Xc, Yc;

	XS = fExpansion * BaseFunctionSecondaryX(fmod(u + _PI * fRotation, 2 * _PI), a, b);
	YS = fExpansion * BaseFunctionSecondaryY(fmod(u + _PI * fRotation, 2 * _PI), a, b);
	Xc = fExpansion * BaseFunctionSecondaryX(u, a, b);
	Yc = fExpansion * BaseFunctionSecondaryY(u, a, b);

	double alpha = 2 * _PI * fRotation;

	XS = cos(alpha) * Xc - sin(alpha) * Yc;
	YS = sin(alpha) * Xc + cos(alpha) * Yc;

	if ((int)a % 2 == 0 && bDivision)
	{
//		double nStruts = (myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor * myDisplay.bDisplayMinor);
		nSideNumber = (int)floor(((double)myDome.nBaseSides * u + (_PI / a)) / (2 * _PI));
		if (nSideNumber % 2)
			Uc = (nSideNumber + 1) * 2 * _PI / a;
		else
//			Uc = 2 * _PI * (nSideNumber * nStruts + 2 * (u * nStruts * (double)myDome.nBaseSides / (2 * _PI) - nSideNumber * nStruts)) / ((double)myDome.nBaseSides * nStruts);
			Uc = (2 * u * (double)myDome.nBaseSides - 2 * _PI * nSideNumber) / (double)myDome.nBaseSides;


		if (w == 0 && dAxe == 17)
		{
			double fraction = fmod((double)myDome.nBaseSides * u, 2 * _PI);

			nn++;
			if (fabs(fabs(fraction - _PI) - _PI) < _PI / (2 * a))
				printf("\n");
			printf("%3d; %8.7f;%8.7f; %8.7f\n", nSideNumber, u, fraction, Uc);
		}

		XR = BaseFunctionX(Uc, a/2, b);
		YR = BaseFunctionY(Uc, a/2, b);
	}
	else
	{
		XR = BaseFunctionX(u, a, b);
		YR = BaseFunctionY(u, a, b);
	}


	X = (dMorphing)*XR + (1 - dMorphing) * XS;
	Y = (dMorphing)*YR + (1 - dMorphing) * YS;

/****************** Gestion de la spiralisation **********************/
	double c = myDome.fBaseSpiral / 4;
	if (c > 0)
	{
		double Ux = atan2(Y, X);
		double Module = sqrt(X * X + Y * Y) * cos(a * c * fmod(u, 2 * _PI / a));

		X = cos(Ux) * Module;
		Y = sin(Ux) * Module;
	}
/********************************************************************/

	if (Round == 0 && Inter == 0)
		return (dAxe == BASE_X ? X : Y);

	if (Round == 0 && Inter != 0)
		return (dAxe == BASE_X ? ((1 - Inter) * X + Inter * BaseFunctionSecondaryX(u, a, b)) :
		                         ((1 - Inter) * Y + Inter * BaseFunctionSecondaryY(u, a, b)));

	v = floor(a * u / (2 * _PI)) * 2 * _PI + _PI * (fmod(a * u, 2 * _PI) - Round) / (_PI - Round);
	v = v / a;
	Alpha = Round * (fmod(a * u, 2 * _PI) - _PI) / (a * (Round - _PI));
	Coefficient = 1 / sqrt(BaseFunctionX(0, a, b) * BaseFunctionX(0, a, b) + BaseFunctionY(0, a, b) * BaseFunctionY(0, a, b));



	if ((fmod(a * u, 2 * _PI) < Round) || (fmod(a * u, 2 * _PI) > 2 * _PI - Round))
		return (dAxe == BASE_X ? BaseFunctionSecondaryX(u, a, b) : BaseFunctionSecondaryY(u, a, b));


	X = BaseFunctionX(v, a, b);
	Y = BaseFunctionY(v, a, b);


	X1 = X * cos(Alpha) - Y * sin(Alpha);;
	Y1 = X * sin(Alpha) + Y * cos(Alpha);

	X1 *= Coefficient;
	Y1 *= Coefficient;
	if (Inter == 0)
		return (dAxe == BASE_X ? X1 : Y1);
	else
		return (dAxe == BASE_X ? ((1 - Inter) * X1 + Inter * BaseFunctionSecondaryX(u, a, b)) :
			((1 - Inter) * Y1 + Inter * BaseFunctionSecondaryY(u, a, b)));

}


void Sort(Tri** Liste, int nbElements)
{
	Tri* sTemp;

	for (int i = 0; i < nbElements - 1; i++)
		for (int j = i + 1; j < nbElements; j++)
			if (strcmp(Liste[i]->Mot, Liste[j]->Mot) > 0)
			{
				sTemp = Liste[i];
				Liste[i] = Liste[j];
				Liste[j] = sTemp;
			}

}
void Sort(char** Liste, int nbElements)
{
	char* sTemp;



	for (int i = 0; i < nbElements - 1; i++)
		for (int j = i + 1; j < nbElements; j++)
			if (strcmp(Liste[i], Liste[j]) > 0)
			{
				sTemp = Liste[i];
				Liste[i] = Liste[j];
				Liste[j] = sTemp;
			}
}

void Calcul(double u, double v, double *dX, double *dY, double *dZ, bool bDerivee)
{
	double dRayon;
	double X, Y, Z;
	double dSinus, dCosinus;
	double fWindowHeight;
	static double fMaxNorme = 0;
	double fNorme;
	double Z1;
	double fLower = myDome.fWindowLower;
//	static double MinU = 100, MaxU = -100, MinV = 100, MaxV = -100;
	 // 0 <= u <= 2pi
	 // 0 <= v <= 1

	if (u < MinU)MinU = u;
	if (u > MaxU)MaxU = u;
	if (v < MinV)MinV = v;
	if (v > MaxV)MaxV = v;

//	X = BaseFunctionX(u, myDome.nBaseSides, myDome.fBaseParameter);
//	Y = BaseFunctionY(u, myDome.nBaseSides, myDome.fBaseParameter);

	X = BaseFunction(u, myDome.nBaseSides, myDome.fBaseParameter, BASE_X, v);
	Y = BaseFunction(u, myDome.nBaseSides, myDome.fBaseParameter, BASE_Y, v);

	if (u == 0 && v == 0)
		fMaxNorme = 0;


	dRayon = sqrt(X * X + Y * Y);

	if (dRayon != 0)
	{
		X *= pow(dRayon, myDome.fBaseDamping - 1);
		Y *= pow(dRayon, myDome.fBaseDamping - 1);
	}

	float R = 1; //fabs(cos((1-v)*_PI/2));
	if (dRayon != 0)
	{
		X *= ((dRayon - R)*myDome.fBaseExpand + R) / dRayon;
		Y *= ((dRayon - R)*myDome.fBaseExpand + R) / dRayon;
	}



	if (dRayon != 0)
	{
		X *= (1 + (dRayon - 1) * (1 - (1 - v) * myDome.fBaseSmoothing)) / dRayon;
		Y *= (1 + (dRayon - 1) * (1 - (1 - v) * myDome.fBaseSmoothing)) / dRayon;
	}

	double pct = myDome.fSlopeDelay;
	double alpha;

	if (myDome.nSlopeTwist != 0)
	{
		if (pct > 0)
			alpha = (1 - v > pct ? (1 - v - pct) / (1 - pct) : 0);
		else
			alpha = (1 - v < 1 + pct ? 1 - v : 1 + pct);

		dCosinus = cos((1.0 - v) * myDome.nSlopeTwist * deg2rad);
		dSinus   = sin((1.0 - v) * myDome.nSlopeTwist * deg2rad);
		dCosinus = cos(alpha * myDome.nSlopeTwist * deg2rad);
		dSinus   = sin(alpha * myDome.nSlopeTwist * deg2rad);

		*dX = dCosinus * X * SlopeFunction(v, myDome.fSlopeParameterXY, SLOPE_XY) - dSinus   * Y * SlopeFunction(v, myDome.fSlopeParameterXY, SLOPE_XY);
		*dY = dSinus   * X * SlopeFunction(v, myDome.fSlopeParameterXY, SLOPE_XY) + dCosinus * Y * SlopeFunction(v, myDome.fSlopeParameterXY, SLOPE_XY);
	}
	else
	{
		if (myDome.nBase != BASE_ROOF)
			*dX = X * SlopeFunction(v, myDome.fSlopeParameterXY, SLOPE_XY);
//			*dX = X * pow(SlopeFunction(v, myDome.fSlopeParameterXY, SLOPE_XY), myDome.fWindowHeight);
		else
			*dX = X;

		*dY = Y * SlopeFunction(v, myDome.fSlopeParameterXY, SLOPE_XY);
//		*dY = Y * pow(SlopeFunction(v, myDome.fSlopeParameterXY, SLOPE_XY), myDome.fWindowHeight);

	}



//	*dX = ;
//	*dY = ;

	fWindowHeight = ComputeWindowCoefficient((float)(u + 2 * myDome.fWindowShift * _PI / myDome.nBaseSides), (float)myDome.nBaseSides, (float)myDome.fWindowPerSide) * myDome.fWindowHeight;

	if (v <= myDome.fWindowDamping)
		fWindowHeight = 0;
	else
		fWindowHeight *= (v - myDome.fWindowDamping) / ((1 - myDome.fWindowDamping) * v);


	if (fLower == 0 && myDome.bWindowBlind)
	{
		if (v == 1 && myDome.bWindowBlind) //£
			fWindowHeight = 0;
	}
	else if (v > 1 - fLower)
		fWindowHeight = 0;

	if (myDome.bWindowInverse)
		Z = (1.0 - SlopeFunction(v, myDome.fSlopeParameterZ, SLOPE_Z) * (1.0 - fWindowHeight * (1.0 - WindowFunctionZ(u + 2 * myDome.fWindowShift * _PI / myDome.nBaseSides, myDome.nBaseSides, myDome.fWindowPerSide, myDome.fWindowPower, false))));
	else
		Z = (1.0 - SlopeFunction(v, myDome.fSlopeParameterZ, SLOPE_Z) * (1.0 - fWindowHeight * WindowFunctionZ(u + 2 * myDome.fWindowShift * _PI / myDome.nBaseSides, myDome.nBaseSides, myDome.fWindowPerSide, myDome.fWindowPower, false)));

	Z1 = abs((1.0 - SlopeFunction(v, myDome.fSlopeParameterZ, SLOPE_Z)) - Z);

	//		*dZ = (1.0 - SlopeFunctionZ(v))*myDome.fDomeScaleZ;

/*	if (myDome.fNoiseAmplitudeR != 0 || myDome.fNoiseAmplitudeZ != 0)
	{
		float   w1, w2;
		double  c1, c2;

		c1 = myDome.fNoiseAmplitudeR * (1 - (1 - v) * myDome.fNoiseProgressiveR);
		c2 = myDome.fNoiseAmplitudeZ * (1 - (1 - v) * myDome.fNoiseProgressiveZ);

		w1 = (float)fabs(sin(myDome.fNoiseFrequencyU * u / 2));
		w2 = (float)mathfmod(myDome.fNoiseFrequencyV * v, 1);

		*dX *= (1 + c1 * PerlinNoiseGet2D(w1, w2));
		*dY *= (1 + c1 * PerlinNoiseGet2D(w1, w2));

		*dZ = Z * (1 + c2 * PerlinNoiseGet2D(w1, w2));

	}
	else*/
		
	*dZ = Z;
		
	
	*dX *= (myDome.fWindowAwning * Z1 + 1);
	*dY *= (myDome.fWindowAwning * Z1 + 1);



	double dThick = myDome.fCrenelThick;
	double dFront = myDome.fCrenelFront;
	double dHeight = myDome.fCrenelHeight;
	double dLargeur = myDome.fCrenelWide;
	int nNombre = myDome.nCrenelPerSide;
	double dPuissance = myDome.fCrenelPower;
	double dDelay = myDome.fCrenelDelay;

	double x;
	double dPush = 0.0;
	if (v >= 1 - (dThick + dFront) && v <= 1 - dFront)
	{
		double w;
		w = (v - (1 - (dThick + dFront))) / dThick;
		dPush = dHeight * pow(4 * (w - w * w), 0.1);
		dPush = dHeight * pow(4 * (w - w * w), dPuissance);
	}

	double uu = nNombre*myDome.nBaseSides*(u+dDelay)/(2*_PI);
	double d, f;
	
	f = modf(uu, &d);
//	printf("\n");

	x = (1 - dLargeur) / 2;

	d = (f - x) / (1 - 2 * x);

	if (nNombre == -1)
		dPush *= 1.0;
	else if (f < x || f > 1 - x)
		dPush = 0;
	else
		dPush *= pow(4 * (d - d * d), dPuissance);

	double phi = atan2(*dY, *dX);
	double psi = myDome.fCrenelAngle * _PI;
	if (dPush != 0)
	{
		double x0, y0, x1, y1;
		double epsilon = 0.000001;
		x0 = BaseFunction(u - epsilon, myDome.nBaseSides, myDome.fBaseParameter, BASE_X, v);
		y0 = BaseFunction(u - epsilon, myDome.nBaseSides, myDome.fBaseParameter, BASE_Y, v);
		x1 = BaseFunction(u + epsilon, myDome.nBaseSides, myDome.fBaseParameter, BASE_X, v);
		y1 = BaseFunction(u + epsilon, myDome.nBaseSides, myDome.fBaseParameter, BASE_Y, v);
		phi = _PI+atan2(x1 - x0, y0 - y1);
	}
	*dX += dPush * cos(phi)*cos(psi);
	*dY += dPush * sin(phi)*cos(psi);
	*dZ += dPush*sin(psi);

	*dZ *= myDome.fDomeScaleZ;

	*dX *= .5;
	*dY *= .5;
	*dZ *= .5;
	if (bDerivee && myDome.fDerivative != 0 )
	{
		double u1, v1;
		double fStep = 0.00001;
		Vector3D C0, C1;
		Vector3D V0, V1, V2;

		C0.x = *dX;
		C0.y = *dY;
		C0.z = *dZ;

		if (u + fStep > 2*_PI)
		{
			u1 = u - fStep;
			Calcul(u1, v, &C1.x, &C1.y, &C1.z, false);
			V0.x = (C0.x - C1.x) / fStep;
			V0.y = (C0.y - C1.y) / fStep;
			V0.z = (C0.z - C1.z) / fStep;
		}
		else
		{
			u1 = u + fStep;
			Calcul(u1, v, &C1.x, &C1.y, &C1.z, false);
			V0.x = (C1.x - C0.x) / fStep;
			V0.y = (C1.y - C0.y) / fStep;
			V0.z = (C1.z - C0.z) / fStep;
		}
		if (v + fStep > 1)
		{
			v1 = v - fStep;

			Calcul(u, v1, &C1.x, &C1.y, &C1.z, false);
			V1.x = (C0.x - C1.x) / fStep;
			V1.y = (C0.y - C1.y) / fStep;
			V1.z = (C0.z - C1.z) / fStep;
		}
		else
		{
			v1 = v + fStep;
			Calcul(u, v1, &C1.x, &C1.y, &C1.z, false);
			V1.x = (C1.x - C0.x) / fStep;
			V1.y = (C1.y - C0.y) / fStep;
			V1.z = (C1.z - C0.z) / fStep;
		}
		/*		if (V0.x*V0.x + V0.y*V0.y + V0.z*V0.z == 0)
			V2 = V1;
		else if (V1.x*V1.x + V1.y*V1.y + V1.z*V1.z == 0)
			V2 = V0;
		else*/
		CrossProduct(&V0, &V1, &V2);
//		printf("u=%lf, v=%lf, V2.x=%lf, V2.y=%lf, V2.z=%lf\n", u, v, V2.x, V2.y, V2.z);
		
		fNorme = sqrt(V2.x*V2.x + V2.y*V2.y + V2.z*V2.z);
		if (fNorme > fMaxNorme && false)
		{
			fMaxNorme = fNorme;
//			printf("Max = %lf (%lf, %lf)\n", fNorme, u, v);
		}
		if (fNorme > 0 && myDome.bDerivativeNorme)
		{
			V2.x /= (fNorme);
			V2.y /= (fNorme);
			V2.z /= (fNorme);
		}
		*dX -= myDome.fDerivative * V2.x;
		*dY -= myDome.fDerivative * V2.y;
		*dZ -= myDome.fDerivative * V2.z;

		double Factor;
		if (myDome.bDerivativeNorme)
			Factor = (2 * myDome.fDerivative) + 1;
		else
			Factor = (_PI * myDome.fDerivative / 4) + 1;
		*dX /= Factor;
		*dY /= Factor;
		*dZ /= Factor;
//		return;
	}

	if (bDerivee && (myDome.fNoiseAmplitudeR != 0 || myDome.fNoiseAmplitudeZ != 0))
	{
		float   w1, w2;
		double  c1, c2;

		c1 = myDome.fNoiseAmplitudeR * (1 - (1 - v) * myDome.fNoiseProgressiveR);
		c2 = myDome.fNoiseAmplitudeZ * (1 - (1 - v) * myDome.fNoiseProgressiveZ);

		w1 = (float)fabs(sin(myDome.fNoiseFrequencyU * u / 2));
		w2 = (float)mathfmod(myDome.fNoiseFrequencyV * (u == 0 ? 0 : v), 1);

		*dX *= (1 + c1 * PerlinNoiseGet2D(w1, w2));
		*dY *= (1 + c1 * PerlinNoiseGet2D(w1, w2));
//		*dZ = Z * (1 + c2 * PerlinNoiseGet2D(w1, w2)) / 2;

		*dZ *= (1 + c2 * PerlinNoiseGet2D(w1, w2));
	}

//	*dZ *= myDome.fDomeScaleZ;


}

float DistancePlan(GLpoint P1, GLpoint P2, GLpoint P3, GLpoint Point)
{
	float a, b, c, d;
	float Norme;

	a = (P2.z - P1.z)*(P3.y - P1.y) - (P2.y - P1.y)*(P3.z - P1.z);
	b = (P3.z - P1.z)*(P2.x - P1.x) - (P3.x - P1.x)*(P2.z - P1.z);
	c = (P2.y - P1.y)*(P3.x - P1.x) - (P2.x - P1.x)*(P3.y - P1.y);
	d =  P1.z * (P2.x - P1.x)*(P3.y - P1.y) + P1.y * (P3.x - P1.x)*(P2.z - P1.z) + P1.x * (P2.y - P1.y)*(P3.z - P1.z) -
		(P1.z * (P3.x - P1.x)*(P2.y - P1.y) + P1.y * (P2.x - P1.x)*(P3.z - P1.z) + P1.x * (P3.y - P1.y)*(P2.z - P1.z));
	Norme = sqrtf(a * a + b * b + c * c);

	if (Norme == 0.0f)
		return 0.0f;

	return (fabsf(a*Point.x + b * Point.y + c * Point.z + d) / Norme);
}
void ChargeBase(int nBase)
{

	BaseFunctionX = BaseFunctionListe[nBase].FunctionX;
	BaseFunctionY = BaseFunctionListe[nBase].FunctionY;
	return;

}
void ChargeTestSlope(int nSlopeXY, int nSlopeZ)
{
	SlopeFunctionXY = TestSlopeFunctionListe[nSlopeXY].Function;
	SlopeFunctionZ = TestSlopeFunctionListe[nSlopeZ].Function;

}
int GetSlopeFunction(int nFunctionXY, int nFunctionZ)
{
	int i = 0;


	while (i < nbSlopes)
	{
		if (SlopeFunctionListe[i].FunctionXY == TestSlopeFunctionListe[nFunctionXY].Function &&
			SlopeFunctionListe[i].FunctionZ == TestSlopeFunctionListe[nFunctionZ].Function)
			return i;
		i++;
	}

	return -1;
}
void ChargeSlope(int nSlope)
{
	Parametric_2 fTemp;

	if (nSlope == CUSTOM_SLOPE)
	{
		fTemp = SlopeFunctionXY;
		SlopeFunctionXY = SlopeFunctionZ;
		SlopeFunctionZ  = fTemp;
	}
	else
	{
		if (myDome.bSlopeInverse)
		{
			SlopeFunctionXY = SlopeFunctionListe[nSlope].FunctionZ;
			SlopeFunctionZ = SlopeFunctionListe[nSlope].FunctionXY;
		}
		else
		{
			SlopeFunctionXY = SlopeFunctionListe[nSlope].FunctionXY;
			SlopeFunctionZ = SlopeFunctionListe[nSlope].FunctionZ;
		}
	}
	return;

}
void ChargeWindow(int nWindow)
{
	WindowFunctionZ = WindowFunctionListe[nWindow].WindowFunctionZ;
	return;
}

void ChargeFrise(int nFrise)
{
	FriseFunction = WindowFunctionListe[nFrise].WindowFunctionZ;
	return;
}


void DefaultSlopeParameters()
{

	myDome.fSlopeParameterXY = SlopeFunctionListe[myDome.nSlope].fSlopeParameterXY;
	myDome.fSlopeParameterZ = SlopeFunctionListe[myDome.nSlope].fSlopeParameterZ;
	if (myDome.nSlope == SLOPE_RUSSIAN_DOME && myDome.fDomeScaleZ == 1.0)
		myDome.fDomeScaleZ = 2.0;

	glParameter->sync_live();
	glSlope->sync_live();

	return;


}


void ChargeShrink(int nShrink)
{
	switch (nShrink)
	{
	case SHRINK_NONE:
		ShrinkCoefficientFunction = ShrinkNoneFunction;
		break;
	case SHRINK_SQUARE_ROOT:
		ShrinkCoefficientFunction = ShrinkSquareRootFunction;
		break;
	case SHRINK_LINEAR:
		ShrinkCoefficientFunction = ShrinkLinearFunction;
		break;
	case SHRINK_SQUARE:
		ShrinkCoefficientFunction = ShrinkSquareFunction;
		break;
	default:
		ShrinkCoefficientFunction = ShrinkNoneFunction;
		break;
	}
}

double CalculDerivee(GLpoint *myVERTICES, int i, int j, Vector3D *V, int nMaxI, int nMaxJ)
{
	Vector3D v0, v1, v2, v3, w0, w1;
	double dNorme;

	v0.x = myVERTICES[(i + 1) * nMaxJ + (j)].x - myVERTICES[i * nMaxJ + j].x;
	v0.y = myVERTICES[(i + 1) * nMaxJ + (j)].y - myVERTICES[i * nMaxJ + j].y;
	v0.z = myVERTICES[(i + 1) * nMaxJ + (j)].z - myVERTICES[i * nMaxJ + j].z;

	v1.x = myVERTICES[(i) * nMaxJ + j-1].x - myVERTICES[i * nMaxJ + j].x;
	v1.y = myVERTICES[(i) * nMaxJ + j-1].y - myVERTICES[i * nMaxJ + j].y;
	v1.z = myVERTICES[(i) * nMaxJ + j-1].z - myVERTICES[i * nMaxJ + j].z;


	v2.x = myVERTICES[((i - 1 + nMaxI - 1) % (nMaxI - 1)) * nMaxJ + j].x - myVERTICES[(i) * nMaxJ + (j)].x;
	v2.y = myVERTICES[((i - 1 + nMaxI - 1) % (nMaxI - 1)) * nMaxJ + j].y - myVERTICES[(i) * nMaxJ + (j)].y;
	v2.z = myVERTICES[((i - 1 + nMaxI - 1) % (nMaxI - 1)) * nMaxJ + j].z - myVERTICES[(i) * nMaxJ + (j)].z;



	if (j < nMaxJ - 1)
	{
		v3.x = myVERTICES[i * nMaxJ + (j + 1)].x - myVERTICES[(i) * nMaxJ + (j)].x;
		v3.y = myVERTICES[i * nMaxJ + (j + 1)].y - myVERTICES[(i) * nMaxJ + (j)].y;
		v3.z = myVERTICES[i * nMaxJ + (j + 1)].z - myVERTICES[(i) * nMaxJ + (j)].z;
	}
	else
	{
		v3.x = 0;
		v3.y = 0;
		v3.z = 0;
	}

	CrossProduct(&v0, &v1, &w0);
	CrossProduct(&v1, &v2, &w1);

	V->x = (w0.x + w1.x);
	V->y = (w0.y + w1.y);
	V->z = (w0.z + w1.z);
	CrossProduct(&v2, &v3, &w0);
	CrossProduct(&v3, &v0, &w1);
	V->x += (w0.x + w1.x);
	V->y += (w0.y + w1.y);
	V->z += (w0.z + w1.z);



	dNorme = sqrt(V->x * V->x + V->y*V->y + V->z*V->z);
	if (dNorme > 0)
	{
		V->x /= dNorme;
		V->y /= dNorme;
		V->z /= dNorme;
	}
	return dNorme;
}





void ChargeMesh()
{
	//	int i,j;
	//	double dX, dY, dZ;
	int nMaxI, nMaxJ;

	ChargeBase(myDome.nBase);
//	ChargeSlope(myDome.nSlope);
	ChargeWindow(myDome.nWindow);

	//printf("Allocation = %d ; vertices = %d",bAllocation, (int)myOBJ_VERTICE );

	nMaxI = (myDome.nMeridianMain + (myDome.nMeridianMain)*myDome.nMeridianMinor*myDisplay.bDisplayMinor) * myDome.nBaseSides + 1;
	nMaxJ = myDome.nParallelMain + (myDome.nParallelMain - 1)*myDome.nParallelMinor*myDisplay.bDisplayMinor;
//	printf("MaxI = % d, MaxJ = %d\n", nMaxI, nMaxJ);
	if (bAllocation)
	{
		if (myOBJ_VERTICE)
			free(myOBJ_VERTICE);


		myOBJ_VERTICE = (GLpoint *)malloc(nMaxI*nMaxJ * sizeof(GLpoint));

		if (myOBJ_VERTICE == 0)
		{
			Message("     Memory Allocation Problem   ");
			return;
		}
	}


	CalculMesh(nMaxI, nMaxJ, myOBJ_VERTICE);
	for (int i = 0; i <  myDome.nLissageNombre; i++)
		Lissage(myOBJ_VERTICE, nMaxI, nMaxJ);
	for (int i = 0; i < myDome.nDiffGeomNombre; i++)
		DifferentialGeometry(myOBJ_VERTICE, nMaxI, nMaxJ);
	for (int i = 0; i < myDome.nNormaleNombre; i++)
		NormaleSharpening(myOBJ_VERTICE, nMaxI, nMaxJ);
	Bumps(nMaxI, nMaxJ, myOBJ_VERTICE, false);
	Marquees(nMaxI, nMaxJ, myOBJ_VERTICE, false);
	Enhance(myOBJ_VERTICE, nMaxI, nMaxJ);
	CalculComb();
	CalculSnake();
	CalculFrise();
	CalculSpirale();

	ComputeBaseLength();
}



void AdjustSphere()
{
	float Radius;


	Radius = myDome.fMeridianMainRadius  > myDome.fParallelMainRadius ? myDome.fMeridianMainRadius : myDome.fParallelMainRadius;
	Radius = Radius > myDome.fDiagonalURRadius ? Radius : myDome.fDiagonalURRadius;
	Radius = Radius > myDome.fDiagonalDRRadius ? Radius : myDome.fDiagonalDRRadius;
	myDome.fSphereMPRadius = Radius;
	if (!strcmp(myDome.sSphereMPMaterial, S_NONE))
		strcpy_s(myDome.sSphereMPMaterial, (Radius == myDome.fMeridianMainRadius ? myDome.sMeridianMainMaterial :
										   (Radius == myDome.fParallelMainRadius ? myDome.sParallelMainMaterial :
										   (Radius == myDome.fDiagonalURRadius   ? myDome.sDiagonalURMaterial : myDome.sDiagonalDRMaterial))));



	Radius = myDome.fMeridianMainRadius  > myDome.fParallelMinorRadius ? myDome.fMeridianMainRadius : myDome.fParallelMinorRadius;
	Radius = Radius > myDome.fDiagonalURRadius ? Radius : myDome.fDiagonalURRadius;
	Radius = Radius > myDome.fDiagonalDRRadius ? Radius : myDome.fDiagonalDRRadius;
	myDome.fSphereMpRadius = Radius;
	if (!strcmp(myDome.sSphereMpMaterial, S_NONE))
		strcpy_s(myDome.sSphereMpMaterial, (Radius == myDome.fMeridianMainRadius  ? myDome.sMeridianMainMaterial  :
									       (Radius == myDome.fParallelMinorRadius ? myDome.sParallelMinorMaterial :
									       (Radius == myDome.fDiagonalURRadius    ? myDome.sDiagonalURMaterial    : myDome.sDiagonalDRMaterial))));


	Radius = myDome.fMeridianMinorRadius  > myDome.fParallelMainRadius ? myDome.fMeridianMinorRadius : myDome.fParallelMainRadius;
	Radius = Radius > myDome.fDiagonalURRadius ? Radius : myDome.fDiagonalURRadius;
	Radius = Radius > myDome.fDiagonalDRRadius ? Radius : myDome.fDiagonalDRRadius;
	myDome.fSpheremPRadius = Radius;
	if (!strcmp(myDome.sSpheremPMaterial, S_NONE))
		strcpy_s(myDome.sSpheremPMaterial, (Radius == myDome.fMeridianMinorRadius ? myDome.sMeridianMinorMaterial :
									       (Radius == myDome.fParallelMainRadius  ? myDome.sParallelMainMaterial  :
									       (Radius == myDome.fDiagonalURRadius    ? myDome.sDiagonalURMaterial    : myDome.sDiagonalDRMaterial))));


	Radius = myDome.fMeridianMinorRadius  > myDome.fParallelMinorRadius ? myDome.fMeridianMinorRadius : myDome.fParallelMinorRadius;
	Radius = Radius > myDome.fDiagonalURRadius ? Radius : myDome.fDiagonalURRadius;
	Radius = Radius > myDome.fDiagonalDRRadius ? Radius : myDome.fDiagonalDRRadius;
	myDome.fSpherempRadius = Radius;
	if (!strcmp(myDome.sSpherempMaterial, S_NONE))
		strcpy_s(myDome.sSpherempMaterial, (Radius == myDome.fMeridianMinorRadius ? myDome.sMeridianMinorMaterial :
									       (Radius == myDome.fParallelMinorRadius ? myDome.sParallelMinorMaterial :
									       (Radius == myDome.fDiagonalURRadius ? myDome.sDiagonalURMaterial : myDome.sDiagonalDRMaterial))));

	glParameter->sync_live();
}

/**************************************** control_cb() *******************/
/* GLUI control callback                                                 */
bool FileValide(char* filename)
{
	char* dot = strrchr(filename, '.');
	//return dot && !strcmp(dot, Pattern);
	return dot && !strcmp(dot, ".DOM");
}

int fill_lbfile(GLUI_Listbox *lb, char *sSearch, char *sFile)
{
	WIN32_FIND_DATA FN;
	HANDLE hFind;
	int dum;
	int i;
	GLUI_String item;
	int nListe = -1;

	i = 0;
	do
	{
		dum = lb->delete_item(i++);
	} while (dum == 1);

	hFind = FindFirstFile(L"*.*", &FN);
	if (hFind != INVALID_HANDLE_VALUE)
	{
		lb->add_item(0, "   ");
		i = 1;
		do
		{
			char filename[MAX_PATH];
			sprintf_s(filename, "%ws", FN.cFileName);
			if ((FN.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) )
			{
				int len = strlen(filename);
				item = '\\';
				item += filename;
				lb->add_item(i, item.c_str());
				i++;
			}
		} while (FindNextFile(hFind, &FN) != 0);

		if (GetLastError() == ERROR_NO_MORE_FILES)
			FindClose(hFind);
		else
			perror("fbreaddir");
	}

	hFind = FindFirstFile(L"*.DOM*", &FN);
	if (hFind != INVALID_HANDLE_VALUE)
	{
		do
		{
			char filename[MAX_PATH];
			sprintf_s(filename, "%ws", FN.cFileName);
			if ( FileValide(filename))
			{
				int len = strlen(filename);
				item = filename;
				lb->add_item(i, item.c_str());

				if (sFile && strcmp(filename, sFile) == 0)
					nListe = i;

				i++;
			}
		} while (FindNextFile(hFind, &FN) != 0);

		if (GetLastError() == ERROR_NO_MORE_FILES)
			FindClose(hFind);
		else
			perror("fbreaddir");
	}
	lb->set_int_val(0);
	if (nListe != -1)
		lb->set_int_val(nListe);
	//£
	//lb->add_item(i, "                                                                                 ");
	return 1;
}
int fill_lbfile_(GLUI_Listbox *lb, char *sSearch)
{
	WIN32_FIND_DATA FN;
	HANDLE hFind;
	int dum;
	int i;
	GLUI_String item;

	hFind = FindFirstFile(L"*.*", &FN);
	i = 0;
	do
	{
		dum = lb->delete_item(i++);
	} while (dum == 1);

	if (hFind != INVALID_HANDLE_VALUE)
	{
		lb->add_item(0, "   ");
		i = 1;
		do
		{
			char filename[MAX_PATH];
			sprintf_s(filename, "%ws", FN.cFileName);
			if ((FN.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) || FileValide(filename))
			{
				int len = strlen(filename);
				if (FN.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
				{
					item = '\\';
					item += filename;
				}
				else
				{
					item = filename;
				}
				lb->add_item(i, item.c_str());
//				printf("File  %3d = %s", i, filename);
				i++;
			}
		} while (FindNextFile(hFind, &FN) != 0);

		if (GetLastError() == ERROR_NO_MORE_FILES)
			FindClose(hFind);
		else
			perror("fbreaddir");
	}
	lb->sort_items();
	return 1;
}
#pragma warning(disable:4996)
static wchar_t* ToWChar(const char* text)
{
	const size_t size = strlen(text) + 1;
	wchar_t* wText = new wchar_t[size];
	mbstowcs(wText, text, size);
	return wText;
}
/*
void control_lb(int control)
{
	GLUI_String sFile;
	char sSelected[100];
	char* p = 0;

	sFile = lbDomeFile->curr_text;
	strcpy_s(sSelected, sFile.c_str());

	printf("File = %s\n", sFile.c_str());
	if (sSelected[0] != '\\')
	{

		p = strstr(sSelected, ".");

		p[0] = 0;

		strcpy_s(myDome.sDomeName, sSelected);
		LoadDome();
		ChargeMesh();

		glParameter->sync_live();

	}
	else
	{
		wchar_t* dir = ToWChar(sSelected + 1);
		int b = SetCurrentDirectory(dir);
		fill_lbfile(lbDomeFile, (char *)"*.DOM", NULL);
	}
}
*/

void ControlSkeleton(int Control)
{
	switch (Control)
	{
	case BT_BASE:
		btBase->disable();
		glBase->show();
		break;

	case BT_SLOPE:
		btSlope->disable();
		glSlope->show();
		break;

	case BT_WINDOW:
		btWindow->disable();
		glWindow->show();
		break;
	case BT_TOP:
		btTop->disable();
		glTop->show();
		break;
	case SP_MERIDIAN_MAIN:
	case SP_PARALLEL_MAIN:
	case SP_MERIDIAN_MINOR:
	case SP_PARALLEL_MINOR:
		bAllocation = 1;
		ChargeMesh();
		break;
	case BT_ADJUST_SCALE:
		myDome.fDomeScaleZ = 1 / (1 - myDome.fDomeSocle);
		glParameter->sync_live();
		ChargeMesh();
		break;
	case BT_SPHERE_ADJUST:
		AdjustSphere();
		break;
	case BT_SPIKE_DEFAULT:
		SpikeDefault();
		glParameter->sync_live();
		break;

	case ET_MERIDIAN_MAIN_MATERIAL:
	case ET_MERIDIAN_MINOR_MATERIAL:
	case ET_PARALLEL_MAIN_MATERIAL:
	case ET_PARALLEL_MINOR_MATERIAL:
	case ET_DIAGONALE_DR_MATERIAL:
	case ET_DIAGONALE_UR_MATERIAL:
	case ET_SPHERE_MP_MATERIAL:
	case ET_SPHERE_Mp_MATERIAL:
	case ET_SPHERE_mP_MATERIAL:
	case ET_SPHERE_mp_MATERIAL:
	case SP_SPHERE_MP_RADIUS:
	case SP_SPHERE_Mp_RADIUS:
	case SP_SPHERE_mP_RADIUS:
	case SP_SPHERE_mp_RADIUS:
	case SP_MERIDIAN_MAIN_RADIUS:
	case SP_MERIDIAN_MINOR_RADIUS:
	case SP_PARALLEL_MAIN_RADIUS:
	case SP_PARALLEL_MINOR_RADIUS:
		break;
	case SP_SPIKE_MP:
	case SP_SPIKE_Mp:
	case SP_SPIKE_mP:
	case SP_SPIKE_mp:
	case SP_SPIKE_ANGLE:
		break;
	}
}

void ControlNoise(int Control)
{
	switch (Control)
	{
	case SP_NOISE_OCTAVE:
	case SP_NOISE_DENSITY:
	case SP_NOISE_FREQUENCY_U:
	case SP_NOISE_FREQUENCY_V:
	case SP_NOISE_AMPLITUDE_R:
	case SP_NOISE_AMPLITUDE_Z:
	case SP_NOISE_PROGRESSIVE_R:
	case SP_NOISE_PROGRESSIVE_Z:
		ChargeMesh();
		break;

	case BT_NOISE_DEFAULT:
		NoiseDefault();
		glPreProcess->sync_live();
		ChargeMesh();
		break;
	case SP_NOISE_SEED:
		myNoise.bStart = true;
		myNoise.nSeed = myDome.nNoiseSeed;
		ChargeMesh();
		break;



	}
}

void ControlFrame(int Control)
{
	std::string sFileName;
	char* p = 0;
	char sLocal[400];

	int nSlopeFunction;

	switch (Control)
	{
	case LB_WINDOW:
		tbWindow->set_text(WindowFunctionListe[myDome.nWindow].sNom);
		myDome.fWindowParameter = WindowFunctionListe[myDome.nWindow].fParameter;
		myDome.fWindowPower = WindowFunctionListe[myDome.nWindow].fPower;
		glParameter->sync_live();

	case SP_WINDOW_ATTACK:
	case SP_WINDOW_DAMPING:
	case SP_WINDOW_DELAY:
	case SP_WINDOW_HEIGHT:
	case SP_WINDOW_PERSIDE:
	case SP_WINDOW_POWER:
	case SP_WINDOW_SHIFT:
	case CK_WINDOW_INVERSE:
	case CK_WINDOW_BLIND:
	case SP_WINDOW_AWNING:
	case SP_WINDOW_PARAMETER:
		glWindow->sync_live();
		ChargeMesh();
		break;
	case SP_WINDOW_LOWER:
		if (myDome.fWindowLower > 0)
		{
			myDome.bWindowBlind = 1;
			ckWindowBlind->disable();
		}
		else
			ckWindowBlind->enable();
		glWindow->sync_live();
		ChargeMesh();
		break;
	case BT_WINDOW_DEFAULT:
		WindowDefault();

		glWindow->sync_live();
		tbWindow->set_text(WindowFunctionListe[myDome.nWindow].sNom);
		glParameter->sync_live();
		ChargeMesh();
		break;
	case BT_WINDOW_QUIT:
		btWindow->enable();
		glWindow->hide();
		break;

	case SP_CRENEL_PERSIDE:
		if(myDome.nCrenelPerSide == -1)
		{
			spCrenelDelay->disable();
			spCrenelWide->disable();
		}
		else
		{ 
			spCrenelDelay->enable();
			spCrenelWide->enable();
		}
	case SP_SLOPE_PARAMETER:
	case SP_SLOPE_PARAMETER_XY:
	case SP_CRENEL_THICK:
	case SP_CRENEL_HEIGHT:
	case SP_CRENEL_POWER:
	case SP_CRENEL_FRONT:
	case SP_CRENEL_ANGLE:
	case SP_CRENEL_WIDE:
	case SP_CRENEL_DELAY:
		ComputeSlopeLength();
		ChargeMesh();
		break;

	case SP_SLOPE_TWIST:
	case RG_SLOPE_SPACED:
	case SP_SLOPE_DELAY:
		ChargeMesh();
		break;

	case BT_SLOPE_DEFAULT:
		SlopeDefault();

		GetFunctions(myDome.nSlope);

		glSlope->sync_live();
		tbSlope->set_text(SlopeFunctionListe[myDome.nSlope].sNom);
		ChargeSlope(myDome.nSlope);

		ComputeSlopeLength();

		glParameter->sync_live();
		ChargeMesh();
		break;
	case BT_SLOPE_QUIT:
		btSlope->enable();
		glSlope->hide();
		break;

	case LB_TEST_SLOPE_XY:
	case LB_TEST_SLOPE_Z:

		ChargeTestSlope(myDome.nTestSlopeXY, myDome.nTestSlopeZ);
		ComputeSlopeLength();
		nSlopeFunction = GetSlopeFunction(myDome.nTestSlopeXY, myDome.nTestSlopeZ);
		if (nSlopeFunction == -1)
			lbSlopeType->set_int_val(nbSlopes - 1);
		else
			lbSlopeType->set_int_val(nSlopeFunction);

		glParameter->sync_live();
		glSlope->sync_live();

		ChargeMesh();

		break;

	case LB_SLOPE:
		DefaultSlopeParameters();
		tbSlope->set_text(SlopeFunctionListe[myDome.nSlope].sNom);
		GetFunctions(myDome.nSlope);
		glParameter->sync_live();
		glSlope->sync_live();

	case CK_SLOPE_INVERSE:
		ChargeSlope(myDome.nSlope);
		ComputeSlopeLength();

		ChargeMesh();
		break;

	case SP_BASE_DAMPING:
	case SP_BASE_EXPAND:
	case SP_BASE_PARAMETER:
	case SP_BASE_SMOOTHING:
	case SP_BASE_ROUND:
	case SP_BASE_BARYCENTRE:
	case SP_BASE_MORPHING:
	case SP_BASE_SECONDARY_ROTATION:
	case SP_BASE_SECONDARY_EXPANSION:
	case CK_BASE_MORPHING_DIVISION:
	case SP_BASE_SPIRAL:
	case RG_BASE_SPACED:
		ChargeMesh();
		break;

	case SP_BASE_SIDE:
		if (myDome.nBase == BASE_FOURIER)
			BaseFourierRho(0, 0, true);
		if (myDome.nBaseSides % 2)
		{
			ckDivision->disable();
			myDome.bBaseMorphingDivision = 0;
			glBase->sync_live();
		}
		else
			ckDivision->enable();
		bAllocation = 1;
		ChargeMesh();
		break;

	case BT_BASE_DEFAULT:
		BaseDefault();

		if (myDome.nBaseSides % 2)
			ckDivision->disable();


		glBase->sync_live();
//		spBaseRound->disable();

		tbBase->set_text(BaseFunctionListe[myDome.nBase].sNom);
		glBase->sync_live();
		tbTop->set_text(TopListe[myDome.nTopType].sNom);
		glTop->sync_live();
		glParameter->sync_live();
		ChargeMesh();
		break;

	case LB_BASE:
		if (myDome.nBase == BASE_FOURIER)
			BaseFourierRho(0, 0, true);

	case LB_BASE_2:
		//		spBaseRound->enable();
//		if (myDome.nBase == BASE_CIRCLE)
//		{
//			myDome.fBaseRound = 0;
//			spBaseRound->disable();
//		}
		if (myDome.nBase == BASE_FOURIER)
			BaseFourierRho(0, 0, true);

		myDome.fBaseParameter = BaseFunctionListe[myDome.nBase].fBaseParameter;
		myDome.nBaseSides = __max(myDome.nBaseSides, BaseFunctionListe[myDome.nBase].iMinSide);
		spSide->set_int_limits(BaseFunctionListe[myDome.nBase].iMinSide, 30);
		if (myDome.nBaseSides % 2)
		{
			ckDivision->disable();
			myDome.bBaseMorphingDivision = 0;
		}



		tbBase->set_text(BaseFunctionListe[myDome.nBase].sNom);
		glParameter->sync_live();
		glBase->sync_live();
		ChargeMesh();
		break;	
	case BT_BASE_QUIT:
		btBase->enable();
		glBase->hide();
		break;
	case BT_TOP_QUIT:
		btTop->enable();
		glTop->hide();
		break;

	case SP_DOME_SOCLE:
		ComputeSlopeLength();
		ChargeMesh();
		break;
	case LB_TOP_TYPE:
		myDome.fTopPart = TopListe[myDome.nTopType].fPart;
		myDome.fTopParameter = TopListe[myDome.nTopType].fParameter;
		tbTop->set_text(TopListe[myDome.nTopType].sNom);
		glTop->sync_live();
		glParameter->sync_live();
//		printf("%s\n", TopListe[myDome.nTopType].sNom);
	case SP_TOP_PARAMETER:
	case SP_DOME_TOP:
		ComputeSlopeLength();
		glTop->sync_live();
		ChargeMesh();
		break;

	case SP_DOME_END:
	case SP_DOME_SCALE_Z:
	case SP_DOME_INFOLDING:
		ChargeMesh();
		break;

	case  BT_DOME_INC:
		PrintPOV();
		break;
	case  BT_DOME_STL:
		PrintSTL();
		break;
	case  BT_DOME_OBJ:
		PrintOBJ();
		break;
	case  BT_DOME_DEFAULT:
		InitParameters();
		if (myDome.bFaceGenerate)
			ckFaceSmooth->enable();
		else
		{
			myDome.bFaceSmooth = 0;
			ckFaceSmooth->disable();
		}
		GetFunctions(myDome.nSlope);

		glParameter->sync_live();
		glDecor->sync_live();
		glPostProcess->sync_live();
		glPreProcess->sync_live();
		glBase->sync_live();
		glSlope->sync_live();
		glWindow->sync_live();
		glTop->sync_live();

		glDisplay->sync_live();
		ControlDisplay(BT_DEFAULT_VIEW);
		break;
	case BT_DOME_LOAD:
		LoadDome();
		ChargeMesh();
	case CK_FACE:
		if (myDome.bFaceGenerate)
			ckFaceSmooth->enable();
		else
		{
			myDome.bFaceSmooth = 0;
			ckFaceSmooth->disable();
		}
		glParameter->sync_live();
		break;
	case BT_DOME_SAVE:
		SaveDome();
		break;

	case LB_DOME_FILE:
		sFileName = lbDomeFile->curr_text;
		strcpy_s(sLocal, sFileName.c_str());

		printf("File = %s\n", sFileName.c_str());
		if (sLocal[0] != '\\')
		{

			p = strstr(sLocal, ".");

			p[0] = 0;

			strcpy_s(myDome.sDomeName, sLocal);
			LoadDome();
			ChargeMesh();

			if (myDome.bFaceGenerate)
				ckFaceSmooth->enable();
			else
			{
				myDome.bFaceSmooth = 0;
				ckFaceSmooth->disable();
			}
			glParameter->sync_live();

		}
		else
		{
			wchar_t* dir = ToWChar(sLocal + 1);
			int b = SetCurrentDirectory(dir);
			fill_lbfile(lbDomeFile, (char *)"*.DOM", NULL);
		}
		char Buffer[501];
		_getcwd(Buffer, 500);
		tbDirectory->set_text(Buffer);

		break;
	case ET_DOME_MATERIAL:
	case ET_DOME_NAME:
		break;


	case ET_FACE_CENTER_MATERIAL:
		if (!strcmp(myDome.sFaceEdgeMaterial, S_NONE))
			strcpy_s(myDome.sFaceEdgeMaterial, myDome.sFaceCenterMaterial);
		glParameter->sync_live();
		break;

	case CK_SMOOTH:
	case SP_MERIDIAN_EDGE_SIZE:
	case SP_PARALLEL_EDGE_SIZE:
	case ET_FACE_EDGE_MATERIAL:
	case CK_WIRE:
	case SP_WIRE_SHRINKING:
	case SP_WIRE_SHRINKING_SPEED:
	case LB_WIRE_SHRINKING_TYPE:
		break;

	}
}

void SetView(int nView)
{
	switch (nView)
	{
	case BT_TOP_VIEW:
		for (int i = 0; i < 16; i++)      // { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 }
		{
			if ((i % 5) == 0)
				DomeRotate[i] = 1.0;
			else
				DomeRotate[i] = 0.0;
		}
		DomeRotate[0] = -1;
		DomeRotate[10] = -1;
		rotationX = 0.0;
		rotationY = 0.0;
		glDisplay->sync_live();

		break;
	case BT_SIDE_VIEW:
		for (int i = 0; i < 16; i++)      // { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 }
		{
			DomeRotate[i] = 0.0;
		}
		DomeRotate[0] = 1;
		DomeRotate[6] = -1;
		DomeRotate[9] = 1;
		DomeRotate[15] = 1;
		rotationX = 0.0;
		rotationY = 0.0;

		glDisplay->sync_live();
		break;
	}
}
void ControlDisplay(int Control)
{
	switch (Control)
	{
	case LB_DISPLAY:
		if (myDisplay.nDisplay == DISPLAY_Base)
		{
			SetView(BT_TOP_VIEW);
			DrawBase();
		}
		if (myDisplay.nDisplay == DISPLAY_Slope)
		{
			SetView(BT_SIDE_VIEW);
			DrawSlope();
		}
		if (myDisplay.nDisplay == DISPLAY_Window)
		{
			SetView(BT_SIDE_VIEW);
			DrawWindow();
		}
		if (myDisplay.nDisplay == DISPLAY_Snake)
		{
			SetView(BT_TOP_VIEW);
			DrawSnake();
		}
		if (myDisplay.nDisplay == DISPLAY_Spirale)
		{
			SetView(BT_TOP_VIEW);
			DrawSpirale();
		}
		if (myDisplay.nDisplay == DISPLAY_Frise)
		{
			SetView(BT_SIDE_VIEW);
			DrawFrise();
		}
		if (myDisplay.nDisplay == DISPLAY_Comb)
		{
			SetView(BT_TOP_VIEW);
			DrawComb();
		}

		break;
	case CK_WIRE_FRAME:
		if (!myDisplay.bWireFrame && !myDisplay.bFlatFace && !myDisplay.bSmooth)
			myDisplay.bWireFrame = 1;
		glDisplay->sync_live();
		break;
	case CK_FLAT_FACE:
		if (myDisplay.bFlatFace)
		{
			myDisplay.bSmooth = 0;
			ckWireFrame->enable();
		}
		else if (!myDisplay.bSmooth)
		{
			myDisplay.bWireFrame = 1;
			ckWireFrame->disable();
		}
		glDisplay->sync_live();
		break;
	case CK_SMOOTH_FACE:
		if (myDisplay.bSmooth)
		{
			myDisplay.bFlatFace = 0;
			ckWireFrame->enable();
		}
		else if (!myDisplay.bFlatFace)
		{
			myDisplay.bWireFrame = 1;
			ckWireFrame->disable();
		}
		glDisplay->sync_live();
		break;
	case CK_DISPLAY_MINOR:
		ChargeMesh();
		break;
	case CK_LIGHT_1:
		if (myDisplay.nLight_1)
			glEnable(GL_LIGHT0);
		else
			glDisable(GL_LIGHT0);
		break;
	case CK_LIGHT_2:
		if (myDisplay.nLight_2)
		{
			glEnable(GL_LIGHT1);
			glEnable(GL_LIGHT2);
		}
		else
		{
			glDisable(GL_LIGHT1);
			glDisable(GL_LIGHT2);
		}
		break;
	case BT_DEFAULT_VIEW:
		DomePosition[0] = 0.0;
		DomePosition[1] = 0.0;
		DomePosition[2] = 0.0;

/*		for (int i = 0; i < 16; i++)      // { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 }
		{
			if ((i % 5) == 0)
				DomeRotate[i] = 1.0;
			else
				DomeRotate[i] = 0.0;
		}
		DomeRotate[0] = -1;
		DomeRotate[10] = -1;*/

		DomeRotate[0] = 1.000000;
		DomeRotate[1] = 0.000000;
		DomeRotate[2] = 0.000000;
		DomeRotate[3] = 0.000000;
		DomeRotate[4] = 0.000000;
		DomeRotate[5] = -0.453990f;
		DomeRotate[6] = 0.891007f;
		DomeRotate[7] = 0.000000;
		DomeRotate[8] = 0.000000;
		DomeRotate[9] = 0.891007f;
		DomeRotate[10] = 0.453990f;
		DomeRotate[11] = 0.000000;
		DomeRotate[12] = 0.000000;
		DomeRotate[13] = 0.000000;
		DomeRotate[14] = 0.000000;
		DomeRotate[15] = 1.000000;

		rotationX = 0.0;
		rotationY = 0.0;

		glDisplay->sync_live();
		myGlutDisplay();
		break;
	case BT_TOP_VIEW:
		for (int i = 0; i < 16; i++)      // { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 }
		{
			if ((i % 5) == 0)
				DomeRotate[i] = 1.0;
			else
				DomeRotate[i] = 0.0;
		}
		DomeRotate[0] = -1;
		DomeRotate[10] = -1;
		rotationX = 0.0;
		rotationY = 0.0;
		glDisplay->sync_live();
		myGlutDisplay();
		break;
	case BT_SIDE_VIEW:
		for (int i = 0; i < 16; i++)      // { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 }
		{
			DomeRotate[i] = 0.0;
		}
		DomeRotate[0] = 1;
		DomeRotate[6] = -1;
		DomeRotate[9] = 1;
		DomeRotate[15] = 1;
		rotationX = 0.0;
		rotationY = 0.0;

		glDisplay->sync_live();
		myGlutDisplay();
		break;
/*	case LB_MATERIAL:
		for (int i = 0; i < 16;i++)
			printf("DomeRotate[%2d] = %lf\n", i, InstantRotate[i]);
		break;*/
	}
}

void ControlPreProcessing(int Control)
{
	switch (Control)
	{
	case BT_PRE_PROCESSING:
		btPreProcessing->disable();
		glPreProcess->show();
		break;
	case BT_PRE_PROCESSING_QUIT:
		btPreProcessing->enable();
		glPreProcess->hide();
		break;
	case SP_DERIVATIVE:
	case CK_DERIVATIVE_NORME:
		ChargeMesh();
		break;

	}
}

void ControlPostProcessing(int Control)
{
	int nFiltre;

	nFiltre = lbFiltreType->get_int_val();
	switch (Control)
	{
	case LB_LISSAGE_FILTRE:
		FilterInput((TYPE_FILTER)nFiltre, true);
	case CK_EVOLUTE_MERIDIEN:
	case CK_EVOLUTE_PARALLEL:
	case CK_NORMALE:
	case CK_LISSAGE:
	
	case SP_LISSAGE_NOMBRE:
	case SP_LISSAGE_PUISSANCE:
	case SP_LISSAGE_TAILLE:
	case SP_LISSAGE_POIDS_DIAGONALES:
	case SP_LISSAGE_DISTANCE_TYPE:
	case SP_LISSAGE_POIDS_SELF:
	case SP_LISSAGE_SIGMA_1:
	case SP_LISSAGE_SIGMA_2:
	case SP_LISSAGE_POIDS_HORIZONTAL:
	case SP_LISSAGE_POIDS_VERTICAL:
	case LB_LISSAGE_DISTANCE_TYPE:
	case LB_LISSAGE_DIRECTION:

		glPostProcess->sync_live();
		glParameter->sync_live();
		ChargeMesh();
		break;

	case SP_NORMALE_LAMBDA:
	case SP_NORMALE_NOMBRE:
		if (myDome.bNormale)
			ChargeMesh();
		break;
	case SP_EVOLUTE_LAMBDA:
	case SP_EVOLUTE_NOMBRE:
	case LB_DIFF_GEOM:
		if (myDome.bDiffGeomMeridien || myDome.bDiffGeomParallel)
			ChargeMesh();
		break;




	case BT_LISSAGE_DEFAULT:
		LissageDefault();
		glPostProcess->sync_live();
		ChargeMesh();
		break;

	case CK_MARQUEES:
		ChargeMesh();
		break;

	case LB_MARQUEE_FUNCTION:
	case SP_MARQUEE_PARAMETER:
	case SP_MARQUEE_LAMBDA:
	case SP_MARQUEE_DAMPING:
	case SP_MARQUEE_ALPHA:
	case SP_MARQUEE_BETA:
	case CK_MARQUEE_NORMALIZE:
	case SP_MARQUEE_TOP:
		if (myDome.bMarquee)
			ChargeMesh();
		break;
	case LB_ENHANCE_TYPE:
	case CK_ENHANCE_MERIDIAN:
	case CK_ENHANCE_PARALLEL:
		ChargeMesh();
		break;


	case LB_ENHANCE_MERGE:
	case LB_ENHANCE_CALCUL:
	case SP_ENHANCE_LAMBDA:
	case SP_ENHANCE_DAMPING:
		if (myDome.bEnhanceMeridian || myDome.bEnhanceParallel)
			ChargeMesh();
		break;

	case SP_POST_MERIDIAN_EVERY:
	case SP_POST_MERIDIAN_DELAY:
	case CK_POST_MERIDIAN_PER_SIDE:
	case SP_POST_PARALLEL_EVOLVE:
	case SP_POST_PARALLEL_FIRST:
	case SP_POST_PARALLEL_EVERY:

	case CK_BUMP:
		ChargeMesh();
		break;

	case SP_BUMP_DAMPING:
	case SP_BUMP_POWER:
	case CK_BUMP_NORMALIZE:
	case LB_BUMP_NORMAL:
	case LB_BUMP_TYPE:
	case SP_BUMP_LAMBDA:
	case SP_BUMP_DISTANCE:
	case LB_BUMP_DISTANCE_TYPE:
	case CK_BUMP_DISTANCE_INVERSE:
		if (myDome.bBump)
			ChargeMesh();
		break;
	case BT_POST_PROCESSING:
		btPostProcessing->disable();
		glPostProcess->show();
		break;

	case BT_POST_PROCESSING_QUIT:
		btPostProcessing->enable();
		glPostProcess->hide();
		break;

	}
}

void ControlDecor(int Control)
{
//	printf("Decor : %d\n", Control);
	switch (Control)
	{
	case LB_FRISE_FUNCTION:
		myDome.fFriseParameter = WindowFunctionListe[myDome.nFriseFunction].fParameter;
		myDome.fFrisePower = WindowFunctionListe[myDome.nFriseFunction].fPower;
	case SP_SNAKE:
	case SP_SNAKE_TOUR:
	case LB_SNAKE_EVOLVE:
	case SP_SNAKE_OFFSET:
	case LB_SNAKE_CHIRALE:
	case SP_SNAKE_STEP:
	case SP_SNAKE_SHIFT:
	case LB_SNAKE_SLOPE:
	case SP_SNAKE_SLOPE_PARAMETER:

	case SP_FRISE:
	case SP_FRISE_HEIGHT:
	case SP_FRISE_PERSIDE:
	case SP_FRISE_POWER:
	case SP_FRISE_SHIFT:
	case SP_FRISE_START:
	case SP_FRISE_STEP:
	case SP_FRISE_OFFSET:
	case SP_FRISE_GAP:
	case SP_FRISE_PARAMETER:
	case CK_FRISE_INVERSE:
	case CK_FRISE_REVERSE:
	case SP_FRISE_SPACING:

	case LB_SPIRALE_CHIRALE:
	case SP_SPIRALE:
	case SP_SPIRALE_DECALAGE:
	case SP_SPIRALE_DELAY:
	case SP_SPIRALE_SPEED:
	case SP_SPIRALE_STEP:
	case SP_SPIRALE_TOUR:
	case LB_SPIRALE_SLOPE:
	case SP_SPIRALE_SLOPE_PARAMETER:

	case SP_COMB:
	case SP_COMB_SHIFT:
	case SP_COMB_STEP:
	case LB_COMB_TYPE:
	case SP_COMB_THREAD:
	case SP_COMB_THREAD_GAP:
	case SP_COMB_THREAD_GAP_EVOLVE:
		glDecor->sync_live();
		ChargeMesh();
		break;
	case BT_FRISE_ADJUST:
		myDome.fFriseSpacing = (myDome.nFrise == 0 ? 0 : 1 / (float)myDome.nFrise);
		glDecor->sync_live();
		ChargeMesh();
		break;

	case BT_DECOR:
		btDecor->disable();
		glDecor->show();
		break;
	case BT_DECOR_QUIT:
		btDecor->enable();
		glDecor->hide();
		break;
	case SP_SNAKE_RADIUS:
	case LB_SNAKE_TOWARD:
	case LB_SNAKE_SHRINK:
	case ET_SNAKE_MATERIAL:

	case SP_FRISE_RADIUS:
	case LB_FRISE_SHRINK:
	case LB_FRISE_TOWARD:
	case ET_FRISE_MATERIAL:
	case ET_SPIRALE_MATERIAL:
	case LB_SPIRALE_SHRINK:
	case LB_SPIRALE_TOWARD:
	case SP_SPIRALE_RADIUS:
	case SP_COMB_RADIUS:
	case LB_COMB_SHRINK:
	case ET_COMB_MATERIAL:
	case LB_COMB_TOWARD:
		break;
		
	}
}


void myGlutKeyboard(unsigned char Key, int x, int y)
{
	switch (Key)
	{
	case 27:
	case 'q':
		exit(0);
		break;
	};

	glutPostRedisplay();
}

void myGlutMenu(int value)
{
	myGlutKeyboard(value, 0, 0);
}

void myGlutIdle(void)
{
	/* According to the GLUT specification, the current window is
	undefined during an idle callback.  So we need to explicitly change
	it if necessary */
	if (glutGetWindow() != MainWindow)
		glutSetWindow(MainWindow);

	/*  GLUI_Master.sync_live_all();  -- not needed - nothing to sync in this
	application  */

	glutPostRedisplay();
}
void myGlutMouse(int button, int button_state, int x, int y)
{
	nButtonState = -1;
//	printf_s("%d %d %s", button, button_state, CRLF);

	if (button == 3)
	{
		myDisplay.fScaleDisplay += (float)0.01;
		glDisplay->sync_live();
	}
	if (button == 4)
	{
		myDisplay.fScaleDisplay -= (float)0.01;
		glDisplay->sync_live();
	}

	//	if (button == GLUT_LEFT_BUTTON && button_state == GLUT_DOWN)
	if (button_state == GLUT_DOWN)
	{
		nButtonState = button;
		last_x = x;
		last_y = y;
	}
}

void myGlutMotion(int x, int y)
{
	if (nButtonState == GLUT_LEFT_BUTTON)
	{
		rotationX -= (float)(y - last_y);
		rotationY += (float)(x - last_x);
	}
	if (nButtonState == GLUT_RIGHT_BUTTON)
	{
		DomePosition[0] += (float)(x - last_x) / 350;
		DomePosition[1] -= (float)(y - last_y) / 350;
	}
	last_x = x;
	last_y = y;

	glutPostRedisplay();
}

/*void myGlutMouse(int button, int button_state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON && button_state == GLUT_DOWN)
	{
		last_x = x;
		last_y = y;
	}
}

void myGlutMotion(int x, int y)
{
	rotationX -= (float)(y - last_y);
	rotationY += (float)(x - last_x);

	last_x = x;
	last_y = y;

	glutPostRedisplay();
}*/

void myGlutReshape(int x, int y)
{
	int tx, ty, tw, th;
	GLUI_Master.get_viewport_area(&tx, &ty, &tw, &th);
	glViewport(tx, ty, tw, th);

	xy_aspect = (float)tw / (float)th;

	glutPostRedisplay();
}

void DrawBase()
{
	int nMaxI = ((myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides + 1);
	int nmi;

	if (myDome.nDomeEnd == 360)
		nmi = nMaxI;
	else
		nmi = (int)((float)(nMaxI * myDome.nDomeEnd) / 360.0) + 1;

	double  u, X, Y, Z;

	glLineWidth(3);

	glBegin(GL_LINE_STRIP);
	for (int ix = 0; ix <= nmi; ix++)
	{
		u = 2*_PI*((double)ix / (double)nMaxI);
		X = BaseFunction(u, myDome.nBaseSides, myDome.fBaseParameter, BASE_X, 1)/2;
		Y = BaseFunction(u, myDome.nBaseSides, myDome.fBaseParameter, BASE_Y, 1)/2;
		Z = 0;
//		Calcul( u, 1, &X, &Y, &Z);
		glVertex3f((GLfloat)X, (GLfloat)Y, (GLfloat)Z);
	}
	glEnd();

	if (myDome.fBaseRound != 0.0 || myDome.fBaseBaryCentre != 0.0 || myDome.fBaseMorphing != 0.0)
	{
		BaseFunctionSecondaryX = BaseFunctionListe[myDome.nBaseSecondary].FunctionX;
		BaseFunctionSecondaryY = BaseFunctionListe[myDome.nBaseSecondary].FunctionY;



//		XS = fExpansion * BaseFunctionSecondaryX(fmod(u + _PI * fRotation, 2 * _PI), a, b);
//		YS = fExpansion * BaseFunctionSecondaryY(fmod(u + _PI * fRotation, 2 * _PI), a, b);
//		Xc = fExpansion * BaseFunctionSecondaryX(u, a, b);
//		Yc = fExpansion * BaseFunctionSecondaryY(u, a, b);

		glBegin(GL_LINE_STRIP);
		for (int ix = 0; ix <= nmi; ix++)
		{
			u = 2 * _PI * ((double)ix / (double)nMaxI);
			X = BaseFunctionSecondaryX(u, myDome.nBaseSides, myDome.fBaseParameter) / 4;
			Y = BaseFunctionSecondaryY(u, myDome.nBaseSides, myDome.fBaseParameter) / 4;
			Z = 0;
			//		Calcul( u, 1, &X, &Y, &Z);
			glVertex3f((GLfloat)X, (GLfloat)Y, (GLfloat)Z);
		}
		glEnd();
	}
	glLineWidth(1);

}
void DrawSlope()
{
	int nMaxJ = myDome.nParallelMain + (myDome.nParallelMain - 1) * myDome.nParallelMinor;
	double  u, X, Y, Z;

	glLineWidth(3);

	glBegin(GL_LINE_STRIP);
	for (int ix = 0; ix <= nMaxJ; ix++)
	{
		u = (double)ix / (double)nMaxJ;
		Calcul(0, u, &X, &Y, &Z, false);
		glVertex3f((GLfloat)X, (GLfloat)Y, (GLfloat)Z);
	}
	glEnd();

	glBegin(GL_LINE_STRIP);
	for (int ix = 0; ix <= nMaxJ; ix++)
	{
		u = (double)ix / (double)nMaxJ;
		Calcul(0, u, &X, &Y, &Z, false);
		glVertex3f((GLfloat)-X, (GLfloat)-Y, (GLfloat)Z);
	}
	glEnd();
	glLineWidth(1);

}
void DrawWindow()
{
	int nMaxI = ((myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides + 1);
	double u, X, Y, Z;
	double fWindowHeight;

	glLineWidth(3);

	glBegin(GL_LINE_STRIP);
	for (int ix = 0; ix <= nMaxI; ix++)
	{
		u = 2*_PI*((double)ix / (double)nMaxI);
		fWindowHeight = ComputeWindowCoefficient((float)(u + 2 * myDome.fWindowShift * _PI / myDome.nBaseSides), (float)myDome.nBaseSides, (float)myDome.fWindowPerSide) * myDome.fWindowHeight;
		X = BaseFunction(u, myDome.nBaseSides, myDome.fBaseParameter, BASE_X, 1);
		Y = BaseFunction(u, myDome.nBaseSides, myDome.fBaseParameter, BASE_Y, 1);

		if (myDome.nBase != BASE_ROOF)
			X *= SlopeFunction(1, myDome.fSlopeParameterXY, SLOPE_XY);
		Y *= SlopeFunction(1, myDome.fSlopeParameterXY, SLOPE_XY);

		if (myDome.bWindowInverse)
			Z = (1.0 - SlopeFunction(1, myDome.fSlopeParameterZ, SLOPE_Z) * (1.0 - fWindowHeight * (1.0 - WindowFunctionZ(u + 2 * myDome.fWindowShift * _PI / myDome.nBaseSides, myDome.nBaseSides, myDome.fWindowPerSide, myDome.fWindowPower, false))));
		else
			Z = (1.0 - SlopeFunction(1, myDome.fSlopeParameterZ, SLOPE_Z) * (1.0 - fWindowHeight * WindowFunctionZ(u + 2 * myDome.fWindowShift * _PI / myDome.nBaseSides, myDome.nBaseSides, myDome.fWindowPerSide, myDome.fWindowPower, false)));
		glVertex3f((GLfloat)X/2, (GLfloat)Y/2, (GLfloat)Z/2);

	}
	glEnd();
	glLineWidth(1);

}


/************************************************** draw_axes() **********/
/* Disables lighting, then draws RGB axes                                */

void draw_axes(float scale)
{
	glDisable(GL_LIGHTING);

	glPushMatrix();
	glScalef(scale, scale, scale);

	glBegin(GL_LINES);

	glColor3f(1.0, 0.0, 0.0);
	glVertex3f(.8f, 0.05f, 0.0);  glVertex3f(1.0, 0.25f, 0.0); /* Letter X */
	glVertex3f(0.8f, .25f, 0.0);  glVertex3f(1.0, 0.05f, 0.0);
	glVertex3f(0.0, 0.0, 0.0);  glVertex3f(1.0, 0.0, 0.0); /* X axis      */

	glColor3f(0.0, 1.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);  glVertex3f(0.0, 1.0, 0.0); /* Y axis      */

	glColor3f(0.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, 0.0);  glVertex3f(0.0, 0.0, 1.0); /* Z axis    */
	glEnd();

	glPopMatrix();

	glEnable(GL_LIGHTING);
}

void myGlutDisplay(void)
{
	GLpoint Normal;
	Vector3D V;
	int nMaxI, nMaxJ, nmi;

	//	glClearColor( .78f, .78f, .78f, 1.0f );
	glClearColor(.5f, .5f, .5f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	GLfloat ambientColor[] = { 0.5f, 0.5f, 0.5f, 1.0f }; //Color (0.2, 0.2, 0.2)
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientColor);

//	if (myDisplay.nDisplayType == RB_FLAT)
	if(myDisplay.bFlatFace)
		glShadeModel(GL_FLAT);

//	if (myDisplay.nDisplayType == RB_SMOOTH)
	if (myDisplay.bSmooth)
		glShadeModel(GL_SMOOTH);



	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, myMaterial[myDisplay.nMaterial].MaterialDiffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, myMaterial[myDisplay.nMaterial].MaterialSpecular);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, myMaterial[myDisplay.nMaterial].MaterialShininess);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, myMaterial[myDisplay.nMaterial].MaterialAmbient);



	glLoadIdentity();
	glPushMatrix();
	glScalef(myDisplay.fScaleDisplay, myDisplay.fScaleDisplay, myDisplay.fScaleDisplay);

	glTranslatef(DomePosition[0], DomePosition[1], -DomePosition[2]);

	//    glRotatef( rotationY, 0.0, 1.0, 0.0 );
	//    glRotatef( -rotationX, 1.0, 0.0, 0.0 );


	for (int i = 0; i < 16; i++)
		InstantRotate[i] = 0;

	InstantRotate[0] = (float)cos(rotationY*deg2rad);
	InstantRotate[2] = (float)sin(rotationY*deg2rad);
	InstantRotate[5] = 1;
	InstantRotate[8] = (float)(-sin(rotationY*deg2rad));
	InstantRotate[10] = (float)cos(rotationY*deg2rad);
	InstantRotate[15] = 1;

	glMultMatrixf(InstantRotate);

	for (int i = 0; i < 16; i++)
		InstantRotate[i] = 0;

	InstantRotate[0] = 1;
	InstantRotate[5] = (float)cos(-rotationX * deg2rad);
	InstantRotate[6] = (float)(-sin(-rotationX * deg2rad));
	InstantRotate[9] = (float)sin(-rotationX * deg2rad);
	InstantRotate[10] = (float)cos(-rotationX * deg2rad);
	InstantRotate[15] = 1;

	glMultMatrixf(InstantRotate);
	glMultMatrixf(DomeRotate);

	glColor3f(0, 1, 1);
	nMaxI = (myDome.nMeridianMain + (myDome.nMeridianMain)*myDome.nMeridianMinor*myDisplay.bDisplayMinor) * myDome.nBaseSides + 1;
	nMaxJ = myDome.nParallelMain + (myDome.nParallelMain - 1)*myDome.nParallelMinor*myDisplay.bDisplayMinor;

	if (myDome.nDomeEnd == 360)
		nmi = nMaxI;
	else
		nmi = (int)((float)(nMaxI * myDome.nDomeEnd) / 360.0) + 1;

//	if (myDisplay.nDisplayType == RB_WIRE || 1)
	if(myDisplay.bWireFrame && myDome.bWireGenerate && myDisplay.nDisplay == DISPLAY_Dome)
	{
		for (int ix = 0; ix < nmi; ix++)
		{
			if (ix % (myDome.nMeridianMinor + 1) == 0 || myDisplay.bDisplayMinor == false)
				glLineWidth(3);

			glBegin(GL_LINE_STRIP);
			for (int jx = 0; jx < nMaxJ; jx++)
			{
				glVertex3fv((float *)&myOBJ_VERTICE[ix * nMaxJ + jx]);
			}
			glEnd();
			glLineWidth(1);
		}

		for (int jx = 0; jx < nMaxJ; jx++)
		{
			if (jx % (myDome.nParallelMinor + 1) == 0 || myDisplay.bDisplayMinor == false)
				glLineWidth(3);

			glBegin(GL_LINE_STRIP);
			for (int ix = 0; ix < nmi; ix++)
			{
				glVertex3fv((float *)&myOBJ_VERTICE[ix * nMaxJ + jx]);
			}
			glEnd();
			glLineWidth(1);
		}


		glLineWidth(3);
		glBegin(GL_LINES);
		glVertex3f(myOBJ_VERTICE[0].x, myOBJ_VERTICE[0].y, myOBJ_VERTICE[0].z);
		V.x = 0;
		V.y = 0;
		V.z = 0;
		SpikeLength(myOBJ_VERTICE, 0, 0, nmi, nMaxJ, &V);
		glVertex3f(myOBJ_VERTICE[0].x + .5f*(float)V.x, myOBJ_VERTICE[0].y + .5f*(float)V.y, myOBJ_VERTICE[0].z + .5f*(float)V.z);
		glEnd();
		 
		for (int ix = 0; ix < nmi -1; ix++)
			for (int jx = 1; jx < nMaxJ; jx++)
			{

				if(((jx % (myDome.nParallelMinor + 1) != 0) && (ix % (myDome.nMeridianMinor + 1) != 0)) && myDisplay.bDisplayMinor)
					glLineWidth(1);
				else if(((jx % (myDome.nParallelMinor + 1) != 0) && (ix % (myDome.nMeridianMinor + 1) == 0)) && myDisplay.bDisplayMinor)
					glLineWidth(2);
				else if (((jx % (myDome.nParallelMinor + 1) == 0) && (ix % (myDome.nMeridianMinor + 1) != 0)) && myDisplay.bDisplayMinor)
				glLineWidth(2);
				else
					glLineWidth(3);
				glBegin(GL_LINES);
				glVertex3f(myOBJ_VERTICE[ix * nMaxJ + jx].x, myOBJ_VERTICE[ix * nMaxJ + jx].y, myOBJ_VERTICE[ix * nMaxJ + jx].z);
				V.x = 0;
				V.y = 0;
				V.z = 0;
				SpikeLength(myOBJ_VERTICE, ix, jx, nmi, nMaxJ, &V);
				glVertex3f(myOBJ_VERTICE[ix * nMaxJ + jx].x + .5f*(float)V.x, myOBJ_VERTICE[ix * nMaxJ + jx].y + .5f*(float)V.y, myOBJ_VERTICE[ix * nMaxJ + jx].z + .5f*(float)V.z);
				glEnd();
			}
	}
//	if ((myDisplay.nDisplayType == RB_FLAT) || (myDisplay.nDisplayType == RB_SMOOTH))
	if((myDisplay.bSmooth || myDisplay.bFlatFace) && myDome.bFaceGenerate && myDisplay.nDisplay == DISPLAY_Dome)
	{
		for (int ix = 0; ix < nmi - 1; ix++)
		{
			glBegin(GL_TRIANGLE_STRIP);
			for (int jx = 0; jx < nMaxJ; jx++) // jx = 0
			{
				Normales(ix, jx % nMaxJ, nMaxI, nMaxJ, myOBJ_VERTICE, &Normal);
				glNormal3f(Normal.x, Normal.y, Normal.z);
				glVertex3f(myOBJ_VERTICE[ix * nMaxJ + (jx % nMaxJ)].x, myOBJ_VERTICE[ix * nMaxJ + (jx % nMaxJ)].y, myOBJ_VERTICE[ix * nMaxJ + (jx % nMaxJ)].z);
				Normales(ix + 1, jx % nMaxJ, nMaxI, nMaxJ, myOBJ_VERTICE, &Normal);
				glNormal3f(Normal.x, Normal.y, Normal.z);
				glVertex3f(myOBJ_VERTICE[(ix + 1) * nMaxJ + (jx % nMaxJ)].x, myOBJ_VERTICE[(ix + 1) * nMaxJ + (jx % nMaxJ)].y, myOBJ_VERTICE[(ix + 1) * nMaxJ + (jx % nMaxJ)].z);
			}
			glEnd();
		}
		glBegin(GL_TRIANGLE_FAN);
		Normales(0, 0, nMaxI, nMaxJ, myOBJ_VERTICE, &Normal);
		glNormal3f(Normal.x, Normal.y, Normal.z);
		glVertex3f(myOBJ_VERTICE[0].x, myOBJ_VERTICE[0].y, myOBJ_VERTICE[0].z);

		for (int ix = 1; ix < nMaxI - 1; ix++)
		{
			Normales(ix, 0, nMaxI, nMaxJ, myOBJ_VERTICE, &Normal);
			glNormal3f(Normal.x, Normal.y, Normal.z);
			glVertex3f(myOBJ_VERTICE[ix * nMaxJ].x, myOBJ_VERTICE[ix * nMaxJ].y, myOBJ_VERTICE[ix * nMaxJ].z);
		}
		glEnd();

		if (myDome.bLid)
		{
			GLpoint C0;
			C0.x = 0;
			C0.y = 0;
			C0.z = myDome.nWindow == WINDOW_NONE ? 0 : myDome.fWindowHeight / 2;
			glBegin(GL_TRIANGLE_FAN);
			glVertex3f(C0.x, C0.y, C0.z);
			for (int ix = 0; ix < nMaxI ; ix++)
				glVertex3f(myOBJ_VERTICE[ix * nMaxJ + nMaxJ - 1].x, myOBJ_VERTICE[ix * nMaxJ + nMaxJ - 1].y, myOBJ_VERTICE[ix * nMaxJ + nMaxJ - 1].z);
			glEnd();

		}



	}

	if (myDisplay.nDisplay == DISPLAY_Base)
		DrawBase();
	if (myDisplay.nDisplay == DISPLAY_Slope)
		DrawSlope();
	if (myDisplay.nDisplay == DISPLAY_Window)
		DrawWindow();
	if (myDisplay.nDisplay == DISPLAY_Snake)
		DrawSnake();
	if (myDisplay.nDisplay == DISPLAY_Spirale)
		DrawSpirale();
	if (myDisplay.nDisplay == DISPLAY_Frise)
		DrawFrise();
	if (myDisplay.nDisplay == DISPLAY_Comb)
		DrawComb();

	glColor3f(0.0, 0.0, 1.0);
	glLineWidth(4);
	if (myDome.nSpirale > 0 && myDisplay.nDisplay == DISPLAY_Dome)
		DrawSpirale();
	if (myDome.nFrise > 0 && myDisplay.nDisplay == DISPLAY_Dome)
		DrawFrise();
	if (myDome.nSnake > 0 && myDisplay.nDisplay == DISPLAY_Dome)
		DrawSnake();
	if (myDome.nComb > 0 && myDisplay.nDisplay == DISPLAY_Dome)
		DrawComb();

	glLineWidth(1);

			draw_axes(1);
	glPopMatrix();

	glutSwapBuffers();

}

void ChargeListeCouleur()
{
	CouleurBase[0] = STL_RGB(11, 0, 0);
	CouleurBase[1] = STL_RGB(24, 0, 0);
	CouleurBase[2] = STL_RGB(0, 11, 0);
	CouleurBase[3] = STL_RGB(26, 12, 31);
	CouleurBase[4] = STL_RGB(0, 0, 11);
	CouleurBase[5] = STL_RGB(0, 0, 28);
	CouleurBase[6] = STL_RGB(8, 8, 8);
	CouleurBase[7] = STL_RGB(11, 24, 0);
	CouleurBase[8] = STL_RGB(0, 11, 11);
	CouleurBase[9] = STL_RGB(0, 24, 24);
	CouleurBase[10] = STL_RGB(2, 31, 0);
	CouleurBase[11] = STL_RGB(8, 8, 0);
	CouleurBase[12] = STL_RGB(3, 0, 3);
	CouleurBase[13] = STL_RGB(31, 31, 0);
	CouleurBase[14] = STL_RGB(11, 11, 0);
	CouleurBase[15] = STL_RGB(5, 12, 24);
	CouleurBase[16] = STL_RGB(31, 31, 31);
	CouleurBase[17] = STL_RGB(24, 3, 24);
	CouleurBase[18] = 0xffff;

}
void InitCouleurs()
{
	int i;
	for (i = 0;CouleurBase[i] != 0xffff;i++)
	{
		stCouleur[i].sCouleur[0] = 0;
		stCouleur[i].bUsed = false;
		stCouleur[i].nCouleur = CouleurBase[i];
	}

	stCouleur[i].sCouleur[0] = 0;
	stCouleur[i].bUsed = false;
	stCouleur[i].nCouleur = CouleurBase[i];
}

unsigned short GetColor(char *sColor)
{

	if (!strncmp(sColor, "&_", strlen("&_")))
		return (unsigned short)atoi(sColor + strlen("&_"));
	for (int i = 0;stCouleur[i].nCouleur != 0xffff ;i++)
	{
		if (!strcmp(sColor, stCouleur[i].sCouleur))
		{
			stCouleur[i].bUsed = true;
			return stCouleur[i].nCouleur;
		}
	}

	for (int i = 0;stCouleur[i].nCouleur != 0xffff;i++)
	{
		if (stCouleur[i].sCouleur[0] == 0)
		{
			strcpy_s(stCouleur[i].sCouleur, sColor);
			stCouleur[i].bUsed = true;
			return stCouleur[i].nCouleur;
		}
	}

	for (int i = 0;stCouleur[i].nCouleur != 0xffff;i++)
	{
		if (stCouleur[i].bUsed == false)
		{
			strcpy_s(stCouleur[i].sCouleur, sColor);
			stCouleur[i].bUsed = true;
			return stCouleur[i].nCouleur;
		}
	}

	return 0;
}
/*
bool Checkxy(int x, int y, int x0, int y0, int x1, int y1) 
// ax +by +c=0  
{
	float z;

	float dx = (float)(x0 - x1);
	float dy = (float)(y0 - y1);

	float ddx = (float)(x - x1);
	float ddy = (float)(y - y1);

	float det = dx * ddy - dy * ddx;
	if (det != 0)
		return false;
	if (dx != 0)
		z = ddx / dx;
	else
		z = ddy / dy;
	if (floor(z) == z)
		return true;
	else
		return false;
}
*/


void lbSortItem(GLUI_Listbox *lb, int nbItem,  char* Structure, int nStructureSize, int nOffset)
{
	Tri** Liste;
	Tri* myTab;

	GLUI_String s;

	myTab = (Tri*)malloc(sizeof(Tri) * nbItem);
	Liste = (Tri**)malloc(sizeof(Tri*) * nbItem);
	if (myTab && Liste)
	{
		for (int i = 0; i < nbItem; i++)
		{
			myTab[i].Index = i;
			myTab[i].Mot = Structure + i * nStructureSize + nOffset;
			Liste[i] = &myTab[i];
			lb->delete_item(i);
		}


		Sort(Liste, nbItem);

		for (int i = 0; i < nbItem; i++)
			lb->add_item(Liste[i]->Index, Liste[i]->Mot);

		free(myTab);
		free(Liste);
	}
}


int main(int argc, char* argv[])
{
	/****************************************/
	/*   Initialize GLUT and create window  */
	/****************************************/
	char Buffer[501];


/*	bool test;
    test = Checkxy(0, 2, 0, 0, 0, 2);
	test = Checkxy(0, 0, 0, 0, 0, 2);
	test = Checkxy(0, 3, 0, 0, 0, 2);
	test = Checkxy(0, 4, 0, 0, 0, 2);
	test = Checkxy(1, 2, 0, 0, 0, 2);

	test = Checkxy(0, 0, 2, 0, 0, 0);
	test = Checkxy(2, 0, 2, 0, 0, 0);
	test = Checkxy(4, 0, 2, 0, 0, 0);
	test = Checkxy(5, 0, 2, 0, 0, 0);
	test = Checkxy(40, 1, 2, 0, 0, 0);

	test = Checkxy(0, 0, 2, 1, 0, 0);
	test = Checkxy(2, 1, 2, 1, 0, 0);
	test = Checkxy(4, 2, 2, 1, 0, 0);
	test = Checkxy(4, 1, 2, 1, 0, 0);
	test = Checkxy(4, 3, 2, 1, 0, 0);
	test = Checkxy(10, 5, 2, 1, 0, 0);
*/

	int   Largeur = 370;
//	Tri** Liste;
//	Tri* myTab;
	int col_x, col_y, col_w, col_h, col_x_off, col_y_off;

	nbBases = (int)(sizeof(BaseFunctionListe) / sizeof(BaseFunctionListe[0]));
	nbSlopes = (int)(sizeof(SlopeFunctionListe) / sizeof(SlopeFunctionListe[0]));
	nbWindows = (int)(sizeof(WindowFunctionListe) / sizeof(WindowFunctionListe[0]));
	nbTestSlopes = (int)(sizeof(TestSlopeFunctionListe) / sizeof(TestSlopeFunctionListe[0]));
	nbFiltres = (int)(sizeof(FiltreListe) / sizeof(FiltreListe[0]));
	nbTops = (int)(sizeof(TopListe) / sizeof(TopListe[0]));

	CUSTOM_SLOPE = nbSlopes - 1;


/*	for (int i = 0; i < nbBases; i++)
		printf("%s, %lf, %d\n", BaseFunctionListe[i].sNom, BaseFunctionListe[i].fBaseParameter, BaseFunctionListe[i].iMinSide);
	printf("\n");
	for (int i = 0; i < nbSlopes; i++)
		printf("%s, %lf, %lf\n", SlopeFunctionListe[i].sNom, SlopeFunctionListe[i].fSlopeParameterXY, SlopeFunctionListe[i].fSlopeParameterZ);
	for (int i = 0; i < nbWindows; i++)
		printf("%s\n", WindowFunctionListe[i].sNom);


	printf("Bases = %d		Slopes = %d		Windows = %d\n", nbBases, nbSlopes, nbWindows);*/

	ChargeListeCouleur();
	InitCouleurs();
	InitParameters();
	ChargeShrink(myDome.nWireShrinkingType);
	if ("Glui")
	{//	int nbColors = ChargeCouleur();
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
		glutInitWindowPosition(120, 10);
		//	glutInitWindowSize(1024, 768);
		glutInitWindowSize(1120, 818);

		MainWindow = glutCreateWindow("Domenator Version 13.0");
		glutDisplayFunc(myGlutDisplay);
		GLUI_Master.set_glutReshapeFunc(myGlutReshape);
		GLUI_Master.set_glutKeyboardFunc(myGlutKeyboard);
		GLUI_Master.set_glutSpecialFunc(NULL);
		GLUI_Master.set_glutMouseFunc(myGlutMouse);
		glutMotionFunc(myGlutMotion);


		//	gluPerspective(65.0, 1024 / 768, 0, 300);
			/****************************************/
			/*       Set up OpenGL lights           */
			/****************************************/

		glEnable(GL_LIGHTING);
		glEnable(GL_NORMALIZE);

		glEnable(GL_LIGHT0);
		glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
		glLightfv(GL_LIGHT0, GL_POSITION, light0_position);

		glEnable(GL_LIGHT1);
		glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
		glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
		glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

		glEnable(GL_LIGHT2);
		glLightfv(GL_LIGHT2, GL_AMBIENT, light2_ambient);
		glLightfv(GL_LIGHT2, GL_DIFFUSE, light2_diffuse);
		glLightfv(GL_LIGHT2, GL_POSITION, light2_position);

		/****************************************/
		/*          Enable z-buferring          */
		/****************************************/

		glEnable(GL_DEPTH_TEST);

		/****************************************/
		/*         Here's the GLUI code         */
		/****************************************/

		printf("GLUI version: %3.2f\n", GLUI_Master.get_version());

		/*** Create the side subwindow ***/

		glMessage = GLUI_Master.create_glui("Message", 0, 400, 50);
		new GLUI_StaticText(glMessage, "                         ");
		stMessage = new GLUI_StaticText(glMessage, "     Message test     ");
		new GLUI_StaticText(glMessage, "                         ");

		GLUI_Panel *pnInvisibleMessage = new GLUI_Panel(glMessage, "", GLUI_PANEL_NONE);

		GLUI_Button *btQuit = new GLUI_Button(pnInvisibleMessage, "Quit", BT_QUIT_MESSAGE, ProcessMessage);
		btQuit->set_w(80);
		glMessage->add_column_to_panel(pnInvisibleMessage, false);

		GLUI_Button *btOK = new GLUI_Button(pnInvisibleMessage, "OK", BT_OK_MESSAGE, ProcessMessage);
		btOK->set_w(80);
		btOK->hidden = true;
		glMessage->hide();


	}
	if ("glPreProcessing")
	{
		glPreProcess = GLUI_Master.create_glui("Pre Processing", 0, 900, 50);
		if ("Derivative")
		{
			GLUI_Panel *pnDerivative = new GLUI_Panel(glPreProcess, "Derivative");
			pnDerivative->set_w(370);
			GLUI_Panel* pnDerivativeInv = new GLUI_Panel(pnDerivative, "", GLUI_PANEL_NONE);
			GLUI_Spinner *spDerivative = new GLUI_Spinner(pnDerivativeInv, "     Lambda:", GLUI_SPINNER_FLOAT, &myDome.fDerivative, SP_DERIVATIVE, ControlPreProcessing);
			spDerivative->set_float_limits(-10, 10);
			spDerivative->set_alignment(GLUI_ALIGN_RIGHT);
			glPreProcess->add_column_to_panel(pnDerivativeInv, true);

			new GLUI_Checkbox(pnDerivativeInv, "Normalized    ", &myDome.bDerivativeNorme, CK_DERIVATIVE_NORME, ControlPreProcessing);
	
			new GLUI_StaticText(pnDerivative, "                                                                            ");

		}
		if ("Noise")
		{
			GLUI_Panel *roNoise = new GLUI_Panel(glPreProcess, "Noise");
			roNoise->set_w(Largeur);
			GLUI_Spinner *spNoiseOctave = new GLUI_Spinner(roNoise, "        Octave:", GLUI_SPINNER_INT, &myDome.nNoiseOctave, SP_NOISE_OCTAVE, ControlNoise);
			spNoiseOctave->set_int_limits(1, 16);
			spNoiseOctave->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spNoiseDensity = new GLUI_Spinner(roNoise, "        Density:", GLUI_SPINNER_FLOAT, &myDome.fNoiseDensity, SP_NOISE_DENSITY, ControlNoise);
			spNoiseDensity->set_float_limits(1, 16);
			spNoiseDensity->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spNoiseSeed = new GLUI_Spinner(roNoise, "           Seed:", GLUI_SPINNER_INT, &myDome.nNoiseSeed, SP_NOISE_SEED, ControlNoise);
			spNoiseSeed->set_int_limits(0, 32000);
			spNoiseSeed->set_alignment(GLUI_ALIGN_RIGHT);
			spNoiseSeed->set_speed((float)(0.1*myDisplay.fSpeedCoefficient));
			GLUI_Spinner *spNoiseFrequencyU = new GLUI_Spinner(roNoise, "Frequency U:", GLUI_SPINNER_FLOAT, &myDome.fNoiseFrequencyU, SP_NOISE_FREQUENCY_U, ControlNoise);
			spNoiseFrequencyU->set_float_limits(0, 16);
			spNoiseFrequencyU->set_alignment(GLUI_ALIGN_RIGHT);
			spNoiseFrequencyU->set_speed((float)(1 * myDisplay.fSpeedCoefficient));
			GLUI_Spinner *spNoiseFrequencyV = new GLUI_Spinner(roNoise, "Frequency V:", GLUI_SPINNER_FLOAT, &myDome.fNoiseFrequencyV, SP_NOISE_FREQUENCY_V, ControlNoise);
			spNoiseFrequencyV->set_float_limits(0, 16);
			spNoiseFrequencyV->set_alignment(GLUI_ALIGN_RIGHT);
			spNoiseFrequencyV->set_speed((float)(1 * myDisplay.fSpeedCoefficient));

			glParameter->add_column_to_panel(roNoise, true);

			GLUI_Spinner *spNoiseAmplitudeR = new GLUI_Spinner(roNoise, "   Amplitude R:", GLUI_SPINNER_FLOAT, &myDome.fNoiseAmplitudeR, SP_NOISE_AMPLITUDE_R, ControlNoise);
			spNoiseAmplitudeR->set_float_limits(-100, 100);
			spNoiseAmplitudeR->set_speed((float)(0.1*myDisplay.fSpeedCoefficient));
			spNoiseAmplitudeR->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spNoiseAmplitudeZ = new GLUI_Spinner(roNoise, "   Amplitude Z:", GLUI_SPINNER_FLOAT, &myDome.fNoiseAmplitudeZ, SP_NOISE_AMPLITUDE_Z, ControlNoise);
			spNoiseAmplitudeZ->set_float_limits(-100, 100);
			spNoiseAmplitudeZ->set_speed((float)(0.1*myDisplay.fSpeedCoefficient));
			spNoiseAmplitudeZ->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spNoiseProgressiveR = new GLUI_Spinner(roNoise, "Progressive R:", GLUI_SPINNER_FLOAT, &myDome.fNoiseProgressiveR, SP_NOISE_PROGRESSIVE_R, ControlNoise);
			spNoiseProgressiveR->set_float_limits(-10, 10);
			spNoiseProgressiveR->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spNoiseProgressiveZ = new GLUI_Spinner(roNoise, "Progressive Z:", GLUI_SPINNER_FLOAT, &myDome.fNoiseProgressiveZ, SP_NOISE_PROGRESSIVE_Z, ControlNoise);
			spNoiseProgressiveZ->set_float_limits(-10, 10);
			spNoiseProgressiveZ->set_alignment(GLUI_ALIGN_RIGHT);



			//  GLUI_Panel *InvisibleNoise = new GLUI_Panel(Noise, "", GLUI_PANEL_NONE);
			new GLUI_Button(roNoise, "Default", BT_NOISE_DEFAULT, ControlNoise);
		}
		GLUI_Panel *pnPreProcessQuit = new GLUI_Panel(glPreProcess, "", GLUI_PANEL_NONE);
		btPreProcessingQuit = new GLUI_Button(pnPreProcessQuit, "Quit", BT_PRE_PROCESSING_QUIT, ControlPreProcessing);
		btPreProcessingQuit->set_alignment(GLUI_ALIGN_CENTER);
		glPreProcess->hide();

	}
	if ("glPostProcessing")
	{
		char sStatic[] = "                                                                                 ";
		glPostProcess = GLUI_Master.create_glui("Post Processing", 0, 700, 50);

		GLUI_Panel* pnPostProcess = new GLUI_Panel(glPostProcess, "Post", GLUI_PANEL_NONE);

		if ("Smoothing")
		{

			GLUI_Panel *pnLissage = new GLUI_Panel(pnPostProcess, "Smoothing");
			pnLissage->set_w(314);
//			new GLUI_StaticText(pnLissage, sStatic);

			GLUI_Checkbox *ckSmoothing = new GLUI_Checkbox(pnLissage, "Smoothing", &myDome.bLissage, CK_LISSAGE, ControlPostProcessing);
			ckSmoothing->set_alignment(GLUI_ALIGN_CENTER);
			GLUI_Panel *pnInvisibleLissage = new GLUI_Panel(pnLissage, "", GLUI_PANEL_NONE);


			lbFiltreType = new GLUI_Listbox(pnInvisibleLissage, "   Filtre:", &myDome.nLissageTypeFiltre, LB_LISSAGE_FILTRE, ControlPostProcessing);
			lbSortItem(lbFiltreType, nbFiltres, (char*)FiltreListe, sizeof(FiltreListe[0]), FiltreListe[0].sNom - (char*)FiltreListe);
			lbFiltreType->set_alignment(GLUI_ALIGN_RIGHT);
			lbFiltreType->set_int_val(FILTER_Gauss);
//£
			GLUI_Spinner* spLissageSize = new GLUI_Spinner(pnInvisibleLissage, "        Size:", GLUI_SPINNER_INT, &myDome.nLissageTaille, SP_LISSAGE_TAILLE, ControlPostProcessing);
			spLissageSize->set_int_limits(1, 16);
			spLissageSize->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spLissagePuissance = new GLUI_Spinner(pnInvisibleLissage, "     Power:", GLUI_SPINNER_FLOAT, &myDome.fLissagePuissance, SP_LISSAGE_PUISSANCE, ControlPostProcessing);
			spLissagePuissance->set_int_limits(0, 20);
			spLissagePuissance->set_alignment(GLUI_ALIGN_RIGHT);
			
			spLissageSigma1 = new GLUI_Spinner(pnInvisibleLissage, "   Sigma 1:", GLUI_SPINNER_FLOAT, &myDome.fLissageSigma1, SP_LISSAGE_SIGMA_1, ControlPostProcessing);
			spLissageSigma1->set_float_limits(0, 20);
			spLissageSigma1->set_alignment(GLUI_ALIGN_RIGHT);

			spPoidsSelf = new GLUI_Spinner(pnInvisibleLissage, "         Self:", GLUI_SPINNER_FLOAT, &myDome.fLissagePoidsSelf, SP_LISSAGE_POIDS_SELF, ControlPostProcessing);
			spPoidsSelf->set_float_limits(-100, 100);
			spPoidsSelf->set_alignment(GLUI_ALIGN_RIGHT);

			spPoidsVertical = new GLUI_Spinner(pnInvisibleLissage, "    Vertical:", GLUI_SPINNER_FLOAT, &myDome.fLissagePoidsVertical, SP_LISSAGE_POIDS_VERTICAL, ControlPostProcessing);
			spPoidsVertical->set_float_limits(-100, 100);
			spPoidsVertical->set_alignment(GLUI_ALIGN_RIGHT);



			glParameter->add_column_to_panel(pnInvisibleLissage, true);

			GLUI_Listbox* lbDistanceType = new GLUI_Listbox(pnInvisibleLissage, "   Distance:", &myDome.nLissageDistanceType, LB_LISSAGE_DISTANCE_TYPE, ControlPostProcessing);
			for (int i = 0; i < 5; i++)
				lbDistanceType->add_item(i, DistanceListe[i]);
			lbDistanceType->set_int_val(DISTANCE_EUCLIDE);
			lbDistanceType->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spLissageNombre = new GLUI_Spinner(pnInvisibleLissage, "      Times", GLUI_SPINNER_INT, &myDome.nLissageNombre, SP_LISSAGE_NOMBRE, ControlPostProcessing);
			spLissageNombre->set_int_limits(1, 16);
			spLissageNombre->set_alignment(GLUI_ALIGN_RIGHT);


			lbLissageDirection = new GLUI_Listbox(pnInvisibleLissage, "  Direction", &myDome.nLissageDirection, LB_LISSAGE_DIRECTION, ControlPostProcessing);
			lbLissageDirection->add_item(DIRECTION_NORTH, "North");
			lbLissageDirection->add_item(DIRECTION_WEST, "West");
			lbLissageDirection->add_item(DIRECTION_SOUTH, "South");
			lbLissageDirection->add_item(DIRECTION_EAST, "East");
			lbLissageDirection->add_item(DIRECTION_NORTH_WEST, "North-West");
			lbLissageDirection->add_item(DIRECTION_NORTH_EAST, "North-East");
			lbLissageDirection->add_item(DIRECTION_SOUTH_WEST, "South-West");
			lbLissageDirection->add_item(DIRECTION_SOUTH_EAST, "South-East");
			lbLissageDirection->set_alignment(GLUI_ALIGN_RIGHT);
			
			
			spLissageSigma2 = new GLUI_Spinner(pnInvisibleLissage, "  Sigma 2:", GLUI_SPINNER_FLOAT, &myDome.fLissageSigma2, SP_LISSAGE_SIGMA_2, ControlPostProcessing);
			spLissageSigma2->set_float_limits(0, 20);
			spLissageSigma2->set_alignment(GLUI_ALIGN_RIGHT);
			spPoidsHorizontal = new GLUI_Spinner(pnInvisibleLissage, "Horizontal:", GLUI_SPINNER_FLOAT, &myDome.fLissagePoidsHorizontal, SP_LISSAGE_POIDS_HORIZONTAL, ControlPostProcessing);
			spPoidsHorizontal->set_float_limits(-100, 100);
			spPoidsHorizontal->set_alignment(GLUI_ALIGN_RIGHT);
			spPoidsDiagonales = new GLUI_Spinner(pnInvisibleLissage, " Diagonal:", GLUI_SPINNER_FLOAT, &myDome.fLissagePoidsDiagonales, SP_LISSAGE_POIDS_DIAGONALES, ControlPostProcessing);
			spPoidsDiagonales->set_float_limits(-100, 100);
			spPoidsDiagonales->set_alignment(GLUI_ALIGN_RIGHT);



			//	new GLUI_Panel(pnInvisibleLissage, "", GLUI_PANEL_NONE);
			
//			new GLUI_Button(pnInvisibleLissage, "Laplacian", BT_LISSAGE_LAPLACIEN, ControlPostProcessing);
//			new GLUI_Button(pnInvisibleLissage, "Gaussian", BT_LISSAGE_GAUSSIEN, ControlPostProcessing);

//			glParameter->add_column_to_panel(pnInvisibleLissage, true);


			//	new GLUI_Panel(pnInvisibleLissage, "", GLUI_PANEL_NONE);
			new GLUI_Button(pnLissage, "Default", BT_LISSAGE_DEFAULT, ControlPostProcessing);
			spLissageSigma1->enable();
			spLissageSigma2->disable();
			lbLissageDirection->disable();
			spPoidsSelf->disable();
			spPoidsVertical->disable();
			spPoidsHorizontal->disable();
			spPoidsDiagonales->disable();


		}

		if ("Differential Geometry")
		{
			GLUI_Panel*pnDiffentialGeometry = new GLUI_Panel(pnPostProcess, "Differential Geometry");
			GLUI_Listbox *lbDifferentialGeometry = new GLUI_Listbox(pnDiffentialGeometry, "   Operator:", &myDome.nDifferentialGeometry, LB_DIFF_GEOM, ControlPostProcessing);
			lbDifferentialGeometry->add_item(DIFF_GEOM_EVOLUTE, "Evolute    ");
			lbDifferentialGeometry->add_item(DIFF_GEOM_PEDAL, "Pedal");
			lbDifferentialGeometry->add_item(DIFF_GEOM_RADIAL, "Radial");
			lbDifferentialGeometry->set_alignment(GLUI_ALIGN_CENTER);

			GLUI_Panel *pnEvolute = new GLUI_Panel(pnDiffentialGeometry, "", GLUI_PANEL_NONE);
			pnEvolute->set_w(314);
//			new GLUI_StaticText(pnEvolute, sStatic);

			GLUI_Panel* pnEvoluteInv = new GLUI_Panel(pnEvolute, "", GLUI_PANEL_NONE);

			GLUI_Checkbox *ckEvoluteMeridien = new GLUI_Checkbox(pnEvoluteInv, "Meridian             ", &myDome.bDiffGeomMeridien, CK_EVOLUTE_MERIDIEN, ControlPostProcessing);
			ckEvoluteMeridien->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spEvoluteLambda = new GLUI_Spinner(pnEvoluteInv, "     Lambda:", GLUI_SPINNER_FLOAT, &myDome.fDiffGeomLambda, SP_EVOLUTE_LAMBDA, ControlPostProcessing);
			spEvoluteLambda->set_float_limits(-5, 5);
			spEvoluteLambda->set_alignment(GLUI_ALIGN_RIGHT);

			glParameter->add_column_to_panel(pnEvoluteInv, true);

			GLUI_Checkbox *ckEvoluteParallel = new GLUI_Checkbox(pnEvoluteInv, "Parallel            ", &myDome.bDiffGeomParallel, CK_EVOLUTE_PARALLEL, ControlPostProcessing);
			ckEvoluteParallel->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spEvoluteNombre = new GLUI_Spinner(pnEvoluteInv, "         Times:", GLUI_SPINNER_INT, &myDome.nDiffGeomNombre, SP_EVOLUTE_NOMBRE, ControlPostProcessing);
			spEvoluteNombre->set_int_limits(1, 16);
			spEvoluteNombre->set_alignment(GLUI_ALIGN_RIGHT);
		}

		if ("Enhance")
		{
			GLUI_Panel* pnEnhance = new GLUI_Panel(pnPostProcess, "Enhance");
			GLUI_Listbox* lbEnhanceType = new GLUI_Listbox(pnEnhance, "Type:", &myDome.nEnhanceType, LB_ENHANCE_TYPE, ControlPostProcessing);
			lbEnhanceType->set_alignment(GLUI_ALIGN_CENTER);
			lbEnhanceType->add_item(ENHANCE_NONE, "        None    ");
			lbEnhanceType->add_item(ENHANCE_RADIAL, "Radial");
			lbEnhanceType->add_item(ENHANCE_VERTICAL, "Vertical");
			lbEnhanceType->add_item(ENHANCE_CYLINDRICAL, "Cylindrical");
			lbEnhanceType->add_item(ENHANCE_OB_OA, "OB_OA");
			lbEnhanceType->add_item(ENHANCE_DB_AC, "DB_AC");
			lbEnhanceType->add_item(ENHANCE_D_B_O, "D_B_O");
			lbEnhanceType->add_item(ENHANCE_A_C_O, "A_C_O");
			lbEnhanceType->add_item(ENHANCE_ABCD_O, "ABCD_O");
			lbEnhanceType->add_item(ENHANCE_OA_OB_OA_OD, "OA_OB_OA_OD");
			lbEnhanceType->set_int_val(ENHANCE_NONE);

			GLUI_Panel* pnEnhanceParameters = new GLUI_Panel(pnEnhance, "", GLUI_PANEL_NONE);

			GLUI_Checkbox* ckEnhanceMeridien = new GLUI_Checkbox(pnEnhanceParameters, "Meridian             ", &myDome.bEnhanceMeridian, CK_ENHANCE_MERIDIAN, ControlPostProcessing);
			ckEnhanceMeridien->set_alignment(GLUI_ALIGN_RIGHT);
			GLUI_Listbox* lbEnhanceMerge = new GLUI_Listbox(pnEnhanceParameters, "Merge:", &myDome.nEnhanceMerge, LB_ENHANCE_MERGE, ControlPostProcessing);
			lbEnhanceMerge->add_item(MERGE_ZERO, "Flat             ");
			lbEnhanceMerge->add_item(MERGE_UN, "Maximum");
			lbEnhanceMerge->add_item(MERGE_SOMME, "Sum");
			lbEnhanceMerge->set_alignment(GLUI_ALIGN_RIGHT);
			lbEnhanceMerge->set_int_val(MERGE_SOMME);


			GLUI_Spinner* spEnhanceLambda = new GLUI_Spinner(pnEnhanceParameters, "Lambda:", GLUI_SPINNER_FLOAT, &myDome.fEnhanceLambda, SP_ENHANCE_LAMBDA, ControlPostProcessing);
			spEnhanceLambda->set_float_limits(-10, 10);
			spEnhanceLambda->set_alignment(GLUI_ALIGN_RIGHT);

			glParameter->add_column_to_panel(pnEnhanceParameters, true);

			GLUI_Checkbox* ckEnhanceParallel = new GLUI_Checkbox(pnEnhanceParameters, "Parallel             ", &myDome.bEnhanceParallel, CK_ENHANCE_PARALLEL, ControlPostProcessing);
			ckEnhanceParallel->set_alignment(GLUI_ALIGN_RIGHT);
			GLUI_Listbox* lbEnhanceCalcul = new GLUI_Listbox(pnEnhanceParameters, "   Calculation:", &myDome.nEnhanceCalcul, LB_ENHANCE_CALCUL, ControlPostProcessing);
			lbEnhanceCalcul->add_item(CALCUL_Somme, "Sum         ");
			lbEnhanceCalcul->add_item(CALCUL_Produit, "Product");
			lbEnhanceCalcul->set_int_val(CALCUL_Produit);
			lbEnhanceCalcul->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spEnhanceDamping = new GLUI_Spinner(pnEnhanceParameters, "   Damping:", GLUI_SPINNER_FLOAT, &myDome.fEnhanceDamping, SP_ENHANCE_DAMPING, ControlPostProcessing);
			spEnhanceDamping->set_float_limits(-1, 1);
			spEnhanceDamping->set_alignment(GLUI_ALIGN_RIGHT);


		}
		glParameter->add_column_to_panel(pnPostProcess, true);

		if ("Normale")
		{
			GLUI_Panel *pnNormale = new GLUI_Panel(pnPostProcess, "Normale");
			pnNormale->set_w(314);
//			new GLUI_StaticText(pnNormale, sStatic);

			GLUI_Checkbox *ckNormale = new GLUI_Checkbox(pnNormale, "Normale", &myDome.bNormale, CK_NORMALE, ControlPostProcessing);
			ckNormale->set_alignment(GLUI_ALIGN_CENTER);

			GLUI_Panel *pnNormaleParameters = new GLUI_Panel(pnNormale, "", GLUI_PANEL_NONE);

			GLUI_Spinner *spNormaleLambda = new GLUI_Spinner(pnNormaleParameters, "   Lambda:", GLUI_SPINNER_FLOAT, &myDome.fNormaleLambda, SP_NORMALE_LAMBDA, ControlPostProcessing);
			spNormaleLambda->set_float_limits(-15, 15);
			spNormaleLambda->set_alignment(GLUI_ALIGN_RIGHT);

			glParameter->add_column_to_panel(pnNormaleParameters, true);

			GLUI_Spinner *spNormaleNombre = new GLUI_Spinner(pnNormaleParameters, "       Times:", GLUI_SPINNER_INT, &myDome.nNormaleNombre, SP_NORMALE_NOMBRE, ControlPostProcessing);
			spNormaleNombre->set_int_limits(1, 16);
			spNormaleNombre->set_alignment(GLUI_ALIGN_RIGHT);

		}

		if ("Bumps")
		{
			GLUI_Panel *pnBumps = new GLUI_Panel(pnPostProcess, "Bumps");
			pnBumps->set_w(365);
//			new GLUI_StaticText(pnBumps, sStatic);

			GLUI_Checkbox *ckBumps = new GLUI_Checkbox(pnBumps, "Bumps", &myDome.bBump, CK_BUMP, ControlPostProcessing);
			ckBumps->set_alignment(GLUI_ALIGN_CENTER);

			GLUI_Listbox* lbBumpType = new GLUI_Listbox(pnBumps, "                     Type:", &myDome.nBumpType, LB_BUMP_TYPE, ControlPostProcessing);
			lbBumpType->add_item(BUMP_FUNNEL, "Funnel             ");
			lbBumpType->add_item(BUMP_WAVE, "Wave");
			lbBumpType->add_item(BUMP_FOURPEAKS, "Four Peaks");
			lbBumpType->add_item(BUMP_ASTROID, "Astroïd");
			lbBumpType->add_item(BUMP_ELLIPSOID, "Ellipse");
			lbBumpType->add_item(BUMP_TRACTOID, "Tractoïd");
			lbBumpType->add_item(BUMP_CONE, "Cone");
			lbBumpType->set_int_val(BUMP_FUNNEL);
			lbBumpType->set_alignment(GLUI_ALIGN_LEFT);

			GLUI_Panel *pnBumpsParameters = new GLUI_Panel(pnBumps, "", GLUI_PANEL_NONE);

			GLUI_Listbox *lbBumpNormal = new GLUI_Listbox(pnBumpsParameters, "Normal:", &myDome.nBumpNormal, LB_BUMP_NORMAL, ControlPostProcessing);
			lbBumpNormal->add_item(BUMP_1VECTEUR, "1 Vector");
			lbBumpNormal->add_item(BUMP_2VECTEUR, "2 Vectors   ");
			lbBumpNormal->add_item(BUMP_MEAN, "Mean");
			lbBumpNormal->set_int_val(BUMP_MEAN);
			lbBumpNormal->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spBumpPower = new GLUI_Spinner(pnBumpsParameters, "     Power:", GLUI_SPINNER_FLOAT, &myDome.fBumpPower, SP_BUMP_POWER, ControlPostProcessing);
			spBumpPower->set_float_limits(-10, 10);
			spBumpPower->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Listbox *lbBumpDistanceType = new GLUI_Listbox(pnBumpsParameters, "Distance:", &myDome.nBumpDistanceType, LB_BUMP_DISTANCE_TYPE, ControlPostProcessing);
			lbBumpDistanceType->add_item(DISTANCE_MANHATTAN, "Manhattan");
			lbBumpDistanceType->add_item(DISTANCE_EUCLIDE, "Euclidian");
			lbBumpDistanceType->add_item(DISTANCE_3, "Minkowski");
			lbBumpDistanceType->add_item(DISTANCE_ULTRA, "Chebyshev");
			lbBumpDistanceType->add_item(DISTANCE_DISCRETE, "Discrete");
			lbBumpDistanceType->set_int_val(DISTANCE_EUCLIDE);
			lbBumpDistanceType->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Checkbox *ckBumpNormalize = new GLUI_Checkbox(pnBumpsParameters, "Normalized", &myDome.nBumpNormalize, CK_BUMP_NORMALIZE, ControlPostProcessing);
			ckBumpNormalize->set_alignment(GLUI_ALIGN_LEFT);

/*			GLUI_Spinner* spBumpParallelEvery = new GLUI_Spinner(pnBumpsParameters, "  P. Every:", GLUI_SPINNER_INT, &myDome.nBumpParallelEvery, SP_BUMP_PARALLEL_EVERY, ControlPostProcessing);
			spBumpParallelEvery->set_float_limits(1, 30);
			spBumpParallelEvery->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spBumpParallelFirst = new GLUI_Spinner(pnBumpsParameters, "    P. First:", GLUI_SPINNER_INT, &myDome.nBumpParallelFirst, SP_BUMP_PARALLEL_FIRST, ControlPostProcessing);
			spBumpParallelFirst->set_float_limits(1, 30);
			spBumpParallelFirst->set_alignment(GLUI_ALIGN_RIGHT);
*/
			glParameter->add_column_to_panel(pnBumpsParameters, true);

			GLUI_Spinner *spBumpLambda = new GLUI_Spinner(pnBumpsParameters, "   Lambda:", GLUI_SPINNER_FLOAT, &myDome.fBumpLambda, SP_BUMP_LAMBDA, ControlPostProcessing);
			spBumpLambda->set_float_limits(-15, 15);
			spBumpLambda->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spBumpDistance = new GLUI_Spinner(pnBumpsParameters, "  Minimum:", GLUI_SPINNER_FLOAT, &myDome.fBumpDistance, SP_BUMP_DISTANCE, ControlPostProcessing);
			spBumpDistance->set_float_limits(0, 1);
			spBumpDistance->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spBumpDamping = new GLUI_Spinner(pnBumpsParameters, "  Damping:", GLUI_SPINNER_FLOAT, &myDome.fBumpDamping, SP_BUMP_DAMPING, ControlPostProcessing);
			spBumpDamping->set_float_limits(-1, 1);
			spBumpDamping->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Checkbox *ckBumpInverse = new GLUI_Checkbox(pnBumpsParameters, "Inverse", &myDome.nBumpDistanceInverse, CK_BUMP_DISTANCE_INVERSE, ControlPostProcessing);
			ckBumpInverse->set_alignment(GLUI_ALIGN_LEFT);

/*			GLUI_Spinner* spBumpMeridianEvery = new GLUI_Spinner(pnBumpsParameters, "  M. Every:", GLUI_SPINNER_INT, &myDome.nBumpMeridianEvery, SP_BUMP_MERIDIAN_EVERY, ControlPostProcessing);
			spBumpMeridianEvery->set_float_limits(1, 30);
			spBumpMeridianEvery->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spBumpMeridianDelay = new GLUI_Spinner(pnBumpsParameters, "  M. Delay:", GLUI_SPINNER_INT, &myDome.nBumpMeridianDelay, SP_BUMP_MERIDIAN_DELAY, ControlPostProcessing);
			spBumpMeridianDelay->set_float_limits(0, 30);
			spBumpMeridianDelay->set_alignment(GLUI_ALIGN_RIGHT);
*/


		}

		if ("Marquees")
		{
			GLUI_Panel* pnMarquees = new GLUI_Panel(pnPostProcess, "Marquees");
			pnMarquees->set_w(365);

			GLUI_Checkbox* ckMarquees = new GLUI_Checkbox(pnMarquees, "Marquees", &myDome.bMarquee, CK_MARQUEES, ControlPostProcessing);
			ckMarquees->set_alignment(GLUI_ALIGN_CENTER);

			GLUI_Panel* pnMarqueesParameters = new GLUI_Panel(pnMarquees, "", GLUI_PANEL_NONE);

			GLUI_Listbox* lbMarqueeFunction = new GLUI_Listbox(pnMarqueesParameters, "Function:", &myDome.nMarqueeFunction, LB_MARQUEE_FUNCTION, ControlPostProcessing);
			lbSortItem(lbMarqueeFunction, nbTestSlopes, (char*)TestSlopeFunctionListe, sizeof(TestSlopeFunctionListe[0]), TestSlopeFunctionListe[0].sNom - (char*)TestSlopeFunctionListe);
			lbMarqueeFunction->set_int_val(SLOPE_STRAIGHT);

			GLUI_Spinner* spMarqueeParameter = new GLUI_Spinner(pnMarqueesParameters, "Parameter:", GLUI_SPINNER_FLOAT, &myDome.fMarqueeParameter, SP_MARQUEE_PARAMETER, ControlPostProcessing);
			spMarqueeParameter->set_float_limits(-5.0, 10.0);
			spMarqueeParameter->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spMarqueeLambda = new GLUI_Spinner(pnMarqueesParameters, "   Lambda:", GLUI_SPINNER_FLOAT, &myDome.fMarqueeLambda, SP_MARQUEE_LAMBDA, ControlPostProcessing);
			spMarqueeLambda->set_float_limits(-15.0, 20.0);
			spMarqueeLambda->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Checkbox* ckMarqueeNormalize = new GLUI_Checkbox(pnMarqueesParameters, "Normalized", &myDome.bMarqueeNormalize, CK_MARQUEE_NORMALIZE, ControlPostProcessing);
			ckMarqueeNormalize->set_alignment(GLUI_ALIGN_LEFT);

			glParameter->add_column_to_panel(pnMarqueesParameters, true);

			GLUI_Spinner* spMarqueeAlpha = new GLUI_Spinner(pnMarqueesParameters, "Alpha:", GLUI_SPINNER_FLOAT, &myDome.fMarqueeAlpha, SP_MARQUEE_ALPHA, ControlPostProcessing);
			spMarqueeAlpha->set_float_limits(-1, 1);
			spMarqueeAlpha->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spMarqueeBeta = new GLUI_Spinner(pnMarqueesParameters, "Beta:", GLUI_SPINNER_FLOAT, &myDome.fMarqueeBeta, SP_MARQUEE_BETA, ControlPostProcessing);
			spMarqueeBeta->set_float_limits(-1, 1);
			spMarqueeBeta->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spMarqueeDamping = new GLUI_Spinner(pnMarqueesParameters, "Damping:", GLUI_SPINNER_FLOAT, &myDome.fMarqueeDamping, SP_MARQUEE_DAMPING, ControlPostProcessing);
			spMarqueeDamping->set_float_limits(-1, 1);
			spMarqueeDamping->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spMarqueeTop = new GLUI_Spinner(pnMarqueesParameters, "Top:", GLUI_SPINNER_FLOAT, &myDome.fMarqueeTop, SP_MARQUEE_TOP, ControlPostProcessing);
			spMarqueeTop->set_float_limits(0, 10);
			spMarqueeTop->set_alignment(GLUI_ALIGN_RIGHT);


			pnMarquees->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("Marquee w = %d\n", col_w);
		}

		if ("Bumps & Marquees")
		{
			GLUI_Panel* pnPost = new GLUI_Panel(pnPostProcess, "Bumps & Marquees");
			pnPost->set_w(365);
			GLUI_Spinner* spPostParallelEvery = new GLUI_Spinner(pnPost, "   P. Every:", GLUI_SPINNER_INT, &myDome.nPostParallelEvery, SP_POST_PARALLEL_EVERY, ControlPostProcessing);
			spPostParallelEvery->set_float_limits(1, 30);
			spPostParallelEvery->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spPostParallelFirst = new GLUI_Spinner(pnPost, "     P. First:", GLUI_SPINNER_INT, &myDome.nPostParallelFirst, SP_POST_PARALLEL_FIRST, ControlPostProcessing);
			spPostParallelFirst->set_float_limits(1, 30);
			spPostParallelFirst->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spPostMeridianEvolve = new GLUI_Spinner(pnPost, " P. Evolve:", GLUI_SPINNER_INT, &myDome.nPostParallelEvolve, SP_POST_PARALLEL_EVOLVE, ControlPostProcessing);
			spPostMeridianEvolve->set_float_limits(0, 30);
			spPostMeridianEvolve->set_w(300);
			spPostMeridianEvolve->set_alignment(GLUI_ALIGN_CENTER);

			glParameter->add_column_to_panel(pnPost, true);

			GLUI_Spinner* spPostMeridianEvery = new GLUI_Spinner(pnPost, "  M. Every:", GLUI_SPINNER_INT, &myDome.nPostMeridianEvery, SP_POST_MERIDIAN_EVERY, ControlPostProcessing);
			spPostMeridianEvery->set_float_limits(1, 30);
			spPostMeridianEvery->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spPostMeridianDelay = new GLUI_Spinner(pnPost, "    M. Shift:", GLUI_SPINNER_INT, &myDome.nPostMeridianShift, SP_POST_MERIDIAN_DELAY, ControlPostProcessing);
			spPostMeridianDelay->set_float_limits(0, 30);
			spPostMeridianDelay->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Checkbox* ckPostPerSide = new GLUI_Checkbox(pnPost, "  M. Per Side             ", &myDome.bPostMeridianPerSide, CK_POST_MERIDIAN_PER_SIDE, ControlPostProcessing);
			ckPostPerSide->set_alignment(GLUI_ALIGN_RIGHT);

		}
		GLUI_Panel *pnPostProcessQuit = new GLUI_Panel(glPostProcess, "", GLUI_PANEL_NONE);
		btPostProcessingQuit = new GLUI_Button(pnPostProcessQuit, "Quit", BT_POST_PROCESSING_QUIT, ControlPostProcessing);
		btPostProcessingQuit->set_alignment(GLUI_ALIGN_CENTER);
		glPostProcess->hide();
	}
	if ("glDecor")
	{
		glDecor = GLUI_Master.create_glui("Decor", 0, 650, 80);
		GLUI_Panel* pnDecor = new GLUI_Panel(glDecor, "Decor", GLUI_PANEL_NONE);

		char sStatic[] = "                                                                                                  ";
		if ("Spirale")
		{
			GLUI_Panel *pnSpirale = new GLUI_Panel(pnDecor, "Spirale");
			pnSpirale->set_w(386);

			GLUI_StaticText(pnSpirale, sStatic);
			GLUI_Spinner *spSpirale = new GLUI_Spinner(pnSpirale, " Number:", GLUI_SPINNER_INT, &myDome.nSpirale, SP_SPIRALE, ControlDecor);
			spSpirale->set_int_limits(0, 16);
			spSpirale->set_alignment(GLUI_ALIGN_CENTER);
				
			GLUI_Panel* pnSpiraleParameters = new GLUI_Panel(pnSpirale, "" , GLUI_PANEL_NONE);
			pnSpiraleParameters->set_w(350);

			GLUI_Listbox* lbSpiraleSlope = new GLUI_Listbox(pnSpiraleParameters, "Function:", &myDome.nSpiraleSlope, LB_SPIRALE_SLOPE, ControlDecor);
			lbSortItem(lbSpiraleSlope, nbTestSlopes, (char*)TestSlopeFunctionListe, sizeof(TestSlopeFunctionListe[0]), TestSlopeFunctionListe[0].sNom - (char*)TestSlopeFunctionListe);
			lbSpiraleSlope->set_int_val(SLOPE_STRAIGHT);
			lbSpiraleSlope->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spSpiraleSlopeParametrer = new GLUI_Spinner(pnSpiraleParameters, "Parameter:", GLUI_SPINNER_FLOAT, &myDome.fSpiraleSlopeParameter, SP_SPIRALE_SLOPE_PARAMETER, ControlDecor);
			spSpiraleSlopeParametrer->set_float_limits(0, 10);
			spSpiraleSlopeParametrer->set_alignment(GLUI_ALIGN_RIGHT);
			spSpiraleSlopeParametrer->set_w(200);
			spSpiraleSlopeParametrer->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("spSpiraleSlopeParametrer : %d, %d, %d, %d, %d, %d\n", col_x, col_y, col_w, col_h, col_x_off, col_y_off);

			GLUI_Spinner *spSpiraleDelay = new GLUI_Spinner(pnSpiraleParameters, "      Delay:", GLUI_SPINNER_FLOAT, &myDome.fSpiraleDelay, SP_SPIRALE_DELAY, ControlDecor);
			spSpiraleDelay->set_float_limits((float)0.0, (float)(2.0*_PI));
			spSpiraleDelay->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Listbox *lbSpiraleChiral = new GLUI_Listbox(pnSpiraleParameters, "Chirality", &myDome.nSpiraleChirale, LB_SPIRALE_CHIRALE, ControlDecor);
			lbSpiraleChiral->add_item(CHIRAL_LEFT, "Clockwise");
			lbSpiraleChiral->add_item(CHIRAL_RIGHT, "Anti-Clockwise");
			lbSpiraleChiral->add_item(CHIRAL_BOTH, "Both");
			lbSpiraleChiral->set_int_val(CHIRAL_LEFT);
			lbSpiraleChiral->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spSpiraleRadius = new GLUI_Spinner(pnSpiraleParameters, "    Radius:", GLUI_SPINNER_FLOAT, &myDome.fSpiraleRadius, SP_SPIRALE_RADIUS, ControlDecor);
			spSpiraleRadius->set_float_limits(0, 1);
			spSpiraleRadius->set_alignment(GLUI_ALIGN_RIGHT);
			spSpiraleRadius->set_w(153);
			spSpiraleRadius->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("spSpiraleRadius : %d, %d, %d, %d, %d, %d\n", col_x, col_y, col_w, col_h, col_x_off, col_y_off);

			GLUI_Listbox *lbSpiraleToward = new GLUI_Listbox(pnSpiraleParameters, "    Toward:", &myDome.nSpiraleToward, LB_SPIRALE_TOWARD, ControlDecor);
			lbSpiraleToward->add_item(TOWARD_TOP, "Top");
			lbSpiraleToward->add_item(TOWARD_BOTTOM, "Bottom        ");
			lbSpiraleToward->set_int_val(TOWARD_TOP);
			lbSpiraleToward->set_alignment(GLUI_ALIGN_RIGHT);

			glDecor->add_column_to_panel(pnSpiraleParameters, true);


			GLUI_Spinner* spSpiraleTour = new GLUI_Spinner(pnSpiraleParameters, "Tours:", GLUI_SPINNER_FLOAT, &myDome.fSpiraleTour, SP_SPIRALE_TOUR, ControlDecor);
			spSpiraleTour->set_float_limits(0, 100);
			spSpiraleTour->set_alignment(GLUI_ALIGN_RIGHT);


			GLUI_Spinner *spSpiraleStep = new GLUI_Spinner(pnSpiraleParameters, " Steps:", GLUI_SPINNER_INT, &myDome.nSpiraleStep, SP_SPIRALE_STEP, ControlDecor);
			spSpiraleStep->set_int_limits(1, 30);
			spSpiraleStep->set_alignment(GLUI_ALIGN_RIGHT);


			GLUI_Spinner *spSpiraleSpeed = new GLUI_Spinner(pnSpiraleParameters, "    Speed:", GLUI_SPINNER_FLOAT, &myDome.fSpiraleSpeed, SP_SPIRALE_SPEED, ControlDecor);
			spSpiraleSpeed->set_float_limits(0, 6);
			spSpiraleSpeed->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_EditText *etSpiraleMaterial = new GLUI_EditText(pnSpiraleParameters, "Material:", GLUI_EDITTEXT_TEXT, myDome.sSpiraleMaterial, ET_SPIRALE_MATERIAL, ControlDecor);
			etSpiraleMaterial->set_alignment(GLUI_ALIGN_RIGHT);
			etSpiraleMaterial->set_w(145);

			GLUI_Spinner *spSpiraleDecalage = new GLUI_Spinner(pnSpiraleParameters, "   Offset:", GLUI_SPINNER_FLOAT, &myDome.fSpiraleDecalage, SP_SPIRALE_DECALAGE, ControlDecor);
			spSpiraleDecalage->set_float_limits(0, 1);
			spSpiraleDecalage->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Listbox *lbSpiraleShrink = new GLUI_Listbox(pnSpiraleParameters, "   Shrink:", &myDome.nSpiraleShrink, LB_SPIRALE_SHRINK, ControlDecor);
			lbSpiraleShrink->add_item(SHRINK_NONE, S_NONE);
			lbSpiraleShrink->add_item(SHRINK_SQUARE_ROOT, "Square Root");
			lbSpiraleShrink->add_item(SHRINK_LINEAR, "Linear");
			lbSpiraleShrink->add_item(SHRINK_SQUARE, "Square");
			lbSpiraleShrink->set_int_val(SHRINK_NONE);
			lbSpiraleShrink->set_alignment(GLUI_ALIGN_RIGHT);
		}
		if ("Snake")
		{
			GLUI_Panel *pnSnake = new GLUI_Panel(pnDecor, "Snake");
			pnSnake->set_w(386);
			pnSnake->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("pnSnake : %d, %d, %d, %d, %d, %d\n", col_x, col_y, col_w, col_h, col_x_off, col_y_off);

			GLUI_StaticText(pnSnake, sStatic);
			GLUI_Spinner *spSnake = new GLUI_Spinner(pnSnake, " Number:", GLUI_SPINNER_INT, &myDome.nSnake, SP_SNAKE, ControlDecor);
			spSnake->set_int_limits(0, 30);
			spSnake->set_alignment(GLUI_ALIGN_CENTER);

			GLUI_Panel *pnSnakeParameters = new GLUI_Panel(pnSnake, "", GLUI_PANEL_NONE);

			GLUI_Listbox *lbSnakeSlope = new GLUI_Listbox(pnSnakeParameters, "Slope:", &myDome.nSnakeSlope, LB_SNAKE_SLOPE, ControlDecor);
			lbSortItem(lbSnakeSlope, nbTestSlopes, (char*)TestSlopeFunctionListe, sizeof(TestSlopeFunctionListe[0]), TestSlopeFunctionListe[0].sNom - (char*)TestSlopeFunctionListe);
			lbSnakeSlope->set_alignment(GLUI_ALIGN_RIGHT);

			lbSnakeSlope->set_int_val(SLOPE_STRAIGHT);

			GLUI_Spinner* spSnakeSlopeParametrer = new GLUI_Spinner(pnSnakeParameters, "Parameter:", GLUI_SPINNER_FLOAT, &myDome.fSnakeSlopeParameter, SP_SNAKE_SLOPE_PARAMETER, ControlDecor);
			spSnakeSlopeParametrer->set_float_limits(0, 10);
			spSnakeSlopeParametrer->set_alignment(GLUI_ALIGN_RIGHT);
			spSnakeSlopeParametrer->set_w(200);


			GLUI_Listbox *lbSnakeEvolve = new GLUI_Listbox(pnSnakeParameters, "             Evolve:", &myDome.nSnakeEvolve, LB_SNAKE_EVOLVE, ControlDecor);
			lbSnakeEvolve->add_item(SHRINK_SQUARE_ROOT, "Square Root");
			lbSnakeEvolve->add_item(SHRINK_LINEAR, "Linear");
			lbSnakeEvolve->add_item(SHRINK_SQUARE, "Square");
			lbSnakeEvolve->set_int_val(SHRINK_LINEAR);
			lbSnakeEvolve->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spSnakeOffset = new GLUI_Spinner(pnSnakeParameters, "   Offset:", GLUI_SPINNER_FLOAT, &myDome.fSnakeOffset, SP_SNAKE_OFFSET, ControlDecor);
			spSnakeOffset->set_float_limits(-2, 2);
			spSnakeOffset->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spSnakeRadius = new GLUI_Spinner(pnSnakeParameters, "   Radius:", GLUI_SPINNER_FLOAT, &myDome.fSnakeRadius, SP_SNAKE_RADIUS, ControlDecor);
			spSnakeRadius->set_float_limits(0, 6);
			spSnakeRadius->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Listbox *lbSnakeChirale = new GLUI_Listbox(pnSnakeParameters, "             Chirality", &myDome.nSnakeChirale, LB_SNAKE_CHIRALE, ControlDecor);
			lbSnakeChirale->add_item(CHIRAL_LEFT, "Clockwise");
			lbSnakeChirale->add_item(CHIRAL_RIGHT, "Anti-Clockwise");
			lbSnakeChirale->add_item(CHIRAL_BOTH, "Both");
			lbSnakeChirale->set_int_val(CHIRAL_LEFT);
			lbSnakeChirale->set_alignment(GLUI_ALIGN_RIGHT);

			glDecor->add_column_to_panel(pnSnakeParameters, true);

			GLUI_Spinner* spSnakeTour = new GLUI_Spinner(pnSnakeParameters, "    Tour:", GLUI_SPINNER_FLOAT, &myDome.fSnakeTour, SP_SNAKE_TOUR, ControlDecor);
			spSnakeTour->set_float_limits(0, 10);
			spSnakeTour->set_alignment(GLUI_ALIGN_RIGHT);
			spSnakeTour->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("spSnakeTour : %d, %d, %d, %d, %d, %d\n", col_x, col_y, col_w, col_h, col_x_off, col_y_off);


			GLUI_Spinner *spSnakeStep = new GLUI_Spinner(pnSnakeParameters, "     Steps:", GLUI_SPINNER_INT, &myDome.nSnakeStep, SP_SNAKE_STEP, ControlDecor);
			spSnakeStep->set_int_limits(1, 30);
			spSnakeStep->set_alignment(GLUI_ALIGN_RIGHT);
			spSnakeStep->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("spSnakeStep : %d, %d, %d, %d, %d, %d\n", col_x, col_y, col_w, col_h, col_x_off, col_y_off);

			GLUI_Spinner *spSnakeShift = new GLUI_Spinner(pnSnakeParameters, "       Shift:", GLUI_SPINNER_FLOAT, &myDome.fSnakeShift, SP_SNAKE_SHIFT, ControlDecor);
			spSnakeShift->set_float_limits(-1, 1);
			spSnakeShift->set_alignment(GLUI_ALIGN_RIGHT);
			spSnakeShift->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("spSnakeShift : %d, %d, %d, %d, %d, %d\n", col_x, col_y, col_w, col_h, col_x_off, col_y_off);

			GLUI_Listbox *lbSnakeToward = new GLUI_Listbox(pnSnakeParameters, "     Toward:", &myDome.nSnakeToward, LB_SNAKE_TOWARD, ControlDecor);
			lbSnakeToward->add_item(TOWARD_TOP, "Top");
			lbSnakeToward->add_item(TOWARD_BOTTOM, "Bottom        ");
			lbSnakeToward->set_int_val(TOWARD_TOP);
			lbSnakeToward->set_alignment(GLUI_ALIGN_RIGHT);
			lbSnakeToward->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("lbSnakeToward : %d, %d, %d, %d, %d, %d\n", col_x, col_y, col_w, col_h, col_x_off, col_y_off);

			GLUI_Listbox *lbSnakeShrink = new GLUI_Listbox(pnSnakeParameters, "          Shrink:", &myDome.nSnakeShrink, LB_SNAKE_SHRINK, ControlDecor);
			lbSnakeShrink->add_item(SHRINK_NONE, S_NONE);
			lbSnakeShrink->add_item(SHRINK_SQUARE_ROOT, "Square Root");
			lbSnakeShrink->add_item(SHRINK_LINEAR, "Linear");
			lbSnakeShrink->add_item(SHRINK_SQUARE, "Square");
			lbSnakeShrink->set_int_val(SHRINK_NONE);
			lbSnakeShrink->set_alignment(GLUI_ALIGN_RIGHT);
			lbSnakeShrink->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("lbSnakeShrink : %d, %d, %d, %d, %d, %d\n", col_x, col_y, col_w, col_h, col_x_off, col_y_off);

			GLUI_EditText *etSnakeMaterial = new GLUI_EditText(pnSnakeParameters, "Material:", GLUI_EDITTEXT_TEXT, myDome.sSnakeMaterial, ET_SNAKE_MATERIAL, ControlDecor);
			etSnakeMaterial->set_alignment(GLUI_ALIGN_RIGHT);
			etSnakeMaterial->set_w(145);
			etSnakeMaterial->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("SnakeMaterial : %d, %d, %d, %d, %d, %d\n", col_x, col_y, col_w, col_h, col_x_off, col_y_off);


		}
		glParameter->add_column_to_panel(pnDecor, true);

		if ("Frise")
		{
			GLUI_Panel *pnFrise = new GLUI_Panel(pnDecor, "Frieze");
			pnFrise->set_w(386);
			pnFrise->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("pnFrise : %d, %d, %d, %d, %d, %d\n", col_x, col_y, col_w, col_h, col_x_off, col_y_off);
			GLUI_StaticText(pnFrise, sStatic);

			GLUI_Spinner *spFrise = new GLUI_Spinner(pnFrise, " Number:", GLUI_SPINNER_INT, &myDome.nFrise, SP_FRISE, ControlDecor);
			spFrise->set_int_limits(0, 16);
			spFrise->set_alignment(GLUI_ALIGN_CENTER);

			GLUI_Panel *pnFriseParameters = new GLUI_Panel(pnFrise, "", GLUI_PANEL_NONE);

			GLUI_Listbox *lbFriseFunction = new GLUI_Listbox(pnFriseParameters, "Decor:", &myDome.nFriseFunction, LB_FRISE_FUNCTION, ControlDecor);
			lbSortItem(lbFriseFunction, nbWindows, (char*)WindowFunctionListe, sizeof(WindowFunctionListe[0]), WindowFunctionListe[0].sNom - (char*)WindowFunctionListe);

			lbFriseFunction->set_int_val(WINDOW_NONE);

			GLUI_Spinner* spFriseParameter = new GLUI_Spinner(pnFriseParameters, "Parameter:", GLUI_SPINNER_FLOAT, &myDome.fFriseParameter, SP_FRISE_PARAMETER, ControlDecor);
			spFriseParameter->set_float_limits(0, 10);
			spFriseParameter->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spFrisePower = new GLUI_Spinner(pnFriseParameters, "     Power:", GLUI_SPINNER_FLOAT, &myDome.fFrisePower, SP_FRISE_POWER, ControlDecor);
			spFrisePower->set_float_limits(-2, 10);
			spFrisePower->set_alignment(GLUI_ALIGN_RIGHT);


			GLUI_Spinner *spFriseHeight = new GLUI_Spinner(pnFriseParameters, "     Height:", GLUI_SPINNER_FLOAT, &myDome.fFriseHeight, SP_FRISE_HEIGHT, ControlDecor);
			spFriseHeight->set_float_limits(-10, 10);
			spFriseHeight->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spFriseStart = new GLUI_Spinner(pnFriseParameters, "        Start:", GLUI_SPINNER_FLOAT, &myDome.fFriseDepart, SP_FRISE_START, ControlDecor);
			spFriseStart->set_float_limits(0, 1);
			spFriseStart->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner *spFriseOffset = new GLUI_Spinner(pnFriseParameters, "      Offset:", GLUI_SPINNER_FLOAT, &myDome.fFriseOffset, SP_FRISE_OFFSET, ControlDecor);
			spFriseOffset->set_float_limits(-1, 1);
			spFriseOffset->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spFriseSpacing = new GLUI_Spinner(pnFriseParameters, "   Spacing:", GLUI_SPINNER_FLOAT, &myDome.fFriseSpacing, SP_FRISE_SPACING, ControlDecor);
			spFriseSpacing->set_float_limits(0, 1);
			spFriseSpacing->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Button *btFriseAdjust = new GLUI_Button(pnFriseParameters, "Adjust Spacung", BT_FRISE_ADJUST, ControlDecor);
			btFriseAdjust->set_alignment(GLUI_ALIGN_RIGHT);
			btFriseAdjust->set_w(150);

			GLUI_Checkbox *ckFriseInverse = new GLUI_Checkbox(pnFriseParameters, "       Inverse        ", &myDome.bFriseInverse, CK_FRISE_INVERSE, ControlDecor);
			ckFriseInverse->set_alignment(GLUI_ALIGN_RIGHT);

			glDecor->add_column_to_panel(pnFriseParameters, true);


			GLUI_Spinner *spFriseStep = new GLUI_Spinner(pnFriseParameters, "        Step:", GLUI_SPINNER_INT, &myDome.nFriseStep, SP_FRISE_STEP, ControlDecor);
			spFriseStep->set_int_limits(1, 30);
			spFriseStep->set_alignment(GLUI_ALIGN_RIGHT);
			spFriseStep->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, & col_y_off);
//			printf("step w = %d\n", col_w);
			spFriseStep->set_w(175);
			spFriseStep->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("step w = %d\n", col_w);

			GLUI_Spinner* spFriseRadius = new GLUI_Spinner(pnFriseParameters, "     Radius:", GLUI_SPINNER_FLOAT, &myDome.fFriseRadius, SP_FRISE_RADIUS, ControlDecor);
			spFriseRadius->set_float_limits(0, 6);
			spFriseRadius->set_alignment(GLUI_ALIGN_RIGHT);


			GLUI_EditText *etFriseMaterial = new GLUI_EditText(pnFriseParameters, "   Material:", GLUI_EDITTEXT_TEXT, myDome.sFriseMaterial, ET_FRISE_MATERIAL, ControlDecor);
			etFriseMaterial->set_alignment(GLUI_ALIGN_RIGHT);
			etFriseMaterial->set_w(145);
			etFriseMaterial->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("material w = %d\n", col_w);

			GLUI_Spinner *spFrisePerSide = new GLUI_Spinner(pnFriseParameters, "  Per Side:", GLUI_SPINNER_FLOAT, &myDome.fFrisePerSide, SP_FRISE_PERSIDE, ControlDecor);
			spFrisePerSide->set_float_limits((float)0.0, (float)(10.0));
			spFrisePerSide->set_alignment(GLUI_ALIGN_RIGHT);
			spFrisePerSide->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("Per side w = %d\n", col_w);



			GLUI_Spinner *spFriseShift = new GLUI_Spinner(pnFriseParameters, "        Shift:", GLUI_SPINNER_FLOAT, &myDome.fFriseShift, SP_FRISE_SHIFT, ControlDecor);
			spFriseShift->set_float_limits(0, 1);
			spFriseShift->set_alignment(GLUI_ALIGN_RIGHT);
			spFriseShift->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("shift w = %d\n", col_w);

			GLUI_Spinner* spFriseGap = new GLUI_Spinner(pnFriseParameters, "        Gap:", GLUI_SPINNER_FLOAT, &myDome.fFriseGap, SP_FRISE_GAP, ControlDecor);
			spFriseGap->set_float_limits(0, 1);
			spFriseGap->set_alignment(GLUI_ALIGN_RIGHT);
			spFriseGap->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("gap w = %d\n", col_w);


			GLUI_Listbox *lbFriseShrink = new GLUI_Listbox(pnFriseParameters, "     Shrink:", &myDome.nFriseShrink, LB_FRISE_SHRINK, ControlDecor);
			lbFriseShrink->add_item(SHRINK_NONE, S_NONE);
			lbFriseShrink->add_item(SHRINK_SQUARE_ROOT, "Square Root");
			lbFriseShrink->add_item(SHRINK_LINEAR, "Linear");
			lbFriseShrink->add_item(SHRINK_SQUARE, "Square");
			lbFriseShrink->set_int_val(SHRINK_NONE);
			lbFriseShrink->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("shrink w = %d\n", col_w);

			GLUI_Listbox *lbFriseToward = new GLUI_Listbox(pnFriseParameters, "    Toward:", &myDome.nFriseToward, LB_FRISE_TOWARD, ControlDecor);
			lbFriseToward->add_item(TOWARD_TOP, "Top");
			lbFriseToward->add_item(TOWARD_BOTTOM, "Bottom      ");
			lbFriseToward->set_int_val(TOWARD_TOP);
			lbFriseToward->set_alignment(GLUI_ALIGN_RIGHT);
			lbFriseToward->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, &col_x_off, &col_y_off);
//			printf("toward w = %d\n", col_w);

			GLUI_Checkbox* ckFriseReverse = new GLUI_Checkbox(pnFriseParameters, "       Reverse        ", &myDome.bFriseReverse, CK_FRISE_REVERSE, ControlDecor);			
			ckFriseReverse->set_alignment(GLUI_ALIGN_RIGHT);

		}
		if ("Comb")
		{
			GLUI_Panel* pnComb = new GLUI_Panel(pnDecor, "Comb");
			pnComb->set_w(386);
			GLUI_StaticText(pnComb, sStatic);
			GLUI_Spinner* spComb = new GLUI_Spinner(pnComb, "Teeth:", GLUI_SPINNER_INT, &myDome.nComb, SP_COMB, ControlDecor);
			spComb->set_int_limits(0, 30);
			spComb->set_w(345);

			spComb->set_alignment(GLUI_ALIGN_CENTER);

			GLUI_Panel* pnCombParameters = new GLUI_Panel(pnComb, "", GLUI_PANEL_NONE);

			GLUI_Listbox* lbCombType = new GLUI_Listbox(pnCombParameters, "          Type:", &myDome.nCombType, LB_COMB_TYPE, ControlDecor);
			lbCombType->add_item(COMB_PARALLEL_EDGE, "Parallel Edge");
			lbCombType->add_item(COMB_UPRIGHT_BASE, "UpRight Base");
			lbCombType->set_int_val(COMB_PARALLEL_EDGE);
			lbCombType->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spCombRadius = new GLUI_Spinner(pnCombParameters, "   Radius:", GLUI_SPINNER_FLOAT, &myDome.fCombRadius, SP_COMB_RADIUS, ControlDecor);
			spCombRadius->set_float_limits(0, 6);
			spCombRadius->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Listbox* lbCombShrink = new GLUI_Listbox(pnCombParameters, "          Shrink:", &myDome.nCombShrink, LB_COMB_SHRINK, ControlDecor);
			lbCombShrink->add_item(SHRINK_NONE, S_NONE);
			lbCombShrink->add_item(SHRINK_SQUARE_ROOT, "Square Root");
			lbCombShrink->add_item(SHRINK_LINEAR, "Linear");
			lbCombShrink->add_item(SHRINK_SQUARE, "Square");
			lbCombShrink->set_int_val(SHRINK_NONE);
			lbCombShrink->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Listbox* lbCombToward = new GLUI_Listbox(pnCombParameters, "     Toward:", &myDome.nCombToward, LB_COMB_TOWARD, ControlDecor);
			lbCombToward->add_item(TOWARD_TOP, "Top");
			lbCombToward->add_item(TOWARD_BOTTOM, "Bottom        ");
			lbCombToward->set_int_val(TOWARD_TOP);
			lbCombToward->set_alignment(GLUI_ALIGN_RIGHT);



			GLUI_Spinner* spCombSpacing = new GLUI_Spinner(pnCombParameters, " Spacing:", GLUI_SPINNER_FLOAT, &myDome.fCombSpacing, SP_COMB_SHIFT, ControlDecor);
//			spCombSpacing->set_float_limits(-1, 1);
			spCombSpacing->set_float_limits(0, 10);
			spCombSpacing->set_alignment(GLUI_ALIGN_RIGHT);

			glParameter->add_column_to_panel(pnCombParameters, true);

			GLUI_Spinner* spCombThread = new GLUI_Spinner(pnCombParameters, " Threads:", GLUI_SPINNER_INT, &myDome.nCombThread, SP_COMB_THREAD, ControlDecor);
			spCombThread->set_int_limits(1, 30);
			spCombThread->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spCombStep = new GLUI_Spinner(pnCombParameters, "     Steps:", GLUI_SPINNER_INT, &myDome.nCombStep, SP_COMB_STEP, ControlDecor);
			spCombStep->set_int_limits(1, 30);
			spCombStep->set_alignment(GLUI_ALIGN_RIGHT);




			GLUI_Spinner* spCombGap = new GLUI_Spinner(pnCombParameters, "     Gap:", GLUI_SPINNER_INT, &myDome.nCombThreadGap, SP_COMB_THREAD_GAP, ControlDecor);
			spCombGap->set_int_limits(2, 30);
			spCombGap->set_alignment(GLUI_ALIGN_RIGHT);

			GLUI_Spinner* spCombOffset = new GLUI_Spinner(pnCombParameters, " Offset:", GLUI_SPINNER_FLOAT, &myDome.fCombOffset, SP_COMB_THREAD_GAP_EVOLVE, ControlDecor);
			spCombOffset->set_float_limits(-1, 1);
			spCombOffset->set_alignment(GLUI_ALIGN_RIGHT);




			GLUI_EditText* etCombMaterial = new GLUI_EditText(pnCombParameters, "Material:", GLUI_EDITTEXT_TEXT, myDome.sCombMaterial, ET_COMB_MATERIAL, ControlDecor);
			etCombMaterial->set_alignment(GLUI_ALIGN_RIGHT);
			etCombMaterial->set_w(145);


		}
		
		GLUI_Panel *pnDecorQuit = new GLUI_Panel(glDecor, "", GLUI_PANEL_NONE);
		btDecorQuit = new GLUI_Button(pnDecorQuit, "Quit", BT_DECOR_QUIT, ControlDecor);

		btDecorQuit->set_alignment(GLUI_ALIGN_CENTER);
		glDecor->hide();
	}
	if ("glBase")
	{
		glBase = GLUI_Master.create_glui("Base", 0, 1250, 15);
		GLUI_Panel* pnBase = new GLUI_Panel(glBase, "Base");

		GLUI_Listbox* lbBaseType = new GLUI_Listbox(pnBase, "", &myDome.nBase, LB_BASE, ControlFrame);
		lbBaseType->set_alignment(GLUI_ALIGN_RIGHT);
		lbSortItem(lbBaseType, nbBases, (char*)BaseFunctionListe, sizeof(BaseFunctionListe[0]), BaseFunctionListe[0].sNom - (char*)BaseFunctionListe);
		lbBaseType->set_int_val(BASE_CIRCLE);


		spSide = new GLUI_Spinner(pnBase, "      Sides:", GLUI_SPINNER_INT, &myDome.nBaseSides, SP_BASE_SIDE, ControlFrame);
		spSide->set_int_limits(1, 30);
		spSide->set_speed((float)(0.1 * myDisplay.fSpeedCoefficient));

		spSide->set_alignment(GLUI_ALIGN_RIGHT);
		GLUI_Spinner* spParameter = new GLUI_Spinner(pnBase, "Parameter:", GLUI_SPINNER_FLOAT, &myDome.fBaseParameter, SP_BASE_PARAMETER, ControlFrame);
		spParameter->set_float_limits(0.0, 100.0);
		spParameter->set_alignment(GLUI_ALIGN_RIGHT);
		spParameter->set_speed((float)(0.2 * myDisplay.fSpeedCoefficient));

		GLUI_Spinner* spSmoothing = new GLUI_Spinner(pnBase, "Smoothing:", GLUI_SPINNER_FLOAT, &myDome.fBaseSmoothing, SP_BASE_SMOOTHING, ControlFrame);
		spSmoothing->set_float_limits(-10.0, 10.0);
		spSmoothing->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spExpand = new GLUI_Spinner(pnBase, "    Expand:", GLUI_SPINNER_FLOAT, &myDome.fBaseExpand, SP_BASE_EXPAND, ControlFrame);
		spExpand->set_float_limits(-10.0, 10.0);
		spExpand->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spBaseDamping = new GLUI_Spinner(pnBase, "  Damping:", GLUI_SPINNER_FLOAT, &myDome.fBaseDamping, SP_BASE_DAMPING, ControlFrame);
		spBaseDamping->set_float_limits(-10.0, 10.0);
		spBaseDamping->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spBaseSpiral = new GLUI_Spinner(pnBase, "      Spiral:", GLUI_SPINNER_FLOAT, &myDome.fBaseSpiral, SP_BASE_SPIRAL, ControlFrame);
		spBaseSpiral->set_float_limits(0.0, 1.0);
		spBaseSpiral->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Panel* pnBaseAlong = new GLUI_Panel(pnBase, "Evenly Spaced");

		GLUI_RadioGroup* rgBaseAlong = new GLUI_RadioGroup(pnBaseAlong, &myDome.nBaseEvenly, RG_BASE_SPACED, ControlFrame);
		GLUI_RadioButton* rbBaseAlongParameter = new GLUI_RadioButton(rgBaseAlong, "Along Parameter      ");
		rbBaseAlongParameter->set_alignment(GLUI_ALIGN_LEFT);
		GLUI_RadioButton* rbBaseAlongLength = new GLUI_RadioButton(rgBaseAlong, "Along Length");
		rbBaseAlongLength->set_alignment(GLUI_ALIGN_LEFT);
		GLUI_RadioButton* rbBaseAlongAngle = new GLUI_RadioButton(rgBaseAlong, "Along Angle");
		rbBaseAlongAngle->set_alignment(GLUI_ALIGN_LEFT);


		GLUI_Panel* pnBase2 = new GLUI_Panel(pnBase, "Secondary");

		GLUI_Listbox* lbBaseType2 = new GLUI_Listbox(pnBase2, "", &myDome.nBaseSecondary, LB_BASE_2, ControlFrame);
		lbBaseType2->set_alignment(GLUI_ALIGN_RIGHT);
		lbSortItem(lbBaseType2, nbBases, (char*)BaseFunctionListe, sizeof(BaseFunctionListe[0]), BaseFunctionListe[0].sNom - (char*)BaseFunctionListe);
		lbBaseType2->set_int_val(BASE_CONVEX_POLYGON);


		spBaseRound = new GLUI_Spinner(pnBase2, "   Rounding:", GLUI_SPINNER_FLOAT, &myDome.fBaseRound, SP_BASE_ROUND, ControlFrame);
		spBaseRound->set_float_limits(-(float)_PI, (float)_PI);
		spBaseRound->set_alignment(GLUI_ALIGN_RIGHT);
//		spBaseRound->disable();

		GLUI_Spinner* spBaseBaryCentre = new GLUI_Spinner(pnBase2, "       Weight:", GLUI_SPINNER_FLOAT, &myDome.fBaseBaryCentre, SP_BASE_BARYCENTRE, ControlFrame);
		spBaseBaryCentre->set_float_limits(-10.0, 10.0);
		spBaseBaryCentre->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spBaseMorphing = new GLUI_Spinner(pnBase2, "   Morphing:", GLUI_SPINNER_FLOAT, &myDome.fBaseMorphing, SP_BASE_MORPHING, ControlFrame);
		spBaseMorphing->set_float_limits(0.0, 1.0);
		spBaseMorphing->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spBaseSecondaryRotation = new GLUI_Spinner(pnBase2, "     Rotation:", GLUI_SPINNER_FLOAT, &myDome.fBaseSecondaryRotation, SP_BASE_SECONDARY_ROTATION, ControlFrame);
		spBaseSecondaryRotation->set_float_limits((float) (- 2.0 * _PI), (float)(2.0 * _PI));
		spBaseSecondaryRotation->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spBaseSecondaryExpansion = new GLUI_Spinner(pnBase2, "  Expansion:", GLUI_SPINNER_FLOAT, &myDome.fBaseSecondaryExpansion, SP_BASE_SECONDARY_EXPANSION, ControlFrame);
		spBaseSecondaryExpansion->set_float_limits(0.0, 10.0);
		spBaseSecondaryExpansion->set_alignment(GLUI_ALIGN_RIGHT);

		ckDivision = new GLUI_Checkbox(pnBase2, "Division           ", &myDome.bBaseMorphingDivision, CK_BASE_MORPHING_DIVISION, ControlFrame);
		ckDivision->set_alignment(GLUI_ALIGN_RIGHT);
		if (myDome.nBaseSides % 2)
		{
			ckDivision->disable();
			myDome.bBaseMorphingDivision = 0;
		}
		new GLUI_Button(pnBase, "Default", BT_BASE_DEFAULT, ControlFrame);
		btBaseQuit = new GLUI_Button(glBase, "Quit", BT_BASE_QUIT, ControlFrame);
		glBase->hide();

	}
	if ("glSlope")
	{
		glSlope = GLUI_Master.create_glui("Slope", 0, 1250, 260);

		GLUI_Panel* pnSlope = new GLUI_Panel(glSlope, "Slope");

		lbSlopeType = new GLUI_Listbox(pnSlope, "", &myDome.nSlope, LB_SLOPE, ControlFrame);
		lbSortItem(lbSlopeType, nbSlopes, (char*)SlopeFunctionListe, sizeof(SlopeFunctionListe[0]), SlopeFunctionListe[0].sNom - (char*)SlopeFunctionListe);
		lbSlopeType->set_int_val(SLOPE_HALF_CIRCLE);
		lbSlopeType->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spSlopeParameter_XY = new GLUI_Spinner(pnSlope, "      XY:", GLUI_SPINNER_FLOAT, &myDome.fSlopeParameterXY, SP_SLOPE_PARAMETER_XY, ControlFrame);
		spSlopeParameter_XY->set_int_limits(-5, 10);
		spSlopeParameter_XY->set_speed((float)(0.5 * myDisplay.fSpeedCoefficient));
		spSlopeParameter_XY->set_w(100);
		spSlopeParameter_XY->set_alignment(GLUI_ALIGN_RIGHT);


		GLUI_Spinner* spSlopeParameter_Z = new GLUI_Spinner(pnSlope, "         Z:", GLUI_SPINNER_FLOAT, &myDome.fSlopeParameterZ, SP_SLOPE_PARAMETER, ControlFrame);
		spSlopeParameter_Z->set_int_limits(-5, 10);
		spSlopeParameter_Z->set_speed((float)(0.5 * myDisplay.fSpeedCoefficient));
		spSlopeParameter_Z->set_w(100);
		spSlopeParameter_Z->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spTwist = new GLUI_Spinner(pnSlope, "     Twist:", GLUI_SPINNER_INT, &myDome.nSlopeTwist, SP_SLOPE_TWIST, ControlFrame);
		spTwist->set_int_limits(-3600, 3600);
		spTwist->set_speed((float)(0.5 * myDisplay.fSpeedCoefficient));
		spTwist->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spTwistDelay = new GLUI_Spinner(pnSlope, "     Delay:", GLUI_SPINNER_FLOAT, &myDome.fSlopeDelay, SP_SLOPE_DELAY, ControlFrame);
		spTwistDelay->set_float_limits(-1, 1);
		spTwistDelay->set_speed((float)(0.5 * myDisplay.fSpeedCoefficient));
		spTwistDelay->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Checkbox* ckSlopeInverse = new GLUI_Checkbox(pnSlope, "Inverse           ", &myDome.bSlopeInverse, CK_SLOPE_INVERSE, ControlFrame);
		//		ckSlopeInverse->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Panel* pnAlong = new GLUI_Panel(pnSlope, "Evenly Spaced");

		GLUI_RadioGroup* rgWindowAlong = new GLUI_RadioGroup(pnAlong, &myDome.nSlopeEvenly, RG_SLOPE_SPACED, ControlFrame);
		GLUI_RadioButton* rbAlongParameter = new GLUI_RadioButton(rgWindowAlong, "Along Parameter");
		rbAlongParameter->set_alignment(GLUI_ALIGN_LEFT);
		GLUI_RadioButton* rbAlongLength = new GLUI_RadioButton(rgWindowAlong, "Along Length");
		rbAlongLength->set_alignment(GLUI_ALIGN_LEFT);
		GLUI_RadioButton* rbAlongHeigth = new GLUI_RadioButton(rgWindowAlong, "Along Height");
		rbAlongHeigth->set_alignment(GLUI_ALIGN_LEFT);

		new GLUI_Button(pnSlope, "Default", BT_SLOPE_DEFAULT, ControlFrame);

//		GLUI_Panel* pnTestSlope = new GLUI_Panel(glSlope, "Test");
		lbTestSlopeXY = new GLUI_Listbox(pnSlope, "XY: ", &myDome.nTestSlopeXY, LB_TEST_SLOPE_XY, ControlFrame);
		lbSortItem(lbTestSlopeXY, nbTestSlopes, (char*)TestSlopeFunctionListe, sizeof(TestSlopeFunctionListe[0]), TestSlopeFunctionListe[0].sNom - (char*)TestSlopeFunctionListe);

		lbTestSlopeZ = new GLUI_Listbox(pnSlope, "  Z: ", &myDome.nTestSlopeZ, LB_TEST_SLOPE_Z, ControlFrame);
		lbSortItem(lbTestSlopeZ, nbTestSlopes, (char*)TestSlopeFunctionListe, sizeof(TestSlopeFunctionListe[0]), TestSlopeFunctionListe[0].sNom - (char*)TestSlopeFunctionListe);
		GetFunctions(myDome.nSlope);



		GLUI_Panel* pnCrenel = new GLUI_Panel(glSlope, "Crenel");
		spCrenelPerSide = new GLUI_Spinner(pnCrenel, "Per Side:", GLUI_SPINNER_INT, &myDome.nCrenelPerSide, SP_CRENEL_PERSIDE, ControlFrame);
		spCrenelPerSide->set_int_limits(-1, 10);
		spCrenelPerSide->set_speed((float)(0.5 * myDisplay.fSpeedCoefficient));
		spCrenelPerSide->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spCrenelFront = new GLUI_Spinner(pnCrenel, "    Front:", GLUI_SPINNER_FLOAT, &myDome.fCrenelFront, SP_CRENEL_FRONT, ControlFrame);
		spCrenelFront->set_float_limits(0, 1);
		spCrenelFront->set_speed((float)(2 * myDisplay.fSpeedCoefficient));
		spCrenelFront->set_w(100);
		spCrenelFront->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spCrenelThick = new GLUI_Spinner(pnCrenel, "    Thick:", GLUI_SPINNER_FLOAT, &myDome.fCrenelThick, SP_CRENEL_THICK, ControlFrame);
		spCrenelThick->set_float_limits(0, 1);
		spCrenelThick->set_speed((float)(2 * myDisplay.fSpeedCoefficient));
		spCrenelThick->set_w(100);
		spCrenelThick->set_alignment(GLUI_ALIGN_RIGHT);
		
		spCrenelWide = new GLUI_Spinner(pnCrenel, "    Width:", GLUI_SPINNER_FLOAT, &myDome.fCrenelWide, SP_CRENEL_WIDE, ControlFrame);
		spCrenelWide->set_float_limits(0, 1);
		spCrenelWide->set_speed((float)(2 * myDisplay.fSpeedCoefficient));
		spCrenelWide->set_w(100);
		spCrenelWide->set_alignment(GLUI_ALIGN_RIGHT);


		GLUI_Spinner* spCrenelHeight = new GLUI_Spinner(pnCrenel, "    Height:", GLUI_SPINNER_FLOAT, &myDome.fCrenelHeight, SP_CRENEL_HEIGHT, ControlFrame);
		spCrenelHeight->set_float_limits(-2, 2);
		spCrenelHeight->set_speed((float)(2 * myDisplay.fSpeedCoefficient));
		spCrenelHeight->set_w(100);
		spCrenelHeight->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spCrenelPower = new GLUI_Spinner(pnCrenel, "    Power:", GLUI_SPINNER_FLOAT, &myDome.fCrenelPower, SP_CRENEL_POWER, ControlFrame);
		spCrenelPower->set_float_limits(0.05f, 10.0f);
		spCrenelPower->set_speed((float)(2 * myDisplay.fSpeedCoefficient));
		spCrenelPower->set_w(100);
		spCrenelPower->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spCrenelAngle = new GLUI_Spinner(pnCrenel, "    Angle:", GLUI_SPINNER_FLOAT, &myDome.fCrenelAngle, SP_CRENEL_ANGLE, ControlFrame);
		spCrenelAngle->set_float_limits(-1.0, 1.0f);
		spCrenelAngle->set_speed((float)(2 * myDisplay.fSpeedCoefficient));
		spCrenelAngle->set_w(100);
		spCrenelAngle->set_alignment(GLUI_ALIGN_RIGHT);

		spCrenelDelay = new GLUI_Spinner(pnCrenel, "    Delay:", GLUI_SPINNER_FLOAT, &myDome.fCrenelDelay, SP_CRENEL_DELAY, ControlFrame);
		spCrenelDelay->set_float_limits(0.0f, 10.0f);
		spCrenelDelay->set_speed((float)(2 * myDisplay.fSpeedCoefficient));
		spCrenelDelay->set_w(100);
		spCrenelDelay->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_StaticText* st1 = new GLUI_StaticText(pnCrenel, "                                            ");


		btSlope = new GLUI_Button(glSlope, "Quit", BT_SLOPE_QUIT, ControlFrame);

		glSlope->hide();

	}
	if ("glWindow")
	{
		glWindow = GLUI_Master.create_glui("Window", 0, 1250, 440);

		GLUI_Panel* pnWindow = new GLUI_Panel(glWindow, "Window");

		GLUI_Listbox* lbWindowType = new GLUI_Listbox(pnWindow, "", &myDome.nWindow, LB_WINDOW, ControlFrame);
		lbSortItem(lbWindowType, nbWindows, (char*)WindowFunctionListe, sizeof(WindowFunctionListe[0]), WindowFunctionListe[0].sNom - (char*)WindowFunctionListe);
		lbWindowType->set_int_val(WINDOW_NONE);

		GLUI_Spinner* spWindowParameter = new GLUI_Spinner(pnWindow, "Parameter:", GLUI_SPINNER_FLOAT, &myDome.fWindowParameter, SP_WINDOW_PARAMETER, ControlFrame);
		spWindowParameter->set_float_limits(0, 10);
		spWindowParameter->set_alignment(GLUI_ALIGN_RIGHT);
		spWindowParameter->set_w(150);

		GLUI_Spinner* spWindowPower = new GLUI_Spinner(pnWindow, "     Power:", GLUI_SPINNER_FLOAT, &myDome.fWindowPower, SP_WINDOW_POWER, ControlFrame);
		spWindowPower->set_float_limits(0, 10);
		spWindowPower->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spHeight = new GLUI_Spinner(pnWindow, "     Height:", GLUI_SPINNER_FLOAT, &myDome.fWindowHeight, SP_WINDOW_HEIGHT, ControlFrame);
		spHeight->set_float_limits(-15, 15);
		spHeight->set_alignment(GLUI_ALIGN_RIGHT);
		spHeight->set_speed(0.2f * myDisplay.fSpeedCoefficient);
		spHeight->set_w(150);

		GLUI_Spinner* spDamping = new GLUI_Spinner(pnWindow, "  Damping:", GLUI_SPINNER_FLOAT, &myDome.fWindowDamping, SP_WINDOW_DAMPING, ControlFrame);
		spDamping->set_float_limits(0, (float)0.95);
		spDamping->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spPerSide = new GLUI_Spinner(pnWindow, "  Per Side:", GLUI_SPINNER_FLOAT, &myDome.fWindowPerSide, SP_WINDOW_PERSIDE, ControlFrame);
		spPerSide->set_float_limits(0, 10);
		spPerSide->set_alignment(GLUI_ALIGN_RIGHT);
		GLUI_Spinner* spShift = new GLUI_Spinner(pnWindow, "        Shift:", GLUI_SPINNER_FLOAT, &myDome.fWindowShift, SP_WINDOW_SHIFT, ControlFrame);
		spShift->set_float_limits(0, 1);
		spShift->set_alignment(GLUI_ALIGN_RIGHT);
		spShift->set_speed(2 * myDisplay.fSpeedCoefficient);

		GLUI_Spinner* spWindowDelay = new GLUI_Spinner(pnWindow, "      Delay:", GLUI_SPINNER_FLOAT, &myDome.fWindowDelay, SP_WINDOW_DELAY, ControlFrame);
		spWindowDelay->set_float_limits(0, 3);
		spWindowDelay->set_alignment(GLUI_ALIGN_RIGHT);
		spWindowDelay->set_speed(2 * myDisplay.fSpeedCoefficient);

		GLUI_Spinner* spWindowAttack = new GLUI_Spinner(pnWindow, "      Attack:", GLUI_SPINNER_FLOAT, &myDome.fWindowAttack, SP_WINDOW_ATTACK, ControlFrame);
		spWindowAttack->set_float_limits(0, 3);
		spWindowAttack->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spWindowAwning = new GLUI_Spinner(pnWindow, "    Awning:", GLUI_SPINNER_FLOAT, &myDome.fWindowAwning, SP_WINDOW_AWNING, ControlFrame);
		spWindowAwning->set_float_limits(-2.0f, 2.0f);
		spWindowAwning->set_speed((float)(2 * myDisplay.fSpeedCoefficient));
		spWindowAwning->set_w(100);
		spWindowAwning->set_alignment(GLUI_ALIGN_RIGHT);

		spWindowLower = new GLUI_Spinner(pnWindow, "      Lower:", GLUI_SPINNER_FLOAT, &myDome.fWindowLower, SP_WINDOW_LOWER, ControlFrame);
		spWindowLower->set_float_limits(0.0f, 1.0f); 
		spWindowLower->set_speed((float)(2 * myDisplay.fSpeedCoefficient));
		spWindowLower->set_w(100);
		spWindowLower->set_alignment(GLUI_ALIGN_RIGHT);


		GLUI_Checkbox* ckWindowInverse = new GLUI_Checkbox(pnWindow, "Inverse", &myDome.bWindowInverse, CK_WINDOW_INVERSE, ControlFrame);
		ckWindowInverse->set_alignment(GLUI_ALIGN_CENTER);
		ckWindowBlind = new GLUI_Checkbox(pnWindow, "Blind   ", &myDome.bWindowBlind, CK_WINDOW_BLIND, ControlFrame);
		ckWindowBlind->set_alignment(GLUI_ALIGN_CENTER);

		new GLUI_Button(pnWindow, "Default", BT_WINDOW_DEFAULT, ControlFrame);
		btWindowQuit = new GLUI_Button(glWindow, "Quit", BT_WINDOW_QUIT, ControlFrame);

		glWindow->hide();
	}
	if ("glTop")
	{
		glTop = GLUI_Master.create_glui("Top", 0, 1250, 440);

		GLUI_Panel* pnTop = new GLUI_Panel(glTop, "Top");

		GLUI_Listbox* lbTopType = new GLUI_Listbox(pnTop, "Type:        ", &myDome.nTopType, LB_TOP_TYPE, ControlFrame);
		lbSortItem(lbTopType, nbTops, (char*)TopListe, sizeof(TopListe[0]), TopListe[0].sNom - (char*)TopListe);
		lbTopType->set_int_val(TOP_NONE);
		lbTopType->set_alignment(GLUI_ALIGN_RIGHT);
		lbTopType->set_w(145);

		GLUI_Spinner* spPart = new GLUI_Spinner(pnTop, "Part (%):   ", GLUI_SPINNER_FLOAT, &myDome.fTopPart, SP_DOME_TOP, ControlFrame);
		spPart->set_float_limits(0, 1);
		spPart->set_speed((float)(0.5 * myDisplay.fSpeedCoefficient));
		spPart->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner* spParam = new GLUI_Spinner(pnTop, "Parameter:", GLUI_SPINNER_FLOAT, &myDome.fTopParameter, SP_TOP_PARAMETER, ControlFrame);
		spParam->set_float_limits(0, 3);
		spParam->set_speed((float)(0.5 * myDisplay.fSpeedCoefficient));
		spParam->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Button* btTopQuit = new GLUI_Button(glTop, "Quit", BT_TOP_QUIT, ControlFrame);

		glTop->hide();

	}
	glParameter = GLUI_Master.create_glui_subwindow(MainWindow, GLUI_SUBWINDOW_RIGHT);
	Largeur = 370;
	GLUI_Panel *pnDome = new GLUI_Panel(glParameter, "Dome");
	if ("Dome")
	{
		GLUI_Panel *pnInvisibleDome2 = new GLUI_Panel(pnDome, "", GLUI_PANEL_NONE);

		GLUI_EditText *etDomeName = new GLUI_EditText(pnInvisibleDome2, "Name:", GLUI_EDITTEXT_TEXT, myDome.sDomeName, ET_DOME_NAME, ControlFrame);
		etDomeName->set_alignment(GLUI_ALIGN_RIGHT);
		etDomeName->set_w(145);



		//	GLUI_Panel *pnInvisibleDome1 = new GLUI_Panel(pnInvisibleDome2, "", GLUI_PANEL_NONE);

		GLUI_Spinner *spScaleZ = new GLUI_Spinner(pnInvisibleDome2, "Scale Z:", GLUI_SPINNER_FLOAT, &myDome.fDomeScaleZ, SP_DOME_SCALE_Z, ControlFrame);
		spScaleZ->set_float_limits(-100, 100);
		spScaleZ->set_speed((float)(0.1*myDisplay.fSpeedCoefficient));
		spScaleZ->set_alignment(GLUI_ALIGN_RIGHT);
		GLUI_Spinner *spEnd = new GLUI_Spinner(pnInvisibleDome2, "     End:", GLUI_SPINNER_INT, &myDome.nDomeEnd, SP_DOME_END, ControlFrame);
		spEnd->set_int_limits(0, 360);
		spEnd->set_speed((float)(0.5*myDisplay.fSpeedCoefficient));
		spEnd->set_alignment(GLUI_ALIGN_RIGHT);
		spEnd->set_w(145);
		//
		GLUI_Spinner *spSocle = new GLUI_Spinner(pnInvisibleDome2, " Socle:", GLUI_SPINNER_FLOAT, &myDome.fDomeSocle, SP_DOME_SOCLE, ControlFrame);
		spSocle->set_float_limits(0, (float)0.95);
		spSocle->set_speed((float)(5*myDisplay.fSpeedCoefficient));
		spSocle->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Button *btScaleAdjust = new GLUI_Button(pnInvisibleDome2, "Adjust Scale", BT_ADJUST_SCALE, ControlSkeleton);
		btScaleAdjust->set_w(150);
		btScaleAdjust->set_alignment(GLUI_ALIGN_CENTER);

		GLUI_Panel* pnTop = new GLUI_Panel(pnInvisibleDome2, "Top");

		tbTop = new GLUI_TextBox(pnTop, false, -1, ControlFrame);
		tbTop->disable();
		tbTop->set_text(TopListe[myDome.nTopType].sNom);
		tbTop->set_h(22);
		tbTop->set_w(150);
		tbTop->set_alignment(GLUI_ALIGN_CENTER);

//		GLUI_Panel* pnInvisibleBase = new GLUI_Panel(pnBase1, "", GLUI_PANEL_NONE);

		GLUI_Panel *pnTopInvisible = new GLUI_Panel(pnTop, "", GLUI_PANEL_NONE);
		btTop = new GLUI_Button(pnTopInvisible, "Top Parameters", BT_TOP, ControlSkeleton);
		btTop->set_w(150);
		btTop->set_alignment(GLUI_ALIGN_CENTER);


/*		GLUI_Spinner* spTop = new GLUI_Spinner(pnInvisibleDome2, " Top:", GLUI_SPINNER_FLOAT, &myDome.fTopPart, SP_DOME_TOP, ControlFrame);
		spTop->set_float_limits(0, (float)0.95);
		spTop->set_speed((float)(5*myDisplay.fSpeedCoefficient));
		spTop->set_alignment(GLUI_ALIGN_RIGHT);*/

		GLUI_EditText *etDomeMaterial = new GLUI_EditText(pnInvisibleDome2, "Material:", GLUI_EDITTEXT_TEXT, myDome.sDomeMaterial, ET_DOME_MATERIAL, ControlFrame);

		etDomeMaterial->set_alignment(GLUI_ALIGN_RIGHT);
		etDomeMaterial->set_w(145);



		//
		int nLargeurBouton = 45;
		GLUI_Panel *pnInvisibleDome = new GLUI_Panel(pnDome, "", GLUI_PANEL_NONE);
		GLUI_Button *btDefaultParameters = new GLUI_Button(pnInvisibleDome, "Default", BT_DOME_DEFAULT, ControlFrame);
		btDefaultParameters->set_w(59);
		GLUI_Button *btLoadDome = new GLUI_Button(pnInvisibleDome, "ReLoad", BT_DOME_LOAD, ControlFrame);
		btLoadDome->set_w(59);
 
		glParameter->add_column_to_panel(pnInvisibleDome, false);
		
		GLUI_Button *btSaveDome = new GLUI_Button(pnInvisibleDome, "Save", BT_DOME_SAVE, ControlFrame);
		btSaveDome->set_w(nLargeurBouton);
		GLUI_Button *btGenerateDome = new GLUI_Button(pnInvisibleDome, "INC", BT_DOME_INC, ControlFrame);
		btGenerateDome->set_w(nLargeurBouton);
		glParameter->add_column_to_panel(pnInvisibleDome, false);

		GLUI_Button *btGenerateStl = new GLUI_Button(pnInvisibleDome, "STL", BT_DOME_STL, ControlFrame);
		btGenerateStl->set_w(nLargeurBouton);
		GLUI_Button *btGenerateObj = new GLUI_Button(pnInvisibleDome, "OBJ", BT_DOME_OBJ, ControlFrame);
		btGenerateObj->set_w(nLargeurBouton);

//		GLUI_Panel* pnInvisibleDome1 = new GLUI_Panel(pnDome, "", GLUI_PANEL_NONE);

	}

	glParameter->add_column_to_panel(pnDome, true);

	if ("Frame")
	{
		GLUI_Panel* pnBase1 = new GLUI_Panel(pnDome, "Base");

		tbBase = new GLUI_TextBox(pnBase1, false, -1, ControlFrame);
		tbBase->disable();
		tbBase->set_text(BaseFunctionListe[myDome.nBase].sNom);
		tbBase->set_h(22);
		tbBase->set_w(150);
		tbBase->set_alignment(GLUI_ALIGN_CENTER);

		GLUI_Panel* pnInvisibleBase = new GLUI_Panel(pnBase1, "", GLUI_PANEL_NONE);

		btBase = new GLUI_Button(pnInvisibleBase, "Base Parameters", BT_BASE, ControlSkeleton);
		btBase->set_w(150);
		btBase->set_alignment(GLUI_ALIGN_CENTER);

		GLUI_Panel* pnSlope1 = new GLUI_Panel(pnDome, "Slope");

		tbSlope = new GLUI_TextBox(pnSlope1, false, -1, ControlFrame);
		tbSlope->disable();
		tbSlope->set_text(SlopeFunctionListe[myDome.nSlope].sNom);
		tbSlope->set_h(22);
		tbSlope->set_w(150);
		tbSlope->set_alignment(GLUI_ALIGN_CENTER);

		GLUI_Panel* pnInvisibleSlope = new GLUI_Panel(pnSlope1, "", GLUI_PANEL_NONE);

		btSlope = new GLUI_Button(pnInvisibleSlope, "Slope Parameters", BT_SLOPE, ControlSkeleton);
		btSlope->set_w(150);
		btSlope->set_alignment(GLUI_ALIGN_CENTER);


		GLUI_Panel* pnWindow1 = new GLUI_Panel(pnDome, "Window");

		tbWindow = new GLUI_TextBox(pnWindow1, false, -1, ControlFrame);
		tbWindow->disable();
		tbWindow->set_text(WindowFunctionListe[myDome.nWindow].sNom);
		tbWindow->set_h(22);
		tbWindow->set_w(150);
		tbWindow->set_alignment(GLUI_ALIGN_CENTER);

		GLUI_Panel* pnInvisibleWindow = new GLUI_Panel(pnWindow1, "", GLUI_PANEL_NONE);

		btWindow = new GLUI_Button(pnInvisibleWindow, "Window Parameters", BT_WINDOW, ControlSkeleton);
		btWindow->set_w(150);
		btWindow->set_alignment(GLUI_ALIGN_CENTER);

		GLUI_Panel* pnInvisibleDome = new GLUI_Panel(pnDome, "", GLUI_PANEL_NONE);

		GLUI_Spinner* spInFolding = new GLUI_Spinner(pnInvisibleDome, " InFolding:", GLUI_SPINNER_FLOAT, &myDome.fDomeInFolding, SP_DOME_INFOLDING, ControlFrame);
		spInFolding->set_float_limits(0, (float)0.95);
		spInFolding->set_speed((float)(5 * myDisplay.fSpeedCoefficient));
		spInFolding->set_alignment(GLUI_ALIGN_RIGHT);

	}
	char sStatic[] = "                                                                                   ";
	
	new GLUI_Panel(glParameter, "", GLUI_PANEL_NONE);

	int nDomeFile;
	lbDomeFile = new GLUI_Listbox(glParameter, "    File: ", &nDomeFile, LB_DOME_FILE, ControlFrame);
	lbDomeFile->set_w(370);
	fill_lbfile(lbDomeFile, (char*)".DOM", NULL);
	_getcwd(Buffer, 500);
//	printf("*** %s ***", Buffer);

	tbDirectory = new GLUI_TextBox(glParameter, false, -1, ControlFrame);
	tbDirectory->disable();
	tbDirectory->set_text(Buffer);
	tbDirectory->set_h(22);
	tbDirectory->set_w(370);
	tbDirectory->set_alignment(GLUI_ALIGN_CENTER);

	new GLUI_Panel(glParameter, "", GLUI_PANEL_NONE);


	if ("Rendering")
	  {
		  GLUI_Rollout *roWireFace = new GLUI_Rollout(glParameter, "Rendering", false);
		  roWireFace->set_w(Largeur);
//		  new GLUI_StaticText(roWireFace, sStatic);
		  GLUI_Panel *pnFace = new GLUI_Panel(roWireFace, "Face");
		  new GLUI_Checkbox(pnFace, "Generate", &myDome.bFaceGenerate, CK_FACE, ControlFrame);
		  ckFaceSmooth = new GLUI_Checkbox(pnFace, "Smooth", &myDome.bFaceSmooth, CK_SMOOTH, ControlFrame);

		  GLUI_Spinner *spMeridianEdgeSize = new GLUI_Spinner(pnFace, "Meridian Edge:", GLUI_SPINNER_INT, &myDome.nMeridianEdgeSize, SP_MERIDIAN_EDGE_SIZE, ControlFrame);
		  spMeridianEdgeSize->set_int_limits(0, 15);
		  spMeridianEdgeSize->set_alignment(GLUI_ALIGN_RIGHT);
		  spMeridianEdgeSize->set_speed((float)(0.5*myDisplay.fSpeedCoefficient));

		  GLUI_Spinner *spParallelEdgeSize = new GLUI_Spinner(pnFace, "  Parallel Edge:", GLUI_SPINNER_INT, &myDome.nParallelEdgeSize, SP_PARALLEL_EDGE_SIZE, ControlFrame);
		  spParallelEdgeSize->set_int_limits(0, 15);
		  spParallelEdgeSize->set_alignment(GLUI_ALIGN_RIGHT);
		  spParallelEdgeSize->set_speed((float)(0.5*myDisplay.fSpeedCoefficient));


		  GLUI_Panel *pnMaterial = new GLUI_Panel(pnFace, "Material");
		  GLUI_EditText *etFaceMaterial = new GLUI_EditText(pnMaterial, "Center:", GLUI_EDITTEXT_TEXT, myDome.sFaceCenterMaterial, ET_FACE_CENTER_MATERIAL, ControlFrame);
		  GLUI_EditText *etFaceEdgeMaterial = new GLUI_EditText(pnMaterial, "Edge:", GLUI_EDITTEXT_TEXT, myDome.sFaceEdgeMaterial, ET_FACE_EDGE_MATERIAL, ControlFrame);

		  glParameter->add_column_to_panel(roWireFace, true);
		  GLUI_Panel *pnWire = new GLUI_Panel(roWireFace, "Wire");
		  new GLUI_Checkbox(pnWire, "Generate", &myDome.bWireGenerate, CK_WIRE, ControlFrame);
		  GLUI_Spinner *spShrinking = new GLUI_Spinner(pnWire, "Shrinking:", GLUI_SPINNER_FLOAT, &myDome.fWireShrinking, SP_WIRE_SHRINKING, ControlFrame);
		  spShrinking->set_float_limits(-10, 10);
		  spShrinking->set_alignment(GLUI_ALIGN_RIGHT);
		  GLUI_Listbox *lbShrinkingType = new GLUI_Listbox(pnWire, "Type:", &myDome.nWireShrinkingType, LB_WIRE_SHRINKING_TYPE, ControlFrame);
		  lbShrinkingType->add_item(SHRINK_NONE, "None");
		  lbShrinkingType->add_item(SHRINK_SQUARE_ROOT, "Square Root");
		  lbShrinkingType->add_item(SHRINK_LINEAR, "Linear");
		  lbShrinkingType->add_item(SHRINK_SQUARE, "Square");
		  lbShrinkingType->set_int_val(SHRINK_NONE);
		  lbShrinkingType->set_alignment(GLUI_ALIGN_RIGHT);
		  GLUI_Spinner *spSpeed = new GLUI_Spinner(pnWire, "Speed:", GLUI_SPINNER_FLOAT, &myDome.fWireShrinkingSpeed, SP_WIRE_SHRINKING_SPEED, ControlFrame);
		  spSpeed->set_float_limits(-10, 10);
		  spSpeed->set_alignment(GLUI_ALIGN_RIGHT);

		  GLUI_Panel *pnLid = new GLUI_Panel(roWireFace, "Lid");
		  new GLUI_Checkbox(pnLid, "Closed", &myDome.bLid, CK_LID, ControlFrame);
		  GLUI_EditText *etLidMaterial = new GLUI_EditText(pnLid, "Material:", GLUI_EDITTEXT_TEXT, myDome.sLidMaterial, ET_LID_MATERIAL, ControlFrame);


	  }
	if ("Struts")
	  {
		  GLUI_Rollout *roMerPar = new GLUI_Rollout(glParameter, "Struts", false);
		  roMerPar->set_w(Largeur);
		  GLUI_Panel *pnMeridian = new GLUI_Panel(roMerPar, "Meridian");
		  GLUI_Spinner *spMainMeridian = new GLUI_Spinner(pnMeridian, "Main:", GLUI_SPINNER_INT, &myDome.nMeridianMain, SP_MERIDIAN_MAIN, ControlSkeleton);
		  spMainMeridian->set_int_limits(1, 30);
		  spMainMeridian->set_alignment(GLUI_ALIGN_RIGHT);
		  spMainMeridian->set_speed((float)(0.5*myDisplay.fSpeedCoefficient));

		  GLUI_Spinner *spMainMeridianRadius = new GLUI_Spinner(pnMeridian, "Radius:", GLUI_SPINNER_FLOAT, &myDome.fMeridianMainRadius, SP_MERIDIAN_MAIN_RADIUS, ControlSkeleton
		  );
		  spMainMeridianRadius->set_float_limits(0, 1);
		  spMainMeridianRadius->set_alignment(GLUI_ALIGN_RIGHT);
		  GLUI_EditText *etMainMeridianMaterial = new GLUI_EditText(pnMeridian, "Material:", GLUI_EDITTEXT_TEXT, myDome.sMeridianMainMaterial, ET_MERIDIAN_MAIN_MATERIAL, ControlSkeleton);
		  etMainMeridianMaterial->set_alignment(GLUI_ALIGN_RIGHT);
		  etMainMeridianMaterial->set_w(145);

		  new GLUI_Separator(pnMeridian);
		  GLUI_Spinner *spSecMeridian = new GLUI_Spinner(pnMeridian, "Minor:", GLUI_SPINNER_INT, &myDome.nMeridianMinor, SP_MERIDIAN_MINOR, ControlSkeleton);
		  spSecMeridian->set_int_limits(0, 50);
		  spSecMeridian->set_alignment(GLUI_ALIGN_RIGHT);
		  spSecMeridian->set_speed((float)(0.5*myDisplay.fSpeedCoefficient));
		  GLUI_Spinner *spSecMeridianRadius = new GLUI_Spinner(pnMeridian, "Radius:", GLUI_SPINNER_FLOAT, &myDome.fMeridianMinorRadius, SP_MERIDIAN_MINOR_RADIUS, ControlSkeleton);
		  spSecMeridianRadius->set_float_limits(0, 1);
		  spSecMeridianRadius->set_alignment(GLUI_ALIGN_RIGHT);
		  GLUI_EditText *etMinorMeridianMaterial = new GLUI_EditText(pnMeridian, "Material:", GLUI_EDITTEXT_TEXT, myDome.sMeridianMinorMaterial, ET_MERIDIAN_MINOR_MATERIAL, ControlSkeleton);
		  etMinorMeridianMaterial->set_alignment(GLUI_ALIGN_RIGHT);
		  etMinorMeridianMaterial->set_w(145);

		  glParameter->add_column_to_panel(roMerPar, true);

		  GLUI_Panel *pnParallel = new GLUI_Panel(roMerPar, "Parallel");
		  GLUI_Spinner *spMainParallel = new GLUI_Spinner(pnParallel, "Main:", GLUI_SPINNER_INT, &myDome.nParallelMain, SP_PARALLEL_MAIN, ControlSkeleton);
		  spMainParallel->set_int_limits(2, 30);
		  spMainParallel->set_alignment(GLUI_ALIGN_RIGHT);
		  spMainParallel->set_speed((float)(0.5*myDisplay.fSpeedCoefficient));
		  GLUI_Spinner *spMainParallelRadius = new GLUI_Spinner(pnParallel, "Radius:", GLUI_SPINNER_FLOAT, &myDome.fParallelMainRadius, SP_PARALLEL_MAIN_RADIUS, ControlSkeleton);
		  spMainParallelRadius->set_float_limits(0, 1);
		  spMainParallelRadius->set_alignment(GLUI_ALIGN_RIGHT);
		  GLUI_EditText *etMainParallelMaterial = new GLUI_EditText(pnParallel, "Material:", GLUI_EDITTEXT_TEXT, myDome.sParallelMainMaterial, ET_PARALLEL_MAIN_MATERIAL, ControlSkeleton);
		  etMainParallelMaterial->set_alignment(GLUI_ALIGN_RIGHT);
		  etMainParallelMaterial->set_w(145);

		  new GLUI_Separator(pnParallel);
		  GLUI_Spinner *spSecParallel = new GLUI_Spinner(pnParallel, "Minor:", GLUI_SPINNER_INT, &myDome.nParallelMinor, SP_PARALLEL_MINOR, ControlSkeleton);
		  spSecParallel->set_int_limits(0, 50);
		  spSecParallel->set_alignment(GLUI_ALIGN_RIGHT);
		  spSecParallel->set_speed((float)(0.5*myDisplay.fSpeedCoefficient));
		  GLUI_Spinner *spSecParallelRadius = new GLUI_Spinner(pnParallel, "Radius:", GLUI_SPINNER_FLOAT, &myDome.fParallelMinorRadius, SP_PARALLEL_MINOR_RADIUS, ControlSkeleton);
		  spSecParallelRadius->set_float_limits(0, 1);
		  spSecParallelRadius->set_alignment(GLUI_ALIGN_RIGHT);
		  GLUI_EditText *etMinorParallelMaterial = new GLUI_EditText(pnParallel, "Material:", GLUI_EDITTEXT_TEXT, myDome.sParallelMinorMaterial, ET_PARALLEL_MINOR_MATERIAL, ControlSkeleton);
		  etMinorParallelMaterial->set_alignment(GLUI_ALIGN_RIGHT);
		  etMinorParallelMaterial->set_w(145);
	  }
	if ("Diagonal")
	{
		GLUI_Rollout *roDiagonal = new GLUI_Rollout(glParameter, "Diagonals", false);
		roDiagonal->set_w(Largeur);
		GLUI_Panel *pnDiagonalUR = new GLUI_Panel(roDiagonal, "Down Left to Up Right");
		GLUI_Spinner *spDiagonalURRadius = new GLUI_Spinner(pnDiagonalUR, "Radius:", GLUI_SPINNER_FLOAT, &myDome.fDiagonalURRadius);
		spDiagonalURRadius->set_float_limits(0, 1);
		spDiagonalURRadius->set_alignment(GLUI_ALIGN_RIGHT);
		GLUI_EditText *etDiagonalURMaterial = new GLUI_EditText(pnDiagonalUR, "Material:", GLUI_EDITTEXT_TEXT, myDome.sDiagonalURMaterial, ET_DIAGONALE_UR_MATERIAL, ControlSkeleton);
		etDiagonalURMaterial->set_alignment(GLUI_ALIGN_RIGHT);
		etDiagonalURMaterial->set_w(145);

		glParameter->add_column_to_panel(roDiagonal, true);
		GLUI_Panel *pnDiagonalDR = new GLUI_Panel(roDiagonal, "Up Left to Down Right");
		GLUI_Spinner *spDiagonalDRRadius = new GLUI_Spinner(pnDiagonalDR, "Radius:", GLUI_SPINNER_FLOAT, &myDome.fDiagonalDRRadius);
		spDiagonalDRRadius->set_float_limits(0, 1);
		spDiagonalDRRadius->set_alignment(GLUI_ALIGN_RIGHT);
		GLUI_EditText *etDiagonalDRMaterial = new GLUI_EditText(pnDiagonalDR, "Material:", GLUI_EDITTEXT_TEXT, myDome.sDiagonalDRMaterial, ET_DIAGONALE_DR_MATERIAL, ControlSkeleton);
		etDiagonalDRMaterial->set_alignment(GLUI_ALIGN_RIGHT);
		etDiagonalDRMaterial->set_w(145);
	}
	if ("Sphere")
	{
		GLUI_Rollout *roSphere = new GLUI_Rollout(glParameter, "Spheres & Spikes", false);

		roSphere->set_w(Largeur);

		GLUI_Panel *pnEspaceBouton = new GLUI_Panel(roSphere, "", GLUI_PANEL_NONE);
		GLUI_Button *btSphere = new GLUI_Button(pnEspaceBouton, "Adjust Spheres", BT_SPHERE_ADJUST, ControlSkeleton);
		btSphere->set_w(120);
		btSphere->set_alignment(GLUI_ALIGN_LEFT);
		glParameter->add_column_to_panel(pnEspaceBouton, true);
		GLUI_Button *btSpike = new GLUI_Button(pnEspaceBouton, "Defaukt Spikes", BT_SPIKE_DEFAULT, ControlSkeleton);
		btSpike->set_w(120);
		btSpike->set_alignment(GLUI_ALIGN_CENTER);

		GLUI_Spinner* spSpike = new GLUI_Spinner(roSphere, "Spike's Angle:", GLUI_SPINNER_FLOAT, &myDome.fSpikeAngle, SP_SPIKE_ANGLE, ControlSkeleton);
		spSpike->set_float_limits(-2, 2);
		spSpike->set_alignment(GLUI_ALIGN_CENTER);

		GLUI_Panel *pnInvisibleSphere = new GLUI_Panel(roSphere, "", GLUI_PANEL_NONE);

		GLUI_Panel *npMP = new GLUI_Panel(pnInvisibleSphere, "Spheres MP");
		GLUI_Spinner *spSphereMPRadius = new GLUI_Spinner(npMP, "Radius:", GLUI_SPINNER_FLOAT, &myDome.fSphereMPRadius, SP_SPHERE_MP_RADIUS, ControlSkeleton);
		spSphereMPRadius->set_float_limits(0, 1);
		spSphereMPRadius->set_alignment(GLUI_ALIGN_RIGHT);
		GLUI_EditText *etSphereMPMaterial = new GLUI_EditText(npMP, "Material:", GLUI_EDITTEXT_TEXT, &myDome.sSphereMPMaterial[0], ET_SPHERE_MP_MATERIAL, ControlSkeleton);
		etSphereMPMaterial->set_alignment(GLUI_ALIGN_RIGHT);
		etSphereMPMaterial->set_w(145);
		GLUI_Spinner *spSpikeMP = new GLUI_Spinner(npMP, "Spikes:", GLUI_SPINNER_FLOAT, &myDome.fSpikeMP, SP_SPIKE_MP, ControlSkeleton);
		spSpikeMP->set_float_limits(0, 5);
		spSpikeMP->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Panel *npMp = new GLUI_Panel(pnInvisibleSphere, "Spheres Mp");
		GLUI_Spinner *spSphereMpRadius = new GLUI_Spinner(npMp, "Radius:", GLUI_SPINNER_FLOAT, &myDome.fSphereMpRadius, SP_SPHERE_Mp_RADIUS, ControlSkeleton);
		spSphereMpRadius->set_float_limits(0, 1);
		spSphereMpRadius->set_alignment(GLUI_ALIGN_RIGHT);
		GLUI_EditText *etSphereMpMaterial = new GLUI_EditText(npMp, "Material:", GLUI_EDITTEXT_TEXT, myDome.sSphereMpMaterial, ET_SPHERE_Mp_MATERIAL, ControlSkeleton);
		etSphereMpMaterial->set_alignment(GLUI_ALIGN_RIGHT);
		etSphereMpMaterial->set_w(145);
		GLUI_Spinner *spSpikeMp = new GLUI_Spinner(npMp, "Spikes:", GLUI_SPINNER_FLOAT, &myDome.fSpikeMp, SP_SPIKE_Mp, ControlSkeleton);
		spSpikeMp->set_float_limits(0, 5);
		spSpikeMp->set_alignment(GLUI_ALIGN_RIGHT);

		glParameter->add_column_to_panel(pnInvisibleSphere, true);

		GLUI_Panel *npmP = new GLUI_Panel(pnInvisibleSphere, "Spheres mP");
		GLUI_Spinner *spSpheremPRadius = new GLUI_Spinner(npmP, "Radius:", GLUI_SPINNER_FLOAT, &myDome.fSpheremPRadius, SP_SPHERE_mP_RADIUS, ControlSkeleton);
		spSpheremPRadius->set_float_limits(0, 1);
		spSpheremPRadius->set_alignment(GLUI_ALIGN_RIGHT);
		GLUI_EditText *etSpheremPMaterial = new GLUI_EditText(npmP, "Material:", GLUI_EDITTEXT_TEXT, myDome.sSpheremPMaterial, ET_SPHERE_mP_MATERIAL, ControlSkeleton);
		etSpheremPMaterial->set_alignment(GLUI_ALIGN_RIGHT);
		etSpheremPMaterial->set_w(145);
		GLUI_Spinner *spSpikemP = new GLUI_Spinner(npmP, "Spikes:", GLUI_SPINNER_FLOAT, &myDome.fSpikemP, SP_SPIKE_mP, ControlSkeleton);
		spSpikemP->set_float_limits(0, 5);
		spSpikemP->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Panel *npmp = new GLUI_Panel(pnInvisibleSphere, "Spheres mp");
		GLUI_Spinner *spSpherempRadius = new GLUI_Spinner(npmp, "Radius:", GLUI_SPINNER_FLOAT, &myDome.fSpherempRadius, SP_SPHERE_mp_RADIUS, ControlSkeleton);
		spSpherempRadius->set_float_limits(0, 1);
		spSpherempRadius->set_alignment(GLUI_ALIGN_RIGHT);
		GLUI_EditText *etSpherempMaterial = new GLUI_EditText(npmp, "Material:", GLUI_EDITTEXT_TEXT, myDome.sSpherempMaterial, ET_SPHERE_mp_MATERIAL, ControlSkeleton);
		etSpherempMaterial->set_alignment(GLUI_ALIGN_RIGHT);
		etSpherempMaterial->set_w(145);
		GLUI_Spinner *spSpikemp = new GLUI_Spinner(npmp, "Spikes:", GLUI_SPINNER_FLOAT, &myDome.fSpikemp, SP_SPIKE_mp, ControlSkeleton);
		spSpikemp->set_float_limits(0, 5);
		spSpikemp->set_alignment(GLUI_ALIGN_RIGHT);
	}
	if ("Noise" && false)

	{
		GLUI_Rollout *roNoise = new GLUI_Rollout(glParameter, "Noise", false);
		roNoise->set_w(Largeur);
		GLUI_Spinner *spNoiseOctave = new GLUI_Spinner(roNoise, "        Octave:", GLUI_SPINNER_INT, &myDome.nNoiseOctave, SP_NOISE_OCTAVE, ControlNoise);
		spNoiseOctave->set_int_limits(1, 16);
		spNoiseOctave->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner *spNoiseDensity = new GLUI_Spinner(roNoise, "        Density:", GLUI_SPINNER_FLOAT, &myDome.fNoiseDensity, SP_NOISE_DENSITY, ControlNoise);
		spNoiseDensity->set_float_limits(1, 16);
		spNoiseDensity->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner *spNoiseSeed = new GLUI_Spinner(roNoise, "           Seed:", GLUI_SPINNER_INT, &myDome.nNoiseSeed, SP_NOISE_SEED, ControlNoise);
		spNoiseSeed->set_int_limits(0, 32000);
		spNoiseSeed->set_alignment(GLUI_ALIGN_RIGHT);
		spNoiseSeed->set_speed((float)(0.1*myDisplay.fSpeedCoefficient));
		GLUI_Spinner *spNoiseFrequencyU = new GLUI_Spinner(roNoise, "Frequency U:", GLUI_SPINNER_FLOAT, &myDome.fNoiseFrequencyU, SP_NOISE_FREQUENCY_U, ControlNoise);
		spNoiseFrequencyU->set_float_limits(0, 16);
		spNoiseFrequencyU->set_alignment(GLUI_ALIGN_RIGHT);
		spNoiseFrequencyU->set_speed((float)(1 * myDisplay.fSpeedCoefficient));
		GLUI_Spinner *spNoiseFrequencyV = new GLUI_Spinner(roNoise, "Frequency V:", GLUI_SPINNER_FLOAT, &myDome.fNoiseFrequencyV, SP_NOISE_FREQUENCY_V, ControlNoise);
		spNoiseFrequencyV->set_float_limits(0, 16);
		spNoiseFrequencyV->set_alignment(GLUI_ALIGN_RIGHT);
		spNoiseFrequencyV->set_speed((float)(1 * myDisplay.fSpeedCoefficient));

		glParameter->add_column_to_panel(roNoise, true);

		GLUI_Spinner *spNoiseAmplitudeR = new GLUI_Spinner(roNoise, "   Amplitude R:", GLUI_SPINNER_FLOAT, &myDome.fNoiseAmplitudeR, SP_NOISE_AMPLITUDE_R, ControlNoise);
		spNoiseAmplitudeR->set_float_limits(-100, 100);
		spNoiseAmplitudeR->set_speed((float)(0.1*myDisplay.fSpeedCoefficient));
		spNoiseAmplitudeR->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner *spNoiseAmplitudeZ = new GLUI_Spinner(roNoise, "   Amplitude Z:", GLUI_SPINNER_FLOAT, &myDome.fNoiseAmplitudeZ, SP_NOISE_AMPLITUDE_Z, ControlNoise);
		spNoiseAmplitudeZ->set_float_limits(-100, 100);
		spNoiseAmplitudeZ->set_speed((float)(0.1*myDisplay.fSpeedCoefficient));
		spNoiseAmplitudeZ->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner *spNoiseProgressiveR = new GLUI_Spinner(roNoise, "Progressive R:", GLUI_SPINNER_FLOAT, &myDome.fNoiseProgressiveR, SP_NOISE_PROGRESSIVE_R, ControlNoise);
		spNoiseProgressiveR->set_float_limits(-10, 10);
		spNoiseProgressiveR->set_alignment(GLUI_ALIGN_RIGHT);

		GLUI_Spinner *spNoiseProgressiveZ = new GLUI_Spinner(roNoise, "Progressive Z:", GLUI_SPINNER_FLOAT, &myDome.fNoiseProgressiveZ, SP_NOISE_PROGRESSIVE_Z, ControlNoise);
		spNoiseProgressiveZ->set_float_limits(-10, 10);
		spNoiseProgressiveZ->set_alignment(GLUI_ALIGN_RIGHT);



		//  GLUI_Panel *InvisibleNoise = new GLUI_Panel(Noise, "", GLUI_PANEL_NONE);
		new GLUI_Button(roNoise, "Default", BT_NOISE_DEFAULT, ControlNoise);
	}
	if ("Buttons")
	{
		Largeur = 115;
		
		GLUI_Panel *InvisibleAutre = new GLUI_Panel(glParameter, "", GLUI_PANEL_NONE);
	
		btPreProcessing = new GLUI_Button(InvisibleAutre, "Pre-Processing", BT_PRE_PROCESSING, ControlPreProcessing);
		btPreProcessing->set_alignment(GLUI_ALIGN_CENTER);
		btPreProcessing->set_w(Largeur);

		glParameter->add_column_to_panel(InvisibleAutre, true);

		btDecor = new GLUI_Button(InvisibleAutre, "Decor", BT_DECOR, ControlDecor);
		btDecor->set_alignment(GLUI_ALIGN_CENTER);
		btDecor->set_w(Largeur);

		glParameter->add_column_to_panel(InvisibleAutre, true);

		btPostProcessing = new GLUI_Button(InvisibleAutre, "Post-Processing", BT_POST_PROCESSING, ControlPostProcessing);
		btPostProcessing->set_alignment(GLUI_ALIGN_CENTER);
		btPostProcessing->set_w(Largeur);

		if ("Display")
		{
			new GLUI_StaticText(glParameter, "                         ");
			GLUI_Listbox* lbDisplayType = new GLUI_Listbox(glParameter, "", &myDisplay.nDisplay, LB_DISPLAY, ControlDisplay);
			lbDisplayType->set_alignment(GLUI_ALIGN_CENTER);
			lbDisplayType->add_item(DISPLAY_Dome, "Dome                                 ");
			lbDisplayType->add_item(DISPLAY_Base, "Base");
			lbDisplayType->add_item(DISPLAY_Slope, "Slope");
			lbDisplayType->add_item(DISPLAY_Window, "Window");
			lbDisplayType->add_item(DISPLAY_Snake, "Snake");
			lbDisplayType->add_item(DISPLAY_Spirale, "Spirale");
			lbDisplayType->add_item(DISPLAY_Frise, "Frise");
			lbDisplayType->add_item(DISPLAY_Comb, "Comb");
			lbDisplayType->set_int_val(DISPLAY_Dome);
			new GLUI_StaticText(glParameter, "                         ");
			new GLUI_StaticText(glParameter, "                         ");

		}

		GLUI_Panel *InvisibleQuit = new GLUI_Panel(glParameter, "", GLUI_PANEL_NONE);


		GLUI_Button *btQuit = new GLUI_Button(InvisibleQuit, "Quit", 0, (GLUI_Update_CB)exit);
		btQuit->set_w(Largeur+30);
	}




	/**** Link windows to GLUI, and register idle callback ******/

	glParameter->set_main_gfx_window(MainWindow);
	glDecor->set_main_gfx_window(MainWindow);
	glPostProcess->set_main_gfx_window(MainWindow);
	glPreProcess->set_main_gfx_window(MainWindow);
	glBase->set_main_gfx_window(MainWindow);
	glSlope->set_main_gfx_window(MainWindow);
	glWindow->set_main_gfx_window(MainWindow);
	glTop->set_main_gfx_window(MainWindow);
	if ("Display")
	{/*** Create the bottom subwindow ***/
		glDisplay = GLUI_Master.create_glui_subwindow(MainWindow, GLUI_SUBWINDOW_BOTTOM);
		glDisplay->set_main_gfx_window(MainWindow);


		//  new GLUI_Column( glDisplay, false );
		GLUI_Rotation *rcDomeRotation = new GLUI_Rotation(glDisplay, "Rotation", DomeRotate);
		rcDomeRotation->set_spin((float)0.2);
		new GLUI_Column(glDisplay, false);
		GLUI_Translation *tcTranslationXY = new GLUI_Translation(glDisplay, "Translation XY", GLUI_TRANSLATION_XY, DomePosition);
		tcTranslationXY->set_speed((float)(0.005*myDisplay.fSpeedCoefficient));
		new GLUI_Column(glDisplay, false);
		GLUI_Translation *tcTranslationZ = new GLUI_Translation(glDisplay, "Translation Z", GLUI_TRANSLATION_Z, &DomePosition[2]);
		tcTranslationZ->set_speed((float)(0.005*myDisplay.fSpeedCoefficient));

		new GLUI_Column(glDisplay, false);


/*		rgDisplayGroup = new GLUI_RadioGroup(glDisplay, &myDisplay.nDisplayType, RG_DISPLAY, control_cb);
		GLUI_RadioButton *DisplayWireFrame = new GLUI_RadioButton(rgDisplayGroup, "WireFrame");
		GLUI_RadioButton *DisplayFace = new GLUI_RadioButton(rgDisplayGroup, "Flat Faces");
		GLUI_RadioButton *DisplaySmooth = new GLUI_RadioButton(rgDisplayGroup, "Smooth Faces");*/

		ckWireFrame  = new GLUI_Checkbox(glDisplay, "Wire Frame", &myDisplay.bWireFrame, CK_WIRE_FRAME, ControlDisplay);
		ckFlatFace   = new GLUI_Checkbox(glDisplay, "Flat Faces", &myDisplay.bFlatFace, CK_FLAT_FACE, ControlDisplay);
		ckSmoothFace = new GLUI_Checkbox(glDisplay, "Smooth Faces", &myDisplay.bSmooth, CK_SMOOTH_FACE, ControlDisplay);
		ckWireFrame->disable();
		GLUI_Listbox *lbMaterial = new GLUI_Listbox(glDisplay, "", &myDisplay.nMaterial, LB_MATERIAL, ControlDisplay);
		for (int i = 0; i < MAT_DERNIER; i++)
			lbMaterial->add_item(i, myMaterial[i].MaterialName);

		lbMaterial->set_int_val(myDisplay.nMaterial);


		new GLUI_Column(glDisplay, false);

		new GLUI_Checkbox(glDisplay, "Minor Struts", &myDisplay.bDisplayMinor, CK_DISPLAY_MINOR, ControlDisplay);
		GLUI_Spinner *spScaleDisplay = new GLUI_Spinner(glDisplay, "Scale:    ", GLUI_SPINNER_FLOAT, &myDisplay.fScaleDisplay);
		spScaleDisplay->set_float_limits((float)0.1, 10);
		spScaleDisplay->set_alignment(GLUI_ALIGN_LEFT);
		GLUI_Panel *Invisible2 = new GLUI_Panel(glDisplay, "", GLUI_PANEL_NONE);
		new GLUI_Checkbox(Invisible2, "Light 1", &myDisplay.nLight_1, CK_LIGHT_1, ControlDisplay);
		new GLUI_Column(Invisible2, false);
		new GLUI_Checkbox(Invisible2, "Light 2", &myDisplay.nLight_2, CK_LIGHT_2, ControlDisplay);

		new GLUI_Column(glDisplay, false);

		new GLUI_Button(glDisplay, "Default View", BT_DEFAULT_VIEW, ControlDisplay);
		new GLUI_Button(glDisplay, "Top View", BT_TOP_VIEW, ControlDisplay);
		new GLUI_Button(glDisplay, "Side View", BT_SIDE_VIEW, ControlDisplay);

	}


#if 0
	/**** We register the idle callback with GLUI, *not* with GLUT ****/
	GLUI_Master.set_glutIdleFunc(myGlutIdle);
#endif

	/**** Regular GLUT main loop ****/
	ControlDisplay(BT_DEFAULT_VIEW);
	glutMainLoop();

	return EXIT_SUCCESS;
}



/**
* glVertex2f(float x, float y).
* The point (0.0, 0.0) represents the middle of the window (not the top left corner).
* The "2f" suffix means 2 values of float type (x and y).
*/

/*
void displayMe(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	glBegin(GL_POLYGON);
	glVertex2f(0.0, 0.0);                    // bottom left
	glVertex2f(0.5, 0.8);                    // bottom right
	glVertex2f(0.5, 0.5);                    // top right
	glVertex2f(0.0, 0.5);                    // top left
	glEnd();
	glFlush();
}

int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE);
	glutInitWindowSize(300, 300);                    // window size
	glutInitWindowPosition(500, 500);                // distance from the top-left screen
	glutCreateWindow("BadproG - Hello world :D");    // message displayed on top bar window
	glutDisplayFunc(displayMe);
	glutMainLoop();
	return 0;
}

*/
