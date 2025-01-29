#include <string.h>
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#include "stdafx.h"
#include <time.h>
//#include "globales.h"
#include "functions.h"
#include "Types.h"

extern DomeParameters myDome;
extern char CRLF[];
extern DomeDisplay myDisplay;
extern GLpoint* myOBJ_VERTICE;
extern Base BaseFunctionListe[];
extern Slope SlopeFunctionListe[];
extern Window WindowFunctionListe[];
extern TestSlope TestSlopeFunctionListe[];
extern const char* ShrinkListe[];
extern const char* DistanceListe[];
extern Filtre FiltreListe[];

int PrintPOV()
{
	FILE* DomeINC;
	errno_t Erreur;
	char  FileName[400];
	char  DomeName[400];
	int   nMaxI;
	int   nMaxJ;

	clock_t StartProcess, EndProcess;
	double duration;

	GLpoint* myVERTICES;
	GLpoint* myNORMALS;
	int nMeridianEdgeSize = myDome.nMeridianEdgeSize;
	int nParallelEdgeSize = myDome.nParallelEdgeSize;

	float fSpikeMP = myDome.fSpikeMP;
	float fSpikeMp = myDome.fSpikeMp;
	float fSpikemP = myDome.fSpikemP;
	float fSpikemp = myDome.fSpikemp;

	bool bSpike = fSpikeMP > 0 || fSpikeMp > 0 || fSpikemP > 0 || fSpikemp > 0;
	bool bCover = true;

	ChargeShrink(myDome.nWireShrinkingType);
	nMaxI = (myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides + 1;
	nMaxJ = myDome.nParallelMain + (myDome.nParallelMain - 1) * myDome.nParallelMinor;


	if (myDisplay.bDisplayMinor)
		myVERTICES = myOBJ_VERTICE;
	else
		myVERTICES = (GLpoint*)malloc(nMaxI * nMaxJ * sizeof(GLpoint));

	if (myVERTICES == NULL)
	{
		printf("myVERTICES ko\n");
		return 0;
	}
	if (myDisplay.bDisplayMinor == 0)
	{
		CalculMesh(nMaxI, nMaxJ, myVERTICES);
		for (int i = 0; i < myDome.nLissageNombre; i++)
			Lissage(myVERTICES, nMaxI, nMaxJ);
		for (int i = 0; i < myDome.nDiffGeomNombre; i++)
			DifferentialGeometry(myVERTICES, nMaxI, nMaxJ);
		for (int i = 0; i < myDome.nNormaleNombre; i++)
			NormaleSharpening(myVERTICES, nMaxI, nMaxJ);
		Bumps(nMaxI, nMaxJ, myVERTICES, true);
		Marquees(nMaxI, nMaxJ, myVERTICES, true);
		Enhance(myVERTICES, nMaxI, nMaxJ);
	}

	myNORMALS = (GLpoint*)malloc(nMaxI * nMaxJ * sizeof(GLpoint));

	if (myDome.bFaceSmooth)
	{
		myNORMALS = (GLpoint*)malloc(nMaxI * nMaxJ * sizeof(GLpoint));
		if (myNORMALS == NULL)
		{
			printf("myNORMALS ko\n");
			return 0;
		}
		CalculNormales(nMaxI, nMaxJ, myVERTICES, myNORMALS);
	}

	strcpy_s(DomeName, myDome.sDomeName);

	//	strcpy_s(FileName, myDome.sDomeName);
	strcpy_s(FileName, "..\\INC\\");


	strcat_s(FileName, myDome.sDomeName);
	strcat_s(FileName, ".INC");

	Erreur = fopen_s(&DomeINC, FileName, "wb");
	if (Erreur == 2)
	{
		printf("Répertoire INC n'existe pas !\n");
		FileName[0] = 0; 
		strcat_s(FileName, myDome.sDomeName);
		strcat_s(FileName, ".INC");
		Erreur = fopen_s(&DomeINC, FileName, "wb");

	}
	if (Erreur != 0)
	{
		Message("     Error Creating POV File   ");
		printf("File : %s, Erreur = %d\n", FileName, Erreur);

		return 0;
	}
	if ("Comments")
	{
		fprintf(DomeINC, "// Generated Dome Object %s", CRLF);
		fprintf(DomeINC, "// Domenator Version %d %s", myDome.nDomeVersion, CRLF);
		fprintf(DomeINC, "%s%s", CRLF, CRLF);

		fprintf(DomeINC, "// nDomeVersion           = %d %s", myDome.nDomeVersion, CRLF);
		fprintf(DomeINC, "// sDomeName              = %s %s", myDome.sDomeName, CRLF);
		fprintf(DomeINC, "// fDomeScaleZ            = %f %s", myDome.fDomeScaleZ, CRLF);
		fprintf(DomeINC, "// sDomeMaterial          = %s %s", myDome.sDomeMaterial, CRLF);
		fprintf(DomeINC, "// nDomeEnd               = %d %s", myDome.nDomeEnd, CRLF);
		fprintf(DomeINC, "// fDomeSocle             = %f %s", myDome.fDomeSocle, CRLF);
		fprintf(DomeINC, "// fDomeInFolding         = %f %s", myDome.fDomeInFolding, CRLF);

		fprintf(DomeINC, "// nTopType               = %d %s", myDome.nTopType, CRLF);
		fprintf(DomeINC, "// fTopPart               = %f %s", myDome.fTopPart, CRLF);
		fprintf(DomeINC, "// fTopParameter          = %f %s", myDome.fTopParameter, CRLF);

		fprintf(DomeINC, "// nBase                  = %d = %s %s", myDome.nBase, BaseFunctionListe[myDome.nBase].sNom, CRLF);
		fprintf(DomeINC, "// nBaseSides             = %d %s", myDome.nBaseSides, CRLF);
		fprintf(DomeINC, "// fBaseParameter         = %f %s", myDome.fBaseParameter, CRLF);
		fprintf(DomeINC, "// fBaseSmoothing         = %f %s", myDome.fBaseSmoothing, CRLF);
		fprintf(DomeINC, "// fBaseExpand            = %f %s", myDome.fBaseExpand, CRLF);
		fprintf(DomeINC, "// fBaseDamping           = %f %s", myDome.fBaseDamping, CRLF);
		fprintf(DomeINC, "// fBaseSpiral            = %f %s", myDome.fBaseSpiral, CRLF);
		fprintf(DomeINC, "// nBaseEvenly            = %d %s", myDome.nBaseEvenly, CRLF);

		fprintf(DomeINC, "// nBaseSecondary         = %d = %s %s", myDome.nBaseSecondary, BaseFunctionListe[myDome.nBaseSecondary].sNom, CRLF);
		fprintf(DomeINC, "// fBaseSecondaryRotation = %f %s", myDome.fBaseSecondaryRotation, CRLF);
		fprintf(DomeINC, "// fBaseSecondaryExpansion= %f %s", myDome.fBaseSecondaryExpansion, CRLF);
		fprintf(DomeINC, "// bBaseMorphingDivision  = %d %s", myDome.bBaseMorphingDivision, CRLF);

		fprintf(DomeINC, "// fBaseRound             = %f %s", myDome.fBaseRound, CRLF);
		fprintf(DomeINC, "// fBaseBaryCentre        = %f %s", myDome.fBaseBaryCentre, CRLF);
		fprintf(DomeINC, "// fBaseMorphing          = %f %s", myDome.fBaseMorphing, CRLF);

		fprintf(DomeINC, "// nSlope                 = %d = %s %s", myDome.nSlope, SlopeFunctionListe[myDome.nSlope].sNom, CRLF);
		fprintf(DomeINC, "// nSlopeTwist            = %d %s", myDome.nSlopeTwist, CRLF);
		fprintf(DomeINC, "// fSlopeDelay            = %f %s", myDome.fSlopeDelay, CRLF);
		fprintf(DomeINC, "// nSlopeEvenly           = %d %s", myDome.nSlopeEvenly, CRLF);
		fprintf(DomeINC, "// fSlopeParameterZ       = %f %s", myDome.fSlopeParameterZ, CRLF);
		fprintf(DomeINC, "// fSlopeParameterXY      = %f %s", myDome.fSlopeParameterXY, CRLF);
		fprintf(DomeINC, "// bSlopeInverse          = %d %s", myDome.bSlopeInverse, CRLF);
		fprintf(DomeINC, "// nTestSlopeXY           = %d = %s %s", myDome.nTestSlopeXY, SlopeFunctionListe[myDome.nTestSlopeXY].sNom, CRLF);
		fprintf(DomeINC, "// nTestSlopeZ            = %d = %s %s", myDome.nTestSlopeZ, SlopeFunctionListe[myDome.nTestSlopeZ].sNom, CRLF);
			

		fprintf(DomeINC, "// nCrenelPerSide         = %d %s", myDome.nCrenelPerSide, CRLF);
		fprintf(DomeINC, "// fCrenelFront           = %f %s", myDome.fCrenelFront, CRLF);
		fprintf(DomeINC, "// fCrenelThick           = %f %s", myDome.fCrenelThick, CRLF);
		fprintf(DomeINC, "// fCrenelWide            = %f %s", myDome.fCrenelWide, CRLF);
		fprintf(DomeINC, "// fCrenelHeight          = %f %s", myDome.fCrenelHeight, CRLF);
		fprintf(DomeINC, "// fCrenelPower           = %f %s", myDome.fCrenelPower, CRLF);
		fprintf(DomeINC, "// fCrenelAngle           = %f %s", myDome.fCrenelAngle, CRLF);
		fprintf(DomeINC, "// fCrenelDelay           = %f %s", myDome.fCrenelDelay, CRLF);


		fprintf(DomeINC, "// nWindow                = %d = %s %s", myDome.nWindow, WindowFunctionListe[myDome.nWindow].sNom, CRLF);
		fprintf(DomeINC, "// fWindowParameter       = %f %s", myDome.fWindowParameter, CRLF);
		fprintf(DomeINC, "// fWindowPower           = %f %s", myDome.fWindowPower, CRLF);
		fprintf(DomeINC, "// fWindowHeight          = %f %s", myDome.fWindowHeight, CRLF);
		fprintf(DomeINC, "// fWindowDamping         = %f %s", myDome.fWindowDamping, CRLF);
		fprintf(DomeINC, "// fWindowPerSide         = %f %s", myDome.fWindowPerSide, CRLF);
		fprintf(DomeINC, "// fWindowShift           = %f %s", myDome.fWindowShift, CRLF);
		fprintf(DomeINC, "// bWindowInverse         = %d %s", myDome.bWindowInverse, CRLF);
		fprintf(DomeINC, "// fWindowDelay           = %f %s", myDome.fWindowDelay, CRLF);
		fprintf(DomeINC, "// fWindowAttack          = %f %s", myDome.fWindowAttack, CRLF);
		fprintf(DomeINC, "// fWindowAwning          = %f %s", myDome.fWindowAwning, CRLF);
		fprintf(DomeINC, "// bWindowBlind           = %d %s", myDome.bWindowBlind, CRLF);
		fprintf(DomeINC, "// fWindowLower           = %f %s", myDome.fWindowLower, CRLF);

		fprintf(DomeINC, "// bFaceGenerate          = %d %s", myDome.bFaceGenerate, CRLF);
		fprintf(DomeINC, "// bFaceSmooth            = %d %s", myDome.bFaceSmooth, CRLF);
		fprintf(DomeINC, "// sFaceCenterMaterial    = %s %s", myDome.sFaceCenterMaterial, CRLF);
		fprintf(DomeINC, "// sFaceEdgeMaterial      = %s %s", myDome.sFaceEdgeMaterial, CRLF);
		fprintf(DomeINC, "// nMeridianEdgeSize      = %d %s", myDome.nMeridianEdgeSize, CRLF);
		fprintf(DomeINC, "// nParallelEdgeSize      = %d %s", myDome.nParallelEdgeSize, CRLF);
		fprintf(DomeINC, "// bWireGenerate          = %d %s", myDome.bWireGenerate, CRLF);
		fprintf(DomeINC, "// fWireShrinking         = %f %s", myDome.fWireShrinking, CRLF);
		fprintf(DomeINC, "// nWireShrinkingType     = %d = %s %s", myDome.nWireShrinkingType, ShrinkListe[myDome.nWireShrinkingType], CRLF);
		fprintf(DomeINC, "// fWireShrinkingSpeed    = %f %s", myDome.fWireShrinkingSpeed, CRLF);
		fprintf(DomeINC, "// nMeridianMain          = %d %s", myDome.nMeridianMain, CRLF);
		fprintf(DomeINC, "// fMeridianMainRadius    = %f %s", myDome.fMeridianMainRadius, CRLF);
		fprintf(DomeINC, "// sMeridianMainMaterial  = %s %s", myDome.sMeridianMainMaterial, CRLF);
		fprintf(DomeINC, "// nMeridianMinor         = %d %s", myDome.nMeridianMinor, CRLF);
		fprintf(DomeINC, "// fMeridianMinorRadius   = %f %s", myDome.fMeridianMinorRadius, CRLF);
		fprintf(DomeINC, "// sMeridianMinorMaterial = %s %s", myDome.sMeridianMinorMaterial, CRLF);
		fprintf(DomeINC, "// nParallelMain          = %d %s", myDome.nParallelMain, CRLF);
		fprintf(DomeINC, "// fParallelMainRadius    = %f %s", myDome.fParallelMainRadius, CRLF);
		fprintf(DomeINC, "// sParallelMainMaterial  = %s %s", myDome.sParallelMainMaterial, CRLF);
		fprintf(DomeINC, "// nParallelMinor         = %d %s", myDome.nParallelMinor, CRLF);
		fprintf(DomeINC, "// fParallelMinorRadius   = %f %s", myDome.fParallelMinorRadius, CRLF);
		fprintf(DomeINC, "// sParallelMinorMaterial = %s %s", myDome.sParallelMinorMaterial, CRLF);
		fprintf(DomeINC, "// fDiagonalURRadius      = %f %s", myDome.fDiagonalURRadius, CRLF);
		fprintf(DomeINC, "// sDiagonalURMaterial    = %s %s", myDome.sDiagonalURMaterial, CRLF);
		fprintf(DomeINC, "// fDiagonalDRRadius      = %f %s", myDome.fDiagonalDRRadius, CRLF);
		fprintf(DomeINC, "// sDiagonalDRMaterial    = %s %s", myDome.sDiagonalDRMaterial, CRLF);
		fprintf(DomeINC, "// fSphereMPRadius        = %f %s", myDome.fSphereMPRadius, CRLF);
		fprintf(DomeINC, "// sSphereMPMaterial      = %s %s", myDome.sSphereMPMaterial, CRLF);
		fprintf(DomeINC, "// fSphereMpRadius        = %f %s", myDome.fSphereMpRadius, CRLF);
		fprintf(DomeINC, "// sSphereMpMaterial      = %s %s", myDome.sSphereMpMaterial, CRLF);
		fprintf(DomeINC, "// fSpheremPRadius        = %f %s", myDome.fSpheremPRadius, CRLF);
		fprintf(DomeINC, "// sSpheremPMaterial      = %s %s", myDome.sSpheremPMaterial, CRLF);
		fprintf(DomeINC, "// fSpherempRadius        = %f %s", myDome.fSpherempRadius, CRLF);
		fprintf(DomeINC, "// sSpherempMaterial      = %s %s", myDome.sSpherempMaterial, CRLF);
		fprintf(DomeINC, "// nNoiseOctave           = %d %s", myDome.nNoiseOctave, CRLF);
		fprintf(DomeINC, "// fNoiseDensity          = %f %s", myDome.fNoiseDensity, CRLF);
		fprintf(DomeINC, "// nNoiseSeed             = %d %s", myDome.nNoiseSeed, CRLF);
		fprintf(DomeINC, "// fNoiseFrequencyU       = %f %s", myDome.fNoiseFrequencyU, CRLF);
		fprintf(DomeINC, "// fNoiseFrequencyV       = %f %s", myDome.fNoiseFrequencyV, CRLF);
		fprintf(DomeINC, "// fNoiseAmplitudeR       = %f %s", myDome.fNoiseAmplitudeR, CRLF);
		fprintf(DomeINC, "// fNoiseAmplitudeZ       = %f %s", myDome.fNoiseAmplitudeZ, CRLF);
		fprintf(DomeINC, "// fNoiseProgressiveR     = %f %s", myDome.fNoiseProgressiveR, CRLF);
		fprintf(DomeINC, "// fNoiseProgressiveZ     = %f %s", myDome.fNoiseProgressiveZ, CRLF);

		fprintf(DomeINC, "// bLissage               = %d %s", myDome.bLissage, CRLF);
		fprintf(DomeINC, "// fLissagePoidsSelf      = %f %s", myDome.fLissagePoidsSelf, CRLF);
		fprintf(DomeINC, "// fLissagePoidsHorizontal = %f %s", myDome.fLissagePoidsHorizontal, CRLF);
		fprintf(DomeINC, "// fLissagePoidsVertical  = %f %s", myDome.fLissagePoidsVertical, CRLF);
		fprintf(DomeINC, "// fLissagePoidsDiagonales= %f %s", myDome.fLissagePoidsDiagonales, CRLF);
		fprintf(DomeINC, "// nLissageNombre         = %d %s", myDome.nLissageNombre, CRLF);
		fprintf(DomeINC, "// nLissageTaille         = %d %s", myDome.nLissageTaille, CRLF);
		fprintf(DomeINC, "// fLissagePuissance      = %f %s", myDome.fLissagePuissance, CRLF);
		fprintf(DomeINC, "// nLissageDistanceType   = %d = %s %s", myDome.nLissageDistanceType, DistanceListe[myDome.nLissageDistanceType], CRLF);
		fprintf(DomeINC, "// fLissageSigma1         = %f %s", myDome.fLissageSigma1, CRLF);
		fprintf(DomeINC, "// fLissageSigma2         = %f %s", myDome.fLissageSigma2, CRLF);
		fprintf(DomeINC, "// nLissageTypeFiltre     = %d %s %s", myDome.nLissageTypeFiltre, FiltreListe[myDome.nLissageTypeFiltre].sNom, CRLF);
		fprintf(DomeINC, "// nLissageDirection      = %d %s", myDome.nLissageDirection, CRLF);

		fprintf(DomeINC, "// nDifferentialGeometry  = %d %s", myDome.nDifferentialGeometry, CRLF);
		fprintf(DomeINC, "// fDiffGeomLambda        = %f %s", myDome.fDiffGeomLambda, CRLF);
		fprintf(DomeINC, "// nDiffGeomNombre        = %d %s", myDome.nDiffGeomNombre, CRLF);
		fprintf(DomeINC, "// bDiffGeomMeridien      = %d %s", myDome.bDiffGeomMeridien, CRLF);
		fprintf(DomeINC, "// bDiffGeomParallel      = %d %s", myDome.bDiffGeomParallel, CRLF);
		fprintf(DomeINC, "// bNormale               = %d %s", myDome.bNormale, CRLF);
		fprintf(DomeINC, "// fNormaleLambda         = %f %s", myDome.fNormaleLambda, CRLF);
		fprintf(DomeINC, "// nNormaleNombre         = %d %s", myDome.nNormaleNombre, CRLF);
		fprintf(DomeINC, "// fDerivative            = %f %s", myDome.fDerivative, CRLF);
		fprintf(DomeINC, "// bDerivativeNorme       = %d %s", myDome.bDerivativeNorme, CRLF);


		fprintf(DomeINC, "// nSpirale               = %d %s", myDome.nSpirale, CRLF);
		fprintf(DomeINC, "// nSpiraleShrink         = %d %s", myDome.nSpiraleShrink, CRLF);
		fprintf(DomeINC, "// fSpiraleRadius         = %f %s", myDome.fSpiraleRadius, CRLF);
		fprintf(DomeINC, "// sSpiraleMaterial       = %s %s", myDome.sSpiraleMaterial, CRLF);
		fprintf(DomeINC, "// fSpiraleTour           = %f %s", myDome.fSpiraleTour, CRLF);
		fprintf(DomeINC, "// fSpiraleSpeed          = %f %s", myDome.fSpiraleSpeed, CRLF);
		fprintf(DomeINC, "// fSpiraleDelay          = %f %s", myDome.fSpiraleDelay, CRLF);
		fprintf(DomeINC, "// fSpiraleDecalage       = %f %s", myDome.fSpiraleDecalage, CRLF);
		fprintf(DomeINC, "// nSpiraleChirale        = %d %s", myDome.nSpiraleChirale, CRLF);
		fprintf(DomeINC, "// nSpiraleStep           = %d %s", myDome.nSpiraleStep, CRLF);
		fprintf(DomeINC, "// nSpiraleToward         = %d %s", myDome.nSpiraleToward, CRLF);
		fprintf(DomeINC, "// nSpiraleSlope          = %d %s %s", myDome.nSpiraleSlope, TestSlopeFunctionListe[myDome.nSpiraleSlope].sNom, CRLF);
		fprintf(DomeINC, "// fSpiraleSlopeParameter = %f %s", myDome.fSpiraleSlopeParameter, CRLF);

		fprintf(DomeINC, "// nFrise                  = %d %s", myDome.nFrise, CRLF);
		fprintf(DomeINC, "// fFriseHeight            = %f %s", myDome.fFriseHeight, CRLF);
		fprintf(DomeINC, "// fFrisePerSide           = %f %s", myDome.fFrisePerSide, CRLF);
		fprintf(DomeINC, "// fFriseDepart            = %f %s", myDome.fFriseDepart, CRLF);
		fprintf(DomeINC, "// fFriseShift             = %f %s", myDome.fFriseShift, CRLF);
		fprintf(DomeINC, "// fFriseRadius            = %f %s", myDome.fFriseRadius, CRLF);
		fprintf(DomeINC, "// nFriseShrink            = %d %s", myDome.nFriseShrink, CRLF);
		fprintf(DomeINC, "// nFriseFunction          = %d %s", myDome.nFriseFunction, CRLF);
		fprintf(DomeINC, "// sFriseMaterial          = %s %s", myDome.sFriseMaterial, CRLF);
		fprintf(DomeINC, "// nFriseToward            = %d %s", myDome.nFriseToward, CRLF);
		fprintf(DomeINC, "// nFriseStep              = %d %s", myDome.nFriseStep, CRLF);
		fprintf(DomeINC, "// fFrisePower             = %f %s", myDome.fFrisePower, CRLF);
		fprintf(DomeINC, "// fFriseOffset            = %f %s", myDome.fFriseOffset, CRLF);
		fprintf(DomeINC, "// fFriseParameter         = %f %s", myDome.fFriseParameter, CRLF);
		fprintf(DomeINC, "// fFriseGap               = %f %s", myDome.fFriseGap, CRLF);
		fprintf(DomeINC, "// fFriseSpacing           = %f %s", myDome.fFriseSpacing, CRLF);
		fprintf(DomeINC, "// bFriseInverse           = %d %s", myDome.bFriseInverse, CRLF);
		fprintf(DomeINC, "// bFriseReverse           = %d %s", myDome.bFriseReverse, CRLF);


		fprintf(DomeINC, "// nSnake                  = %d %s", myDome.nSnake, CRLF);
		fprintf(DomeINC, "// fSnakeTour              = %f %s", myDome.fSnakeTour, CRLF);
		fprintf(DomeINC, "// nSnakeEvolve            = %d %s", myDome.nSnakeEvolve, CRLF);
		fprintf(DomeINC, "// fSnakeShift             = %f %s", myDome.fSnakeShift, CRLF);
		fprintf(DomeINC, "// fSnakeRadius            = %f %s", myDome.fSnakeRadius, CRLF);
		fprintf(DomeINC, "// nSnakeShrink            = %d %s", myDome.nSnakeShrink, CRLF);
		fprintf(DomeINC, "// sSnakeMaterial          = %s %s", myDome.sSnakeMaterial, CRLF);
		fprintf(DomeINC, "// nSnakeChirale           = %d %s", myDome.nSnakeChirale, CRLF);
		fprintf(DomeINC, "// nSnakeStep              = %d %s", myDome.nSnakeStep, CRLF);
		fprintf(DomeINC, "// nSnakeToward            = %d %s", myDome.nSnakeToward, CRLF);
		fprintf(DomeINC, "// fSnakeOffset            = %f %s", myDome.fSnakeOffset, CRLF);
		fprintf(DomeINC, "// nSnakeSlope             = %d %s %s", myDome.nSnakeSlope, TestSlopeFunctionListe[myDome.nSnakeSlope].sNom, CRLF);
		fprintf(DomeINC, "// fSnakeSlopeParameter    = %f %s", myDome.fSnakeSlopeParameter, CRLF);

		fprintf(DomeINC, "// nComb                   = %d %s", myDome.nComb, CRLF);
		fprintf(DomeINC, "// nCombType               = %d %s", myDome.nCombType, CRLF);
		fprintf(DomeINC, "// fCombRadius             = %f %s", myDome.fCombRadius, CRLF);
		fprintf(DomeINC, "// nCombShrink             = %d %s", myDome.nCombShrink, CRLF);
		fprintf(DomeINC, "// sCombMaterial           = %s %s", myDome.sCombMaterial, CRLF);
		fprintf(DomeINC, "// fCombSpacing            = %f %s", myDome.fCombSpacing, CRLF);
		fprintf(DomeINC, "// nCombStep               = %d %s", myDome.nCombStep, CRLF);
		fprintf(DomeINC, "// nCombToward             = %d %s", myDome.nCombToward, CRLF);
		fprintf(DomeINC, "// nCombThread             = %d %s", myDome.nCombThread, CRLF);
		fprintf(DomeINC, "// nCombThreadGap          = %d %s", myDome.nCombThreadGap, CRLF);
		fprintf(DomeINC, "// fCombOffset             = %f %s", myDome.fCombOffset, CRLF);


		fprintf(DomeINC, "// bBump                   = %d %s", myDome.bBump, CRLF);
		fprintf(DomeINC, "// nBumpType               = %d %s", myDome.nBumpType, CRLF);
		fprintf(DomeINC, "// fBumpDamping            = %f %s", myDome.fBumpDamping, CRLF);
		fprintf(DomeINC, "// nBumpNormalize          = %d %s", myDome.nBumpNormalize, CRLF);
		fprintf(DomeINC, "// fBumpPower              = %f %s", myDome.fBumpPower, CRLF);
		fprintf(DomeINC, "// fBumpLambda             = %f %s", myDome.fBumpLambda, CRLF);
		fprintf(DomeINC, "// nBumpNormal             = %d %s", myDome.nBumpNormal, CRLF);
		fprintf(DomeINC, "// nBumpDistanceType       = %d %s", myDome.nBumpDistanceType, CRLF);
		fprintf(DomeINC, "// fBumpDistance           = %f %s", myDome.fBumpDistance, CRLF);
		fprintf(DomeINC, "// nBumpDistanceInverse    = %d %s", myDome.nBumpDistanceInverse, CRLF);

		fprintf(DomeINC, "// bMarquee                = %d %s", myDome.bMarquee, CRLF);
		fprintf(DomeINC, "// bMarqueeNormalize       = %d %s", myDome.bMarqueeNormalize, CRLF);
		fprintf(DomeINC, "// fMarqueeAlpha           = %f %s", myDome.fMarqueeAlpha, CRLF);
		fprintf(DomeINC, "// fMarqueeBeta            = %f %s", myDome.fMarqueeBeta, CRLF);
		fprintf(DomeINC, "// fMarqueeDamping         = %f %s", myDome.fMarqueeDamping, CRLF);
		fprintf(DomeINC, "// fMarqueeLambda          = %f %s", myDome.fMarqueeLambda, CRLF);
		fprintf(DomeINC, "// fMarqueeParameter       = %f %s", myDome.fMarqueeParameter, CRLF);
		fprintf(DomeINC, "// fMarqueeTop             = %f %s", myDome.fMarqueeTop, CRLF);
		fprintf(DomeINC, "// nMarqueeFunction        = %d %s %s", myDome.nMarqueeFunction, TestSlopeFunctionListe[myDome.nMarqueeFunction].sNom, CRLF);



		fprintf(DomeINC, "// nPostMeridianEvery      = %d %s", myDome.nPostMeridianEvery, CRLF);
		fprintf(DomeINC, "// nPostMeridianShift      = %d %s", myDome.nPostMeridianShift, CRLF);
		fprintf(DomeINC, "// nPostParallelEvolve     = %d %s", myDome.nPostParallelEvolve, CRLF);
		fprintf(DomeINC, "// nPostParallelEvery      = %d %s", myDome.nPostParallelEvery, CRLF);
		fprintf(DomeINC, "// nPostParallelFirst      = %d %s", myDome.nPostParallelFirst, CRLF);
		fprintf(DomeINC, "// bPostMeridianPerSide    = %d %s", myDome.bPostMeridianPerSide, CRLF);


		fprintf(DomeINC, "// fSpikeMP                = %f %s", myDome.fSpikeMP, CRLF);
		fprintf(DomeINC, "// fSpikeMp                = %f %s", myDome.fSpikeMp, CRLF);
		fprintf(DomeINC, "// fSpikemP                = %f %s", myDome.fSpikemP, CRLF);
		fprintf(DomeINC, "// fSpikemp                = %f %s", myDome.fSpikemp, CRLF);
		fprintf(DomeINC, "// fSpikeAngle             = %f %s", myDome.fSpikeAngle, CRLF);

		fprintf(DomeINC, "// bLid                    = %d %s", myDome.bLid, CRLF);
		fprintf(DomeINC, "// sLidMaterial            = %s %s", myDome.sLidMaterial, CRLF);

		fprintf(DomeINC, "// nEnhanceType            = %d %s", myDome.nEnhanceType, CRLF);
		fprintf(DomeINC, "// nEnhanceCalcul          = %d %s", myDome.nEnhanceCalcul, CRLF);
		fprintf(DomeINC, "// nEnhanceMerge           = %d %s", myDome.nEnhanceMerge, CRLF);
		fprintf(DomeINC, "// bEnhanceMeridian        = %d %s", myDome.bEnhanceMeridian, CRLF);
		fprintf(DomeINC, "// bEnhanceParallel        = %d %s", myDome.bEnhanceParallel, CRLF);
		fprintf(DomeINC, "// fEnhanceDamping         = %f %s", myDome.fEnhanceDamping, CRLF);
		fprintf(DomeINC, "// fEnhanceLambda          = %f %s", myDome.fEnhanceLambda, CRLF);


	}


	fprintf(DomeINC, "%s%s", CRLF, CRLF);

	fprintf(DomeINC, "// MaxI = %d ; MaxJ = %d %s", nMaxI, nMaxJ, CRLF);

	int nmi;
	if (myDome.nDomeEnd == 360)
		nmi = nMaxI;
	else
		nmi = (int)((float)(nMaxI * myDome.nDomeEnd) / 360.0) + 1;

	bool bTwoMeshes = (myDome.nParallelEdgeSize > 0 || myDome.nMeridianEdgeSize > 0) &&
		(myDome.nParallelMinor + 1 > 2 * myDome.nParallelEdgeSize && myDome.nMeridianMinor + 1 > 2 * myDome.nMeridianEdgeSize);


	//	fprintf(DomeINC, "#declare %s = union {%s", strrep(DomeName), CRLF);

	if ("Print")
	{
		if (myDome.bFaceGenerate || myDome.bWireGenerate || myDome.nSpirale > 0 || myDome.nFrise > 0 || myDome.nSnake > 0 || myDome.nComb > 0 || bSpike || bCover)
		{
			fprintf(DomeINC, "#declare %s = ", strrep(DomeName));
			if (myDome.bWireGenerate || myDome.nSpirale > 0 || myDome.nFrise > 0 || myDome.nSnake > 0 || myDome.nComb > 0 || bTwoMeshes || bSpike || bCover)
				fprintf(DomeINC, "union {%s", CRLF);
			else
				fprintf(DomeINC, "%s", CRLF);
		}

		printf("POV File\nDur%ce : Sphere      Cones       Faces       Spikes      Spirales      Frises      Snakes       Combs\n", 130);

		StartProcess = clock();

		if (myDome.bWireGenerate)
			PrintStrutsSphere(myVERTICES, DomeINC, nmi, nMaxJ, 0, PRINT_POV);
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("      %8.3lf", duration);
		StartProcess = clock();

		if (myDome.bWireGenerate)
			PrintStrutsCone(myVERTICES, DomeINC, nmi, nMaxJ, 0, PRINT_POV);
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("   %8.3lf", duration);
		StartProcess = clock();


		if (myDome.bFaceGenerate)
			PrintPOVMesh(myVERTICES, myNORMALS, DomeINC, nmi, nMaxJ, strrep(DomeName));
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("    %8.3lf", duration);
		StartProcess = clock();

		if (bSpike)
			PrintSpike(myVERTICES, DomeINC, nmi, nMaxJ, 0, PRINT_POV);

		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("     %8.3lf", duration);
		StartProcess = clock();
		if (myDome.nSpirale > 0)
			PrintSpirale(DomeINC, 0, PRINT_POV);
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("      %8.3lf", duration);
		StartProcess = clock();

		if (myDome.nFrise > 0)
			PrintFrise(DomeINC, 0, PRINT_POV);
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("    %8.3lf", duration);
		StartProcess = clock();

		if (myDome.nSnake > 0)
			PrintSnake(DomeINC, 0, PRINT_POV);
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("    %8.3lf", duration);
		StartProcess = clock();

		if (myDome.nComb > 0)
			PrintComb(DomeINC, 0, PRINT_POV);
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("    %8.3lf\n", duration);
		StartProcess = clock();

		PrintCloture(DomeINC, nmi, nMaxJ, myVERTICES, 0, PRINT_POV);
	}

	if (myDome.bWireGenerate || myDome.nSpirale > 0 || myDome.nFrise > 0 || myDome.nSnake > 0 || myDome.nComb > 0 || bTwoMeshes || bSpike || bCover)
	{
		if (strcmp(myDome.sDomeMaterial, S_NONE))
			fprintf(DomeINC, "  material { %s }} %s", myDome.sDomeMaterial, CRLF);
		else
			fprintf(DomeINC, "}%s", CRLF);
	}


	fflush(DomeINC);
	fclose(DomeINC);

	if (myDisplay.bDisplayMinor == 0)
		free(myVERTICES);

	if (myDome.bFaceSmooth)
		free(myNORMALS);
	return 1;
}

void PrintPOVRoundedCone(GLpoint Centre0, GLpoint Centre1, float fRadius0, float fRadius1, FILE* DomeINC, char* sMaterial)
{
	GLpoint AB, Centre;

	float Longueur, j, Rayon;
	if ((fRadius0 == 0 && fRadius1 == 0) || (Centre0.x == Centre1.x && Centre0.y == Centre1.y && Centre0.z == Centre1.z))
		return;

	AB.x = Centre1.x - Centre0.x;
	AB.y = Centre1.y - Centre0.y;
	AB.z = Centre1.z - Centre0.z;
	Longueur = sqrtf(AB.x * AB.x + AB.y * AB.y + AB.z * AB.z);

	if (Longueur == 0)
		return;

	if (fRadius1 <= fRadius0)
	{
		j = (fRadius0 - fRadius1) * fRadius1 / Longueur;
		Rayon = sqrtf(fRadius1 * fRadius1 + j * j);
		Centre.x = Centre1.x - j * AB.x / Longueur;
		Centre.y = Centre1.y - j * AB.y / Longueur;
		Centre.z = Centre1.z - j * AB.z / Longueur;
	}
	else
	{
		j = (fRadius1 - fRadius0) * fRadius1 / Longueur;
		Rayon = sqrtf(fRadius1 * fRadius1 + j * j);
		Centre.x = Centre1.x + j * AB.x / Longueur;
		Centre.y = Centre1.y + j * AB.y / Longueur;
		Centre.z = Centre1.z + j * AB.z / Longueur;
	}

	PrintPOVSphere(Centre, Rayon, DomeINC, sMaterial);
	PrintPOVCone(Centre0, Centre1, abs(fRadius0), abs(fRadius1), DomeINC, sMaterial);
	return;


	fprintf(DomeINC, "cone { < %lf, %lf, %lf >, %lf, < %lf, %lf, %lf >, %lf ",
		Centre0.x, Centre0.z, Centre0.y, abs(fRadius0),
		Centre1.x, Centre1.z, Centre1.y, abs(fRadius1));

	if (strcmp(sMaterial, S_NONE))
		fprintf(DomeINC, " material{%s}} %s", sMaterial, CRLF);
	else
		fprintf(DomeINC, "}%s", CRLF);

}
void PrintPOVCone(GLpoint Centre0, GLpoint Centre1, float fRadius0, float fRadius1, FILE* DomeINC, char* sMaterial)
{
	if ((fRadius0 == 0 && fRadius1 == 0) || (Centre0.x == Centre1.x && Centre0.y == Centre1.y && Centre0.z == Centre1.z))
		return;

	fprintf(DomeINC, "cone { < %9.8lf, %9.8lf, %9.8lf >, %lf, < %9.8lf, %9.8lf, %9.8lf >, %lf ",
		Centre0.x, Centre0.z, Centre0.y, abs(fRadius0),
		Centre1.x, Centre1.z, Centre1.y, abs(fRadius1));
	if (strcmp(sMaterial, S_NONE))
		fprintf(DomeINC, " material{%s}} %s", sMaterial, CRLF);
	else
		fprintf(DomeINC, "}%s", CRLF);
}
int PrintPOVMesh(GLpoint* Vertice, GLpoint* Normale, FILE* DomeINC, int nMaxI, int nMaxJ, char* DomeName)
{
	bool bEdge;
	int nMeridianEdgeSize = myDome.nMeridianEdgeSize;
	int nParallelEdgeSize = myDome.nParallelEdgeSize;

	if (myDome.nParallelEdgeSize > 0 || myDome.nMeridianEdgeSize > 0)
	{
		fprintf(DomeINC, "mesh {  // %s %s", DomeName, CRLF);

		for (int i = 0; i < nMaxI - 1; i++)
		{
			for (int j = 1; j < nMaxJ; j++)
			{
				bEdge = (j % (myDome.nParallelMinor + 1) >= 1 && j % (myDome.nParallelMinor + 1) <= nParallelEdgeSize)
					|| ((j % (myDome.nParallelMinor + 1) == 0 && nParallelEdgeSize != 0) || j % (myDome.nParallelMinor + 1) > myDome.nParallelMinor + 1 - nParallelEdgeSize)
					|| (i % (myDome.nMeridianMinor + 1) < nMeridianEdgeSize) || ((i % (myDome.nMeridianMinor + 1) >= myDome.nMeridianMinor + 1 - nMeridianEdgeSize));
				if (bEdge)
				{
					//					if ((i == 0)) fprintf(DomeINC, "//C'est là (i = %d ; j = %d) %s", i, j, CRLF);

					if (Colineaire(Vertice[i * nMaxJ + j], Vertice[(i + 1) * nMaxJ + j], Vertice[(i + 1) * nMaxJ + j - 1]) != 0)
						if (myDome.bFaceSmooth && (fabs(Normale[i * nMaxJ + j].x) > 0.00001 || fabs(Normale[i * nMaxJ + j].y) > 0.00001 || fabs(Normale[i * nMaxJ + j].z) > 0.00001)
							&& (fabs(Normale[(i + 1) * nMaxJ + j].x) > 0.00001 || fabs(Normale[(i + 1) * nMaxJ + j].y) > 0.00001 || fabs(Normale[(i + 1) * nMaxJ + j].z) > 0.00001)
							&& (fabs(Normale[(i + 1) * nMaxJ + j - 1].x) > 0.00001 || fabs(Normale[(i + 1) * nMaxJ + j - 1].y) > 0.00001 || fabs(Normale[(i + 1) * nMaxJ + j - 1].z) > 0.00001))
							fprintf(DomeINC, "smooth_triangle { <%lf, %lf, %lf>, <%lf, %lf, %lf>	<%lf, %lf, %lf>, <%lf, %lf, %lf>	<%lf, %lf, %lf>, <%lf, %lf, %lf> } %s",
								2 * Vertice[i * nMaxJ + j].x, 2 * Vertice[i * nMaxJ + j].z, 2 * Vertice[i * nMaxJ + j].y,
								Normale[i * nMaxJ + j].x, Normale[i * nMaxJ + j].z, Normale[i * nMaxJ + j].y,
								2 * Vertice[(i + 1) * nMaxJ + j].x, 2 * Vertice[(i + 1) * nMaxJ + j].z, 2 * Vertice[(i + 1) * nMaxJ + j].y,
								Normale[(i + 1) * nMaxJ + j].x, Normale[(i + 1) * nMaxJ + j].z, Normale[(i + 1) * nMaxJ + j].y,
								2 * Vertice[(i + 1) * nMaxJ + j - 1].x, 2 * Vertice[(i + 1) * nMaxJ + j - 1].z, 2 * Vertice[(i + 1) * nMaxJ + j - 1].y,
								Normale[(i + 1) * nMaxJ + j - 1].x, Normale[(i + 1) * nMaxJ + j - 1].z, Normale[(i + 1) * nMaxJ + j - 1].y,
								CRLF);
						else
							fprintf(DomeINC, "triangle { <%lf, %lf, %lf>, <%lf, %lf, %lf>, <%lf, %lf, %lf> }%s",
								2 * Vertice[i * nMaxJ + j].x, 2 * Vertice[i * nMaxJ + j].z, 2 * Vertice[i * nMaxJ + j].y,
								2 * Vertice[(i + 1) * nMaxJ + j].x, 2 * Vertice[(i + 1) * nMaxJ + j].z, 2 * Vertice[(i + 1) * nMaxJ + j].y,
								2 * Vertice[(i + 1) * nMaxJ + j - 1].x, 2 * Vertice[(i + 1) * nMaxJ + j - 1].z, 2 * Vertice[(i + 1) * nMaxJ + j - 1].y,
								CRLF);
					if (Colineaire(Vertice[i * nMaxJ + j], Vertice[i * nMaxJ + j - 1], Vertice[(i + 1) * nMaxJ + j - 1]) != 0)
						if (myDome.bFaceSmooth && (fabs(Normale[i * nMaxJ + j].x) > 0.00001 || fabs(Normale[i * nMaxJ + j].y) > 0.00001 || fabs(Normale[i * nMaxJ + j].z) > 0.00001)
							&& (fabs(Normale[i * nMaxJ + j - 1].x) > 0.00001 || fabs(Normale[i * nMaxJ + j - 1].y) > 0.00001 || fabs(Normale[i * nMaxJ + j - 1].z) > 0.00001)
							&& (fabs(Normale[(i + 1) * nMaxJ + j - 1].x) > 0.00001 || fabs(Normale[(i + 1) * nMaxJ + j - 1].y) > 0.00001 || fabs(Normale[(i + 1) * nMaxJ + j - 1].z) > 0.00001))
							fprintf(DomeINC, "smooth_triangle { <%lf, %lf, %lf>, <%lf, %lf, %lf>	<%lf, %lf, %lf>, <%lf, %lf, %lf>	<%lf, %lf, %lf>, <%lf, %lf, %lf> } %s",
								2 * Vertice[i * nMaxJ + j].x, 2 * Vertice[i * nMaxJ + j].z, 2 * Vertice[i * nMaxJ + j].y,
								Normale[i * nMaxJ + j].x, Normale[i * nMaxJ + j].z, Normale[i * nMaxJ + j].y,
								2 * Vertice[i * nMaxJ + j - 1].x, 2 * Vertice[i * nMaxJ + j - 1].z, 2 * Vertice[i * nMaxJ + j - 1].y,
								Normale[i * nMaxJ + j - 1].x, Normale[i * nMaxJ + j - 1].z, Normale[i * nMaxJ + j - 1].y,
								2 * Vertice[(i + 1) * nMaxJ + j - 1].x, 2 * Vertice[(i + 1) * nMaxJ + j - 1].z, 2 * Vertice[(i + 1) * nMaxJ + j - 1].y,
								Normale[(i + 1) * nMaxJ + j - 1].x, Normale[(i + 1) * nMaxJ + j - 1].z, Normale[(i + 1) * nMaxJ + j - 1].y,
								CRLF);
						else
							fprintf(DomeINC, "triangle { <%lf, %lf, %lf>, <%lf, %lf, %lf>, <%lf, %lf, %lf> }%s",
								2 * Vertice[i * nMaxJ + j].x, 2 * Vertice[i * nMaxJ + j].z, 2 * Vertice[i * nMaxJ + j].y,
								2 * Vertice[i * nMaxJ + j - 1].x, 2 * Vertice[i * nMaxJ + j - 1].z, 2 * Vertice[i * nMaxJ + j - 1].y,
								2 * Vertice[(i + 1) * nMaxJ + j - 1].x, 2 * Vertice[(i + 1) * nMaxJ + j - 1].z, 2 * Vertice[(i + 1) * nMaxJ + j - 1].y,
								CRLF);
				}
				//			printf("i = %d, j = %d, Meridian = %d, Parallele = %d, nMeridian = %d, nParallele = %d\n", i, j, myDome.nMeridianMain * myDome.nBaseSides, myDome.nParallelMain, i / (myDome.nMeridianMinor + 1), j / (myDome.nParallelMinor + 1));

			}
		}



		if (strcmp(myDome.sFaceEdgeMaterial, S_NONE))
			fprintf(DomeINC, " material{%s}} %s", myDome.sFaceEdgeMaterial, CRLF);
		else
			fprintf(DomeINC, "}%s", CRLF);
	}
	if (myDome.nParallelMinor + 1 > 2 * myDome.nParallelEdgeSize && myDome.nMeridianMinor + 1 > 2 * myDome.nMeridianEdgeSize) //Center
	{
		fprintf(DomeINC, "mesh {  // %s %s", strrep(DomeName), CRLF);

		for (int i = 0; i < nMaxI - 1; i++)
		{
			for (int j = 1; j < nMaxJ; j++)
			{
				bEdge = (j % (myDome.nParallelMinor + 1) >= 1 && j % (myDome.nParallelMinor + 1) <= nParallelEdgeSize)
					|| ((j % (myDome.nParallelMinor + 1) == 0 && nParallelEdgeSize != 0) || j % (myDome.nParallelMinor + 1) > myDome.nParallelMinor + 1 - nParallelEdgeSize)
					|| (i % (myDome.nMeridianMinor + 1) < nMeridianEdgeSize) || ((i % (myDome.nMeridianMinor + 1) >= myDome.nMeridianMinor + 1 - nMeridianEdgeSize));
				if (!bEdge)
				{
					//					if ((i == 0)) fprintf(DomeINC, "//C'est là (i = %d ; j = %d) %s", i, j, CRLF);

					if (Colineaire(Vertice[i * nMaxJ + j], Vertice[(i + 1) * nMaxJ + j], Vertice[(i + 1) * nMaxJ + j - 1]) != 0)
						if (myDome.bFaceSmooth && (fabs(Normale[i * nMaxJ + j].x) > 0.00001 || fabs(Normale[i * nMaxJ + j].y) > 0.00001 || fabs(Normale[i * nMaxJ + j].z) > 0.00001)
							&& (fabs(Normale[(i + 1) * nMaxJ + j].x) > 0.00001 || fabs(Normale[(i + 1) * nMaxJ + j].y) > 0.00001 || fabs(Normale[(i + 1) * nMaxJ + j].z) > 0.00001)
							&& (fabs(Normale[(i + 1) * nMaxJ + j - 1].x) > 0.00001 || fabs(Normale[(i + 1) * nMaxJ + j - 1].y) > 0.00001 || fabs(Normale[(i + 1) * nMaxJ + j - 1].z) > 0.00001))
							fprintf(DomeINC, "smooth_triangle { <%lf, %lf, %lf>, <%lf, %lf, %lf>	<%lf, %lf, %lf>, <%lf, %lf, %lf>	<%lf, %lf, %lf>, <%lf, %lf, %lf> } %s",
								2 * Vertice[i * nMaxJ + j].x, 2 * Vertice[i * nMaxJ + j].z, 2 * Vertice[i * nMaxJ + j].y,
								Normale[i * nMaxJ + j].x, Normale[i * nMaxJ + j].z, Normale[i * nMaxJ + j].y,
								2 * Vertice[(i + 1) * nMaxJ + j].x, 2 * Vertice[(i + 1) * nMaxJ + j].z, 2 * Vertice[(i + 1) * nMaxJ + j].y,
								Normale[(i + 1) * nMaxJ + j].x, Normale[(i + 1) * nMaxJ + j].z, Normale[(i + 1) * nMaxJ + j].y,
								2 * Vertice[(i + 1) * nMaxJ + j - 1].x, 2 * Vertice[(i + 1) * nMaxJ + j - 1].z, 2 * Vertice[(i + 1) * nMaxJ + j - 1].y,
								Normale[(i + 1) * nMaxJ + j - 1].x, Normale[(i + 1) * nMaxJ + j - 1].z, Normale[(i + 1) * nMaxJ + j - 1].y,
								CRLF);
						else
							fprintf(DomeINC, "triangle { <%lf, %lf, %lf>, <%lf, %lf, %lf>, <%lf, %lf, %lf> }%s",
								2 * Vertice[i * nMaxJ + j].x, 2 * Vertice[i * nMaxJ + j].z, 2 * Vertice[i * nMaxJ + j].y,
								2 * Vertice[(i + 1) * nMaxJ + j].x, 2 * Vertice[(i + 1) * nMaxJ + j].z, 2 * Vertice[(i + 1) * nMaxJ + j].y,
								2 * Vertice[(i + 1) * nMaxJ + j - 1].x, 2 * Vertice[(i + 1) * nMaxJ + j - 1].z, 2 * Vertice[(i + 1) * nMaxJ + j - 1].y,
								CRLF);
					if (Colineaire(Vertice[i * nMaxJ + j], Vertice[i * nMaxJ + j - 1], Vertice[(i + 1) * nMaxJ + j - 1]) != 0)
						if (myDome.bFaceSmooth && (fabs(Normale[i * nMaxJ + j].x) > 0.00001 || fabs(Normale[i * nMaxJ + j].y) > 0.00001 || fabs(Normale[i * nMaxJ + j].z) > 0.00001)
							&& (fabs(Normale[i * nMaxJ + j - 1].x) > 0.00001 || fabs(Normale[i * nMaxJ + j - 1].y) > 0.00001 || fabs(Normale[i * nMaxJ + j - 1].z) > 0.00001)
							&& (fabs(Normale[(i + 1) * nMaxJ + j - 1].x) > 0.00001 || fabs(Normale[(i + 1) * nMaxJ + j - 1].y) > 0.00001 || fabs(Normale[(i + 1) * nMaxJ + j - 1].z) > 0.00001))
							fprintf(DomeINC, "smooth_triangle { <%lf, %lf, %lf>, <%lf, %lf, %lf>	<%lf, %lf, %lf>, <%lf, %lf, %lf>	<%lf, %lf, %lf>, <%lf, %lf, %lf> } %s",
								2 * Vertice[i * nMaxJ + j].x, 2 * Vertice[i * nMaxJ + j].z, 2 * Vertice[i * nMaxJ + j].y,
								Normale[i * nMaxJ + j].x, Normale[i * nMaxJ + j].z, Normale[i * nMaxJ + j].y,
								2 * Vertice[i * nMaxJ + j - 1].x, 2 * Vertice[i * nMaxJ + j - 1].z, 2 * Vertice[i * nMaxJ + j - 1].y,
								Normale[i * nMaxJ + j - 1].x, Normale[i * nMaxJ + j - 1].z, Normale[i * nMaxJ + j - 1].y,
								2 * Vertice[(i + 1) * nMaxJ + j - 1].x, 2 * Vertice[(i + 1) * nMaxJ + j - 1].z, 2 * Vertice[(i + 1) * nMaxJ + j - 1].y,
								Normale[(i + 1) * nMaxJ + j - 1].x, Normale[(i + 1) * nMaxJ + j - 1].z, Normale[(i + 1) * nMaxJ + j - 1].y,
								CRLF);
						else
							fprintf(DomeINC, "triangle { <%lf, %lf, %lf>, <%lf, %lf, %lf>, <%lf, %lf, %lf> }%s",
								2 * Vertice[i * nMaxJ + j].x, 2 * Vertice[i * nMaxJ + j].z, 2 * Vertice[i * nMaxJ + j].y,
								2 * Vertice[i * nMaxJ + j - 1].x, 2 * Vertice[i * nMaxJ + j - 1].z, 2 * Vertice[i * nMaxJ + j - 1].y,
								2 * Vertice[(i + 1) * nMaxJ + j - 1].x, 2 * Vertice[(i + 1) * nMaxJ + j - 1].z, 2 * Vertice[(i + 1) * nMaxJ + j - 1].y,
								CRLF);
				}
				//				printf("i = %d, j = %d, Meridian = %d, Parallele = %d, nMeridian = %d, nParallele = %d\n", i, j, myDome.nMeridianMain* myDome.nBaseSides, myDome.nParallelMain, i / (myDome.nMeridianMinor + 1), j / (myDome.nParallelMinor + 1));

			}
		}

		if (strcmp(myDome.sFaceCenterMaterial, S_NONE))
			fprintf(DomeINC, " material{%s}} %s", myDome.sFaceCenterMaterial, CRLF);
		else
			fprintf(DomeINC, "}%s", CRLF);

	}
	return 0;
}
void PrintPOVSphere(GLpoint Centre, float fRadius, FILE* DomeINC, char* sMaterial)
{
	if (fRadius != 0.0)
	{
		fprintf(DomeINC, "sphere   { < %lf, %lf, %lf >, %lf ", Centre.x, Centre.z, Centre.y, abs(fRadius));
		if (strcmp(sMaterial, S_NONE))
			fprintf(DomeINC, " material{%s}} %s", sMaterial, CRLF);
		else
			fprintf(DomeINC, "}%s", CRLF);
	}
}
