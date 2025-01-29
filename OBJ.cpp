#include <string.h>
#include <stdio.h>
#include "stdafx.h"
#include <time.h>
//#include "globales.h"
#include "functions.h"
#include "Types.h"

extern DomeParameters myDome;
extern char CRLF[];
extern DomeDisplay myDisplay;
extern GLpoint* myOBJ_VERTICE;
extern int I_Sphere;
extern int J_Sphere;
extern int I_Cone;
extern int J_Cone;


int PrintOBJ()
{
	FILE* DomeOBJ;
	errno_t Erreur;
	char  FileName[400];
	char  DomeName[400];
	int   nMaxI;
	int   nMaxJ;
	GLpoint* myVERTICES;
	int nmi;

	int Start = 0;

	clock_t StartProcess, EndProcess;
	double duration;

	float fSpikeMP = myDome.fSpikeMP;
	float fSpikeMp = myDome.fSpikeMp;
	float fSpikemP = myDome.fSpikemP;
	float fSpikemp = myDome.fSpikemp;

	bool bSpike = fSpikeMP > 0 || fSpikeMp > 0 || fSpikemP > 0 || fSpikemp > 0;


	nMaxI = (myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides + 1;
	nMaxJ = myDome.nParallelMain + (myDome.nParallelMain - 1) * myDome.nParallelMinor;

	if (myDome.nDomeEnd != 360)
		nMaxI = (int)((float)(nMaxI * myDome.nDomeEnd) / 360.0) + 1;


	strcpy_s(DomeName, myDome.sDomeName);

	strcpy_s(FileName, "..\\");
	strcat_s(FileName, myDome.sDomeName);
	strcat_s(FileName, ".OBJ");

	Erreur = fopen_s(&DomeOBJ, FileName, "wb");
	if (Erreur != 0)
	{
		Message("     Error Creating OBJ File   ");
		printf("Erreur = %d\n", Erreur);
		return 0;
	}


	fprintf(DomeOBJ, "# %s %s", strrep(DomeName), CRLF);
	fprintf(DomeOBJ, "# nMaxI = %d ; nMaxJ = %d %s", nMaxI, nMaxJ, CRLF);
	fprintf(DomeOBJ, "mtllib Colors.mtl%s", CRLF);

	////////////////////////////////////////////////////////////////

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
	if (myDome.nDomeEnd == 360)
		nmi = nMaxI;
	else
		nmi = (int)((float)(nMaxI * myDome.nDomeEnd) / 360.0) + 1;

	if ("Print")
	{
		Start = 0;

		printf("OBJ File\nDur%ce : Sphere      Cones       Faces       Spikes      Spirales      Frises      Snakes       Combs\n", 130);

		StartProcess = clock();
		if (myDome.bWireGenerate)
			Start = PrintStrutsSphere(myVERTICES, DomeOBJ, nmi, nMaxJ, Start, PRINT_OBJ);
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("      %8.3lf", duration);

		StartProcess = clock();
		if (myDome.bWireGenerate)
			Start = PrintStrutsCone(myVERTICES, DomeOBJ, nmi, nMaxJ, Start, PRINT_OBJ);
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("   %8.3lf", duration);

		StartProcess = clock();
		if (myDome.bFaceGenerate)
			Start = PrintOBJMesh(DomeOBJ, myVERTICES, Start, nmi, nMaxJ);
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("       %5.3lf", duration);

		StartProcess = clock();
		if (bSpike)
			Start = PrintSpike(myVERTICES, DomeOBJ, nmi, nMaxJ, Start, PRINT_OBJ);
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("     %8.3lf", duration);

		StartProcess = clock();
		if (myDome.nSpirale > 0)
			Start = PrintSpirale(DomeOBJ, Start, PRINT_OBJ);
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("      %8.3lf", duration);

		StartProcess = clock();
		if (myDome.nFrise > 0)
			Start = PrintFrise(DomeOBJ, Start, PRINT_OBJ);
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("    %8.3lf", duration);

		StartProcess = clock();
		if (myDome.nSnake > 0)
			Start = PrintSnake(DomeOBJ, Start, PRINT_OBJ);
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("    %8.3lf", duration);

		StartProcess = clock();
		if (myDome.nComb > 0)
			Start = PrintComb(DomeOBJ, Start, PRINT_OBJ);
		EndProcess = clock();
		duration = (double)(EndProcess - StartProcess) / CLOCKS_PER_SEC;
		printf("    %8.3lf\n", duration);

		Start = PrintCloture(DomeOBJ, nMaxI, nMaxJ, myVERTICES, Start, PRINT_OBJ);


	}

	fflush(DomeOBJ);
	fclose(DomeOBJ);

	return 1;
}
int PrintOBJCone(GLpoint Centre1, GLpoint Centre2, double Rayon1, double Rayon2, FILE* DomeOBJ, int Start, const char* sTitre, char* sMaterial, bool bStraight)
{
	GLpoint retour;

	if ((Rayon1 == 0 && Rayon2 == 0) || (Centre1.x == Centre2.x && Centre1.y == Centre2.y && Centre1.z == Centre2.z))
		return 0;

	fprintf(DomeOBJ, "o %s %s", sTitre, CRLF);
	fprintf(DomeOBJ, "usemtl color_%s%s", sMaterial, CRLF);

	fprintf(DomeOBJ, "# Centre 1 : %10lf %10lf %10lf %s", Centre1.x, Centre1.y, Centre1.z, CRLF);
	fprintf(DomeOBJ, "# Centre 2 : %10lf %10lf %10lf %s", Centre2.x, Centre2.y, Centre2.z, CRLF);
	fprintf(DomeOBJ, "# Rayon 1  : %10lf  %s", Rayon1, CRLF);
	fprintf(DomeOBJ, "# Rayon 2  : %10lf  %s", Rayon2, CRLF);


	for (int j = 0; j < J_Cone; j++)
		for (int i = 0; i < I_Cone; i++)
		{
			retour = ConeFormula(Centre1, Centre2, Rayon1, Rayon2, (double)((double)i * 2.0 * _PI / (double)I_Cone), (double)((double)j / (double)(J_Cone - 1)), bStraight);
			fprintf(DomeOBJ, "v %10lf %10lf %10lf %s", 2 * retour.x, 2 * retour.y, 2 * retour.z, CRLF);

		}

	for (int j = 0; j < J_Cone - 1; j++)
		for (int i = 0; i < I_Cone; i++)
		{
			fprintf(DomeOBJ, "f %d %d %d %s", j * I_Cone + i + 1 + Start, (j + 1) * I_Cone + i + 1 + Start, (j + 1) * I_Cone + (i + 1) % I_Cone + (1) + Start, CRLF);
			fprintf(DomeOBJ, "f %d %d %d %s", j * I_Cone + i + 1 + Start, j * I_Cone + 1 + (i + I_Cone - 1) % I_Cone + Start, (j + 1) * I_Cone + i + 1 + Start, CRLF);
		}
	return I_Cone * J_Cone;
}


int PrintOBJMesh(FILE* DomeOBJ, GLpoint* myVERTICES, int Start, int nMaxI, int nMaxJ)
{
	int nMeridianEdgeSize = myDome.nMeridianEdgeSize;
	int nParallelEdgeSize = myDome.nParallelEdgeSize;
	bool Edge;

	nMaxI--;//£

	fprintf(DomeOBJ, "v %lf %lf %lf %s", 4 * myVERTICES[0].x, 4 * myVERTICES[0].y, 4 * myVERTICES[0].z, CRLF);

	for (int j = 1; j < nMaxJ; j++)
	{
		for (int i = 0; i < nMaxI; i++)
		{
			fprintf(DomeOBJ, "v %lf %lf %lf %s", 4 * myVERTICES[i * nMaxJ + j].x, 4 * myVERTICES[i * nMaxJ + j].y, 4 * myVERTICES[i * nMaxJ + j].z, CRLF);

		}
	}

	if (myDome.nParallelEdgeSize > 0 || myDome.nMeridianEdgeSize > 0)
	{
		fprintf(DomeOBJ, "o Mesh %s", CRLF);
		fprintf(DomeOBJ, "usemtl color_%s%s", myDome.sFaceEdgeMaterial, CRLF);
		for (int i = 0; i < nMaxI - (myDome.nDomeEnd != 360); i++)
		{
			Edge = myDome.nParallelEdgeSize > 0 || (i % (myDome.nMeridianMinor + 1) == 0 || i % (myDome.nMeridianMinor + 1) == myDome.nMeridianMinor);
			Edge = myDome.nParallelEdgeSize > 0 || (i % (myDome.nMeridianMinor + 1) < nMeridianEdgeSize || i % (myDome.nMeridianMinor + 1) > myDome.nMeridianMinor - nMeridianEdgeSize);
			if (Edge)
				fprintf(DomeOBJ, "f %d %d %d %s", Start + 2 + i, Start + 2 + (i + 1) % nMaxI, Start + 1, CRLF);
		}


		for (int j = 2; j < nMaxJ; j++)
			for (int i = 0; i < nMaxI - (myDome.nDomeEnd != 360); i++)
			{
				Edge = (j % (myDome.nParallelMinor + 1) >= 1 && j % (myDome.nParallelMinor + 1) <= nParallelEdgeSize)
					|| ((j % (myDome.nParallelMinor + 1) == 0 && nParallelEdgeSize != 0) || j % (myDome.nParallelMinor + 1) > myDome.nParallelMinor + 1 - nParallelEdgeSize)
					|| (i % (myDome.nMeridianMinor + 1) < nMeridianEdgeSize) || ((i % (myDome.nMeridianMinor + 1) >= myDome.nMeridianMinor + 1 - nMeridianEdgeSize));
				if (Edge)
				{
					fprintf(DomeOBJ, "f %d %d %d %s", Start + 2 + (j - 1) * nMaxI + i, Start + 2 + (j - 1) * nMaxI + (i + 1) % nMaxI, Start + 2 + (j - 2) * nMaxI + i, CRLF);
					fprintf(DomeOBJ, "f %d %d %d %s", Start + 2 + (j - 1) * nMaxI + (i + 1) % nMaxI, Start + 2 + (j - 2) * nMaxI + (i + 1) % nMaxI, Start + 2 + (j - 2) * nMaxI + i, CRLF);
				}
			}
	}
	if (myDome.nParallelMinor + 1 > 2 * nParallelEdgeSize && myDome.nParallelMinor + 1 > 2 * nMeridianEdgeSize)
	{
		fprintf(DomeOBJ, "o Mesh %s", CRLF);
		fprintf(DomeOBJ, "usemtl color_%s%s", myDome.sFaceCenterMaterial, CRLF);

		for (int i = 0; i < nMaxI - (myDome.nDomeEnd != 360); i++)
		{
			Edge = myDome.nParallelEdgeSize > 0 || (i % (myDome.nMeridianMinor + 1) == 0 || i % (myDome.nMeridianMinor + 1) == myDome.nMeridianMinor);
			Edge = myDome.nParallelEdgeSize > 0 || (i % (myDome.nMeridianMinor + 1) < nMeridianEdgeSize || i % (myDome.nMeridianMinor + 1) > myDome.nMeridianMinor - nMeridianEdgeSize);
			if (!Edge)
				fprintf(DomeOBJ, "f %d %d %d %s", Start + 2 + i, Start + 2 + (i + 1) % nMaxI, Start + 1, CRLF);
		}


		for (int j = 2; j < nMaxJ; j++)
			for (int i = 0; i < nMaxI - (myDome.nDomeEnd != 360); i++)
			{
				Edge = (j % (myDome.nParallelMinor + 1) >= 1 && j % (myDome.nParallelMinor + 1) <= nParallelEdgeSize)
					|| ((j % (myDome.nParallelMinor + 1) == 0 && nParallelEdgeSize != 0) || j % (myDome.nParallelMinor + 1) > myDome.nParallelMinor + 1 - nParallelEdgeSize)
					|| (i % (myDome.nMeridianMinor + 1) < nMeridianEdgeSize) || ((i % (myDome.nMeridianMinor + 1) >= myDome.nMeridianMinor + 1 - nMeridianEdgeSize));
				if (!Edge)
				{
					fprintf(DomeOBJ, "f %d %d %d %s", Start + 2 + (j - 1) * nMaxI + i, Start + 2 + (j - 1) * nMaxI + (i + 1) % nMaxI, Start + 2 + (j - 2) * nMaxI + i, CRLF);
					fprintf(DomeOBJ, "f %d %d %d %s", Start + 2 + (j - 1) * nMaxI + (i + 1) % nMaxI, Start + 2 + (j - 2) * nMaxI + (i + 1) % nMaxI, Start + 2 + (j - 2) * nMaxI + i, CRLF);
				}
			}

	}


	return Start + nMaxI * (nMaxJ - 1) + 1;

}

int PrintOBJSphere(GLpoint Centre, double Rayon, FILE* DomeOBJ, int Start, const char* sTitre, char* sMaterial)
{
	GLpoint myV;

	if (Rayon == 0)
		return 0;

	fprintf(DomeOBJ, "o %s %s", sTitre, CRLF);
	fprintf(DomeOBJ, "usemtl color_%s %s", sMaterial, CRLF);


	fprintf(DomeOBJ, "# Centre : %10lf %10lf %10lf %s", Centre.x, Centre.y, Centre.z, CRLF);
	fprintf(DomeOBJ, "# Rayon  : %10lf %s", Rayon, CRLF);


	myV = SphereFormula(Centre, Rayon, 0, 0);
	fprintf(DomeOBJ, "v %10lf %10lf %10lf %s", 2 * myV.x, 2 * myV.y, 2 * myV.z, CRLF);


	for (int j = 1; j < J_Sphere - 1; j++)
	{
		for (int i = 0; i < I_Sphere; i++)
		{
			myV = SphereFormula(Centre, Rayon, 2 * i * _PI / (I_Sphere), j * _PI / (J_Sphere - 1));
			fprintf(DomeOBJ, "v %10lf %10lf %10lf %s", 2 * myV.x, 2 * myV.y, 2 * myV.z, CRLF);
		}
	}
	myV = SphereFormula(Centre, Rayon, 2 * _PI, _PI);
	fprintf(DomeOBJ, "v %10lf %10lf %10lf %s", 2 * myV.x, 2 * myV.y, 2 * myV.z, CRLF);


	for (int i = 0; i < I_Sphere; i++)
		fprintf(DomeOBJ, "f %d %d %d %s", 1 + Start, ((i + 1) % I_Sphere) + 2 + Start, (i + 2) + Start, CRLF);

	for (int j = 1; j < J_Sphere - 2; j++)
		for (int i = 0; i < I_Sphere; i++)
		{
			fprintf(DomeOBJ, "f %d %d %d %s", (j - 1) * I_Sphere + i + 2 + Start, j * I_Sphere + (i + 1) % I_Sphere + 2 + Start, j * I_Sphere + i + 2 + Start, CRLF);
			fprintf(DomeOBJ, "f %d %d %d %s", (j - 1) * I_Sphere + i + 2 + Start, j * I_Sphere + i + 2 + Start, (j - 1) * I_Sphere + (i + I_Sphere - 1) % I_Sphere + 2 + Start, CRLF);

		}

	for (int i = 0; i < I_Sphere; i++)
	{
		fprintf(DomeOBJ, "f %d %d %d %s", I_Sphere * (J_Sphere - 2) + 2 + Start, I_Sphere * (J_Sphere - 3) + 2 + i + Start, I_Sphere * (J_Sphere - 3) + 2 + (i + 1) % I_Sphere + Start, CRLF);

	}
	return I_Sphere * (J_Sphere - 2) + 2;
}

