#include <string.h>
#include <stdio.h>

#include "stdafx.h"
#include <time.h>
#pragma once
#include "Types.h"
#include "Dome.h"
#include <GL\glui.h>
#include "Functions.h"

//#include "Globales.h"
extern Tube mySpiraleData;
extern DomeParameters myDome;
extern DomeDisplay myDisplay;
extern GLpoint* myOBJ_VERTICE;
extern Couleur stCouleur[MAX_COULEUR];

extern int I_Sphere;
extern int J_Sphere;
extern int I_Cone;
extern int J_Cone;
extern char CRLF[];
/*
int PrintStrutsSphere(GLpoint* Vertices, FILE* Dome, int nMaxI, int nMaxJ, int nbObjects, PRINT_TYPE Type);
int PrintStrutsCone(GLpoint* Vertices, FILE* Dome, int nMaxI, int nMaxJ, int nbObjects, PRINT_TYPE Type);
int PrintSpike(GLpoint* Vertice, FILE* Dome, int nMaxI, int nMaxJ, int nbObjects, PRINT_TYPE Type);
int PrintCloture(FILE* Dome, int nMaxI, int nMaxJ, GLpoint* Poly, int nbObjects, PRINT_TYPE Type);
void  ChargeShrink(int nShrink);
void  CalculMesh(int nMaxI, int nMaxJ, GLpoint* Vertices);
void Lissage(GLpoint* Vertices, int nMaxI, int nMaxJ);
void Evolute(GLpoint* Vertices, int nMaxI, int nMaxJ);
void NormaleSharpening(GLpoint* Vertices, int nMaxI, int nMaxJ);
int PrintSpirale(FILE* Dome, int nbObjects, PRINT_TYPE Type);
int PrintSnake(FILE* Dome, int nbObjects, PRINT_TYPE Type);
int PrintFrise(FILE* Dome, int nbObjects, PRINT_TYPE Type);
void Bumps(int nMaxI, int nMaxJ, GLpoint* Vertices, bool bDraw);
int PrintSTLMesh(FILE* DomeSTL, GLpoint* myVERTICES, int nmi, int nMaxJ);
unsigned long PrintColor(FILE* DomeSTL);
void  Message(const char* sMessage, int Type = 0);
void  CrossProduct(GLpoint* V, GLpoint* W, GLpoint* Resultat);
void CrossProduct(Vector3D* V, Vector3D* W, Vector3D* Resultat);
GLpoint SphereFormula(GLpoint Centre, double Rayon, double u, double v);
GLpoint ConeFormula(GLpoint Centre1, GLpoint Centre2, double Rayon1, double Rayon2, double u, double v, bool Straight = true);
unsigned short GetColor(char* sColor);*/

int PrintSTL()
{
	FILE* DomeSTL;
	errno_t Erreur;
	char  FileName[400];
	char  DomeName[400];
	int   nMaxI;
	int   nMaxJ;


	GLpoint* myVERTICES;


	float fSpikeMP = myDome.fSpikeMP;
	float fSpikeMp = myDome.fSpikeMp;
	float fSpikemP = myDome.fSpikemP;
	float fSpikemp = myDome.fSpikemp;


	bool bSpike = fSpikeMP > 0 || fSpikeMp > 0 || fSpikemP > 0 || fSpikemp > 0;

	clock_t start, finish;
	double duration;

	STLHeader stHeader;
	long nbTriangles = 0;

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

	strcpy_s(DomeName, myDome.sDomeName);

	strcpy_s(FileName, "..\\");
	strcat_s(FileName, myDome.sDomeName);
	strcat_s(FileName, ".STL");

	Erreur = fopen_s(&DomeSTL, FileName, "wb");
	if (Erreur != 0)
	{
		Message("     Error Creating STL File   ");
		printf("Erreur = %d\n", Erreur);
		return 0;
	}

	strcpy_s(stHeader.sHeader, 80, (const char*)FileName);
	strcat_s(stHeader.sHeader, 80, "    UNITS=cm");
	stHeader.nbTriangles = 0;


	int nmi;

	if (myDome.nDomeEnd == 360)
		nmi = nMaxI;
	else
		nmi = (int)((float)(nMaxI * myDome.nDomeEnd) / 360.0) + 1;
	for (int i = 0; i < MAX_COULEUR; i++)
		stCouleur[i].bUsed = false;

	fwrite(&stHeader, 1, sizeof(STLHeader), DomeSTL);

	if (strncmp(myDome.sDomeName, "Color", 5) == 0)
	{
		nbTriangles = PrintColor(DomeSTL);
		stHeader.nbTriangles = nbTriangles;
		fseek(DomeSTL, sizeof(stHeader.sHeader), SEEK_SET);
		fwrite(&nbTriangles, 1, sizeof(long), DomeSTL);

		fflush(DomeSTL);
		fclose(DomeSTL);

		return 0;

	}

	if ("Print")
	{
		printf("STL File\nDur%ce : Sphere      Cones       Faces       Spikes      Spirales      Frises      Snakes       Combs\n", 130);

		start = clock();
		nbTriangles = 0;
		if (myDome.bWireGenerate)
			nbTriangles += PrintStrutsSphere(myVERTICES, DomeSTL, nmi, nMaxJ, nbTriangles, PRINT_STL_BINARY);
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("      %8.3lf", duration);

		start = clock();
		if (myDome.bWireGenerate)
			nbTriangles += PrintStrutsCone(myVERTICES, DomeSTL, nmi, nMaxJ, nbTriangles, PRINT_STL_BINARY);
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("   %8.3lf", duration);

		start = clock();
		if (myDome.bFaceGenerate)
			nbTriangles += PrintSTLMesh(DomeSTL, myVERTICES, nmi, nMaxJ);
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("       %5.3lf", duration);

		start = clock();
		if (bSpike)
			nbTriangles += PrintSpike(myVERTICES, DomeSTL, nmi, nMaxJ, 0, PRINT_STL_BINARY);
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("     %8.3lf", duration);

		start = clock();

		if (myDome.nSpirale > 0)
			nbTriangles += PrintSpirale(DomeSTL, 0, PRINT_STL_BINARY);
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("      %8.3lf", duration);

		start = clock();
		if (myDome.nFrise > 0)
			nbTriangles += PrintFrise(DomeSTL, 0, PRINT_STL_BINARY);
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("    %8.3lf", duration);

		start = clock();
		if (myDome.nSnake > 0)
			nbTriangles += PrintSnake(DomeSTL, 0, PRINT_STL_BINARY);
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("    %8.3lf", duration);

		start = clock();
		if (myDome.nComb > 0)
			nbTriangles += PrintComb(DomeSTL, 0, PRINT_STL_BINARY);
		finish = clock();
		duration = (double)(finish - start) / CLOCKS_PER_SEC;
		printf("    %8.3lf\n", duration);


		nbTriangles += PrintCloture(DomeSTL, nMaxI, nMaxJ, myVERTICES, nbTriangles, PRINT_STL_BINARY);
		stHeader.nbTriangles = nbTriangles;
		fseek(DomeSTL, sizeof(stHeader.sHeader), SEEK_SET);
		fwrite(&nbTriangles, 1, sizeof(long), DomeSTL);
	}

	fflush(DomeSTL);
	fclose(DomeSTL);

	if (myDisplay.bDisplayMinor == 0)
		free(myVERTICES);

	return 1;
}

bool PrintSTLTriangle(FILE* Dome, GLpoint P1, GLpoint P2, GLpoint P3, int nCouleur)
{
	GLpoint N, V1, V2;
	STLBinary stBinary;

	V1.x = P2.x - P1.x;
	V1.y = P2.y - P1.y;
	V1.z = P2.z - P1.z;
	V2.x = P3.x - P2.x;
	V2.y = P3.y - P2.y;
	V2.z = P3.z - P2.z;

	CrossProduct(&V1, &V2, &N);

	stBinary.Normal.x = (float)N.x;
	stBinary.Normal.y = (float)N.y;
	stBinary.Normal.z = (float)N.z;

	stBinary.Vertex1.x = 2 * P1.x;
	stBinary.Vertex1.y = 2 * P1.y;
	stBinary.Vertex1.z = 2 * P1.z;
	stBinary.Vertex2.x = 2 * P2.x;
	stBinary.Vertex2.y = 2 * P2.y;
	stBinary.Vertex2.z = 2 * P2.z;
	stBinary.Vertex3.x = 2 * P3.x;
	stBinary.Vertex3.y = 2 * P3.y;
	stBinary.Vertex3.z = 2 * P3.z;

	stBinary.Attribut = nCouleur;

	fwrite(&stBinary, 1, sizeof(STLBinary), Dome);

	return true;

}

int PrintSTLSphere(GLpoint Centre, float Rayon, FILE* DomeSTL)
{
	GLpoint P1, P2, P3;
	Vector3D N, V1, V2;

	if (Rayon == 0)
		return 0;

	P1 = SphereFormula(Centre, Rayon, 0, 0);

	for (int i = 0; i < I_Sphere; i++)
	{
		P2 = SphereFormula(Centre, Rayon, (float)(2.0f * _PI * (float)i) / (float)(I_Sphere), _PI / (float)(J_Sphere - 1));
		P3 = SphereFormula(Centre, Rayon, (float)(2.0f * _PI * (float)((i + 1) % I_Sphere)) / (float)(I_Sphere), _PI / (float)(J_Sphere - 1));

		V1.x = P2.x - P1.x;
		V1.y = P2.y - P1.y;
		V1.z = P2.z - P1.z;
		V2.x = P3.x - P2.x;
		V2.y = P3.y - P2.y;
		V2.z = P3.z - P2.z;

		CrossProduct(&V1, &V2, &N);


		fprintf(DomeSTL, "facet normal %16.12lf %16.12lf %16.12lf %s", N.x, N.y, N.z, CRLF);
		fprintf(DomeSTL, "outer loop %s", CRLF);
		fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P1.x, P1.y, P1.z, CRLF);
		fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P2.x, P2.y, P2.z, CRLF);
		fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P3.x, P3.y, P3.z, CRLF);
		fprintf(DomeSTL, "endloop %s", CRLF);
		fprintf(DomeSTL, "endfacet %s", CRLF);

	}

	for (int i = 0; i < I_Sphere; i++)
		for (int j = 1; j < J_Sphere - 2; j++)
		{
			P1 = SphereFormula(Centre, Rayon, (float)i * 2.0f * _PI / (float)(I_Sphere), (float)(_PI * (float)j) / (float)(J_Sphere - 1));
			P2 = SphereFormula(Centre, Rayon, (float)(i + 1) * 2.0f * _PI / (float)(I_Sphere), (float)(_PI * (float)j) / (float)(J_Sphere - 1));
			P3 = SphereFormula(Centre, Rayon, (float)(i + 1) * 2.0f * _PI / (float)(I_Sphere), (float)(_PI * (float)(j + 1)) / (float)(J_Sphere - 1));
			V1.x = P2.x - P1.x;
			V1.y = P2.y - P1.y;
			V1.z = P2.z - P1.z;
			V2.x = P3.x - P2.x;
			V2.y = P3.y - P2.y;
			V2.z = P3.z - P2.z;

			CrossProduct(&V1, &V2, &N);


			fprintf(DomeSTL, "facet normal %16.12lf %16.12lf %16.12lf %s", N.x, N.y, N.z, CRLF);
			fprintf(DomeSTL, "outer loop %s", CRLF);
			fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P1.x, P1.y, P1.z, CRLF);
			fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P2.x, P2.y, P2.z, CRLF);
			fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P3.x, P3.y, P3.z, CRLF);
			fprintf(DomeSTL, "endloop %s", CRLF);
			fprintf(DomeSTL, "endfacet %s", CRLF);

			P2 = SphereFormula(Centre, Rayon, (float)(i) * 2.0f * _PI / (float)(I_Sphere), (float)(_PI * (float)(j + 1)) / (float)(J_Sphere - 1));

			V1.x = P2.x - P1.x;
			V1.y = P2.y - P1.y;
			V1.z = P2.z - P1.z;
			V2.x = P3.x - P2.x;
			V2.y = P3.y - P2.y;
			V2.z = P3.z - P2.z;

			CrossProduct(&V1, &V2, &N);


			fprintf(DomeSTL, "facet normal %16.12lf %16.12lf %16.12lf %s", N.x, N.y, N.z, CRLF);
			fprintf(DomeSTL, "outer loop %s", CRLF);
			fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P1.x, P1.y, P1.z, CRLF);
			fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P2.x, P2.y, P2.z, CRLF);
			fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P3.x, P3.y, P3.z, CRLF);
			fprintf(DomeSTL, "endloop %s", CRLF);
			fprintf(DomeSTL, "endfacet %s", CRLF);


		}

	P1 = SphereFormula(Centre, Rayon, 0, _PI);
	for (int i = 0; i < I_Sphere; i++)
	{
		P2 = SphereFormula(Centre, Rayon, (float)(2.0f * _PI * (float)i) / (float)(I_Sphere), _PI * (J_Sphere - 2) / (float)(J_Sphere - 1));
		P3 = SphereFormula(Centre, Rayon, (float)(2.0f * _PI * (float)((i + 1) % I_Sphere)) / (float)(I_Sphere), _PI * (J_Sphere - 2) / (float)(J_Sphere - 1));

		V1.x = P2.x - P1.x;
		V1.y = P2.y - P1.y;
		V1.z = P2.z - P1.z;
		V2.x = P3.x - P2.x;
		V2.y = P3.y - P2.y;
		V2.z = P3.z - P2.z;

		CrossProduct(&V1, &V2, &N);


		fprintf(DomeSTL, "facet normal %16.12lf %16.12lf %16.12lf %s", N.x, N.y, N.z, CRLF);
		fprintf(DomeSTL, "outer loop %s", CRLF);
		fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P1.x, P1.y, P1.z, CRLF);
		fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P2.x, P2.y, P2.z, CRLF);
		fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P3.x, P3.y, P3.z, CRLF);
		fprintf(DomeSTL, "endloop %s", CRLF);
		fprintf(DomeSTL, "endfacet %s", CRLF);

	}

	return 2 * I_Sphere * (J_Sphere - 2);

}

int PrintSTLConeBinary(GLpoint Centre1, GLpoint Centre2, float Rayon1, float Rayon2, unsigned short nCouleur, FILE* DomeSTL, bool bStraight)
{
	GLpoint P1, P2, P3;
	Vector3D N, V1, V2;
	STLBinary stBinary;


	if ((Rayon1 == 0 && Rayon2 == 0) || (Centre1.x == Centre2.x && Centre1.y == Centre2.y && Centre1.z == Centre2.z))
		return 0;
	for (int i = 0; i < I_Cone; i++)
		for (int j = 0; j < J_Cone - 1; j++)
		{
			P1 = ConeFormula(Centre1, Centre2, Rayon1, Rayon2, (float)i * 2.0f * _PI / (float)(I_Cone), (float)(j) / (float)(J_Cone - 1), bStraight);
			P2 = ConeFormula(Centre1, Centre2, Rayon1, Rayon2, (float)((i + 1) % I_Cone) * 2.0f * _PI / (float)(I_Cone), (float)(j) / (float)(J_Cone - 1), bStraight);
			P3 = ConeFormula(Centre1, Centre2, Rayon1, Rayon2, (float)((i + 1) % I_Cone) * 2.0f * _PI / (float)(I_Cone), (float)((j + 1)) / (float)(J_Cone - 1), bStraight);
			V1.x = P2.x - P1.x;
			V1.y = P2.y - P1.y;
			V1.z = P2.z - P1.z;
			V2.x = P3.x - P2.x;
			V2.y = P3.y - P2.y;
			V2.z = P3.z - P2.z;

			CrossProduct(&V1, &V2, &N);


			stBinary.Normal.x = (float)N.x;
			stBinary.Normal.y = (float)N.y;
			stBinary.Normal.z = (float)N.z;

			stBinary.Vertex1.x = P1.x;
			stBinary.Vertex1.y = P1.y;
			stBinary.Vertex1.z = P1.z;
			stBinary.Vertex2.x = P2.x;
			stBinary.Vertex2.y = P2.y;
			stBinary.Vertex2.z = P2.z;
			stBinary.Vertex3.x = P3.x;
			stBinary.Vertex3.y = P3.y;
			stBinary.Vertex3.z = P3.z;

			stBinary.Attribut = nCouleur;
			fwrite(&stBinary, 1, sizeof(STLBinary), DomeSTL);

			P2 = ConeFormula(Centre1, Centre2, Rayon1, Rayon2, (float)(i) * 2.0f * _PI / (float)(I_Cone), (float)(j + 1) / (float)(J_Cone - 1), bStraight);

			V1.x = P2.x - P1.x;
			V1.y = P2.y - P1.y;
			V1.z = P2.z - P1.z;
			V2.x = P3.x - P2.x;
			V2.y = P3.y - P2.y;
			V2.z = P3.z - P2.z;

			CrossProduct(&V1, &V2, &N);


			stBinary.Normal.x = (float)N.x;
			stBinary.Normal.y = (float)N.y;
			stBinary.Normal.z = (float)N.z;

			stBinary.Vertex1.x = P1.x;
			stBinary.Vertex1.y = P1.y;
			stBinary.Vertex1.z = P1.z;
			stBinary.Vertex2.x = P2.x;
			stBinary.Vertex2.y = P2.y;
			stBinary.Vertex2.z = P2.z;
			stBinary.Vertex3.x = P3.x;
			stBinary.Vertex3.y = P3.y;
			stBinary.Vertex3.z = P3.z;

			stBinary.Attribut = nCouleur;
			fwrite(&stBinary, 1, sizeof(STLBinary), DomeSTL);

		}

	return 2 * I_Cone * (J_Cone - 1);

}

int PrintSTLCone(GLpoint Centre1, GLpoint Centre2, float Rayon1, float Rayon2, FILE* DomeSTL)
{
	GLpoint P1, P2, P3;
	Vector3D N, V1, V2;

	if ((Rayon1 == 0 && Rayon2 == 0) || (Centre1.x == Centre2.x && Centre1.y == Centre2.y && Centre1.z == Centre2.z))
		return 0;
	for (int i = 0; i < I_Cone; i++)
		for (int j = 0; j < J_Cone - 1; j++)
		{
			P1 = ConeFormula(Centre1, Centre2, Rayon1, Rayon2, (float)i * 2.0f * _PI / (float)(I_Cone), (float)(j) / (float)(J_Cone - 1));
			P2 = ConeFormula(Centre1, Centre2, Rayon1, Rayon2, (float)((i + 1) % I_Cone) * 2.0f * _PI / (float)(I_Cone), (float)(j) / (float)(J_Cone - 1));
			P3 = ConeFormula(Centre1, Centre2, Rayon1, Rayon2, (float)((i + 1) % I_Cone) * 2.0f * _PI / (float)(I_Cone), (float)((j + 1)) / (float)(J_Cone - 1));
			V1.x = P2.x - P1.x;
			V1.y = P2.y - P1.y;
			V1.z = P2.z - P1.z;
			V2.x = P3.x - P2.x;
			V2.y = P3.y - P2.y;
			V2.z = P3.z - P2.z;

			CrossProduct(&V1, &V2, &N);


			fprintf(DomeSTL, "facet normal %16.12lf %16.12lf %16.12lf %s", N.x, N.y, N.z, CRLF);
			fprintf(DomeSTL, "outer loop %s", CRLF);
			fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P1.x, P1.y, P1.z, CRLF);
			fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P2.x, P2.y, P2.z, CRLF);
			fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P3.x, P3.y, P3.z, CRLF);
			fprintf(DomeSTL, "endloop %s", CRLF);
			fprintf(DomeSTL, "endfacet %s", CRLF);

			P2 = ConeFormula(Centre1, Centre2, Rayon1, Rayon2, (float)(i) * 2.0f * _PI / (float)(I_Cone), (float)(j + 1) / (float)(J_Cone - 1));

			V1.x = P2.x - P1.x;
			V1.y = P2.y - P1.y;
			V1.z = P2.z - P1.z;
			V2.x = P3.x - P2.x;
			V2.y = P3.y - P2.y;
			V2.z = P3.z - P2.z;

			CrossProduct(&V1, &V2, &N);


			fprintf(DomeSTL, "facet normal %16.12lf %16.12lf %16.12lf %s", N.x, N.y, N.z, CRLF);
			fprintf(DomeSTL, "outer loop %s", CRLF);
			fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P1.x, P1.y, P1.z, CRLF);
			fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P2.x, P2.y, P2.z, CRLF);
			fprintf(DomeSTL, "vertex %16.12lf %16.12lf %16.12lf %s", P3.x, P3.y, P3.z, CRLF);
			fprintf(DomeSTL, "endloop %s", CRLF);
			fprintf(DomeSTL, "endfacet %s", CRLF);


		}

	return 2 * I_Cone * (J_Cone - 1);

}

int PrintSTLMesh(FILE* DomeSTL, GLpoint* myVERTICES, int nMaxI, int nMaxJ)
{
	bool  bEdge;
	GLpoint P1, P2, P3;
	Vector3D N, V1, V2;
	STLBinary stBinary;

	int nMeridianEdgeSize = myDome.nMeridianEdgeSize;
	int nParallelEdgeSize = myDome.nParallelEdgeSize;

	int nbTriangles = 0;
	for (int i = 0; i < nMaxI - 1; i++)
	{
		for (int j = 1; j < nMaxJ; j++)
		{
			bEdge = (j % (myDome.nParallelMinor + 1) >= 1 && j % (myDome.nParallelMinor + 1) <= nParallelEdgeSize)
				|| ((j % (myDome.nParallelMinor + 1) == 0 && nParallelEdgeSize != 0) || j % (myDome.nParallelMinor + 1) > myDome.nParallelMinor + 1 - nParallelEdgeSize)
				|| (i % (myDome.nMeridianMinor + 1) < nMeridianEdgeSize)
				|| ((i % (myDome.nMeridianMinor + 1) >= myDome.nMeridianMinor + 1 - nMeridianEdgeSize));
			P1.x = 2 * myVERTICES[i * nMaxJ + j].x;
			P1.y = 2 * myVERTICES[i * nMaxJ + j].y;
			P1.z = 2 * myVERTICES[i * nMaxJ + j].z;
			P2.x = 2 * myVERTICES[(i + 1) * nMaxJ + j].x;
			P2.y = 2 * myVERTICES[(i + 1) * nMaxJ + j].y;
			P2.z = 2 * myVERTICES[(i + 1) * nMaxJ + j].z;
			P3.x = 2 * myVERTICES[(i + 1) * nMaxJ + j - 1].x;
			P3.y = 2 * myVERTICES[(i + 1) * nMaxJ + j - 1].y;
			P3.z = 2 * myVERTICES[(i + 1) * nMaxJ + j - 1].z;

			V1.x = P2.x - P1.x;
			V1.y = P2.y - P1.y;
			V1.z = P2.z - P1.z;
			V2.x = P3.x - P2.x;
			V2.y = P3.y - P2.y;
			V2.z = P3.z - P2.z;

			CrossProduct(&V1, &V2, &N);


			stBinary.Normal.x = (float)N.x;
			stBinary.Normal.y = (float)N.y;
			stBinary.Normal.z = (float)N.z;

			stBinary.Vertex1.x = P1.x;
			stBinary.Vertex1.y = P1.y;
			stBinary.Vertex1.z = P1.z;
			stBinary.Vertex2.x = P2.x;
			stBinary.Vertex2.y = P2.y;
			stBinary.Vertex2.z = P2.z;
			stBinary.Vertex3.x = P3.x;
			stBinary.Vertex3.y = P3.y;
			stBinary.Vertex3.z = P3.z;

			if (bEdge)
				stBinary.Attribut = GetColor(myDome.sFaceEdgeMaterial);
			else
				stBinary.Attribut = GetColor(myDome.sFaceCenterMaterial);
			fwrite(&stBinary, 1, sizeof(STLBinary), DomeSTL);
			nbTriangles++;

			if (j != 1)
			{
				P2.x = 2 * myVERTICES[(i)*nMaxJ + j - 1].x;
				P2.y = 2 * myVERTICES[(i)*nMaxJ + j - 1].y;
				P2.z = 2 * myVERTICES[(i)*nMaxJ + j - 1].z;

				V1.x = P2.x - P1.x;
				V1.y = P2.y - P1.y;
				V1.z = P2.z - P1.z;
				V2.x = P3.x - P2.x;
				V2.y = P3.y - P2.y;
				V2.z = P3.z - P2.z;

				CrossProduct(&V1, &V2, &N);


				stBinary.Normal.x = (float)N.x;
				stBinary.Normal.y = (float)N.y;
				stBinary.Normal.z = (float)N.z;

				stBinary.Vertex1.x = P1.x;
				stBinary.Vertex1.y = P1.y;
				stBinary.Vertex1.z = P1.z;
				stBinary.Vertex2.x = P2.x;
				stBinary.Vertex2.y = P2.y;
				stBinary.Vertex2.z = P2.z;
				stBinary.Vertex3.x = P3.x;
				stBinary.Vertex3.y = P3.y;
				stBinary.Vertex3.z = P3.z;

				if (bEdge)
					stBinary.Attribut = GetColor(myDome.sFaceEdgeMaterial);
				else
					stBinary.Attribut = GetColor(myDome.sFaceCenterMaterial);
				fwrite(&stBinary, 1, sizeof(STLBinary), DomeSTL);
				nbTriangles++;
			}
		}
	}

	return nbTriangles;
}

int PrintSTLSphereBinary(GLpoint Centre, float Rayon, unsigned short nCouleur, FILE* DomeSTL)
{
	GLpoint P1, P2, P3;
	Vector3D N, V1, V2;
	STLBinary stBinary;

	if (Rayon == 0)
		return 0;

	P1 = SphereFormula(Centre, Rayon, 0, 0);

	for (int i = 0; i < I_Sphere; i++)
	{
		P2 = SphereFormula(Centre, Rayon, (float)(2.0f * _PI * (float)i) / (float)(I_Sphere), _PI / (float)(J_Sphere - 1));
		P3 = SphereFormula(Centre, Rayon, (float)(2.0f * _PI * (float)((i + 1) % I_Sphere)) / (float)(I_Sphere), _PI / (float)(J_Sphere - 1));

		V1.x = P2.x - P1.x;
		V1.y = P2.y - P1.y;
		V1.z = P2.z - P1.z;
		V2.x = P3.x - P2.x;
		V2.y = P3.y - P2.y;
		V2.z = P3.z - P2.z;

		CrossProduct(&V1, &V2, &N);

		stBinary.Normal.x = (float)N.x;
		stBinary.Normal.y = (float)N.y;
		stBinary.Normal.z = (float)N.z;

		stBinary.Vertex1.x = P1.x;
		stBinary.Vertex1.y = P1.y;
		stBinary.Vertex1.z = P1.z;
		stBinary.Vertex2.x = P2.x;
		stBinary.Vertex2.y = P2.y;
		stBinary.Vertex2.z = P2.z;
		stBinary.Vertex3.x = P3.x;
		stBinary.Vertex3.y = P3.y;
		stBinary.Vertex3.z = P3.z;

		stBinary.Attribut = nCouleur;

		fwrite(&stBinary, 1, sizeof(STLBinary), DomeSTL);

	}

	for (int i = 0; i < I_Sphere; i++)
		for (int j = 1; j < J_Sphere - 2; j++)
		{
			P1 = SphereFormula(Centre, Rayon, (float)i * 2.0f * _PI / (float)(I_Sphere), (float)(_PI * (float)j) / (float)(J_Sphere - 1));
			P2 = SphereFormula(Centre, Rayon, (float)(i + 1) * 2.0f * _PI / (float)(I_Sphere), (float)(_PI * (float)j) / (float)(J_Sphere - 1));
			P3 = SphereFormula(Centre, Rayon, (float)(i + 1) * 2.0f * _PI / (float)(I_Sphere), (float)(_PI * (float)(j + 1)) / (float)(J_Sphere - 1));
			V1.x = P2.x - P1.x;
			V1.y = P2.y - P1.y;
			V1.z = P2.z - P1.z;
			V2.x = P3.x - P2.x;
			V2.y = P3.y - P2.y;
			V2.z = P3.z - P2.z;

			CrossProduct(&V1, &V2, &N);

			stBinary.Normal.x = (float)N.x;
			stBinary.Normal.y = (float)N.y;
			stBinary.Normal.z = (float)N.z;

			stBinary.Vertex1.x = P1.x;
			stBinary.Vertex1.y = P1.y;
			stBinary.Vertex1.z = P1.z;
			stBinary.Vertex2.x = P2.x;
			stBinary.Vertex2.y = P2.y;
			stBinary.Vertex2.z = P2.z;
			stBinary.Vertex3.x = P3.x;
			stBinary.Vertex3.y = P3.y;
			stBinary.Vertex3.z = P3.z;

			stBinary.Attribut = nCouleur;

			fwrite(&stBinary, 1, sizeof(STLBinary), DomeSTL);


			P2 = SphereFormula(Centre, Rayon, (float)(i) * 2.0f * _PI / (float)(I_Sphere), (float)(_PI * (float)(j + 1)) / (float)(J_Sphere - 1));

			V1.x = P2.x - P1.x;
			V1.y = P2.y - P1.y;
			V1.z = P2.z - P1.z;
			V2.x = P3.x - P2.x;
			V2.y = P3.y - P2.y;
			V2.z = P3.z - P2.z;

			CrossProduct(&V1, &V2, &N);


			stBinary.Normal.x = (float)N.x;
			stBinary.Normal.y = (float)N.y;
			stBinary.Normal.z = (float)N.z;

			stBinary.Vertex1.x = P1.x;
			stBinary.Vertex1.y = P1.y;
			stBinary.Vertex1.z = P1.z;
			stBinary.Vertex2.x = P2.x;
			stBinary.Vertex2.y = P2.y;
			stBinary.Vertex2.z = P2.z;
			stBinary.Vertex3.x = P3.x;
			stBinary.Vertex3.y = P3.y;
			stBinary.Vertex3.z = P3.z;

			stBinary.Attribut = nCouleur;

			fwrite(&stBinary, 1, sizeof(STLBinary), DomeSTL);
		}

	P1 = SphereFormula(Centre, Rayon, 0, _PI);
	for (int i = 0; i < I_Sphere; i++)
	{
		P2 = SphereFormula(Centre, Rayon, (float)(2.0f * _PI * (float)i) / (float)(I_Sphere), _PI * (J_Sphere - 2) / (float)(J_Sphere - 1));
		P3 = SphereFormula(Centre, Rayon, (float)(2.0f * _PI * (float)((i + 1) % I_Sphere)) / (float)(I_Sphere), _PI * (J_Sphere - 2) / (float)(J_Sphere - 1));

		V1.x = P2.x - P1.x;
		V1.y = P2.y - P1.y;
		V1.z = P2.z - P1.z;
		V2.x = P3.x - P2.x;
		V2.y = P3.y - P2.y;
		V2.z = P3.z - P2.z;

		CrossProduct(&V1, &V2, &N);


		stBinary.Normal.x = (float)N.x;
		stBinary.Normal.y = (float)N.y;
		stBinary.Normal.z = (float)N.z;

		stBinary.Vertex1.x = P1.x;
		stBinary.Vertex1.y = P1.y;
		stBinary.Vertex1.z = P1.z;
		stBinary.Vertex2.x = P2.x;
		stBinary.Vertex2.y = P2.y;
		stBinary.Vertex2.z = P2.z;
		stBinary.Vertex3.x = P3.x;
		stBinary.Vertex3.y = P3.y;
		stBinary.Vertex3.z = P3.z;

		stBinary.Attribut = nCouleur;

		fwrite(&stBinary, 1, sizeof(STLBinary), DomeSTL);
	}

	return 2 * I_Sphere * (J_Sphere - 2);

}
