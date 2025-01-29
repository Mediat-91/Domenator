#pragma once
#include "Types.h"
#include "Dome.h"
#include "Functions.h"


Tube *mySpiraleData;
Tube *myFriseData;
Tube *mySnakeData;
Tube* myCombData;

float xy_aspect;
int   MainWindow;

char ListeCouleur[150][40];

float scale = 1.0;
float DomeRotate[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
float InstantRotate[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
float DomePosition[] = { 0.0, 0.0, 0.0 };

int   last_x, last_y;
float rotationX = 0.0, rotationY = 0.0;

char  CRLF[] = { 0x0d, 0x0a, 0 };
std::string text;

GLpoint *myOBJ_VERTICE;

Parametric_3 BaseFunctionX;
Parametric_3 BaseFunctionY;
Parametric_2 SlopeFunctionXY;
Parametric_2 SlopeFunctionZ;
Parametric_4 WindowFunctionZ;
Parametric_4 FriseFunction;

Parametric_3 BaseFunctionSecondaryX;
Parametric_3 BaseFunctionSecondaryY;

Parametric_2 ShrinkFunction;
Parametric_3 ShrinkCoefficientFunction;



DomeParameters myDome;
DomeDisplay myDisplay;

int bAllocation = 1;


/** Pointers to the windows and some of the controls we'll create **/
GLUI			 *glParameter;
GLUI             *glDisplay;
GLUI             *glMessage;
GLUI			 *glPreProcess;
GLUI			 *glPostProcess;
GLUI			 *glDecor;
GLUI			 *glBase;
GLUI			 *glSlope;
GLUI			 *glWindow;
GLUI			 *glTop;
GLUI_TextBox* tbDirectory;

GLUI_StaticText  *stMessage;
GLUI_Spinner *spSide;

GLUI_RadioGroup  *rgDisplayGroup;

GLUI_TextBox* tbBase;
GLUI_TextBox* tbSlope;
GLUI_TextBox* tbWindow;
GLUI_TextBox* tbTop;

GLUI_Checkbox* ckDivision;
GLUI_Checkbox *ckFaceSmooth;
GLUI_Checkbox *ckWireFrame;
GLUI_Checkbox *ckFlatFace;
GLUI_Checkbox *ckSmoothFace;
GLUI_Button *btDecor, *btDecorQuit;
GLUI_Button *btPostProcessing, *btPostProcessingQuit;
GLUI_Button *btPreProcessing, *btPreProcessingQuit;
GLUI_Button *btBase, *btBaseQuit;
GLUI_Button *btSlope, *btSlopeQuit;
GLUI_Button *btWindow, *btWindowQuit;
GLUI_Button *btTop, *btTopQuit;

GLUI_Spinner *spBaseRound;
GLUI_Listbox *lbTestSlopeXY;
GLUI_Listbox *lbTestSlopeZ;
GLUI_Listbox* lbFiltreType;
GLUI_Spinner* spLissageSigma1;
GLUI_Spinner* spPoidsSelf;
GLUI_Spinner* spPoidsVertical;
GLUI_Listbox* lbLissageDirection;
GLUI_Spinner* spLissageSigma2;
GLUI_Spinner* spPoidsHorizontal;
GLUI_Spinner* spPoidsDiagonales;
GLUI_Listbox* lbSlopeType;

float fSlopeLength;


GLUI_FileBrowser *fbDomeFile;
GLUI_Listbox *lbDomeFile = NULL;

/********** Miscellaneous global variables **********/

//GLfloat light0_ambient[] =  {0.2f, 0.2f, 0.6f, 1.0f};
//GLfloat light0_diffuse[] =  {.6f, .6f, 1.0f, 1.0f};
//GLfloat light0_position[] = {.5f, .5f, 1.0f, 0.0f};

GLfloat light0_ambient[] = { 0.7f, 0.7f, 0.7f, 1.0f };
GLfloat light0_diffuse[] = { .6f, .6f, 0.6f, 1.0f };
//GLfloat light0_position[] = {.5f, 2.5f, 0.5f, 0.0f};
GLfloat light0_position[] = { 0.5f, 2.5f, 5.0f, 0.0f };

//GLfloat light1_ambient[] =  {0.1f, 0.1f, 0.3f, 1.0f};
//GLfloat light1_diffuse[] =  {.9f, .6f, 0.0f, 1.0f};
//GLfloat light1_position[] = {-1.0f, -1.0f, 1.0f, 0.0f};

GLfloat light1_ambient[] = { 0.7f, 0.3f, 0.9f, 1.0f };
GLfloat light1_diffuse[] = { .9f, .5f, 0.0f, 1.0f };
GLfloat light1_position[] = { -1.0f, 1.0f, 1.0f, 0.0f };

//GLfloat light2_ambient[] = { 1.0f, 1.0f, 1.0f, 1.0f };
//GLfloat light2_diffuse[] = { .9f, 1.0f, 1.0f, 1.0f };
//GLfloat light2_position[] = { -1.0f, 1.0f, -5.0f, 0.0f };

GLfloat light2_ambient[] = { 1.0f, 1.0f, 1.0f, 1.0f };
GLfloat light2_diffuse[] = { .9f, 1.0f, 1.0f, 1.0f };
GLfloat light2_position[] = { -3.0f, 1.0f, 5.0f, -5.0f };

PerlinNoise myNoise;

Material myMaterial[] = { {"Black Plastic", { (float)0.0f, (float)0.0f, (float)0.0f, (float)1.0f}, { (float)0.01f, (float)0.01f, (float)0.01f, (float)1.0f}, { (float)0.5f, (float)0.5f, (float)0.5f, (float)1.0f}, (float)32.0f},
						 {"Cyan Plastic", { (float)0.0f, (float)0.1f, (float)0.06f, (float)1.0f}, { (float)0.0f, (float)0.50980392f, (float)0.50980392f, (float)1.0f}, { (float)0.50196078f, (float)0.50196078f, (float)0.50196078f, (float)1.0f}, (float)32.0f},
						 {"Green Plastic", {0.0f, 0.0f, 0.0f, 1.0f}, {0.1f, 0.35f, 0.1f, 1.0f}, {0.45f, 0.55f, 0.45f, 1.0f}, 32.0f},
						 {"Red Plastic", {0.0f, 0.0f, 0.0f, 1.0f}, {0.5f, 0.0f, 0.0f, 1.0f}, {0.7f, 0.6f, 0.6f, 1.0f}, 32.0f},
						 {"White Plastic", {0.0f, 0.0f, 0.0f, 1.0f}, {0.55f, 0.55f, 0.55f, 1.0f}, {0.70f, 0.70f, 0.70f, 1.0f}, 32.0f},
						 {"Yellow Plastic", {0.0f, 0.0f, 0.0f, 1.0f}, {0.5f, 0.5f, 0.0f, 1.0f}, {0.60f, 0.60f, 0.50f, 1.0f}, 32.0f},
						 {"Black Rubber", {0.02f, 0.02f, 0.02f, 1.0f}, {0.01f, 0.01f, 0.01f, 1.0f}, {0.4f, 0.4f, 0.4f, 1.0f}, 10.0f},
						 {"Cyan Rubber", {0.0f, 0.05f, 0.05f, 1.0f}, {0.4f, 0.5f, 0.5f, 1.0f}, {0.04f, 0.7f, 0.7f, 1.0f}, 10.0f},
						 {"Green Rubber", {0.0f, 0.05f, 0.0f, 1.0f}, {0.4f, 0.5f, 0.4f, 1.0f}, {0.04f, 0.7f, 0.04f, 1.0f}, 10.0f},
						 {"Red Rubber", {0.05f, 0.0f, 0.0f, 1.0f}, {0.5f, 0.4f, 0.4f, 1.0f}, {0.7f, 0.04f, 0.04f, 1.0f}, 10.0f},
						 {"White Rubber", {0.05f, 0.05f, 0.05f, 1.0f}, {0.5f, 0.5f, 0.5f, 1.0f}, {0.7f, 0.7f, 0.7f, 1.0f}, 10.0f},
						 {"Yellow Rubber", {0.05f, 0.05f, 0.0f, 1.0f}, {0.5f, 0.5f, 0.4f, 1.0f}, {0.7f, 0.7f, 0.04f, 1.0f}, 10.0f},
						 {"Glass", {0.19225f, 0.19225f, 0.19225f, 0.6f}, {0.50754f, 0.50754f, 0.50754f, 0.6f}, {0.508273f, 0.508273f, 0.508273f, 0.6f}, 51.2f} ,
						 {"Brass", {0.329412f, 0.223529f, 0.027451f, 1.0f}, {0.780392f, 0.568627f, 0.113725f, 1.0f}, {0.992157f, 0.941176f, 0.807843f, 1.0f}, 27.8974f},
						 {"Bronze", {0.2125f, 0.1275f, 0.054f, 1.0f}, {0.714f, 0.4284f, 0.18144f, 1.0f}, {0.393548f, 0.271906f, 0.166721f, 1.0f}, 25.6f},
						 {"Chrome", {0.25f, 0.25f, 0.25f, 1.0f}, {0.4f, 0.4f, 0.4f, 1.0f}, {0.774597f, 0.774597f, 0.774597f, 1.0f}, 76.8f},
						 {"Copper", {0.19125f, 0.0735f, 0.0225f, 1.0f}, {0.7038f, 0.27048f, 0.0828f, 1.0f}, {0.256777f, 0.137622f, 0.086014f, 1.0f}, 12.8f},
						 {"Gold", {0.24725f, 0.1995f, 0.0745f, 1.0f}, {0.75164f, 0.60648f, 0.22648f, 1.0f}, {0.628281f, 0.555802f, 0.366065f, 1.0f}, 51.2f},
						 {"Pewter", {0.10588f, 0.058824f, 0.113725f, 1.0f}, {0.427451f, 0.470588f, 0.541176f, 1.0f}, {0.3333f, 0.3333f, 0.521569f, 1.0f}, 9.84615f},
						 {"Silver", {0.19225f, 0.19225f, 0.19225f, 1.0f}, {0.50754f, 0.50754f, 0.50754f, 1.0f}, {0.508273f, 0.508273f, 0.508273f, 1.0f}, 51.2f},
						 {"Polished Silver", {0.23125f, 0.23125f, 0.23125f, 1.0f}, {0.2775f, 0.2775f, 0.2775f, 1.0f}, {0.773911f, 0.773911f, 0.773911f, 1.0f}, 89.6f},
						 {"Obsidian", {0.05375f, 0.05f, 0.06625f, 1.0f}, {0.18275f, 0.17f, 0.22525f, 1.0f}, {0.332741f, 0.328634f, 0.346435f, 1.0f}, 38.4f},
						 {"Emerald", {0.0215f, 0.1745f, 0.0215f, 1.0f}, {0.07568f, 0.61424f, 0.07568f, 1.0f}, {0.633f, 0.727811f, 0.633f, 1.0f}, 76.8f},
						 {"Jade", {0.135f, 0.2225f, 0.1575f, 1.0f}, {0.54f, 0.89f, 0.63f, 1.0f}, {0.316228f, 0.316228f, 0.316228f, 1.0f}, 12.8f},
						 {"Pearl", {0.25f, 0.20725f, 0.20725f, 1.0f}, {1.0f, 0.829f, 0.829f, 1.0f}, {0.296648f, 0.296648f, 0.296648f, 1.0f}, 10.24f},
						 {"Turquoise", {0.1f, 0.18725f, 0.1745f, 1.0f}, {0.396f, 0.74151f, 0.63102f, 1.0f}, {0.297254f, 0.30829f, 0.306678f, 1.0f}, 12.8f},
						 {"Ruby", {0.1745f, 0.01175f, 0.01175f, 1.0f}, {0.61424f, 0.04136f, 0.04136f, 1.0f}, {0.727811f, 0.626959f, 0.626959f, 1.0f}, 76.8f},
						 {"Purple", {(GLfloat)0.2706f, (GLfloat)0.1922f, (GLfloat)0.2745, (GLfloat)1.0f}, {(GLfloat)0.6863, (GLfloat)0.4706, (GLfloat)0.6745, (GLfloat)1}, {(GLfloat)0.898, (GLfloat)0.898, (GLfloat)0.898, (GLfloat)1}, (GLfloat)14.0}
};

int   PerlinNoise_p[SAMPLE_SIZE + SAMPLE_SIZE + 2];
float PerlinNoise_g3[SAMPLE_SIZE + SAMPLE_SIZE + 2][3];
float PerlinNoise_g2[SAMPLE_SIZE + SAMPLE_SIZE + 2][2];
float PerlinNoise_g1[SAMPLE_SIZE + SAMPLE_SIZE + 2];
int   nButtonState;
int   nbBases;
int   nbSlopes;
int	  nbWindows;
int   nbTestSlopes;
int   nbFiltres;
int   nbTops;


double deg2rad = _PI / 180;


const char *BaseListe[] = { "Clover", "EpiCycloid", "HypoCycloid", "Shuriken", "Tablemat", "Circle",
					  "Convex Polygon", "Hook", "Star", "Arrow", "Bloody Cranesbill", "Lotus",
					  "Windmill", "Maltese Cross", "Basileus", "Truncated Star", "Saint Andrew's Cross",
					  "Star Polygon", "Saw Toooth", "Gear", "Snow Flake", "Roof", "Cat's ear", "Daisy", 
					  "Fourier-gon", "Savior"};

const char *SlopeListe[] = { "Kazack Hat", "Half Circle", "Triangle", "Russian Dome",
					  "Yurt", "Mongolian Tent", "Sabanci Dome", "Javanese Hat", "Cloche Hat",
					  "Saint Sophia Dome", "Persian Dome", "Tower", "Cylinder",
					  "Bun", "NeMatollah Vali Dome", "Sponge Cake", "Emperor Cap", "Staircase",
					  "Byzantine Dome", "Empress Hat", "Nasi Tumpeng", "Overhang", "Mushroom",
					  "Hut", "Rounded Hut", "Karen Hat", "Mary Magdalena Dome",
					  "Tri-Foiled Dome", "Ispahan Dome", "Tent", "Rim", "Tri-Cuspid", "Hindi", "BollyWood",
					  "Saint Pierre", "Cuve", "Pine Tree", "Casablanca"};

const char *WindowListe[] = { "M", "Cycloid", "Saw Tooth", "Sine", "Straight", "Hat", "Step", "None                         ", "Podium", "Staircase",
					   "Moscow", "Pointed Steps", "Cathedral", "Pointed Square", "Middle Age", "Oriental", "Curly", "TriFoil", "Square", "Mexican Hat"};

const char *ShrinkListe[] = { "Dummy", "Square Root", "Linear", "Square" };


//float fSlopeParameterXY;
bool  bDomeModifie;

Couleur stCouleur[MAX_COULEUR];
unsigned short CouleurBase[MAX_COULEUR];
int I_Sphere = 20;
int J_Sphere = 20;
int I_Cone = 20;
int J_Cone = 20;

Base BaseFunctionListe[] = {
	{BaseClover_X, BaseClover_Y, "Clover", 0.4f, 1},
	{BaseEpiCycloid_X, BaseEpiCycloid_Y, "EpiCycloid", 1, 1},
	{BaseHypoCycloid_X, BaseHypoCycloid_Y, "HypoCycloid", 1, 3},
	{BaseShuriken_X, BaseShuriken_Y, "Shuriken", 1, 2},
	{BaseTablemat_X, BaseTablemat_Y, "Tablemat", 80, 1},
	{BaseCircle_X, BaseCircle_Y, "Circle",1, 1},
	{BasePolygon_X, BasePolygon_Y, "Convex Polygon", 1, 3},
	{BaseHook_X, BaseHook_Y, "Hook", 1, 1},
	{BaseStar_X, BaseStar_Y, "Star", 1, 1},
	{BaseArrow_X, BaseArrow_Y, "Arrow", 1, 1},
	{BaseBloodyCranesbill_X, BaseBloodyCranesbill_Y, "Bloody Cranesbill", 1, 1},
	{BaseLotus_X, BaseLotus_Y, "Lotus", 1, 2},
	{BaseWindmill_X, BaseWindmill_Y, "Windmill", 1, 2},
	{BaseMalteseCross_X, BaseMalteseCross_Y, "Maltese Cross", 1, 2},
	{BaseBasileus_X, BaseBasileus_Y, "Basileus", 1, 1},
	{BaseTruncatedStar_X, BaseTruncatedStar_Y, "Truncated Star", 1, 2},
	{BaseSaintAndrew_X, BaseSaintAndrew_Y, "Saint Andrew's Cross", 1, 2},
	{BaseStarPolygon_X, BaseStarPolygon_Y, "Star Polygon", 0.5f, 2},
	{BaseSawToooth_X, BaseSawToooth_Y, "Saw Toooth", .5, 3},
	{BaseGear_X, BaseGear_Y, "Gear", 1, 3 },
	{BaseSnowFlake_X, BaseSnowFlake_Y, "Snow Flake", 4.0, 3},
	{BaseRoof_X, BaseRoof_Y, "Roof", 1, 1},
	{BaseCatsEar_X, BaseCatsEar_Y, "Cat's ear", 65.0, 3},
	{BaseDaisy_X, BaseDaisy_Y, "Daisy", 2.0, 1},
	{BaseMandala_X, BaseMandala_Y, "Mandala", 45, 1},
	{BaseRoissy_X, BaseRoissy_Y, "Roissy", 2, 2},
	{BaseBrokenClover_X, BaseBrokenClover_Y, "Broken Clover", 90, 1},
	{BasePropeller_X, BasePropeller_Y, "Propeller", 4, 1},
	{BaseBrokenCircle_X, BaseBrokenCircle_Y, "Broken Circle", 30, 1},
	{BaseSpiral_X, BaseSpiral_Y, "Spiral", 0.25, 1},
	{BasePointedPolygon_X, BasePointedPolygon_Y, "Pointed Polygon", 10, 5},
	{BaseFourier_X, BaseFourier_Y, "Fourier Polygon", 20, 3},
	{BaseSavior_X, BaseSavior_Y, "Savior", 2, 3},
	{BaseRoundedPolygon_X, BaseRoundedPolygon_Y, "Rounded Polygon",40, 3 }
};


Slope SlopeFunctionListe[] = {
	{SlopeKazackHatXY, SlopeKazackHatZ, "Kazack Hat", 8.0f, 7.0f},
	{SlopeSinus, SlopeCosinus, "Half Circle", 1, 2},
	{SlopeStraight, SlopeStraight, "Triangle", 1, 1},
	{SlopeIranian, SlopeStraight, "Russian Dome", 1, 1},
	{SlopeSinus, SlopeInflected, "Yurt", 1, 2 },
	{SlopeSinus, SlopeStraight, "Mongolian Tent", 1, 1},
	{SlopeCosinus, SlopeSinus, "Sabanci Dome", 2, 5},
	{SlopeStraight, SlopeSinus, "Javanese Hat", 2, 2},
	{SlopeCubic2, SlopeStraight, "Cloche Hat", 1, 2},
	{SlopeInflected, SlopeSinus, "Saint Sophia Dome", 3, 0.5f},
	{SlopeIranian, SlopeInflected, "Persian Dome", 1, 3},
	{SlopeSinus, SlopeFloor, "Tower", 1, 0},
	{SlopeCylinderXY, SlopeCylinderZ, "Cylinder", 0.5f, 0.5f},
	{SlopeRoman, SlopeCosinus, "Bun", 1, 2},
	{SlopeRoman, SlopeStraight, "NeMatollah Vali Dome", 1, 1},
	{SlopeRoman, SlopeSinus, "Sponge Cake", 1, 5},
	{SlopeRoman, SlopeInflected,"Emperor Cap", 0.5f, 3},
	{SlopeStraight, SlopeFloor, "Staircase", 1, 3},
	{SlopeInflected, SlopeStraight, "Byzantine Dome", 4, 0.5f},
	{SlopeStraight, SlopeInflected, "Empress Hat", 1, 4},
	{SlopeInflected, SlopeInflected, "Nasi Tumpeng", 2, 4},
	{SlopeExponential, SlopeInflected, "Overhang", 0.3f, 4},
	{SlopeExponential, SlopeStraight, "Mushroom", 0.3f, 1},
	{SlopeSqrtPower4, SlopeInflected, "Hut", 0.5f, 3},
	{SlopeSqrtPower4, SlopeRoundedPower4, "Rounded Hut", 1, 3},
	{SlopeEigthSquare, SlopeStraight, "Karen Hat", 0.125f, 1},
	{SlopeIranian3, SlopeStraight, "Mary Magdalena Dome", 1, 1},
	{SlopeStraight, SlopeTriFoil, "Tri-Foiled Dome", 1, 5},
	{SlopeIranian2, SlopeStraight, "Ispahan Dome", 1, 1},
	{SlopeStraight, SlopeTente, "Tent", 1, 5},
	{SlopeStraight, SlopeRim,"Rim", 1, 7},
	{SlopeStraight, SlopeTriCusp, "Tri-Cuspid", 1, 5},
	{SlopeHindi, SlopeStraight, "Hindi", 2, 1},
	{SlopeHindi, SlopeMongol, "BollyWood", 2, 1},
	{SlopeSaintPierre, SlopeStraight, "Saint Pierre", 3, 1},
	{SlopeCuve, SlopeStraight, "Cuve", 7.5f, 1},
	{SlopePineTreeXY, SlopePineTreeZ, "Pine Tree", 0, 0},
	{SlopeCasablanca, SlopeJoliCubic, "Casablanca", 1, 1},
	{SlopeAngleX, SlopeAngleY, "Obtuse", 8, 8},
	{SlopeMaison_XY, SlopeMaison_Z, "House", 1.0f, 1.33f},
	{SlopeKazackHatXY, SlopeRoman, "_Customized", 1.0f, 2.0f},
	{SlopeBulbXY, SlopeBulbZ, "Bulb", 3.0f, 5.0f},
	{SlopeMoroccoXY, SlopeMoroccoZ, "Morocco", 3.5f, 3.5f},
	{SlopeDamaoHatXY, SlopeDamaoHatZ, "Damao Hat", 3.0f, 0.13f},
	{SlopeDamaoHatXY, SlopeDamaoHatXY, "   Custom Slope", 3.0f, 3.0f}
};

int CUSTOM_SLOPE;

Window WindowFunctionListe[] = {
	{WindowM, "M Shape", 8, 1},
	{WindowCycloid, "Cycloid", 0, 1},
	{WindowSawTooth, "Saw Tooth", 8, 2},
	{WindowSinus, "Sine", 0, 3},
	{WindowStraight, "Straight", 1, 1},
	{WindowHat, "Hat", 1, 1},
	{WindowStep, "Step", 3, 0.75f},
	{WindowNone, "         None              ", 0, 1},
	{WindowPodium, "Podium", 1.5f, 1.7f},
	{WindowStaircase, "Staircase", 4, 1.5f},
	{WindowMoscow, "Moscow", 2, 1},
	{WindowPointedSteps, "Pointed Steps", 3, 1},
	{WindowCathedral, "Cathedral", 7, 0.4f},
	{WindowPointedSquare, "Pointed Square", 4, 0.8f},
	{WindowMiddleAge, "Middle Age", 5, 0.75f},
	{WindowOriental, "Oriental", 1, 0.2f},
	{WindowCurly, "Curly", 6, 0.8f},
	{WindowTriFoil, "TriFoil", 2, 1.5f},
	{WindowSquare, "Square", 8, 6},
	{WindowMexicanHat, "Mexican Hat", 5, 10},
	{WindowAngkorWat, "Angkor Wat", 0, 0.35f}
};


TestSlope TestSlopeFunctionListe[] =
{
	{ SlopeInflected,"Inflected"},
	{ SlopeRim,"Rim"},
	{ SlopeTente,"Tent"},
	{ SlopeTriFoil,"Tri-Foil"},
	{ SlopeCylinderXY,"Cylinder XY"},
	{ SlopeCylinderZ,"Cylinder Z"},
	{ SlopeCuve,"Cuve"},
	{ SlopeSaintPierre,"Saint-Pierre"},
	{ SlopeTriCusp,"Tri-Cusp"},
	{ SlopeStraight,"Straight"},
	{ SlopeMongol,"Mongol"},
	{ SlopeHindi,"Hindi"},
	{ SlopeAngle,"Angle"},
	{ SlopeAngleX,"Angle X"},
	{ SlopeAngleY,"Angle Y"},
	{ SlopeJoliCubic,"Nice Cubic"},
	{ SlopeCasablanca,"Casablanca"},
	{ SlopeSinus,"Sin"},
	{ SlopeSquareRootSinus,"Square Root Sine"},
	{ SlopeCosinus,"Cosine"},
	{ SlopeKazackHatXY,"Kazack Hat XY"},
	{ SlopeKazackHatZ,"Kazack Hat Z"},
	{ SlopePineTreeXY,"Pine Tree XY"},
	{ SlopePineTreeZ,"Pine Tree Z"},
	{ SlopeIranian,"Iranian"},
	{ SlopeFloor,"Floor"},
	{ SlopeStaircase,"Staircase"},
	{ SlopeRoman,"Roman"},
	{ SlopeFloorAside,"Floor-Aside"},
	{ SlopeExponential,"Exponential"},
	{ SlopePower4,"Power 4"},
	{ SlopeSqrtPower4,"Sqrt Power 4"},
	{ SlopeRoundedPower4,"Rounded Power 4"},
	{ SlopeEigthSquare,"Eigth Square"},
	{ SlopeIranian2,"Iranian 2"},
	{ SlopeIranian3,"Iranian 3"},
	{ SlopeCubic2,"Cubic"},
	{ SlopeMaison_XY,"House XY"},
	{ SlopeMaison_Z,"House Z"},
	{ SlopeMoroccoXY, "Morocco XY"},
	{ SlopeMoroccoZ, "Morocco Z"},
	{ SlopeBulbXY, "Bulb XY"},
	{ SlopeBulbZ, "Bulb Z"},
	{ SlopeDamaoHatXY, "Damao Hat XY"},
	{ SlopeDamaoHatZ, "Damao Hat Z"},
	{ SlopeVibrate, "Vibrate"}
};


const char* DistanceListe[] = {
			"Manhattan",
			"Euclidian",
			"3 Minkowski",
			"Chebyshev",
			"Discrete"
};

Filtre FiltreListe[] =
{
	{FILTER_Mean, Mean, "Mean", DISTANCE_EUCLIDE, 1, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, DIRECTION_NONE, S_NONE },
	{FILTER_Perso, Perso, "Personnal", DISTANCE_EUCLIDE, 1, 1, "Self", 0, "Horizontal", 0, "Vertical", 0, "Diagonal", 0, S_NONE, 0, S_NONE, DIRECTION_NONE, S_NONE},
	{FILTER_Gauss, Gauss, "Gauss" , DISTANCE_EUCLIDE, 2, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 1, "Sigma1", 0, S_NONE, DIRECTION_NONE, S_NONE},
	{FILTER_Binomial, Binomial, "Binomial" , DISTANCE_EUCLIDE, 0, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, DIRECTION_NONE, S_NONE},
	{FILTER_Pyramidal, Pyramidal, "Pyramidal" , DISTANCE_EUCLIDE, 0, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, DIRECTION_NONE, S_NONE},
	{FILTER_Circular, Circular, "Circular" , DISTANCE_EUCLIDE, 1, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 3, "Rayon", 0, S_NONE, DIRECTION_NONE, S_NONE},
	{FILTER_LoG, LoG, "LoG" , DISTANCE_EUCLIDE, 2, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0.5, "Sigma1", 0, S_NONE, DIRECTION_NONE, S_NONE},
	{FILTER_Exponential, Exponential, "Exponential" , DISTANCE_EUCLIDE, 1, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 1, "Gamma", 0, S_NONE, DIRECTION_NONE, S_NONE},
	{FILTER_Laplacian, Laplacian, "Laplacian" , DISTANCE_EUCLIDE, 1, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, DIRECTION_NONE, S_NONE},
	{FILTER_Sobel, Sobel, "Sobel" , DISTANCE_EUCLIDE, 1, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, DIRECTION_NONE, "Direction"},
	{FILTER_MDIF, MDIF, "MDIF" , DISTANCE_EUCLIDE, 2, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0.5, "Sigma1", 0, S_NONE, DIRECTION_NONE, "Direction"},
	{FILTER_Conic, Conic, "Conic" , DISTANCE_EUCLIDE, 1, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, DIRECTION_NONE, S_NONE},
	{FILTER_HighPass, HighPass, "HighPass" , DISTANCE_EUCLIDE, 1, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 1, "Gamma", 0, S_NONE, DIRECTION_NONE, S_NONE},
	{FILTER_Prewitt, Prewitt, "Prewitt" , DISTANCE_EUCLIDE, 0, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, DIRECTION_NONE, "Direction"},
	{FILTER_LowPass, LowPass, "LowPass" , DISTANCE_EUCLIDE, 1, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 1, "Alpha", 0, S_NONE, DIRECTION_NONE, S_NONE},
	{FILTER_DoG, DoG, "DoG" , DISTANCE_EUCLIDE, 2, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 3, "Sigma1", 1, "Sigma2", DIRECTION_NONE, S_NONE},
	{FILTER_Enhancer, Enhancer, "Enhancer", DISTANCE_EUCLIDE, 1, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, 0, S_NONE, DIRECTION_NONE, S_NONE }
};

Top TopListe[] =
{
	{"    None    ", 0, 0},
	{"Lantern", .5, 1},
	{"Arrow", .5, 2},
	{"Meringue", .5, .5},
	{"Foot", .5, 1}
};
