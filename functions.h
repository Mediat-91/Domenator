#include "Types.h"
//#include "Globales.h"
#include "Dome.h"

double mathfmod(double x, double y);


int PrintStrutsCone(GLpoint* Vertices, FILE* Dome, int nMaxI, int nMaxJ, int nbObjects, PRINT_TYPE Type);
int PrintStrutsSphere(GLpoint* Vertices, FILE* Dome, int nMaxI, int nMaxJ, int nbObjects, PRINT_TYPE Type);
int PrintSpike(GLpoint* Vertice, FILE* Dome, int nMaxI, int nMaxJ, int nbObjects, PRINT_TYPE Type);
int PrintCloture(FILE* Dome, int nMaxI, int nMaxJ, GLpoint* Poly, int nbObjects, PRINT_TYPE Type);

double BaseEpiCycloid_X(double, double, double);
double BaseEpiCycloid_Y(double, double, double);
double BaseHypoCycloid_X(double, double, double);
double BaseHypoCycloid_Y(double, double, double);
double BaseMalteseCross_X(double, double, double);
double BaseMalteseCross_Y(double, double, double);
double BaseTablemat_X(double, double, double);
double BaseTablemat_Y(double, double, double);
double BaseCircle_X(double, double, double);
double BaseCircle_Y(double, double, double);
double BaseStar_X(double, double, double);
double BaseStar_Y(double, double, double);
double BaseArrow_X(double, double, double);
double BaseArrow_Y(double, double, double);
double BaseLotus_X(double u, double a, double b);
double BaseLotus_Y(double u, double a, double b);
double BaseBasileus_X(double u, double a, double b);
double BaseBasileus_Y(double u, double a, double b);
double BaseSawToooth_X(double u, double a, double b);
double BaseSawToooth_Y(double u, double a, double b);
double BaseGear_X(double u, double a, double b);
double BaseGear_Y(double u, double a, double b);
double BaseSaintAndrew_X(double u, double a, double b);
double BaseSaintAndrew_Y(double u, double a, double b);
double BaseSnowFlake_X(double u, double a, double b);
double BaseSnowFlake_Y(double u, double a, double b);
double BaseDaisy_X(double u, double a, double b);
double BaseDaisy_Y(double u, double a, double b);
double BaseClover_X(double u, double a, double b);
double BaseClover_Y(double u, double a, double b);
double BasePolygon_X(double u, double a, double b);
double BasePolygon_Y(double u, double a, double b);
double BaseFourierRho(double v, double fParam, bool bCalcul);
double BaseRoundedPolygon_X(double u, double a, double b);
double BaseRoundedPolygon_Y(double u, double a, double b);

double BaseStarPolygon_X(double u, double a, double b);
double BaseStarPolygon_Y(double u, double a, double b);
double BaseBloodyCranesbill_X(double, double, double);
double BaseBloodyCranesbill_Y(double, double, double);
double BaseTruncatedStar_X(double u, double a, double b);
double BaseTruncatedStar_Y(double u, double a, double b);
double BaseShuriken_X(double u, double a, double b);
double BaseShuriken_Y(double u, double a, double b);
double BaseConvexPolygon_X(double, double, double);
double BaseConvexPolygon_Y(double, double, double);
double BaseHook_X(double, double, double);
double BaseHook_Y(double, double, double);
double BaseWindmill_X(double u, double a, double b);
double BaseWindmill_Y(double u, double a, double b);
double BaseCatsEar_X(double u, double a, double b);
double BaseCatsEar_Y(double u, double a, double b);
double BaseMandala_X(double u, double a, double b);
double BaseMandala_Y(double u, double a, double b);
double BaseRoissy_X(double u, double a, double b);
double BaseRoissy_Y(double u, double a, double b);
double BaseBrokenClover_X(double u, double a, double b);
double BaseBrokenClover_Y(double u, double a, double b);
double BasePropeller_X(double u, double a, double b);
double BasePropeller_Y(double u, double a, double b);
double BaseBrokenCircle_X(double u, double a, double b);
double BaseBrokenCircle_Y(double u, double a, double b);
double BaseSpiral_X(double u, double a, double b);
double BaseSpiral_Y(double u, double a, double b);
double BasePointedPolygon_X(double u, double a, double b);
double BasePointedPolygon_Y(double u, double a, double b);
double BaseFourier_X(double u, double a, double b);
double BaseFourier_Y(double u, double a, double b);
double BaseSavior_X(double u, double a, double b);
double BaseSavior_Y(double u, double a, double b);



double BaseRoof_X(double u, double a, double b);
double BaseRoof_Y(double u, double a, double b);

double SlopeFunction(double, double, double);

double SlopeSinus(double, double);
double SlopeCosinus(double, double);
double SlopeStraight(double, double);
double SlopeIranian(double, double);
double SlopeSquareRootSinus(double, double);
double SlopeStaircase(double, double);
double SlopeRoman(double, double);
double SlopeFloor(double, double);
double SlopeFloorAside(double, double);
double SlopeInflected(double, double);
double SlopeExponential(double, double);
//double SlopeSquareRoot4(double, double);
double SlopePower4(double, double);
double SlopeSqrtPower4(double, double);
double SlopeRoundedPower4(double, double);
double SlopeEigthSquare(double, double);
double SlopeIranian2(double, double);
double SlopeIranian3(double, double);
double SlopeCubic2(double, double);
double SlopeHindi(double, double);
double SlopeMongol(double, double);
double SlopeKazackHatXY(double v, double fParam);
double SlopeKazackHatZ(double v, double fParam);
double SlopeCylinderXY(double v, double fParam);
double SlopeCylinderZ(double v, double fParam);
double SlopeTriFoil(double v, double fParam);
double SlopeTente(double v, double fParam);
double SlopeRim(double v, double fParam);
double SlopeVibrate(double v, double fParam);
double SlopeTriCusp(double v, double fParam);
double SlopeSaintPierre(double v, double fParam);
double SlopeCuve(double v, double fParam);
double SlopePineTreeXY(double v, double fParam);
double SlopePineTreeZ(double v, double fParam);
double SlopeCasablanca(double v, double fParam);
double SlopeJoliCubic(double v, double fParam);
double SlopeAngleX(double v, double fParam);
double SlopeAngleY(double v, double fParam);
double SlopeMaison_XY(double v, double fParam);
double SlopeMaison_Z(double v, double fParam);
double SlopeAngle(double v, double fParam);
double SlopeBulbXY(double v, double fParam);
double SlopeBulbZ(double v, double fParam);
double SlopeMoroccoXY(double v, double fParam);
double SlopeMoroccoZ(double v, double fParam);
double SlopeDamaoHatZ(double v, double fParam);
double SlopeDamaoHatXY(double v, double fParam);

double WindowCycloid(double, double, double, double, bool bFrise);
double WindowSawTooth(double, double, double, double, bool bFrise);
double WindowSinus(double, double, double, double, bool bFrise);
double WindowStraight(double, double, double, double, bool bFrise);
double WindowHat(double, double, double, double, bool bFrise);
double WindowStep(double, double, double, double, bool bFrise);
double WindowSquare(double, double, double, double, bool bFrise);
double WindowNone(double, double, double, double, bool bFrise);
double WindowPodium(double, double, double, double, bool bFrise);
double WindowStaircase(double, double, double, double, bool bFrise);
double WindowMoscow(double, double, double, double, bool bFrise);
double WindowPointedSteps(double u, double a, double d, double p, bool bFrise);
double WindowCathedral(double u, double a, double d, double p, bool bFrise);
double WindowOriental(double u, double a, double d, double p, bool bFrise);
double WindowPointedSquare(double u, double a, double d, double p, bool bFrise);
double WindowMiddleAge(double u, double a, double d, double p, bool bFrise);
double WindowCurly(double u, double a, double d, double p, bool bFrise);
double WindowM(double u, double a, double d, double p, bool bFrise);
double WindowTriFoil(double u, double a, double d, double p, bool bFrise);
double WindowMexicanHat(double u, double a, double d, double p, bool bFrise);
double WindowAngkorWat(double u, double a, double d, double p, bool bFrise);


void Lissage(GLpoint *Vertices, int nMaxI, int nMaxJ);
void Evolute(GLpoint *Vertices, int nMaxI, int nMaxJ);
void DifferentialGeometry(GLpoint* Vertices, int nMaxI, int nMaxJ);
void Enhance(GLpoint* Vertices, int nMaxI, int nMaxJ);
void NormaleSharpening(GLpoint *Vertices, int nMaxI, int nMaxJ);
int PrintOBJCone(GLpoint Centre1, GLpoint Centre2, double Rayon1, double Rayon2, FILE *DomeOBJ, int Start, const char *sTitre, char *sMaterial, bool bStraight = true);
int PrintSpirale(FILE *Dome, int nbObjects, PRINT_TYPE Type);
int PrintSnake(FILE* Dome, int nbObjects, PRINT_TYPE Type);
int PrintComb(FILE* Dome, int nbObjects, PRINT_TYPE Type);
int PrintFrise(FILE *Dome, int nbObjects, PRINT_TYPE Type);
int PrintSTLConeBinary(GLpoint Centre1, GLpoint Centre2, float Rayon1, float Rayon2, unsigned short nCouleur, FILE *DomeSTL, bool bStraight = true);
int PrintOBJSphere(GLpoint Centre, double Rayon, FILE *DomeOBJ, int Start, const char *sTitre, char * sMaterial);

int CalculSpirale();
int CalculFrise();
int CalculSnake();
void ChargeFrise(int nFrise);
//void PrintPOVSnake(FILE *DomeINC);
void PrintPOVCone(GLpoint Centre0, GLpoint Centre1, float fRadius0, float fRadius1, FILE *DomeINC, char *sMaterial);
//double ShrinkCoefficient(int j, int MaxJ);
bool PrintSTLTriangle(FILE *Dome, GLpoint P1, GLpoint P2, GLpoint P3, int nCouleur);

int main(int argc, char* argv[]);

void ChargeSlope(int nSlope);

void  ChargeWindow(int nWindow);
void  Calcul(double u, double v, double *dX, double *dY, double *dZ, bool bDerivee = true);
void  ChargeMesh();
void  InitParameters();
void  DefaultSlopeParameters();
void  ProcessMessage(int Control);
void  Message(const char *sMessage, int Type = 0);
void  CrossProduct(GLpoint *V, GLpoint *W, GLpoint *Resultat);
void CrossProduct(Vector3D* V, Vector3D* W, Vector3D* Resultat);
float Colineaire(GLpoint O, GLpoint M, GLpoint N);
float  ComputeSlopeLength();
float ComputeWindowCoefficient(float u, float a, float d);
void  LoadDome();
void  SaveDome();
void  CalculMesh(int nMaxI, int nMaxJ, GLpoint *Vertices);
void  Normales(int i, int j, int nMaxI, int nMaxJ, GLpoint *Vertices, GLpoint *Normals);
void  CalculNormales(int nMaxI, int nMaxJ, GLpoint *Vertices, GLpoint *Normals);
//void  PovCone(FILE *INCFile, GLpoint *VERTICES, int i, int j, int nType, float fRadius, char *sMaterial, int nMax);
int   PrintPOV();
int   PrintSTL();
int   PrintOBJ();
void  ChargeShrink(int nShrink);
void  AdjustSphere();
int   fill_lbfile(GLUI_Listbox *lb, char *sSearch, char *sFile);
GLpoint SphereFormula(GLpoint Centre, double Rayon, double u, double v);
GLpoint ConeFormula(GLpoint Centre1, GLpoint Centre2, double Rayon1, double Rayon2, double u, double v, bool Straight = true);
void ChargeTestSlope(int nSlopeXY, int nSlopeZ);
void ControlDisplay(int Control);
void Bumps(int nMaxI, int nMaxJ, GLpoint* Vertices, bool bDraw);
void Marquees(int nMaxI, int nMaxJ, GLpoint* Vertices, bool bDraw);
void Coord(int i, int j, int a, int b, int nMaxI, int nMaxJ, int* n_I, int* n_J);
float Poids(int a, int b, float p0, float p1, float p2, float p3, int puissance);
void Normalize(Vector3D* V);

void DrawBase();
void DrawSlope();
void DrawWindow();

void myGlutKeyboard(unsigned char Key, int x, int y);
void myGlutMenu(int value);
void myGlutIdle(void);
void myGlutMouse(int button, int button_state, int x, int y);
void myGlutMotion(int x, int y);
void myGlutReshape(int x, int y);
void myGlutDisplay(void);

void draw_axes(float scale);
double CalculDerivee(GLpoint *Vertices, int i, int j, Vector3D *V, int nMaxI, int nMaxJ);
void ChargeListeCouleur();
void InitCouleurs();
unsigned short GetColor(char *sColor);

void  PerlinNoiseNormalize2(float v[2]);
void  PerlinNoiseInit(void);
float PerlinNoise2(float vec[2]);
float PerlinNoise2D(float vec[2]);
float PerlinNoiseGet2D(float x, float y);

int PrintSTLSphereBinary(GLpoint Centre, float Rayon, unsigned short nCouleur, FILE *DomeSTL);
int PrintSTLMesh(FILE *DomeSTL, GLpoint *myVERTICES, int nmi, int nMaxJ);
float DistancePlan(GLpoint P1, GLpoint P2, GLpoint P3, GLpoint Point);

int PrintOBJMesh(FILE *Dome, GLpoint *Vertice, int Start, int nMaxI, int nMaxJ);
void ChargeShrink(int nShrink);

int PrintPOVMesh(GLpoint* Vertice, GLpoint* Normale, FILE* DomeINC, int nMaxI, int nMaxJ, char* DomeName);
void PrintPOVSphere(GLpoint Centre, float fRadius, FILE* DomeINC, char* sMaterial);
void PrintPOVRoundedCone(GLpoint Centre0, GLpoint Centre1, float fRadius0, float fRadius1, FILE* DomeINC, char* sMaterial);

void GetFunctions(int nSlope);
unsigned short STL_RGB(unsigned short red, unsigned short green, unsigned short blue);

char* strrep(char* a);

double ShrinkNoneFunction(double u, double fShrinking, double fSpeed);
double ShrinkLinearFunction(double u, double fShrinking, double fSpeed);
double ShrinkSquareRootFunction(double u, double fShrinking, double fSpeed);
double ShrinkSquareFunction(double u, double fShrinking, double fSpeed);



float ShrinkStruts(int j, int MaxJ);
double ShrinkSpirale(int j, int MaxJ);
double ShrinkSnake(int j, int MaxJ);
double ShrinkComb(int j, int MaxJ);
double ShrinkFrise(int j, int MaxJ);

void DrawFrise();
void DrawSnake();
void DrawComb();
void DrawSpirale();

int CalculComb();
int CalculSnake();
int CalculFrise();
int CalculSpirale();

void FilterInput(TYPE_FILTER nFiltre, bool bDefault);
unsigned long PrintColor(FILE* DomeSTL);
float myMod(float x, float y);
double myMod(double x, double y);
int myMod(int x, int y);

float Perso(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float Gauss(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float Mean(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float Binomial(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float Pyramidal(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float Circular(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float LoG(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float Exponential(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float Sobel(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float Laplacian(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float MDIF(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float Conic(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float HighPass(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float Prewitt(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float LowPass(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float DoG(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);
float Enhancer(float* Matrix, int nCote, DISTANCE_Type nDistance, float fPower, FilterParameter* Parameters);

void ControlPostProcessing(int Control);

double BaseFunction(double u, double a, double b, double dAxe, double w);
