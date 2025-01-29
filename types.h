#pragma once

#include <GL\glui.h>

typedef double(*Parametric_1) (double);
typedef double(*Parametric_2) (double, double);
typedef double(*Parametric_3) (double, double, double);
typedef double(*Parametric_4) (double, double, double, double, bool);

typedef struct _DomeParameters
{
	int    nDomeVersion;
	char   sDomeName[400];
	float  fDomeScaleZ;
	char   sDomeMaterial[400];


	int    nBase;
	int    nBaseSides;
	float  fBaseParameter;
	float  fBaseSmoothing;

	int    nSlope;
	int    nSlopeTwist;
	int    nSlopeEvenly;

	int    nWindow;
	float  fWindowHeight;
	float  fWindowPerSide;
	float  fWindowShift;
	int    bWindowInverse;
	float  fWindowDelay;
	float  fWindowAttack;

	int    bFaceGenerate;
	int    bFaceSmooth;
	char   sFaceCenterMaterial[400];

	int    bWireGenerate;
	float  fWireShrinking;
	int    nWireShrinkingType;
	float  fWireShrinkingSpeed;

	int    nMeridianMain;
	float  fMeridianMainRadius;
	char   sMeridianMainMaterial[400];

	int    nMeridianMinor;
	float  fMeridianMinorRadius;
	char   sMeridianMinorMaterial[400];

	int    nParallelMain;
	float  fParallelMainRadius;
	char   sParallelMainMaterial[400];

	int    nParallelMinor;
	float  fParallelMinorRadius;
	char   sParallelMinorMaterial[400];

	float  fDiagonalURRadius;
	char   sDiagonalURMaterial[400];

	float  fDiagonalDRRadius;
	char   sDiagonalDRMaterial[400];

	float  fSphereMPRadius;
	char   sSphereMPMaterial[400];

	float  fSphereMpRadius;
	char   sSphereMpMaterial[400];

	float  fSpheremPRadius;
	char   sSpheremPMaterial[400];

	float  fSpherempRadius;
	char   sSpherempMaterial[400];

	int    nNoiseOctave;
	float  fNoiseDensity;
	int    nNoiseSeed;
	float  fNoiseFrequencyU;
	float  fNoiseFrequencyV;
	float  fNoiseAmplitudeR;
	float  fNoiseAmplitudeZ;
	float  fNoiseProgressiveR;
	float  fNoiseProgressiveZ;

	float  fWindowDamping;
	float  fBaseExpand;
	float  fSlopeParameterZ;
	float  fBaseDamping;      // printpov
	float  fSlopeParameterXY; //Attention load and save
	float  fWindowPower;
	int    nDomeEnd;
	float  fDomeSocle;

	int    bLissage;
	float  fLissagePoidsSelf;
	float  fLissagePoidsHorizontal;
	float  fLissagePoidsDiagonales;
	int    nLissageNombre;
	int    nLissageTaille;
	float  fLissagePuissance;
	int    nLissageDistanceType;

	float  fDiffGeomLambda;
	int    nDiffGeomNombre;
	int    bDiffGeomMeridien;
	int    bDiffGeomParallel;

	int    bNormale;
	float  fNormaleLambda;
	int    nNormaleNombre;

	int   nSpirale;
	int   nSpiraleShrink;
	float fSpiraleRadius;
	char  sSpiraleMaterial[400];
	float fSpiraleTour;
	float fSpiraleSpeed;
	float fSpiraleDelay;
	float fSpiraleDecalage;
	int   nSpiraleChirale;

	int   nFrise;
	float fFriseHeight;
	float fFrisePerSide;
	float fFriseDepart;
	float fFriseShift;
	float fFriseRadius;
	int   nFriseShrink;
	int   nFriseFunction;
	char  sFriseMaterial[400];

	int   nSnake;
	float fSnakeTour;
	int   nSnakeEvolve;
	float fSnakeShift;
	float fSnakeRadius;
	int   nSnakeShrink;
	char  sSnakeMaterial[400];
	int   nSnakeChirale;
	int   nSnakeStep;
	int   nSnakeToward;
	float fSnakeOffset;

	char  sFaceEdgeMaterial[400];
	int   nMeridianEdgeSize;
	int   nParallelEdgeSize;

	int   nSpiraleStep;
	int   nSpiraleToward;

	int   nFriseToward;
	int   nFriseStep;
	float fFrisePower;

	float fDerivative;
	int   bDerivativeNorme;

	float fFriseOffset;

	float fBumpDamping;
	float fBumpPower;
	int   nBumpNormalize;
	int   nBumpNormal;
	float fBumpLambda;
	int   bBump;
	float fBumpDistance;
	int   nBumpDistanceType;
	int   nBumpDistanceInverse;

	float fSpikeMP;
	float fSpikeMp;
	float fSpikemP;
	float fSpikemp;

	int   bLid;
	char  sLidMaterial[400];

	float fTopPart;
	int   bSlopeInverse;
	float fBaseRound;

	int   nBaseSecondary;
	float fBaseBaryCentre;

	int	  nCrenelPerSide;
	float fCrenelFront;
	float fCrenelThick;
	float fCrenelWide;
	float fCrenelHeight;
	float fCrenelPower;
	float fWindowAwning;
	float fWindowParameter;
	float fSlopeDelay;
	float fBaseMorphing;
	float fBaseSecondaryRotation;
	float fBaseSecondaryExpansion;
	int   bBaseMorphingDivision;
	int   bWindowBlind;
	float fFriseParameter;
	float fFriseGap;
	int   nTestSlopeXY;
	int   nTestSlopeZ;
	float fBaseSpiral;
	int   nSnakeSlope;
	float fSnakeSlopeParameter;
	int   nSpiraleSlope;
	float fSpiraleSlopeParameter;

	float  fLissagePoidsVertical;
	int nLissageTypeFiltre;
	int nLissageDirection;
	float fLissageSigma1;
	float fLissageSigma2;

	int   nTopType;
	float fTopParameter;

	float fDomeInFolding;
	int   nBumpType;

	int nPostMeridianEvery;
	int nPostMeridianShift;
	int nPostParallelEvolve;
	int nPostParallelFirst;
	int nPostParallelEvery;

	float fCrenelAngle;
	float fCrenelDelay;

	int nBaseEvenly;

	float fSpikeAngle;

	int bFriseInverse;
	int bFriseReverse;

	int   bMarquee;
	int   nMarqueeFunction;
	float fMarqueeParameter;
	float fMarqueeLambda;
	float fMarqueeDamping;
	float fMarqueeAlpha;
	float fMarqueeBeta;
	int   bMarqueeNormalize;
	float fMarqueeTop;
	int   bPostMeridianPerSide;
	int   nDifferentialGeometry;

	int   nEnhanceType;
	int   bEnhanceParallel;
	int   bEnhanceMeridian;
	int   nEnhanceMerge;
	int   nEnhanceCalcul;
	float fEnhanceLambda;
	float fEnhanceDamping;

	float fWindowLower;

	int   nComb;
	float fCombRadius;
	int   nCombShrink;
	char  sCombMaterial[400];
	float fCombSpacing;
	int   nCombStep;
	int   nCombToward;

	float fFriseSpacing;
	int   nCombType;
	int   nCombThread;
	int   nCombThreadGap;
	float fCombOffset;
}DomeParameters;

typedef struct _Base
{
	Parametric_3 FunctionX;
	Parametric_3 FunctionY;
	char         sNom[50];
	float        fBaseParameter;
	int          iMinSide;
}Base;
 

typedef struct _Slope
{
	Parametric_2 FunctionXY;
	Parametric_2 FunctionZ;
	char         sNom[50];
	float        fSlopeParameterXY;
	float        fSlopeParameterZ;
}Slope;

typedef struct _Window
{
	Parametric_4 WindowFunctionZ;
	char         sNom[50];
	float        fParameter;
	float        fPower;
}Window;

typedef struct _Top
{
	char		sNom[50];
	float		fPart;
	float		fParameter;
}Top;

typedef struct _TestSlope
{
	Parametric_2 Function;
	char		 sNom[50];
}TestSlope;


typedef struct _DomeDisplay
{
	int   nDisplayType;
	int   bDisplayMinor;
	float fScaleDisplay;
	int   nLight_1;
	int   nLight_2;
	float fSpeedCoefficient;
	int   nMaterial;
	int   bWireFrame;
	int   bFlatFace;
	int   bSmooth;
	int   nDisplay;
}DomeDisplay;


typedef struct _Tri
{
	int Index;
	char* Mot;
}Tri;

typedef struct _Vector3D
{
	double x;
	double y;
	double z;
}Vector3D;

struct GLpoint {
	GLfloat x, y, z;
};

typedef struct _STLHeader
{
	char sHeader[80];
	unsigned long  nbTriangles;
}STLHeader;

#pragma pack(2)
typedef struct _STLBinary
{
	GLpoint Normal;
	GLpoint Vertex1;
	GLpoint Vertex2;
	GLpoint Vertex3;
	unsigned short  Attribut;
}STLBinary;


typedef struct _PerlinNoise
{
	int   nOctaves;
	float fFrequency;
	float fAmplitude;
	int   nSeed;
	bool  bStart;
}PerlinNoise;

typedef struct _Material
{
	char  MaterialName[20];
	GLfloat MaterialAmbient[4];
	GLfloat MaterialDiffuse[4];
	GLfloat MaterialSpecular[4];
	GLfloat MaterialShininess;
}Material;

typedef struct _Tube
{
	GLpoint Centre;
	float   Radius;
}Tube;

typedef struct _GLFace
{
	unsigned int v1, v2, v3;
}GLface;

typedef struct _Couleur
{
	char sCouleur[400];
	bool bUsed;
	unsigned short nCouleur;
}Couleur;

enum TopType
{
	TOP_NONE = 0,
	TOP_LANTERN,
	TOP_ARROW,
	TOP_MERINGUE,
	TOP_FOOT
};
enum DISTANCE_Type
{
	DISTANCE_MANHATTAN = 0,
	DISTANCE_EUCLIDE,
	DISTANCE_3,
	DISTANCE_ULTRA,
	DISTANCE_DISCRETE
};
enum FILTER_Sens
{
	SENS_Vertical = 0,
	SENS_Horizontal
};
enum DIRECTION_Type
{
	DIRECTION_NONE = -1,
	DIRECTION_NORTH = SENS_Vertical,
	DIRECTION_WEST,
	DIRECTION_SOUTH,
	DIRECTION_EAST,
	DIRECTION_NORTH_WEST,
	DIRECTION_NORTH_EAST,
	DIRECTION_SOUTH_WEST,
	DIRECTION_SOUTH_EAST
};
enum ENHANCE_Type
{
	ENHANCE_NONE = 0,
	ENHANCE_RADIAL,
	ENHANCE_VERTICAL,
	ENHANCE_CYLINDRICAL,
	ENHANCE_OB_OA,
	ENHANCE_DB_AC,
	ENHANCE_D_B_O,
	ENHANCE_A_C_O,
	ENHANCE_ABCD_O,
	ENHANCE_OA_OB_OA_OD
};

enum ENHANCE_Merge
{
	MERGE_ZERO = 0,
	MERGE_UN,
	MERGE_SOMME
};
enum ENHANCE_Calcul
{
	CALCUL_Somme = 0,
	CALCUL_Produit
};
enum SLOPE_Function
{
	SLOPE_STRAIGHT = 9
};

enum COMB_Type
{
	COMB_PARALLEL_EDGE = 0,
	COMB_UPRIGHT_BASE
};

typedef struct _FilterParameter
{
	float fSelf;
	float fHorizontal;
	float fVertical;
	float fDiagonal;
	float fSigma1;
	float fSigma2;
	DIRECTION_Type nDirection;
}FilterParameter;

 
typedef float(*FilterFunction)(float *, int, DISTANCE_Type, float, FilterParameter*); 

enum TYPE_FILTER
{
	FILTER_Mean = 0,
	FILTER_Perso,
	FILTER_Gauss,
	FILTER_Binomial,
	FILTER_Pyramidal,
	FILTER_Circular,
	FILTER_LoG,
	FILTER_Exponential,
	FILTER_Laplacian,
	FILTER_Sobel,
	FILTER_MDIF,
	FILTER_Conic,
	FILTER_HighPass,
	FILTER_Prewitt,
	FILTER_LowPass,
	FILTER_DoG,
	FILTER_Enhancer
};


typedef struct _Filtre
{
	TYPE_FILTER nType;
	FilterFunction Function;
	char sNom[50];
	DISTANCE_Type nDistance;
	float fPower;
	float fSelf;
	char sSelf[20];
	float fHorizontal;
	char sHorizontal[20];
	float fVertical;
	char sVertical[20];
	float fDiagonal;
	char sDiagonal[20];
	float fSigma1;
	char sSigma1[20];
	float fSigma2;
	char sSigma2[20];
	DIRECTION_Type nDirection;
	char sDirection[20];
}Filtre;




/////////////////////////           ENUM    //////////////////////////////

enum BASE_Type
{
	BASE_CLOVER = 0, 
	BASE_EPICYCLOID,
	BASE_HYPOCYCLOID,
	BASE_SHURIKEN,
	BASE_TABLE_MAT,
	BASE_CIRCLE,
	BASE_CONVEX_POLYGON,
	BASE_HOOK,
	BASE_STAR,
	BASE_ARROW,
	BASE_BLOODY_CRANESBILL,
	BASE_LOTUS,
	BASE_WINDMILL,
	BASE_MALTESE_CROSS,
	BASE_BASILEUS,
	BASE_TRUNCATED_STAR,
	BASE_ST_ANDREW,
	BASE_STAR_POLYGON,
	BASE_SAW_TOOTH,
	BASE_GEAR,
	BASE_SNOW_FLAKE,
	BASE_ROOF,
	BASE_CATSEAR,
	BASE_DAISY,
	BASE_MANDALA,
	BASE_ROISSY, 
	BASE_BROKEN_CLOVER,
	BASE_PROPELLER,
	BASE_BROKEN_CIRCLE,
	BASE_SPIRAL,
	BASE_POINTED_POLYGON,
	BASE_FOURIER,
	BASE_SAVIOR,
	BASE_ROUNDED_POLYGON

};

enum SLOPE_Type
{
	SLOPE_KAZACK_HAT = 0,
	SLOPE_HALF_CIRCLE,
	SLOPE_TRIANGLE,
	SLOPE_RUSSIAN_DOME,
	SLOPE_YURT,
	SLOPE_MONGOLIAN_TENT,
	SLOPE_SABANCI_DOME,
	SLOPE_JAVANESE_HAT,
	SLOPE_CLOCHE_HAT,
	SLOPE_SAINT_SOPHIA_DOME,
	SLOPE_PERSIAN_DOME,
	SLOPE_TOWER,
	SLOPE_CYLINDER,
	SLOPE_BUN,
	SLOPE_NEMATOLLAH_VALI_DOME,
	SLOPE_SPONGE_CAKE,
	SLOPE_EMPEROR_CAP,
	SLOPE_STAIRCASE,
	SLOPE_BYZANTINE_DOME,
	SLOPE_EMPRESS_HAT,
	SLOPE_NASI_TUMPEG,
	SLOPE_OVERHANG,
	SLOPE_MUSHROOM,
	SLOPE_HUT,
	SLOPE_ROUNDED_HUT,
	SLOPE_KAREN_HAT,
	SLOPE_MARY_MAGDALENA_DOME,
	SLOPE_TRIFOIL,
	SLOPE_ISPAHAN_DOME,
	SLOPE_TENTE,
	SLOPE_RIM,
	SLOPE_TRICUSP,
	SLOPE_HINDI,
	SLOPE_BOLLYWOOD,
	SLOPE_SAINT_PIERRE,
	SLOPE_CUVE,
	SLOPE_PINE_TREE,
	SLOPE_CASABLANCA
	

};

enum WINDOW_Type
{
	WINDOW_M =0,
	WINDOW_CYCLOID,
	WINDOW_SAW_TOOTH,
	WINDOW_SINE,
	WINDOW_STRAIGHT,
	WINDOW_HAT,
	WINDOW_STEP,
	WINDOW_NONE,
	WINDOW_PODIUM,
	WINDOW_STAIRCASE,
	WINDOW_MOSCOW,
	WINDOW_POINTED_STEPS,
	WINDOW_CATHEDRAL,
	WINDOW_POINTED_SQUARE,
	WINDOW_MIDDLE_AGE,
	WINDOW_ORIENTAL,
	WINDOW_CURLY,
	WINDOW_TRIFOIL,
	WINDOW_SQUARE
};

enum SHRINK_Type
{
	SHRINK_NONE = 0,
	SHRINK_SQUARE_ROOT ,
	SHRINK_LINEAR,
	SHRINK_SQUARE,
};

enum SLOPE_Axe
{
	SLOPE_XY = 0,
	SLOPE_Z
};

enum BASE_Axe
{
	BASE_X = 0,
	BASE_Y
};

enum SPIRALE_CHIRALITE
{
	CHIRAL_LEFT = 0,
	CHIRAL_RIGHT,
	CHIRAL_BOTH
};

enum BUMP_Normal
{
	BUMP_1VECTEUR = 0,
	BUMP_2VECTEUR,
	BUMP_MEAN
};

enum BUMP_Type
{
	BUMP_FUNNEL = 0,
	BUMP_WAVE,
	BUMP_FOURPEAKS,
	BUMP_ASTROID,
	BUMP_ELLIPSOID,
	BUMP_TRACTOID,
	BUMP_CONE
};

enum SPHERE_Type
{
	SPHERE_TOP = 0,
	SPHERE_MP,
	SPHERE_Mp,
	SPHERE_mP,
	SPHERE_mp
};

enum TOWARD
{
	TOWARD_BOTTOM = 0,
	TOWARD_TOP
};

enum MATERIAL_Type
{
	MAT_BLACK_PLASTIC = 0,
	MAT_CYAN_PLASTIC,
	MAT_GREEN_PLASTIC,
	MAT_RED_PLASTIC,
	MAT_WHITE_PLASTIC,
	MAT_YELLOW_PLASTIC,
	MAT_BLACK_RUBBER,
	MAT_CYAN_RUBBER,
	MAT_GREEN_RUBBER,
	MAT_RED_RUBBER,
	MAT_WHITE_RUBBER,
	MAT_YELLOW_RUBBER,
	MAT_GLASS,
	MAT_BRASS,
	MAT_BRONZE,
	MAT_CHROME,
	MAT_COPPER,
	MAT_GOLD,
	MAT_PEWTER,
	MAT_SILVER,
	MAT_POLISHED_SILVER,
	MAT_OBSIDIAN,
	MAT_EMERALD,
	MAT_JADE,
	MAT_PEARL,
	MAT_TURQUOISE,
	MAT_RUBY,
	MAT_PURPLE,
	MAT_DERNIER
};

enum PRINT_TYPE
{
	PRINT_POV = 0,
	PRINT_OBJ,
	PRINT_STL_ASCII,
	PRINT_STL_BINARY
};

enum MATRIX_PRINT
{
	TYPE_Brut = 0,
	TYPE_Norme,
	TYPE_Integer,
	TYPE_Round,
	TYPE_Standard
};

enum DISPLAY_TYPE
{
	DISPLAY_Dome = 0,
	DISPLAY_Base,
	DISPLAY_Slope,
	DISPLAY_Window,
	DISPLAY_Snake,
	DISPLAY_Spirale,
	DISPLAY_Frise,
	DISPLAY_Comb
};

enum DIFF_GEOM
{
	DIFF_GEOM_EVOLUTE = 0,
	DIFF_GEOM_RADIAL,
	DIFF_GEOM_PEDAL
};
