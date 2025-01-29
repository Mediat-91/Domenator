#include <stdio.h>
#include "stdafx.h"
#include <math.h>
#include "Types.h"

//#include "Functions.h"
//#include "Globales.h"
#include "dome.h"


extern DomeParameters myDome;
extern double deg2rad;
extern Parametric_3 ShrinkCoefficientFunction;
extern int nbTestSlopes;
extern Slope SlopeFunctionListe[];
extern TestSlope TestSlopeFunctionListe[];
extern GLUI_Listbox* lbTestSlopeXY;
extern GLUI_Listbox* lbTestSlopeZ;
extern double myMod(double x, double y);

double mathfmod(double x, double y)
{
	double h = x - floor(x / y) * y;
	return	x - floor(x / y) * y;
}
char* strrep(char* a)
{
	char* b = a;
	while (*a != '\0')
	{
		if (*a == ' ')
			*a = '_';
		a++;
	}
	return b;
}
float ShrinkStruts(int j, int MaxJ)
{
	return (float)abs(ShrinkCoefficientFunction((float)j / (float)(MaxJ - 1), myDome.fWireShrinking, myDome.fWireShrinkingSpeed));
}
double ShrinkSpirale(int j, int MaxJ)
{
	int nSens = (myDome.nSpiraleToward == TOWARD_TOP ? 1 : -1);
	return abs(ShrinkCoefficientFunction((float)j / (float)(MaxJ - 1), 0, nSens));
}
double ShrinkSnake(int j, int MaxJ)
{
	return abs(ShrinkCoefficientFunction((float)j / (float)(MaxJ - 1), 0, (myDome.nSnakeToward == TOWARD_TOP ? -1 : 1)));
}
double ShrinkFrise(int j, int MaxJ)
{
	return abs(ShrinkCoefficientFunction(abs(sin(myDome.nBaseSides * myDome.fFrisePerSide * _PI * (float)j / (float)(MaxJ - 1))), 0.1, (myDome.nFriseToward == TOWARD_TOP ? -1 : 1)));
}
double ShrinkComb(int j, int MaxJ)
{
	return abs(ShrinkCoefficientFunction((float)j / (float)(MaxJ - 1), 0, (myDome.nCombToward == TOWARD_TOP ? 1 : -1)));
}


/*Fonction de Base*/
double BaseLotus_X(double u, double a, double b)
{
	return (cos(a * u) + (2 * b)) * (((2 * a - 1) * cos(u) + cos(u * (2 * a - 1))) / (2 * a * (1 + (2 * b))));
}
double BaseLotus_Y(double u, double a, double b)
{
	return (cos(a * u) + (2 * b)) * (((2 * a - 1) * sin(u) - sin(u * (2 * a - 1))) / (2 * a * (1 + (2 * b))));
}
double BaseLotus(double u, double a, double b, BASE_Axe Axe)
{
	double x = (cos(a * u) + (2 * b)) / (2 * a * (1 + (2 * b)));

	return Axe == BASE_X ? x * ((2 * a - 1) * cos(u) + cos(u * (2 * a - 1))) : x * ((2 * a - 1) * sin(u) - sin(u * (2 * a - 1)));
}
double BaseEpiCycloid_X(double u, double a, double b)
{
	double X = (((a + 1) * cos(u + _PI / a) - b * cos((u + _PI / a) * (a + 1))) / (a + 1 - b));
	double Y = (((a + 1) * sin(u + _PI / a) - b * sin((u + _PI / a) * (a + 1))) / (a + 1 - b));
	double Xa = (((a + 1) * cos(_PI / a) - b * cos(_PI / a * (a + 1))) / (a + 1 - b));
	double Ya = (((a + 1) * sin(_PI / a) - b * sin(_PI / a * (a + 1))) / (a + 1 - b));

	double X1 = cos(_PI / a) * X + sin(_PI / a) * Y;
	double X1a = cos(_PI / a) * Xa + sin(_PI / a) * Ya;

	return X1 / X1a;

	return	(5.0 / 7.0) * (((a + 1) * cos(u) - b * cos(u * (a + 1))) / (a + 1 - b));
}
double BaseEpiCycloid_Y(double u, double a, double b)
{
	double X = (((a + 1) * cos(u + _PI / a) - b * cos((u + _PI / a) * (a + 1))) / (a + 1 - b));
	double Y = (((a + 1) * sin(u + _PI / a) - b * sin((u + _PI / a) * (a + 1))) / (a + 1 - b));
	double Xa = (((a + 1) * cos(_PI / a) - b * cos(_PI / a * (a + 1))) / (a + 1 - b));
	double Ya = (((a + 1) * sin(_PI / a) - b * sin(_PI / a * (a + 1))) / (a + 1 - b));

	double Y1 = sin(_PI / a) * X - cos(_PI / a) * Y;
	double X1a = cos(_PI / a) * Xa + sin(_PI / a) * Ya;

	return -Y1 / X1a;
	return	(5.0 / 7.0) * (((a + 1) * sin(u) - b * sin(u * (a + 1))) / (a + 1 - b));
}
double BaseHypoCycloid_X(double u, double a, double b)
{
	return	(((a - 1) * cos(u) + b * cos(u * (a - 1))) / (a - 1 + b));
}
double BaseHypoCycloid_Y(double u, double a, double b)
{
	return	(((a - 1) * sin(u) - b * sin(u * (a - 1))) / (a - 1 + b));
}
double BaseShuriken_X(double u, double a, double b)
{
	return (((a - 1) * cos(u) + b * fabs(cos(a * u)) * cos(u * (a - 1))) / (a - 1 + b));
}
double BaseShuriken_Y(double u, double a, double b)
{
	return (((a - 1) * sin(u) - b * fabs(cos(a * u)) * sin(u * (a - 1))) / (a - 1 + b));
}
double BaseTablemat_XOld(double u, double a, double b)
{
	double X, Y, X1, Y1;

	X = (cos(u) * (1 + (.2 * b) * cos(a * u)) / (1 + (.2 * b)));
	Y = (sin(u) * (1 + (.2 * b) * cos(a * u)) / (1 + (.2 * b)));

	X1 = -Y - cos(u) * (.2 * a * b * sin(a * u)) / (1 + (.2 * b));
	Y1 = X - sin(u) * (.2 * a * b * sin(a * u)) / (1 + (.2 * b));

	return X + .2 * Y1;
	return (cos(u) * (1 + (.2 * b) * cos(a * u)) / (1 + (.2 * b)));
}
double BaseTablemat_YOld(double u, double a, double b)
{
	double X, Y, X1, Y1;

	X = (cos(u) * (1 + (.2 * b) * cos(a * u)) / (1 + (.2 * b)));
	Y = (sin(u) * (1 + (.2 * b) * cos(a * u)) / (1 + (.2 * b)));

	X1 = -Y - cos(u) * (.2 * a * b * sin(a * u)) / (1 + (.2 * b));
	Y1 = X - sin(u) * (.2 * a * b * sin(a * u)) / (1 + (.2 * b));

	return Y - .2 * X1;

	return (sin(u) * (1 + (.2 * b) * cos(a * u)) / (1 + (.2 * b)));
}
double BaseTablemat_X(double u, double a, double b)
{
	b /= 100;
	double r = (1 + b) / 2 + (1 - b) * sin(a * u) / 2;

	return r * sin(u);
}
double BaseTablemat_Y(double u, double a, double b)
{
	b /= 100;
	double r = (1 + b) / 2 + (1 - b) * sin(a * u) / 2;

	return r * cos(u);
}
double BaseCircle_X(double u, double a, double b)
{
	return cos(u);
}
double BaseCircle_Y(double u, double a, double b)
{
	return sin(u);
}
double BasePolygon_X(double u, double a, double b)
{
	double x = fmod(round(100 + a * u / _PI), 2);
	b -= 1;
	return (1 - fabs(sin(a * u)) * b * x / 100) * (cos(u) * cos(_PI / a) / cos(_PI / a - fmod((200 * _PI / a) + u, 2 * _PI / a)));
}
double BasePolygon_Y(double u, double a, double b)
{
	double x = fmod(round(100 + a * u / _PI), 2);
	b -= 1;
	return (1 - fabs(sin(a * u)) * b * x / 100) * (sin(u) * cos(_PI / a) / cos(_PI / a - fmod((200 * _PI / a) + u, 2 * _PI / a)));
}
double BaseCatsEar_X(double u, double a, double b)
{
	double x = (1 - fmod(100 + round(a * u / _PI), 2));
	b = 50 - b;
	return (1 - fabs(pow(sin(2 * a * u), 7)) * b * x / 100) * (cos(u) * cos(_PI / a) / cos(_PI / a - fmod((200 * _PI / a) + u, 2 * _PI / a)));
}
double BaseCatsEar_Y(double u, double a, double b)
{
	double x = 1 - fmod(100 + round(a * u / _PI), 2);
	b = 50 - b;
	return (1 - fabs(pow(sin(2 * a * u), 7)) * b * x / 100) * (sin(u) * cos(_PI / a) / cos(_PI / a - fmod((200 * _PI / a) + u, 2 * _PI / a)));
}
double BaseConvexPolygon_X(double u, double a, double b)
{
	return (cos(u) * cos(_PI / a) / cos(_PI / a - b * mathfmod(u, 2 * _PI / a)));
}
double BaseConvexPolygon_Y(double u, double a, double b)
{
	return (sin(u) * cos(_PI / a) / cos(_PI / a - b * mathfmod(u, 2 * _PI / a)));
}
double BaseHook_X(double u, double a, double b)
{
	return (1 / ((1 + b / 10) * (1 + b / 10))) * (cos(u) * (1 + (b / 10) * pow(cos(a * u), 7)) * (1 + (b / 10) * pow(cos(a * u), 7)) * (1 + (b / 10) * pow(cos(a * u), 7)) / (1 + b / 10));


	return (0.5 / 0.605) * (cos(u) * (1 + (b / 10) * pow(cos(a * u), 7)) * (1 + (b / 10) * pow(cos(a * u), 7)) * (1 + (b / 10) * pow(cos(a * u), 7)) / (1 + b / 10));
}
double BaseHook_Y(double u, double a, double b)
{
	return (1 / ((1 + b / 10) * (1 + b / 10))) * (sin(u) * (1 + (b / 10) * pow(cos(a * u), 7)) * (1 + (b / 10) * pow(cos(a * u), 7)) * (1 + (b / 10) * pow(cos(a * u), 7)) / (1 + b / 10));


	return (0.5 / 0.605) * (sin(u) * (1 + (b / 10) * pow(cos(a * u), 7)) * (1 + (b / 10) * pow(cos(a * u), 7)) * (1 + (b / 10) * pow(cos(a * u), 7)) / (1 + b / 10));
}
double BaseStar_X(double u, double a, double b)
{
	return (cos(u) * (1 + (.8 * b) * pow(cos(a * u / 2), 14)) / (1 + (.8 * b)));
}
double BaseStar_Y(double u, double a, double b)
{
	return (sin(u) * (1 + (.8 * b) * pow(cos(a * u / 2), 14)) / (1 + (.8 * b)));
}
double BaseArrow_X(double u, double a, double b)
{
	double alpha = (1 + (.2 * b) * (sqrt(2))) / (1 + (.4 * b));
	return (1 / alpha) * (cos(u) * (1 + (.4 * b) * (cos(a * u) - .8 * fabs(cos(2 * a * u)))) / (1 + (.4 * b)));

	return (0.5 / 0.458158) * (cos(u) * (1 + (.4 * b) * (cos(a * u) - .8 * fabs(cos(2 * a * u)))) / (1 + (.4 * b)));

}
double BaseArrow_Y(double u, double a, double b)
{
	double alpha = (1 + (.2 * b) * (sqrt(2))) / (1 + (.4 * b));
	return (1 / alpha) * (sin(u) * (1 + (.4 * b) * (cos(a * u) - .8 * fabs(cos(2 * a * u)))) / (1 + (.4 * b)));

	return (0.5 / 0.458158) * (sin(u) * (1 + (.4 * b) * (cos(a * u) - .8 * fabs(cos(2 * a * u)))) / (1 + (.4 * b)));
}
double BaseBloodyCranesbill_X(double u, double a, double b)
{
	double alpha = sqrt(3) / 2;
	double X = (cos(u) * (1 + (.35 * b) * (cos(a * u) - alpha * cos(2 * a * u))));
	double beta = acos(alpha / 3) / a;
	double R = (1 + (.35 * b) * (cos(a * beta) - alpha * cos(2 * a * beta)));

	return X / R;
	return (0.5 / 0.493731) * (cos(u) * (1 + (.35 * b) * (cos(a * u) - .8 * cos(2 * a * u))) / (1 + (.35 * b)));
}
double BaseBloodyCranesbill_Y(double u, double a, double b)
{
	double alpha = sqrt(3) / 2;
	double Y = (sin(u) * (1 + (.35 * b) * (cos(a * u) - alpha * cos(2 * a * u))));
	double beta = acos(alpha / 3) / a;
	double R = (1 + (.35 * b) * (cos(a * beta) - alpha * cos(2 * a * beta)));
	return Y / R;

	return (0.5 / 0.493731) * (sin(u) * (1 + (.35 * b) * (cos(a * u) - .8 * cos(2 * a * u))) / (1 + (.35 * b)));
}
double BaseWindmill_X(double u, double a, double b)
{
	double alpha = _PI / a;
	double X1 = (floor(mathfmod((u + alpha + _PI / (2 * a)) * a / _PI, 2.0)) + (1.5 * b)) * (((2 * a - 1) * cos(u + alpha) + cos((u + alpha) * (2 * a - 1))) / (2 * a * (1 + (1.5 * b))));
	double Y1 = (floor(mathfmod((u + alpha + _PI / (2 * a)) * a / _PI, 2.0)) + (1.5 * b)) * (((2 * a - 1) * sin(u + alpha) - sin((u + alpha) * (2 * a - 1))) / (2 * a * (1 + (1.5 * b))));

	return cos(alpha) * X1 + sin(alpha) * Y1;
	return (floor(mathfmod((u + _PI / (2 * a)) * a / _PI, 2.0)) + (1.5 * b)) * (((2 * a - 1) * cos(u) + cos(u * (2 * a - 1))) / (2 * a * (1 + (1.5 * b))));
}
double BaseWindmill_Y(double u, double a, double b)
{
	double alpha = _PI / a;
	double X1 = (floor(mathfmod((u + alpha + _PI / (2 * a)) * a / _PI, 2.0)) + (1.5 * b)) * (((2 * a - 1) * cos(u + alpha) + cos((u + alpha) * (2 * a - 1))) / (2 * a * (1 + (1.5 * b))));
	double Y1 = (floor(mathfmod((u + alpha + _PI / (2 * a)) * a / _PI, 2.0)) + (1.5 * b)) * (((2 * a - 1) * sin(u + alpha) - sin((u + alpha) * (2 * a - 1))) / (2 * a * (1 + (1.5 * b))));

	return -(sin(alpha) * X1 - cos(alpha) * Y1);
	double _X1 = (floor(mathfmod((u + _PI / (2 * a)) * a / _PI, 2.0)) + (1.5 * b)) * (((2 * a - 1) * cos(u) + cos(u * (2 * a - 1))) / (2 * a * (1 + (1.5 * b))));
	double _Y1 = (floor(mathfmod((u + _PI / (2 * a)) * a / _PI, 2.0)) + (1.5 * b)) * (((2 * a - 1) * sin(u) - sin(u * (2 * a - 1))) / (2 * a * (1 + (1.5 * b))));

	return sin(_PI / a) * X1 - cos(_PI / a) * Y1;
	return (floor(mathfmod((u + _PI / (2 * a)) * a / _PI, 2.0)) + (1.5 * b)) * (((2 * a - 1) * sin(u) - sin(u * (2 * a - 1))) / (2 * a * (1 + (1.5 * b))));
}
double BaseMalteseCross_X(double u, double a, double b)
{
	double alpha = 3 * _PI / (2 * a);
	double X1 = (1 - .5 * fabs(sin(a * (u + alpha)) + (.6 * b) * fabs(sin(a * (u + alpha))))) * ((((2 * a - 1) * cos(u + alpha) + cos((u + alpha) * (2 * a - 1))) / (2 * a)));
	double Y1 = (1 - .5 * fabs(sin(a * (u + alpha)) + (.6 * b) * fabs(sin(a * (u + alpha))))) * ((((2 * a - 1) * sin(u + alpha) - sin((u + alpha) * (2 * a - 1))) / (2 * a)));

	return cos(alpha) * X1 + sin(alpha) * Y1;

	return (1 - .5 * fabs(sin(a * (u)) + (.6 * b) * fabs(sin(a * u)))) * ((((2 * a - 1) * cos(u) + cos(u * (2 * a - 1))) / (2 * a)));
}
double BaseMalteseCross_Y(double u, double a, double b)
{
	double alpha = 3 * _PI / (2 * a);
	double X1 = (1 - .5 * fabs(sin(a * (u + alpha)) + (.6 * b) * fabs(sin(a * (u + alpha))))) * ((((2 * a - 1) * cos(u + alpha) + cos((u + alpha) * (2 * a - 1))) / (2 * a)));
	double Y1 = (1 - .5 * fabs(sin(a * (u + alpha)) + (.6 * b) * fabs(sin(a * (u + alpha))))) * ((((2 * a - 1) * sin(u + alpha) - sin((u + alpha) * (2 * a - 1))) / (2 * a)));

	return -(sin(alpha) * X1 - cos(alpha) * Y1);

	return (1 - .5 * fabs(sin(a * (u)) + (.6 * b) * fabs(sin(a * u)))) * ((((2 * a - 1) * sin(u) - sin(u * (2 * a - 1))) / (2 * a)));
}
double BaseBasileus_X(double u, double a, double b)
{

	double alpha = 3 * _PI / (2 * a);

	double X = (1 - .175 * fabs(sin(a * (u + alpha)) + b * fabs(sin(a * (u + alpha))))) * (((2 * a + 1) * cos(u + alpha) - cos((u + alpha) * (2 * a + 1))) / (2 * a));
	double Y = (1 - .175 * fabs(sin(a * (u + alpha)) + b * fabs(sin(a * (u + alpha))))) * (((2 * a + 1) * sin(u + alpha) - sin((u + alpha) * (2 * a + 1))) / (2 * a));

	double X1 = cos(alpha) * X + sin(alpha) * Y;

	double X2 = (1 - .175 * fabs(sin(a * (alpha)) + b * fabs(sin(a * (alpha))))) * (((2 * a + 1) * cos(alpha) - cos((alpha) * (2 * a + 1))) / (2 * a));
	double Y2 = (1 - .175 * fabs(sin(a * (alpha)) + b * fabs(sin(a * (alpha))))) * (((2 * a + 1) * sin(alpha) - sin((alpha) * (2 * a + 1))) / (2 * a));
	double X3 = cos(alpha) * X2 + sin(alpha) * Y2;

	return X1 / X3;

	double _X1 = (5.0 / 6.0) * (1 - .175 * fabs(sin(a * (u)) + b * fabs(sin(a * u)))) * (((2 * a + 1) * cos(u) - cos(u * (2 * a + 1))) / (2 * a));
	double _Y1 = (5.0 / 6.0) * (1 - .175 * fabs(sin(a * (u)) + b * fabs(sin(a * u)))) * (((2 * a + 1) * sin(u) - sin(u * (2 * a + 1))) / (2 * a));

	double _X3 = (5.0 / 6.0) * (1 - .175 * fabs(sin(a * (3 * _PI) / (2 * a)) + b * fabs(sin(a * (3 * _PI) / (2 * a))))) * (((2 * a + 1) * cos((3 * _PI) / (2 * a)) - cos((3 * _PI) / (2 * a) * (2 * a + 1))) / (2 * a));
	double _Y3 = (5.0 / 6.0) * (1 - .175 * fabs(sin(a * (3 * _PI) / (2 * a)) + b * fabs(sin(a * (3 * _PI) / (2 * a))))) * (((2 * a + 1) * sin((3 * _PI) / (2 * a)) - sin((3 * _PI) / (2 * a) * (2 * a + 1))) / (2 * a));

	double _X2 = cos(3 * _PI / (2 * a)) * _X1 + sin(3 * _PI / (2 * a)) * _Y1;

	return _X2 / (cos(3 * _PI / (2 * a)) * _X3 + sin(3 * _PI / (2 * a)) * _Y3);

	return	(5.0 / 6.0) * (1 - .175 * fabs(sin(a * (u)) + b * fabs(sin(a * u)))) * (((2 * a + 1) * cos(u) - cos(u * (2 * a + 1))) / (2 * a));
}
double BaseBasileus_Y(double u, double a, double b)
{
	double alpha = 3 * _PI / (2 * a);

	double X = (1 - .175 * fabs(sin(a * (u + alpha)) + b * fabs(sin(a * (u + alpha))))) * (((2 * a + 1) * cos(u + alpha) - cos((u + alpha) * (2 * a + 1))) / (2 * a));
	double Y = (1 - .175 * fabs(sin(a * (u + alpha)) + b * fabs(sin(a * (u + alpha))))) * (((2 * a + 1) * sin(u + alpha) - sin((u + alpha) * (2 * a + 1))) / (2 * a));

	double Y1 = sin(alpha) * X - cos(alpha) * Y;

	double X2 = (1 - .175 * fabs(sin(a * (alpha)) + b * fabs(sin(a * (alpha))))) * (((2 * a + 1) * cos(alpha) - cos((alpha) * (2 * a + 1))) / (2 * a));
	double Y2 = (1 - .175 * fabs(sin(a * (alpha)) + b * fabs(sin(a * (alpha))))) * (((2 * a + 1) * sin(alpha) - sin((alpha) * (2 * a + 1))) / (2 * a));
	double X3 = cos(alpha) * X2 + sin(alpha) * Y2;

	return -Y1 / X3;


	double _X1 = (5.0 / 6.0) * (1 - .175 * fabs(sin(a * (u)) + b * fabs(sin(a * u)))) * (((2 * a + 1) * cos(u) - cos(u * (2 * a + 1))) / (2 * a));
	double _Y1 = (5.0 / 6.0) * (1 - .175 * fabs(sin(a * (u)) + b * fabs(sin(a * u)))) * (((2 * a + 1) * sin(u) - sin(u * (2 * a + 1))) / (2 * a));

	double _X3 = (5.0 / 6.0) * (1 - .175 * fabs(sin(a * (3 * _PI) / (2 * a)) + b * fabs(sin(a * (3 * _PI) / (2 * a))))) * (((2 * a + 1) * cos((3 * _PI) / (2 * a)) - cos((3 * _PI) / (2 * a) * (2 * a + 1))) / (2 * a));
	double _Y3 = (5.0 / 6.0) * (1 - .175 * fabs(sin(a * (3 * _PI) / (2 * a)) + b * fabs(sin(a * (3 * _PI) / (2 * a))))) * (((2 * a + 1) * sin((3 * _PI) / (2 * a)) - sin((3 * _PI) / (2 * a) * (2 * a + 1))) / (2 * a));

	double _Y2 = sin(3 * _PI / (2 * a)) * _X1 - cos(3 * _PI / (2 * a)) * _Y1;

	return _Y2 / (cos(3 * _PI / (2 * a)) * _X3 + sin(3 * _PI / (2 * a)) * _Y3);

	return	(5.0 / 6.0) * (1 - .175 * fabs(sin(a * (u)) + b * fabs(sin(a * u)))) * (((2 * a + 1) * sin(u) - sin(u * (2 * a + 1))) / (2 * a));
}
double BaseTruncatedStar_X(double u, double a, double b)
{
	double alpha = 3 * _PI / (2 * a);
	double X1 = (1 - .13 * fabs(sin(a * (u + alpha)) + b * fabs(sin(a * (u + alpha))))) * (cos(u + alpha) * cos(_PI / (2 * a)) / cos(_PI / (2 * a) - fmod((100 * _PI / a) + u + alpha, _PI / a)));
	double Y1 = (1 - .13 * fabs(sin(a * (u + alpha)) + b * fabs(sin(a * (u + alpha))))) * (sin(u + alpha) * cos(_PI / (2 * a)) / cos(_PI / (2 * a) - fmod((100 * _PI / a) + u + alpha, _PI / a)));

	return cos(alpha) * X1 + sin(alpha) * Y1;


	double _X1 = (1 - .13 * fabs(sin(a * (u)) + b * fabs(sin(a * u)))) * (cos(u) * cos(_PI / (2 * a)) / cos(_PI / (2 * a) - fmod((100 * _PI / a) + u, _PI / a)));
	double _Y1 = (1 - .13 * fabs(sin(a * (u)) + b * fabs(sin(a * u)))) * (sin(u) * cos(_PI / (2 * a)) / cos(_PI / (2 * a) - fmod((100 * _PI / a) + u, _PI / a)));

	return cos(3 * _PI / (2 * a)) * X1 + sin(3 * _PI / (2 * a)) * Y1;

	return (1 - .13 * fabs(sin(a * (u)) + b * fabs(sin(a * u)))) * (cos(u) * cos(_PI / (2 * a)) / cos(_PI / (2 * a) - fmod((100 * _PI / a) + u, _PI / a)));
}
double BaseTruncatedStar_Y(double u, double a, double b)
{

	double alpha = 3 * _PI / (2 * a);
	double X1 = (1 - .13 * fabs(sin(a * (u + alpha)) + b * fabs(sin(a * (u + alpha))))) * (cos(u + alpha) * cos(_PI / (2 * a)) / cos(_PI / (2 * a) - fmod((100 * _PI / a) + u + alpha, _PI / a)));
	double Y1 = (1 - .13 * fabs(sin(a * (u + alpha)) + b * fabs(sin(a * (u + alpha))))) * (sin(u + alpha) * cos(_PI / (2 * a)) / cos(_PI / (2 * a) - fmod((100 * _PI / a) + u + alpha, _PI / a)));


	return -(sin(alpha) * X1 - cos(alpha) * Y1);

	double _X1 = (1 - .13 * fabs(sin(a * (u)) + b * fabs(sin(a * u)))) * (cos(u) * cos(_PI / (2 * a)) / cos(_PI / (2 * a) - fmod((100 * _PI / a) + u, _PI / a)));
	double _Y1 = (1 - .13 * fabs(sin(a * (u)) + b * fabs(sin(a * u)))) * (sin(u) * cos(_PI / (2 * a)) / cos(_PI / (2 * a) - fmod((100 * _PI / a) + u, _PI / a)));

	return sin(3 * _PI / (2 * a)) * X1 - cos(3 * _PI / (2 * a)) * Y1;

	return (1 - .13 * fabs(sin(a * (u)) + b * fabs(sin(a * u)))) * (sin(u) * cos(_PI / (2 * a)) / cos(_PI / (2 * a) - fmod((100 * _PI / a) + u, _PI / a)));
}
double BaseSaintAndrew_X1(double u, double a, double b)
{
	return (0.5 / 0.612847) * (sin(a * u) + 2 * b) * (((2 * a - 1) * cos(u) + cos(+u * (2 * a - 1))) / (2 * a * (1 + b)));
}
double BaseSaintAndrew_Y1(double u, double a, double b)
{
	return (0.5 / 0.612847) * (sin(a * u) + 2 * b) * (((2 * a - 1) * sin(u) - sin(u * (2 * a - 1))) / (2 * a * (1 + b)));
}
double BaseSaintAndrew_X(double u, double a, double b)
{
	double X = (0.5 / 0.612847) * (sin(a * u + _PI / 2) + 2 * b) * (((2 * a - 1) * sin(u + _PI / 2) - sin(_PI / 2 + u * (2 * a - 1))) / (2 * a * (1 + b)));

	double X0 = (0.5 / 0.612847) * (1 + 2 * b) * (((2 * a - 1) - 1) / (2 * a * (1 + b)));
	return X / X0;
	return (0.5 / 0.612847) * (sin(a * u + _PI / 2) + 2 * b) * (((2 * a - 1) * cos(u + _PI / 2) + cos(_PI / 2 + u * (2 * a - 1))) / (2 * a * (1 + b)));
}
double BaseSaintAndrew_Y(double u, double a, double b)
{
	double Y = (0.5 / 0.612847) * (sin(a * u + _PI / 2) + 2 * b) * (((2 * a - 1) * cos(u + _PI / 2) + cos(_PI / 2 + u * (2 * a - 1))) / (2 * a * (1 + b)));
	double X0 = (0.5 / 0.612847) * (1 + 2 * b) * (((2 * a - 1) - 1) / (2 * a * (1 + b)));

	return -Y / X0;
	return (0.5 / 0.612847) * (sin(a * u + _PI / 2) + 2 * b) * (((2 * a - 1) * sin(u + _PI / 2) - sin(_PI / 2 + u * (2 * a - 1))) / (2 * a * (1 + b)));
}
double BaseStarPolygon_X(double u, double a, double b)
{
	double alpha = _PI / a;
	double X1 = 0.5 * (cos(u + alpha) * (1 + b) * sin(_PI / a) / ((1 + b) * sin(fabs(_PI / a - fmod((200 * _PI / a) + u + alpha, 2 * _PI / a))) + sin(_PI / a - fabs(_PI / a - fmod((200 * _PI / a) + u + alpha, 2 * _PI / a)))));
	double Y1 = 0.5 * (sin(u + alpha) * (1 + b) * sin(_PI / a) / ((1 + b) * sin(fabs(_PI / a - fmod((200 * _PI / a) + u + alpha, 2 * _PI / a))) + sin(_PI / a - fabs(_PI / a - fmod((200 * _PI / a) + u + alpha, 2 * _PI / a)))));

	double X2 = (cos(alpha) * X1 + sin(alpha) * Y1);
	double Y2 = (sin(alpha) * X1 - cos(alpha) * Y1);

	return 2 * X2 / (b + 1);

	return (cos(alpha) * X1 + sin(alpha) * Y1) * 2 / (b + 1);

	return 0.5 * (cos(u) * (1 + b) * sin(_PI / a) / ((1 + b) * sin(fabs(_PI / a - fmod((200 * _PI / a) + u, 2 * _PI / a))) + sin(_PI / a - fabs(_PI / a - fmod((200 * _PI / a) + u, 2 * _PI / a)))));
}
double BaseStarPolygon_Y(double u, double a, double b)
{
	double alpha = _PI / a;
	double X1 = 0.5 * (cos(u + alpha) * (1 + b) * sin(_PI / a) / ((1 + b) * sin(fabs(_PI / a - fmod((200 * _PI / a) + u + alpha, 2 * _PI / a))) + sin(_PI / a - fabs(_PI / a - fmod((200 * _PI / a) + u + alpha, 2 * _PI / a)))));
	double Y1 = 0.5 * (sin(u + alpha) * (1 + b) * sin(_PI / a) / ((1 + b) * sin(fabs(_PI / a - fmod((200 * _PI / a) + u + alpha, 2 * _PI / a))) + sin(_PI / a - fabs(_PI / a - fmod((200 * _PI / a) + u + alpha, 2 * _PI / a)))));

	double X2 = (cos(alpha) * X1 + sin(alpha) * Y1);
	double Y2 = (sin(alpha) * X1 - cos(alpha) * Y1);

	return -2 * Y2 / (b + 1);



	double _X1 = 0.5 * (cos(u) * (1 + b) * sin(_PI / a) / ((1 + b) * sin(fabs(_PI / a - fmod((200 * _PI / a) + u, 2 * _PI / a))) + sin(_PI / a - fabs(_PI / a - fmod((200 * _PI / a) + u, 2 * _PI / a)))));
	double _Y1 = 0.5 * (sin(u) * (1 + b) * sin(_PI / a) / ((1 + b) * sin(fabs(_PI / a - fmod((200 * _PI / a) + u, 2 * _PI / a))) + sin(_PI / a - fabs(_PI / a - fmod((200 * _PI / a) + u, 2 * _PI / a)))));

	return (sin(_PI / (a)) * X1 - cos(_PI / (a)) * Y1) * 2 / (b + 1);

	return 0.5 * (sin(u) * (1 + b) * sin(_PI / a) / ((1 + b) * sin(fabs(_PI / a - fmod((200 * _PI / a) + u, 2 * _PI / a))) + sin(_PI / a - fabs(_PI / a - fmod((200 * _PI / a) + u, 2 * _PI / a)))));
}
double BaseSawToooth_X(double u, double a, double b)
{
	double fMaxi = (double)((myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides + 1);

	double alpha = 2 * _PI / a;
	double Xb = b * cos(alpha);
	double Yb = b * sin(alpha);
	double k = Yb / (Xb - 1);
	double x = fmod(u, alpha);
	double y = fmod(a * u, 2 * _PI);

	if (abs(abs(x - alpha / 2) - alpha / 2) < alpha / (2 * fMaxi))
		x = 0;

	double r = k / ((k - tan(x)) * cos(x));
	return r * cos(u);

}
double BaseSawToooth_Y(double u, double a, double b)
{
	double fMaxi = (double)((myDome.nMeridianMain + (myDome.nMeridianMain) * myDome.nMeridianMinor) * myDome.nBaseSides + 1);
	double alpha = 2 * _PI / a;
	double Xb = b * cos(alpha);
	double Yb = b * sin(alpha);
	double k = Yb / (Xb - 1);
	double x = fmod(u, alpha);
	double y = fmod(a * u, 2 * _PI);

	if (abs(abs(x - alpha / 2) - alpha / 2) < alpha / (2 * fMaxi))
		x = 0;

	double r = k / ((k - tan(x)) * cos(x));
	return r * sin(u);


}
double BaseGear_X(double u, double a, double b)
{
	double alpha = _PI / a;

	double X1 = ((cos(u + alpha) * cos(alpha) / cos(alpha - fmod((200 * alpha) + u + alpha, 2 * alpha))) + cos(u + alpha) * ((1.7 * b) / 2 - fabs((1.7 * b) / 2 * (2 * (fmod((200 * alpha) + u + alpha, 2 * alpha) * a / (2 * _PI)) - 1))));
	double Y1 = ((sin(u + alpha) * cos(alpha) / cos(alpha - fmod((200 * alpha) + u + alpha, 2 * alpha))) + sin(u + alpha) * ((1.7 * b) / 2 - fabs((1.7 * b) / 2 * (2 * (fmod((200 * alpha) + u + alpha, 2 * alpha) * a / (2 * _PI)) - 1))));

	double X3 = ((cos(alpha) * cos(alpha) / cos(alpha - fmod((200 * alpha) + (alpha), 2 * alpha))) + cos((alpha)) * ((1.7 * b) / 2 - fabs((1.7 * b) / 2 * (2 * (fmod((200 * alpha) + (alpha), 2 * alpha) * a / (2 * _PI)) - 1))));
	double Y3 = ((sin(alpha) * cos(alpha) / cos(alpha - fmod((200 * alpha) + (alpha), 2 * alpha))) + sin((alpha)) * ((1.7 * b) / 2 - fabs((1.7 * b) / 2 * (2 * (fmod((200 * alpha) + (alpha), 2 * alpha) * a / (2 * _PI)) - 1))));

	return (cos(alpha) * X1 + sin(alpha) * Y1) / (cos(alpha) * X3 + sin(alpha) * Y3);
	return (0.5 / 0.829508) * ((cos(u) * cos(_PI / a) / cos(_PI / a - fmod((200 * alpha) + u, 2 * _PI / a))) + cos(u) * ((1.7 * b) / 2 - fabs((1.7 * b) / 2 * (2 * (fmod((200 * alpha) + u, 2 * _PI / a) * a / (2 * _PI)) - 1))));
}
double BaseGear_Y(double u, double a, double b)
{
	double alpha = _PI / a;

	double X1 = ((cos(u + alpha) * cos(alpha) / cos(alpha - fmod((200 * alpha) + u + alpha, 2 * alpha))) + cos(u + alpha) * ((1.7 * b) / 2 - fabs((1.7 * b) / 2 * (2 * (fmod((200 * alpha) + u + alpha, 2 * alpha) * a / (2 * _PI)) - 1))));
	double Y1 = ((sin(u + alpha) * cos(alpha) / cos(alpha - fmod((200 * alpha) + u + alpha, 2 * alpha))) + sin(u + alpha) * ((1.7 * b) / 2 - fabs((1.7 * b) / 2 * (2 * (fmod((200 * alpha) + u + alpha, 2 * alpha) * a / (2 * _PI)) - 1))));

	double X3 = ((cos(alpha) * cos(alpha) / cos(alpha - fmod((200 * alpha) + (alpha), 2 * alpha))) + cos((alpha)) * ((1.7 * b) / 2 - fabs((1.7 * b) / 2 * (2 * (fmod((200 * alpha) + (alpha), 2 * alpha) * a / (2 * _PI)) - 1))));
	double Y3 = ((sin(alpha) * cos(alpha) / cos(alpha - fmod((200 * alpha) + (alpha), 2 * alpha))) + sin((alpha)) * ((1.7 * b) / 2 - fabs((1.7 * b) / 2 * (2 * (fmod((200 * alpha) + (alpha), 2 * alpha) * a / (2 * _PI)) - 1))));

	return -(sin(alpha) * X1 - cos(alpha) * Y1) / (cos(alpha) * X3 + sin(alpha) * Y3);
	return (0.5 / 0.829508) * ((sin(u) * cos(_PI / a) / cos(_PI / a - fmod((200 * alpha) + u, 2 * _PI / a))) + sin(u) * ((1.7 * b) / 2 - fabs((1.7 * b) / 2 * (2 * (fmod((200 * alpha) + u, 2 * _PI / a) * a / (2 * _PI)) - 1))));
}
double BaseSnowFlake_X(double u, double a, double b)
{
	b = b / 10;

	double Z = (1.5 - (sqrt(fabs(sin(a * u / 2))) + b * (cos(a * u)))) / (1 + b);
	double Z0 = (1.5 - b) / (1 + b);

	return cos(u) * Z / Z0;
	return ((cos(u) * cos(_PI / a) / cos(_PI / a - b * fmod((200 * _PI / a) + u, 2 * _PI / a))) * (1 + .2 * b * cos(2 * a * u)) / (1 + .2 * b));
}
double BaseSnowFlake_Y(double u, double a, double b)
{
	b = b / 10;

	double Z = (1.5 - (sqrt(fabs(sin(a * u / 2))) + b * (cos(a * u)))) / (1 + b);
	double Z0 = (1.5 - b) / (1 + b);

	return sin(u) * Z / Z0;
	return ((sin(u) * cos(_PI / a) / cos(_PI / a - b * fmod((200 * _PI / a) + u, 2 * _PI / a))) * (1 + .2 * b * cos(2 * a * u)) / (1 + .2 * b));
}
double BaseDaisy_X(double u, double a, double b)
{
	return ((cos(u) * (1 + b * abs(cos(a * u / 2))) / (1 + b)));
}
double BaseDaisy_Y(double u, double a, double b)
{
	return ((sin(u) * (1 + b * abs(cos(a * u / 2))) / (1 + b)));
}
double BaseClover_X(double u, double a, double b)
{
	double alpha = ((1 + b * (cos(2 * _PI / 5) - .8 * cos(4 * _PI / 5))) / (1 + b));

	return (((cos(u) * (1 + b * (cos(a * u) - .8 * cos(2 * a * u))) / (1 + b)))) / alpha;
}
double BaseClover_Y(double u, double a, double b)
{
	double alpha = ((1 + b * (cos(2 * _PI / 5) - .8 * cos(4 * _PI / 5))) / (1 + b));

	return (((sin(u) * (1 + b * (cos(a * u) - .8 * cos(2 * a * u))) / (1 + b)))) / alpha;
}
double BaseMandala_X(double u, double a, double b)
{
	b = -.3 + 1.3 * b / 100;
	double r = ((2.5 + (((fabs(cos(a * u / 2))) + (b - (fabs(cos(a * u / 2 + _PI / 2)))) * 2) / (2 + fabs(cos(2 * a * u + _PI / 2)) * 8)))) / (3 + b);
	return r * cos(u);
}
double BaseMandala_Y(double u, double a, double b)
{
	b = -.3 + 1.3 * b / 100;
	double r = ((2.5 + (((fabs(cos(a * u / 2))) + (b - (fabs(cos(a * u / 2 + _PI / 2)))) * 2) / (2 + fabs(cos(2 * a * u + _PI / 2)) * 8)))) / (3 + b);
	return r * sin(u);

}
double BaseRoissy_X(double u, double a, double b)
{
	double r = (1 + .1 * ceil((pow(fabs(cos(a * u / 2)), b) * 3))) / 1.3;

	return r * cos(u);
}
double BaseRoissy_Y(double u, double a, double b)
{
	double r = (1 + .1 * ceil((pow(fabs(cos(a * u / 2)), b) * 3))) / 1.3;

	return r * sin(u);

}
double BaseBrokenClover_X(double u, double a, double b)
{
	b = b / 500;
	double r = (1 - b * (ceil(fabs(cos(a * (u + _PI / a) + _PI / 2) * 2)) + cos(a * (u + _PI / a))));
	return r * cos(u);
}
double BaseBrokenClover_Y(double u, double a, double b)
{
	b = b / 500;
	double r = (1 - b * (ceil(fabs(cos(a * (u + _PI / a) + _PI / 2) * 2)) + cos(a * (u + _PI / a))));
	return r * sin(u);
}
double BasePropeller_X(double u, double a, double b)
{
	b = 7 * b / 100;
	double r = (b * ceil(fabs(2 * cos(a * u / 2) * cos(a * u / 2))) + fmin(pow(fabs(cos(a * u / 2)), 5), .4)) / (2 * b + .4);
	return r * cos(u);
}
double BasePropeller_Y(double u, double a, double b)
{
	b = 7 * b / 100;
	double r = (b * ceil(fabs(2 * cos(a * u / 2) * cos(a * u / 2))) + fmin(pow(fabs(cos(a * u / 2)), 5), .4)) / (2 * b + .4);
	return r * sin(u);
}
double BaseBrokenCircle_X(double u, double a, double b)

{
	b = b / 10;
	double r = (2 + ceil(b * pow(fabs(cos(a * u / 2)), 2))) / (2 + ceil(b));
	return r * cos(u);
}
double BaseBrokenCircle_Y(double u, double a, double b)

{
	b = b / 10;
	double r = (2 + ceil(b * pow(fabs(cos(a * u / 2)), 2))) / (2 + ceil(b));
	return r * sin(u);
}
double BaseSpiral_X(double u, double a, double b)
{
	b = b / 10;
	double r = (2 * _PI - a * b * mathfmod(u, 2 * _PI / a)) / (2 * _PI);
	return r * cos(u);
}
double BaseSpiral_Y(double u, double a, double b)
{
	b = b / 10;
	double r = (2 * _PI - a * b * mathfmod(u, 2 * _PI / a)) / (2 * _PI);
	return r * sin(u);
}
double BasePointedPolygon_X(double u, double a, double b)
{
	b = b / 10;
	double r;

	r = (int)a % 2 ? (cos(_PI / a) / cos(_PI / a - b * mathfmod(u, 2 * _PI / a))) : (cos(2 * _PI / a) / cos(2 * _PI / a - b * mathfmod(u, 4 * _PI / a)));
	return r * cos(u);
}
double BasePointedPolygon_Y(double u, double a, double b)
{
	b = b / 10;
	double r;

	r = (int)a % 2 ? (cos(_PI / a) / cos(_PI / a - b * mathfmod(u, 2 * _PI / a))) : (cos(2 * _PI / a) / cos(2 * _PI / a - b * mathfmod(u, 4 * _PI / a)));
	return r * sin(u);
}


double BaseRoof_X(double u, double a, double b)
{
	if (u < _PI)
		return b * u / (_PI);
	else
		return (2 * _PI - u) * b / (_PI);
}
double BaseRoof_Y(double u, double a, double b)
{
	if (u < _PI)
		return -1;
	else
		return 1;
}

double CoeffFourrier(int nSide, int nDegree)
{
	double dRetour = 0;
	for (double t = -_PI / nSide; t <= _PI / nSide; t += _PI / 10000)
	{
		dRetour += cos(nSide * nDegree * t) / cos(t);
	}
	return dRetour;
}

double BaseFourierRho(double v, double fParam, bool bCalcul)
{
	static int snSide = -1;
	static double dCoefficient[21];
	int nSide = myDome.nBaseSides;
	int nMaxDegree = (int)(myDome.fBaseParameter / 5) + 1;
	double dRetour = .5;
	if (bCalcul)
	{
		snSide = nSide;
		double dDen = CoeffFourrier(nSide, 0);

		for (int i = 0; i < 21; i++)
			dCoefficient[i] = CoeffFourrier(nSide, i + 1)/dDen;
	}
	else
	{
		if (nSide != snSide)
		{
			printf("Problème Rho : %d vs %d\n", nSide, snSide);
		}
		for (int i = 0; i < nMaxDegree; i++)
			dRetour += dCoefficient[i] * cos((i + 1) * nSide * v);
	}
	return dRetour;
}

double BaseFourier_X(double u, double a, double b)
{
	return cos(u) * BaseFourierRho(u+_PI/a, b, false)/ BaseFourierRho(_PI/a, b, false);
}
double BaseFourier_Y(double u, double a, double b) 
{
	return sin(u) * BaseFourierRho(u+_PI/a, b, false)/ BaseFourierRho(_PI/a, b, false);

}
double BaseSavior_X(double u, double a, double b)
{
	double p = 2*round(b);
	double dRho = (1 + pow(cos(a * u / 2), p)) / (2 * cos(u - 2 * _PI * floor((a * u + _PI) / (2 * _PI)) / a));
	return cos(u) * dRho;
}

double BaseSavior_Y(double u, double a, double b)
{
	double p = 2 * round(b);
	double dRho = (1 + pow(cos(a * u / 2), p)) / (2 * cos(u - 2 * _PI * floor((a * u + _PI) / (2 * _PI)) / a));
	return sin(u) * dRho;

}

double BaseRoundedPolygonRho(double u, double a, double b)
{
	double alpha = (a - 2) * _PI / (2 * a);
	double alpha1 = 2 * _PI / a;
	double r = cos(alpha1 / 2) * b / 100;
	double X1 = 1 - r / sin(alpha);
	double beta = atan(((1 - X1) * tan(alpha)) / (tan(alpha) * tan(alpha) + X1));
	double Y;

	if (fmod(u + beta, alpha1) <= 2 * beta)
	{
		Y = (X1 * cos(fmod(u + beta, alpha1) - beta)) + sqrt(r*r - X1 * X1 * sin(fmod(u + beta, alpha1) - beta) * sin(fmod(u + beta, alpha1) - beta));

	}
	else
	{
		Y = tan(alpha)/(tan(alpha)*cos(fmod(u, alpha1)) + sin(fmod(u,alpha1)) );

	}

	return Y / (X1 + r);

	
}
double BaseRoundedPolygon_X(double u, double a, double b)
{
	return cos(u) * BaseRoundedPolygonRho(u, a, b);
}

double BaseRoundedPolygon_Y(double u, double a, double b)
{
	return sin(u) * BaseRoundedPolygonRho(u, a, b);
}



double SlopeInflected(double v, double fParam)
{

	return (2 * fParam - 2) * v * v * v + (3 - 3 * fParam) * v * v + fParam * v;

}
double SlopeVibrate(double v, double fParam)
{
	double alpha = 3 * fParam / 10;

	if (v == 0)
		return 0;

	return (v + 6 * sin(alpha * 10 * _PI * v) / 40)/(1 + 6 * sin(alpha * 10 * _PI) / 40);

}
double SlopeRim(double v, double fParam)
{
	double alpha = 3 * fParam / 40;

	if (v == 0)
		return 0;


	return ((4 * v * v * v * v * v - 5 * alpha * v * v * v * v) / (4 - 5 * alpha));

}
double SlopeTente(double v, double fParam)
{

	if (v == 0)
		return 0;
	double dTemp;

	if (v < fParam / 10.0)
		dTemp = fParam / 10.0 * v;
	else
		dTemp = (fParam / 10.0 + 1) * (v - 1) + 1;


	return dTemp;

}
double SlopeTriFoil(double v, double fParam)
{

	if (v == 0)
		return 0;

	double dTemp;


	if (v < fParam / 10.0)
		dTemp = fParam / 10.0 - sqrt(fParam * fParam / 100.0 - v * v);
	else
		dTemp = 1 - sqrt((1 - fParam / 10) * (1 - fParam / 10) - (v - fParam / 10) * (v - fParam / 10));

	return dTemp;

}
double SlopeCylinderXY(double v, double fParam)
{
	if (fParam < 0.0)
		fParam = 0.0;
	if (fParam > 1)
		fParam = 1;
	if (v == 0)
		return 0;
	if (v > fParam)
		return 1;
	return v / fParam;
}
double SlopeCylinderZ(double v, double fParam)
{
	if (fParam < 0.0)
		fParam = 0.0;
	if (fParam > 1)
		fParam = 1;
	if (v == 0)
		return 0;
	if (v < fParam)
		return 0;
	return (v - fParam) / (1 - fParam);
}
double SlopeCuve(double v, double fParam)
{
	double alpha = fParam / 10.0;

	if (v == 0)
		return 0;
	double dTemp;


	if (v < .25)
		dTemp = .5;
	else if (v < .5)
		dTemp = alpha;
	else dTemp = 1;

	return dTemp;

}
double SlopeSaintPierre(double v, double fParam)
{
	double alpha = fParam / 20.0;

	if (v == 0)
		return 0;

	double dTemp;

	if (v < .25)
		dTemp = sqrt(2 * v) * sqrt(2 * v) * sqrt(2 * v) * alpha;
	else if (v < .5)
		dTemp = alpha;
	else
		dTemp = (1 - alpha) * sqrt(2 * v - 1) + alpha;



	return dTemp;



}
double SlopeTriCusp(double v, double fParam)
{
	if (v == 0)
		return 0;

	double dTemp;

	if (v < fParam / 10.0)
		dTemp = sqrt(fParam * fParam / 100.0 - (v - fParam / 10) * (v - fParam / 10));
	else
		dTemp = fParam / 10.0 + sqrt((1 - fParam / 10) * (1 - fParam / 10) - (v - 1) * (v - 1));


	return dTemp;
}
double SlopeStraight(double v, double fParam)
{

	if (v == 0)
		return 0.0;

	return pow(v, fParam);


}
double SlopeMongol(double v, double fParam)
{

	if (v == 0)
		return 0.0;

	return pow((-(5.0 * pow(.5 - v, 3) + (.5 - v) - 9.0 / 8.0) * 4.0 / 9.0), fParam);

}
double SlopeHindi(double v, double fParam)
{
	double alpha = fParam * .03;

	if (v == 0)
		return 0.0;


	return  (v + alpha * sin(9 * v * _PI / 2)) / (1 + alpha);
}
double SlopeAngle(double v, double fParam)
{
	if (v == 0)
		return 0.0;
	if (v < fParam / (fParam + 10))
		return  10 * v / fParam;
	else
		return 1 - fParam * (1 - v) / 10;
	return pow(((pow(2 * v - 1, 3) + 1) / 2), fParam);
}
double SlopeAngleX(double v, double fParam)
{
	if (v == 0)
		return 0.0;

	double a = tan(v * _PI / 2);

	if (a < 1)
		return   1 - fParam / ((10 - fParam) * a + fParam);
	else
		return 1 - fParam / (a * fParam + 10 - fParam);


}
double SlopeAngleY(double v, double fParam)
{
	if (v == 0)
		return 0.0;


	double a = tan(v * _PI / 2);

	if (a < 1)
		return a * fParam / ((10 - fParam) * a + fParam);
	else
		return a * fParam / (a * fParam + 10 - fParam);
}
double SlopeJoliCubic(double v, double fParam)
{
	if (v == 0)
		return 0.0;

	return pow(((pow(2 * v - 1, 3) + 1) / 2), fParam);
}
double SlopeCasablanca(double v, double fParam)
{
	if (v == 0)
		return 0.0;
	double x = pow(sin((1 - v) * _PI / 2), 5);
	double dTemp = 1 - pow(((pow(2 * x - 1, 3) + 1) / 2), fParam);

	return dTemp;



}
double SlopeSinus(double v, double fParam)
{
	if (v == 0)
		return 0.0;

	return pow(sin(v * _PI / 2), fParam);

	//	return sin(v * _PI / 2);
}
double SlopeSquareRootSinus(double v, double fParam)
{
	if (v == 0)
		return 0.0;
	return pow(sin(v * _PI / 2), fParam);

	//	return sqrt(sin(v*_PI/2));
}
double SlopeCosinus(double v, double fParam)
{

	if (v == 0)
		return 0.0;

	return (1 - cos(v * _PI / 2)) * pow(sin(v * _PI / 2), (fParam - 2));

	//	return (1 - cos(v * _PI / 2));
}
double SlopeKazackHatXY(double v, double fParam)
{
	double a = fParam * _PI / 10;

	if (v == 0)
		return 0;

	if (v > .5)
		return 1 - .25 + .25 * cos(2 * a * (1 - v));
	else
	{
		double w = 1 - .25 + .25 * cos(a);

		return 2 * w * v;
	}
}
double SlopeKazackHatZ(double v, double fParam)
{
	double b = fParam * _PI / 10;

	if (v == 0)
		return 0;
	if (v > .5)
		return  1 - .25 * sin(2 * b * (1 - v));
	else
	{
		double g = .25 * sin(b);


		return 2 * (1 - g) * v;
	}
}
double SlopePineTreeXY(double v, double fParam)
{

	if (v == 0)
		return 0.0;

	return pow(v, fParam) * (0.5 + ((2 * v - 1) * (pow(2 * v - 1, 2) - 0.6 * 0.6) / (2 * (1 - 0.6 * 0.6))));
}
double SlopePineTreeZ(double v, double fParam)
{

	if (v == 0)
		return 0.0;


	return pow(v, fParam) * ((pow(2 * v - 1, 3) + 1) / 2);
}
double SlopeIranian(double v, double fParam)
{
	if (v == 0)
		return 0.0;
	/*	double y = (sin(53 * deg2rad) + sin(106 * deg2rad) - sin((170 * v + 53) * deg2rad) - sin(2 * (170 * v + 53) * deg2rad)) / (sin(53 * deg2rad) + sin(106 * deg2rad) - sin(223 * deg2rad) - sin(446 * deg2rad));
		double x = (pow(y, fParam));
		if (isnan(x))
			return 0;
		else
			return x;*/
	return (pow((sin(53 * deg2rad) + sin(106 * deg2rad) - sin((170 * v + 53) * deg2rad) - sin(2 * (170 * v + 53) * deg2rad)) / (sin(53 * deg2rad) + sin(106 * deg2rad) - sin(223 * deg2rad) - sin(446 * deg2rad)), fParam));
}
double SlopeFloor(double v, double fParam)
{
	if (v == 1)
		return 1.0;
	if (fParam == -1)
		return 0;

	return (floor((fParam + 1) * v) / (fParam + 1));

	//	return floor(4*v) / 4;
}
double SlopeStaircase(double v, double fParam)
{
	if (v == 0)
		return 0.0;

	if (fParam == -1)
		return 0;


	return (v < fParam / 10 ? fabs(fParam / 10) : 1.0);
	//	return floor(v +.5);
}
double SlopeRoman(double v, double fParam)
{
	if (v == 0)
		return 0.0;

	return  2 * pow(v, fParam) - v * v * v;
}
double SlopeFloorAside(double v, double fParam)
{
	return (floor(4 * v - .5) + 1) / 4;
}
double SlopeExponential(double v, double fParam)
{
	if (v == 0)
		return 0.0;

	return pow(v, fParam) * exp(1 - v);
}
double SlopePower4(double v, double fParam)
{
	if (v == 0)
		return 0;


	return pow(-(9 * v * v * v * v - 24 * v * v * v + 22 * v * v - 8 * v), fParam);
}
double SlopeSqrtPower4(double v, double fParam)
{
	if (v == 0)
		return 0;
	return pow(-(9 * v * v * v * v - 24 * v * v * v + 22 * v * v - 8 * v), fParam);

}
double SlopeRoundedPower4(double v, double fParam)
{


	return (((2 * fParam - 2) * v * v * v + (3 - 3 * fParam) * v * v + fParam * v) + sin(_PI * v / 2)) / 2;

}
double SlopeEigthSquare(double v, double fParam)
{
	if (v == 1)
		return 1;


	return (1 - pow(1 - v, fParam));
}
double SlopeIranian2(double v, double fParam)
{
	double dTemp = (pow((sin(53 * deg2rad) + sin(106 * deg2rad) - sin((170 * v + 53) * deg2rad) - sin(2 * (170 * v + 53) * deg2rad)) / (sin(53 * deg2rad) + sin(106 * deg2rad) - sin(223 * deg2rad) - sin(446 * deg2rad)), fParam));

	return dTemp * v;
	//	return SlopeIranian(v, fParam, dAxe)*v;
}
double SlopeIranian3(double v, double fParam)
{
	if (v == 0)
		return 0;

	return pow(1 - (sin((240 * (1 - v) - 150) * deg2rad) + sin(150 * deg2rad)) / (1 + sin(150 * deg2rad)), fParam);
}
double SlopeCubic2(double v, double fParam)
{
	if (v == 0)
		return 0;
	return pow(4 * v * v * v - 7 * v * v + 4 * v, fParam);
}
double SlopeMaison_XY(double v, double fParam)
{
	double f1 = myDome.fSlopeParameterXY / 10;
	double f2 = -1 + myDome.fSlopeParameterZ;
	double x;
	double S = sqrt((1 + fabs(f2)) * (1 + fabs(f2)) + (1 - f1) * (1 - f1));
	double L = f1 + fabs(f2) + S;


	if (v < S / L)
		x = v * (1 + f2) * L / S;
	else if (v < (S + fabs(f2)) / L)
		x = -v * L * f2 / fabs(f2) + 1 + f2 + S * f2 / fabs(f2);
	else
		x = 1;

	return x;
}
double SlopeMaison_Z(double v, double fParam)
{
	double f1 = myDome.fSlopeParameterXY / 10;
	double f2 = -1 + myDome.fSlopeParameterZ;
	double x;
	double S = sqrt((1 + fabs(f2)) * (1 + fabs(f2)) + (1 - f1) * (1 - f1));
	double L = f1 + fabs(f2) + S;


	if (v < S / L)
		x = 2 * (f1 - 1) * L * L / (S * S) * (v * v / 2 - S * v / L);
	else if (v < (S + fabs(f2)) / L)
		x = 1 - f1;
	else
		x = 1 - L * f1 * (v - 1) / (S + fabs(f2) - L);
	return x;
}
double SlopeBulbXY(double v, double fParam)
{
	double f1 = myDome.fSlopeParameterXY / 10;
	double f2 = myDome.fSlopeParameterZ / 10;

	if (v > f1)
		return sin(v * _PI / 2);
	else
		return v*sin(f1*_PI/2)/f1;

}
double SlopeBulbZ(double v, double fParam)
{
	double f1 = myDome.fSlopeParameterXY / 10;
	double f2 = myDome.fSlopeParameterZ / 10;

	if (v > f1)
		return 1 - f2 * cos(v * _PI / 2);
	else
		return v * (1 - f2 * cos(v * _PI / 2)) / f1;
}
double SlopeMoroccoZ(double v, double fParam)
{
	double a = myDome.fSlopeParameterZ/10;
	double b = myDome.fSlopeParameterXY/10;
	double alpha = 1 - a;
	double beta = 1 - b;
	double h = (alpha - beta) / (alpha + beta);
	h = h * h;
	double p = _PI * (alpha + beta) * (1 + 3 * h / (10 + sqrt(4 - 3 * h)));
	double l = p / 4 + a + b;
	double t1 = a / l;
	double t2 = (a + b) / l;


	t1 = p / (4 * l);
	t2 = (p + 4 * a) / (4 * l);

	if (v < t1)
		return 1-(b + beta * cos(v * _PI / (2 * t1)));
	else if (v < t2)
		return 1-b;
	else
		return 1-(v - 1) * b / (t2 - 1);

	if (v < t1)
		return 1;
	else if (v < t2)
		return ((1 - beta) * v / (t1 - t2) + (beta * t1 - t2) / (t1 - t2));
	else
		return beta * cos(_PI * (v - t2) / (2 * (1 - t2)));
}
double SlopeMoroccoXY(double v, double fParam)
{
	double a = myDome.fSlopeParameterZ/10;
	double b = myDome.fSlopeParameterXY/10;
	double alpha = 1 - a;
	double beta = 1 - b;
	double h = (alpha - beta) / (alpha + beta);
	h = h * h;
	double p = _PI * (alpha + beta) * (1 + 3 * h / (10 + sqrt(4 - 3 * h)));
	double l = p / 4 + a + b;
	double t1 = a / l;
	double t2 = (a + b) / l;

	t1 = p / (4 * l);
	t2 = (p + 4 * a) / (4 * l);

	if (v < t1)
		return  alpha * sin(v * _PI / (2 * t1));
	else if (v < t2)
		return (alpha - 1) * v / (t1 - t2) + (t1 - alpha * t2) / (t1 - t2);
	else
		return 1;

	if (v < t1)
		return a*v/t1;
	else if (v < t2)
		return a;
	else
		return alpha * sin(_PI * (v - t2) / (2 * (1 - t2)));
}
double SlopeDamaoHatXY(double v, double fParam)
{
	return v;
}
double SlopeDamaoHatZ(double v, double fParam)
{

	if (v == 0)
		return 0.0;
	double alpha = myDome.fSlopeParameterZ;
	double beta = myDome.fSlopeParameterXY/10;

	if (v < beta)
		return (1 - alpha) * v * v * v / (beta * beta * beta);
	else
		return alpha * v / (1 - beta) + (1 - alpha - beta) / (1 - beta);

	return pow(v, fParam) * ((pow(2 * v - 1, 3) + 1) / 2);
}


/** Fonctions Window  */

double WindowCycloid(double u, double a, double d, double p, bool bFrise = false)  //Window 4
{
	double v = fabs(_PI - mathfmod(a * d * u, 2 * _PI)) - _PI;
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);
	//0<b<5
	double b = .5 * Parameter;

	return pow(fabs(sin(v / 2) * cos(b * (v - _PI))), p);


	return pow(fabs(sin(a * d * u / 2)), abs(p));
}
double WindowSawTooth(double u, double a, double d, double p, bool bFrise = false)
{

	double w = (mathfmod(u, 2.0 * _PI / (a * d)) * a * d / (2.0 * _PI));

	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);

	//0<b<10
	double b = Parameter / 20;

	return pow(w * (w + b * sin(w * _PI)), p);

	return pow((w + b) / (b + 1), p);

	return pow((mathfmod(u, 2.0 * _PI / (a * d)) * a * d / (2.0 * _PI)), abs(p));
}
double WindowSinus(double u, double a, double d, double p, bool bFrise = false)
{
	double v = mathfmod(a * d * u, 2 * _PI);
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);

	//0<b<20 (entier)
	double b = round(2 * Parameter);

	double k = (_PI + 2 * _PI * floor(b / 2)) / (b + 1);

	return pow(fabs(sin(v / 2) + sin((2 * b + 1) * v / 2) / (2 * b + 1)), p) / pow(fabs(sin(k / 2) + sin((2 * b + 1) * k / 2) / (2 * b + 1)), p);

	return pow((1.0 - cos(a * d * u)) / 2.0, abs(p));
}
double WindowStraight(double u, double a, double d, double p, bool bFrise = false)
{

	double w = (1.0 - fabs(mathfmod(u, 2.0 * _PI / (a * d)) - _PI / (a * d)) * a * d / _PI);
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);

	//-1<b<10	
	double b = 1.1 * Parameter - 1;

	if (b < -.89)
		b = -.89;
	return pow(fmin(abs((w + b) / (b + 1)), 1), p);

	return pow((1.0 - fabs(mathfmod(u, 2.0 * _PI / (a * d)) - _PI / (a * d)) * a * d / _PI), abs(p));
}
double WindowHat(double u, double a, double d, double p, bool bFrise = false)
{
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);

	double b = Parameter;

	return pow((1 + b - cos(a * d * u) - b * cos(2.0 * a * d * u)) / (1 + 2 * b + 1 / 8 * b), p);


	return pow((2.0 - cos(a * d * u) - cos(2.0 * a * d * u)) * 8.0 / 25.0, abs(p));
}
double WindowStep(double u, double a, double d, double p, bool bFrise = false)
{
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);

	double b = 1 + Parameter;

	double v = mathfmod(a * d * u, 2 * _PI);

	return pow(floor(abs(sin(v / 2)) * b) / b, p);


	double x = 1 - fabs((mathfmod(u * a * (d + 1.0), 2 * _PI) - _PI) / _PI);
	double y = abs(p) + 1.0000012354;

	return floor(x * y) / y;
}
double WindowSquare(double u, double a, double d, double p, bool bFrise = false)
{
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);

	double b = Parameter / 10;

	if (fabs(sin(a * d * u / 2)) < b)
		return 0;
	else
		return 1 - pow(abs(cos(a * d * u / 2)), 10 - p);


	double alpha = (p == 0.0 ? 0.0001 : p);

	return floor(alpha * mathfmod((u + (_PI / 2) / (a * d)) * a * d / _PI, 2)) / floor(1.99999999 * alpha);
}
double WindowNone(double u, double a, double d, double p, bool bFrise = false)
{
	return 0.0;
}
double WindowStaircase(double u, double a, double d, double p, bool bFrise = false) // desmos windows3
{
	double v = mathfmod(a * d * u, 2 * _PI);
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);

	double b = 1 + Parameter;

	return pow(floor(sin(v / 2) * (b + 0.012354)) / b, p) / pow(floor(b + 0.012354) / b, p);

	return pow(floor(sin(v / 2) * 4.2) / 4, abs(p));
	//	return (1 - cos (v - fmod(v, _PI/4.0)))/2.0;
}
double WindowTriFoil(double u, double a, double d, double p, bool bFrise = false)
{

	//0<b<10
	double v = fabs((2 * mathfmod(a * d * u, 2 * _PI) / _PI) - 2);
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);

	double b = 0.2 * Parameter;


	if (v < b)
		return pow(4 * (0.75 - 0.25 * (b - sqrt(b * b - v * v))) / 3, p);
	else if (v < 1)
		return pow(4 * (0.75 - 0.25 * (1 - sqrt((1 - b) * (1 - b) - (v - b) * (v - b)))) / 3, p);
	else if (v < 2 - b)
		return pow((2 / 3) * (1 - 2 * sqrt(1.0 / 4.0 - (v - 1.5) * (v - 1.5))), p);
	else
		return 0;

	double alpha = 0.5;

	if (v < 0.5)
		return pow(4 * (0.75 - 0.25 * (alpha - sqrt(alpha * alpha - v * v))) / 3, p);
	else if (v < 1)
		return pow(4 * (0.75 - 0.25 * (1 - sqrt((1 - alpha) * (1 - alpha) - (v - alpha) * (v - alpha)))) / 3, p);
	else if (v < 1.5)
		return pow(2.0 / 3.0, p) * (-2.0 * sqrt(1.0 / 4.0 - (v - 1.5) * (v - 1.5)) + 1);
	else return 0.0;
}
double WindowPodium(double u, double a, double d, double p, bool bFrise = false)
{
	double v = mathfmod(a * d * u, 2 * _PI);
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);

	double b = round(3 * Parameter) + 1.2;

	double x = floor(b * fabs(sin(v / 2))) / floor(b);

	return pow(x, abs(p));

	return pow((1.0 - cos((u * a * d + _PI / 6.0) - fmod(100 * _PI + (u * a * d), _PI / 3.0))) / 2.0, abs(p));
}
double WindowMoscow(double u, double a, double d, double p, bool bFrise = false)
{
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);

	double b = Parameter;

	if (b < .0001)
		b = 0.0001;

	return pow(-(pow(1.0 - b * abs(sin(a * d * u / 2.0)), 3) - 1.0) / 2.0, p) / pow((-(pow(1.0 - b, 3) - 1.0) / 2.0), p);



	return pow((-(pow(1.0 - 2.0 * fabs(sin(a * d * u / 2.0)), 3.0) - 1.0) / 2.0), abs(p));
}
double WindowPointedSteps(double u, double a, double d, double p, bool bFrise = false)
{
	double x, y, z, v;
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);

	double b = Parameter;

	v = mathfmod(a * d * u, 2 * _PI);

	if (v < _PI)
		x = b * v / _PI;
	else
		x = 2 * b - b * v / _PI;

	if (x - floor(x) < 0.5)
		y = floor(x);
	else
		y = (2 * x - floor(x) - 1);

	if (b - floor(b) < 0.5)
		z = floor(b);
	else
		z = (2 * b - floor(b) - 1);

	return pow(y / z, p);

	x = mathfmod(a * d * u, 2 * _PI);

	if (x < _PI)
		y = 3 * x / _PI;
	else
		y = 6 - 3 * x / _PI;
	z = floor(y);

	if (y - z < 0.5)
		return pow(z / 3, abs(p));
	else
		return pow((2 * y - z - 1) / 3, abs(p));
}
double WindowCathedral(double u, double a, double d, double p, bool bFrise = false)
{
	double x = mathfmod(u * a * d, 2 * _PI);
	double y;
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);

	//-pi/2<c<pi
	double b = 3 * _PI * Parameter / 20 - _PI / 2;



	if (x < b)
		y = x * x * pow(sqrt(2) / 2, p / 4) / (b * b);
	else if (x > 2 * _PI - b)
		y = (pow(sqrt(2) / 2, p / 4) / (b * b)) * (2 * _PI - x) * (2 * _PI - x);
	else
		y = pow(fabs(sin(x / 2)), p / 4) - pow(fabs(sin(b / 2)), p / 4) + pow(sqrt(2) / 2, p / 4);

	return y / 1 - pow(fabs(sin(b / 2)), p / 4) + pow(sqrt(2) / 2, p / 4);


	double c = _PI / 4;

	if (x < c)
		return (pow(sqrt(2) / 2, abs(p) / 4) / (c * c)) * x * x;
	if (x > 2 * _PI - c)
		return (pow(sqrt(2) / 2, abs(p) / 4) / (c * c)) * (2 * _PI - x) * (2 * _PI - x);
	else
		return pow(sqrt(sqrt(fabs(sin(u * a * d / 2)))), abs(p));
}
double WindowOriental(double u, double a, double d, double p, bool bFrise = false)
{
	double v = mathfmod(a * d * u, 2 * _PI);

	//	return pow((2/_PI)*sqrt(fabs(_PI*_PI/4 - (v-_PI)*(v-_PI)))*floor(fabs((sin(v/2)))*sqrt((float)2)), abs(p));
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);


	//0<b<10
	double b = Parameter;

	if (b < .32)
		b = 0.32;


	return pow(sin(sqrt(fmax(_PI - b * (_PI - v) * (_PI - v), 0))), p);

	//1<=b<=100000
	return pow(sqrt((abs(_PI * _PI - b * (v - _PI) * (v - _PI)))), p) / pow(fmax(b * abs(v - _PI), _PI), p);

}
double WindowMiddleAge(double u, double a, double d, double p, bool bFrise = false)
{

	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);

	double b = 1.0 + Parameter / 10;
	if (b > 1.99999)
		b = 1.99999;




	//	return pow(fabs(sin(a*d*u/2))*floor(fabs((sin(u*a*d/2)))*sqrt((float)3)), p);


		//1 < b < 1.9999
	return pow(fabs(sin(a * d * u / 2)) * floor(fabs((sin(u * a * d / 2))) * b), abs(p));
}
double WindowCurly(double u, double a, double d, double p, bool bFrise = false)
{
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);

	double b = -10.0 + 2 * Parameter;

	double y;
	double x = pow(fabs(4 * fabs(sin(a * d * u / 2)) - 2 * b * (1 - cos(a * d * u))), p);
	if (b > .5)
		y = fmax(fabs(4 - 4 * b), 1 / b);
	else
		y = fabs(4 - 4 * b);

	return x / pow(y, p);

	return pow((4 * fabs(sin(a * d * u / 2)) - 2 * (1 - cos(a * d * u))), abs(p));
}
double WindowPointedSquare(double u, double a, double d, double p, bool bFrise = false)
{
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);
	double b = 1.51 + Parameter / 10;


	return pow(floor(fabs(sin(u * a * d / 2)) * b) * (1 - abs(mathfmod(u, 2 * _PI / (a * d)) - _PI / (a * d)) * a * d / _PI), p);
	return pow(floor(mathfmod((u + _PI / (2 * a * d)) * a * d / _PI, 2)) * (1 - fabs(mathfmod(u, 2 * _PI / (a * d)) - _PI / (a * d)) * a * d / _PI), abs(p));
}
double WindowM(double u, double a, double d, double p, bool bFrise = false) //windowM
{
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);
	double b = Parameter / 10;
	if (b < .0001)
		b = 0.0001;
	double v = abs((b - fabs(fabs(sin(a * d * u / 2)) - b) / b));
	double y = fmax(1, fabs((b - fabs(1 - b) / b)));

	return pow(v / y, p);

	return (pow(((.9 - fabs(fabs(sin(a * d * u / 2)) - .9)) / .9), fabs(p)));

	//	p=3, a=0.9
}
double WindowMexicanHat(double u, double a, double d, double p, bool bFrise = false)
{
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);
	double b = 2 * Parameter / 10;

	double v = (double)myMod((float)(a * d * u), 2.0 * _PI);

	return fabs(sin(a * d * u / 2)) * pow(sqrt(fmax(0, (_PI * _PI - b * (v - _PI) * (v - _PI)))), p) / pow(fmax(b * fabs(v - _PI), _PI), p);

}
double WindowAngkorWat(double u, double a, double d, double p, bool bFrise = false)
{
	double Parameter = (bFrise ? myDome.fFriseParameter : myDome.fWindowParameter);
	double b = 1 + 9 * Parameter / 10;
	double x = mathfmod(a * d * u, 2 * _PI);
	if (x < _PI / 2 || x>3 * _PI / 2)
		return 0;
	else
	{
		double y = 2 * (2.5 + (((fabs(cos(x))) + (b - (fabs(cos(x + _PI / 2)))) * 2) / (2 + fabs(cos(4 * x + _PI / 2)) * 8))) / (3 + b) - 1;
		return pow(y, p);
	}
}

double ShrinkNoneFunction(double u, double fShrinking, double fSpeed)
{
	return 1.0;
}
double ShrinkLinearFunction(double u, double fShrinking, double fSpeed)
{
	double a = abs(fSpeed);

	if (fSpeed >= 0)
		u = 1 - u;
	return abs(a * u - 1) + (fShrinking - abs(a - 1)) * u;
}
double ShrinkSquareRootFunction(double u, double fShrinking, double fSpeed)
{
	double a = abs(fSpeed);
	if (fSpeed >= 0)
		u = 1 - u;
	double k = a * sqrt(u) + (1 - a) * u;

	return abs((fShrinking - 1) * k + 1);

}
double ShrinkSquareFunction(double u, double fShrinking, double fSpeed)
{
	double a = abs(fSpeed);
	if (fSpeed >= 0)
		u = 1 - u;
	double k = a * u * u + (1 - a) * u;

	return abs((fShrinking - 1) * k + 1);

}





unsigned short STL_RGB(unsigned short red, unsigned short green, unsigned short blue)
{
	return red + 32 * green + 1024 * blue;
}

void GetFunctions(int nSlope)
{
	for (int i = 0; i < nbTestSlopes; i++)
	{
		if (SlopeFunctionListe[nSlope].FunctionXY == TestSlopeFunctionListe[i].Function)
			lbTestSlopeXY->set_int_val(i);
		if (SlopeFunctionListe[nSlope].FunctionZ == TestSlopeFunctionListe[i].Function)
			lbTestSlopeZ->set_int_val(i);

	}
}
