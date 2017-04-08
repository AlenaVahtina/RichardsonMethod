#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include "plots.h"

#include "matrix/normalmatrix.h"
#include "matrix/crsmatrix.h"

using namespace std;

#ifndef M_PI
    const double M_PI=3,1415926535897932384626;
#endif

class Common
{
public:

//    const double EPSELON_SOLUTION=0.00000000000001;
//    const double EPSELON_ERROU=0.03;


    //функция расчета гамма1* гамма1** и гамма2
    static void gammacalculation(double& gamma11,double &gamma12, double & gamma2,BaseMatrix *SLAU, double p, double q, int nAmountPoints);

    //функция расчета гамма1 и гамма2 для модуля без конкурирующих процессов
    static void gammacalculation1(double& gamma1, double & gamma2,BaseMatrix *SLAU, int nAmountPoints);


    //вспомогательная функция расчета значения у для конкурирующх процессов
    static void supportingComputeResultVector (vector<double> &y, BaseMatrix *SLAU,vector<double> f,int fold, vector<double> &deltak,
                                                int iterationNomber,double gamma1,double gamma2,int nAmountPoints,double step,int ya,int yb);


    //функция сортировки
    static void Rsort(vector<int> &index);


    //сделать массив лямбда
    static void flambda (vector <int> &index, vector <double> &lambda, double iterationNomber);


    //функция итерационной ошибки
    static double IterError(vector<double> &y, vector<double> &oldy);


    //является ли s степенью 2
    static bool fold2(int s);

    //критерий для физики
    static bool appliedCriteria(vector<double>y, vector<double>&y1, double relativE, double absolutlyE);

    //проверка на неверное количество итераций и неверные гамма
    static bool error(double gamma1, double gamma2, int iterationNomber);

    // фун. для проверки на пересечение
    static int cross(vector<double> deltak1, vector<double> deltak2);


    static  bool criterion27(vector<double> deltak1, vector<double> deltak2, int istop, double eps);

    static bool issolution(vector<double> deltak, double eps);
};


#endif // FUNCTIONS_H
