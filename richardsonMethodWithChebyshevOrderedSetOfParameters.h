#ifndef RICHARDSONMETHOD_H
#define RICHARDSONMETHOD_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>

#include "plots.h"
#include "matrix/normalmatrix.h"
#include "matrix/crsmatrix.h"
#include "common.h"

using namespace std;

#ifndef M_PI
    const double M_PI=3,1415926535897932384626;
#endif

class RichardsonMethod
{
    RichardsonMethod(){};
public:

    //функция расчета итогового значения вектора у для еденичной матрицы и без конкурирующих процессов
    void computeResultVectorForE (vector<double> &y, BaseMatrix *SLAU, vector<double> f,int fold);

    void computeResultVectorForC (vector<double> &y, BaseMatrix *SLAU, vector<double> f,int fold);


private:
    double a,b,ya,yb,step;
    int nAmountPoints; //чтение колличества ячеек
    double iterationNomber;//колличесство итераций
    vector<double> f;//функция, правая часть
    vector<double> y;//искомое расспределение темпиратур, вектор неизвестных
    vector<vector<double> > Matrix;//матрица оператора A, которую заполняет пользователь
    vector<double> MultVector;//результат умножения A' на y
    vector<int> index;
    vector<double> lambda;
    vector<double> tao;
    vector<double> tao1;
    vector<double> tao2;
    vector<double> deltak;
    vector<double> deltak1;
    vector<double> deltak2;

    double gamma1, gamma2, p0, tao0, tao01, tao02, p, q, gamma12, gamma11, p01, p02;
};

#endif // RICHARDSONMETHOD_H
