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

const double EPSELON_ERROU=0.5;
const double EPSELON_SOLUTION=0.000001;

#ifndef M_PI
    const double M_PI=3,1415926535897932384626;
#endif

class RichardsonMethod
{
public:
    RichardsonMethod(){};
    void setQ(double _q){q=_q;}
    void setP(double _p){p=_p;}
    void setA(double _a){a=_a;}
    void setB(double _b){b=_b;}
    void setS(double _iterationNomber){iterationNomber=_iterationNomber;}
    vector<double> getErrors(){return deltak;}

    void computeResultVectorForE (vector<double> &y, BaseMatrix *SLAU, vector<double> f,int fold, double gamma1=8, double gamma2=40000);

    void computeResultVectorForC (vector<double> &y, BaseMatrix *SLAU, vector<double> f,int fold, double gamma1=8, double gamma2=40000);

    void computeResultVectorForEWithRivalProcess(vector<double> &y, BaseMatrix *SLAU,vector<double> f,int fold, double gamma11=4000, double gamma12=800, double gamma2=40000);

    void computeResultVectorForNotEWithRivalProcess(vector<double> &y, BaseMatrix *SLAU,vector<double> f,int fold, double gamma11=4000, double gamma12=800, double gamma2=40000);

    void  calculate (vector<double> &y, BaseMatrix *SLAU, vector<double> f,int fold, double gamma1, double gamma2, vector<double> &deltak, bool matrixType=true, bool processType=false, int iterationNomberFrom=0, int iterationNomberTo=0);


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
    double divider;
};

#endif // RICHARDSONMETHOD_H
