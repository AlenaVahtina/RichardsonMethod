#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include "plots.h"

#include "normalmatrix.h"
#include "crsmatrix.h"

using namespace std;

#ifndef M_PI
    const double M_PI=3,1415926535897932384626;
#endif

class Richardson
{
public:
    Richardson();
    const double EPSELON_SOLUTION=0.00000000000001;
    const double EPSELON_ERROU=0.3;
    void setA(double _a){a=_a;}
    void setB(double _b){b=_b;}
    void setS(double _iterationNomber){iterationNomber=_iterationNomber;}
    vector<double> getErrors(){return deltak;}


    //функция расчета итогового значения вектора у для еденичной матрицы и без конкурирующих процессов
    void computeResultVectorForE (vector<double> &y, RichardsonSLAU *SLAU, vector<double> f,int fold);

    //функция расчета гамма1* гамма1** и гамма2
    void gammacalculation(double& gamma11,double &gamma12, double & gamma2,RichardsonSLAU *SLAU);


    //вспомогательная функция расчета значения у для конкурирующх процессов
    void supportingComputeResultVector (vector<double> &y, RichardsonSLAU *SLAU,vector<double> f,int fold, vector<double> &deltak,int iterationNomber,double gamma1,double gamma2,int nAmountPoints,double step,int ya,int yb);


    //функция расчета итогового значенияя вектора с канкурирующими процессами для единичной матрицы
    void computeResultVectorForEWithRivalProcess(vector<double> &y, RichardsonSLAU *SLAU,vector<double> f,int fold);


    //функция расчета итогового значенияя вектора с канкурирующими процессами для не единичной матрицы
    void computeResultVectorForNotEWithRivalProcess(vector<double> &y, RichardsonSLAU *SLAU,vector<double> f,int fold);


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


    //функция сортировки
    void Rsort(vector<int> &index);


    //сделать массив лямбда
    void flambda (vector <int> &index, vector <double> &lambda);


    //функция итерационной ошибки
    double IterError(vector<double> &y, vector<double> &oldy);

};


#endif // FUNCTIONS_H
