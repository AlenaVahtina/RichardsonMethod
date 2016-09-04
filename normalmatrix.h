#ifndef NORMALMATRIX_H
#define NORMALMATRIX_H

#include <iostream>
using namespace std;
#include <math.h>
#include <stdio.h>
#include <vector>

#include "richardsonslau.h"

class RichardsonSLAU;

class normalmatrix : public RichardsonSLAU
{
public:
    vector <vector <double> > Matrix;
    normalmatrix(vector <vector <double> > Matrix_);
    normalmatrix(RichardsonSLAU * another);

    int getType() {
        return NORMALMATRIX;
    }

    //чтение матрицы с клавиатуры
    void ReadMatrix(int nAmountPoints);


    //создание матрицы B из A
    RichardsonSLAU *CreateB();

    //преобразование матрицы^(-1/2)
    RichardsonSLAU * ReB();

    //функция умножения матриц
    RichardsonSLAU *MatrixMatrix (RichardsonSLAU * secondM);


    //функция создания матрицы С=B^(-1/2)*A*B^(-1/2)
    RichardsonSLAU *CreateC();


    //функция обращения матрицы(-A)
    RichardsonSLAU * MinesMatrex();


    //функция умножения матрицы на вектор
    vector <double> MultMatrixVector( vector <double> y);


    //функция вывода матрицы
    void WriteMatrix ();
};

#endif // NORMALMATRIX_H
