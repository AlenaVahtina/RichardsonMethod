#ifndef CRSMATRIX_H
#define CRSMATRIX_H

#include <iostream>
using namespace std;
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include "richardsonslau.h"



class crsmatrix : public RichardsonSLAU
{
public:
    crsmatrix(vector <int> pointer_, vector <int> cols_, vector <double> values_);
    crsmatrix(RichardsonSLAU * another);

    vector <int> pointer;
    vector <int> cols;
    vector <double> values;

    int getType() {
        return CRSMATRIX;
    }

    //чтение матрицы с клавиатуры в Йельском формате
    void ReadMatrix (int length, int nomberNotNullElemet);


    //создание матрицы B из А
    RichardsonSLAU *CreateB();


    //преобразование матрицы(-1/2)
    RichardsonSLAU * ReB();


    //фунция транспонирования матрицы
    RichardsonSLAU * TCRSMaatrix();


    //функция умножение матриц
    RichardsonSLAU *MatrixMatrix (RichardsonSLAU * secondM);


    //функция создания матрицы С=B^(-1/2)*A*B^(-1/2)
    RichardsonSLAU *CreateC();

    //функция обращения матрицы (-А)
    RichardsonSLAU * MinesMatrex();


    //функция умножения матрицы на вектор
    vector <double> MultMatrixVector(vector <double> y);


    //функция вывода матрицы в Йельском формате
    void WriteMatrix ();


};

#endif // CRSMATRIX_H
