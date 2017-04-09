#ifndef NORMALMATRIX_H
#define NORMALMATRIX_H

#include <iostream>
using namespace std;
#include <math.h>
#include <stdio.h>
#include <vector>

#include "basematrix.h"

class BaseMatrix;

class NormalMatrix : public BaseMatrix
{
public:
    vector <vector <double> > Matrix;
    NormalMatrix(vector <vector <double> > Matrix_);
    NormalMatrix(BaseMatrix * another);

    int getType() {
        return NORMALMATRIX;
    }

    //чтение матрицы с клавиатуры
    void readMatrix(int nAmountPoints);


    //чтение матрицы из файла
    void readMatrixFile(int nAmountPoints, string fileName);


    //создание матрицы B из A
    BaseMatrix *createB();

    //преобразование матрицы^(-1/2)
    BaseMatrix * reB();

    //функция умножения матриц
    BaseMatrix *matrixMatrix (BaseMatrix * secondM);


    //функция создания матрицы С=B^(-1/2)*A*B^(-1/2)
    BaseMatrix *createC();


    //функция обращения матрицы(-A)
    BaseMatrix * minesMatrex();


    //функция умножения матрицы на вектор
    vector <double> multMatrixVector( vector <double> y);


    //функция вывода матрицы
    void writeMatrix ();

    //взять элемент
    double getElement (int row, int column) const;
};

#endif // NORMALMATRIX_H
