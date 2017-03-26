#ifndef CRSMATRIX_H
#define CRSMATRIX_H

#include <iostream>
using namespace std;
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include "basematrix.h"



class CrsMatrix : public BaseMatrix
{
public:
    CrsMatrix(vector <int> pointer_, vector <int> cols_, vector <double> values_);
    CrsMatrix(BaseMatrix * another);

    vector <int> pointer;
    vector <int> cols;
    vector <double> values;

    int getType() {
        return CRSMATRIX;
    }

    //чтение матрицы с клавиатуры в Йельском формате
    void readMatrix (int length, int nomberNotNullElemet);


    //чтение матрицы из файла в Йельском формате
    void readMatrixFile(int length, int nomberNotNullElemet);


    //создание матрицы B из А
    BaseMatrix *createB();


    //преобразование матрицы(-1/2)
    BaseMatrix * reB();


    //фунция транспонирования матрицы
    BaseMatrix * TCRSMaatrix();


    //функция умножение матриц
    BaseMatrix *matrixMatrix (BaseMatrix * secondM);


    //функция создания матрицы С=B^(-1/2)*A*B^(-1/2)
    BaseMatrix *createC();

    //функция обращения матрицы (-А)
    BaseMatrix * minesMatrex();


    //функция умножения матрицы на вектор
    vector <double> multMatrixVector(vector <double> y);


    //функция вывода матрицы в Йельском формате
    void writeMatrix ();


    //взять элемент
    double getElement (int row, int column) const;
};

#endif // CRSMATRIX_H
