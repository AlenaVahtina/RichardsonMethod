#ifndef BASEMATRIX_H
#define BASEMATRIX_H

#include <iostream>
using namespace std;
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>

class BaseMatrix
{
public:
    BaseMatrix();
    enum TYPE{NO_TYPE, NORMALMATRIX, CRSMATRIX};
    virtual int getType(){
        return NO_TYPE;
    }


    //чтение матрицы с клавиатуры
    virtual void readMatrix(){ }


    //создание матрицы B из A
    virtual BaseMatrix * createB(){ }


    //преобразование матрицы B^(-1/2)
    virtual BaseMatrix * reB(){ }


    //функция умножения матриц (создание матрицы С)
    virtual BaseMatrix * matrixMatrix (BaseMatrix * secondM) { }


    //функция обращения матрицы(-A)
    virtual BaseMatrix * minesMatrex(){ }

    //создание матрицы С
    virtual BaseMatrix * createC(){ }


    //функция умножения матрицы на вектор
    virtual vector <double> multMatrixVector(vector <double> y){ }


    //функция вывода матрицы
    virtual void writeMatrix (){ }

    //взять элемент
    virtual double getElement(int row, int column) const{ }
};

#endif // BASEMATRIX_H
