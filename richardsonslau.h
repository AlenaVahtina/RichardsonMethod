#ifndef RICHARDSONSLAU_H
#define RICHARDSONSLAU_H

#include <iostream>
using namespace std;
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>

class RichardsonSLAU
{
public:
    RichardsonSLAU();
    enum TYPE{NO_TYPE, NORMALMATRIX, CRSMATRIX};
    virtual int getType(){
        return NO_TYPE;
    }


    //чтение матрицы с клавиатуры
    virtual void ReadMatrix(){ }


    //создание матрицы B из A
    virtual RichardsonSLAU * CreateB() { }


    //преобразование матрицы B^(-1/2)
    virtual RichardsonSLAU * ReB(){ }


    //функция умножения матриц (создание матрицы С)
    virtual RichardsonSLAU * MatrixMatrix (RichardsonSLAU * secondM) { }


    //функция обращения матрицы(-A)
    virtual RichardsonSLAU * MinesMatrex(){ }

    //создание матрицы С
    virtual RichardsonSLAU * CreateC() { }


    //функция умножения матрицы на вектор
    virtual vector <double> MultMatrixVector(vector <double> y){ }


    //функция вывода матрицы
    virtual void WriteMatrix (){ }
};

#endif // RICHARDSONSLAU_H
