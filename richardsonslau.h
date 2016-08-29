#ifndef RICHARDSONSLAU_H
#define RICHARDSONSLAU_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>

using namespace std;
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


 private:
    //перевод из Йельского формата в обычный
    void TranslateCRSNormal(vector<vector<double> > &Matrix, vector <double> values,vector <int> cols, vector <int>  pointer )
    {
        int N=pointer.size();
        int colon;
        for (int i=1; i<=N; i++){
            for (int j=pointer[i-1]; j<pointer[i]; j++){
                colon=cols[j];
                Matrix[i-1][colon]=values[j];
            }
        }
    }


    //перевод из обычного формата в Йельский
    void TranslateNormalCRS(vector<vector<double> > Matrix, vector <double>  &values,vector <int> &cols, vector <int>  &pointer )
    {
        int k=0;
        int numb;
        int N=Matrix[0].size();
        pointer.resize(N+1);
        pointer[0]=0;
        for (int i=0; i<N;i++){
            if (i!=0) {pointer[i]=pointer[i-1]+numb;}
            numb=0;
            for (int j=0; j<N; j++){
                if (Matrix[i][j]!=0){
                    k++;
                    numb++;
                    values.resize(k);
                    values[k-1]=Matrix[i][j];
                    cols.resize(k);
                    cols[k-1]=j;
                }
            }
        }
        pointer[N]=pointer[N-1]+numb;
      }

};

#endif // RICHARDSONSLAU_H
