#ifndef NORMALMATRIX_H
#define NORMALMATRIX_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include "richardsonslau.h"

using namespace std;

class normalmatrix : public RichardsonSLAU
{
public:
    vector <vector <double> > Matrix;
    normalmatrix(vector <vector <double> > Matrix_);

    int getType() {
        return NORMALMATRIX;
    }

    //чтение матрицы с клавиатуры
    void ReadMatrix(int nAmountPoints){
        vector <vector <double> > Matrix;
        Matrix.resize(nAmountPoints);
        for (int i=0;i<nAmountPoints; i++){
            for (int j=0;j<nAmountPoints; j++){
                cin>>Matrix[i][j];
            }
        }
    }


    //создание матрицы B из A
    RichardsonSLAU CreateB() {
        double max;
        vector <vector <double> > MatrixB;
        for (int i=0; i<MatrixB.size();i++){
            max=Matrix[i][0];
            for (int j=0; j<Matrix.size(); j++){
                if (Matrix[i][j]!=0) {
                    if (Matrix[i][j]>max){max=Matrix[i][j];}
                }
            }
            if (Matrix[i].empty()) {std::cout<<"Error. In your matrix you have line with all 0 "<<endl;}
            for (int j=0; j<MatrixB.size(); j++){
                if ((i==j) && (Matrix[i][j]!=0) && (!Matrix[i].empty())) {
                    MatrixB[i][j]=Matrix[i][j];
                }
                if ((i==j) && (Matrix[i][j]==0) && (!Matrix[i].empty())) {
                    MatrixB[i][j]=max;
                }
            }
        }

        return normalmatrix (MatrixB);
    }


    //преобразование матрицы^(-1/2)
    void ReB(){
        for (int i=0; i<Matrix.size();i++){
            for (int j=0; j<Matrix.size(); j++){
                if (i==j) {Matrix[i][j]=pow(Matrix[i][j],-0.5);}
                }
            }
    }


    //функция умножения матриц
    void MatrixMatrix (vector<vector<double> > &Matrix, vector<vector<double> > & MatrixB, vector<vector<double> > &MatrixC){
        int nAmountPoints=Matrix.size();
        for (int i=0; i<nAmountPoints;i++){
            for (int j=0; j<nAmountPoints; j++){
                MatrixC[i][j]=0;
                for (int l=0; l<nAmountPoints; l++){
               MatrixC[i][j]+=MatrixB[i][l]*Matrix[l][j];
                }
            }
         }
    }


    //функция создание матрицы С
    void CreatC (vector<vector<double> > &Matrix, vector<vector<double> > & MatrixB, vector<vector<double> > &MatrixC) {
        int nAmountPoints=Matrix.size();
        vector<vector<double> > MatrixAB;
        MatrixAB.resize(nAmountPoints);
        for (int i=0; i<nAmountPoints; i++){
             MatrixAB[i].resize(nAmountPoints);
        }

        MatrixMatrix(MatrixB, Matrix, MatrixAB);
        MatrixMatrix(MatrixAB, MatrixB, MatrixC);

    }


    //функция обращения матрицы(-A)
    void MinesMatrex(){
        for (int i=0; i<Matrix.size(); i++){
            for (int j=0; j<Matrix.size(); j++){
                Matrix[i][j]*=-1;
            }
        }
    }


    //функция умножения матрицы на вектор
    void MultMatrixVector(vector<double> &y, vector<vector<double> > &Matrix, vector <double> &MultVector){
        for (int i=0; i<y.size(); i++){
            MultVector[i]=0;
            for (int j=0; j<y.size(); j++){
                MultVector[i]+=Matrix[i][j]*y[j];
            }
        }
    }


    //функция вывода матрицы
    void WriteMatrix (vector<vector<double> > &Matrix){
        for (int i=0; i<Matrix.size(); i++){
            for (int j=0; j<Matrix.size(); j++){
                cout<<Matrix[i][j]<<"  ";
            }
            cout<<endl;
        }
    }
};

#endif // NORMALMATRIX_H
