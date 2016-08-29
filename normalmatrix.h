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
    RichardsonSLAU *CreateB() {
        double max;
        vector <vector <double> > MatrixB;
        MatrixB.resize(Matrix.size());
        for (int i=0; i<MatrixB.size();i++){
            MatrixB[i].resize(Matrix.size());
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

        return  new normalmatrix (MatrixB);
    }


    //преобразование матрицы^(-1/2)
    RichardsonSLAU * ReB(){
        for (int i=0; i<Matrix.size();i++){
            for (int j=0; j<Matrix.size(); j++){
                if (i==j) {Matrix[i][j]=pow(Matrix[i][j],-0.5);}
                }
            }
        return this;
    }


    //функция умножения матриц
    RichardsonSLAU *MatrixMatrix (RichardsonSLAU * secondM)
    {
        int nAmountPoints=Matrix.size();
        vector <vector <double> > resultMatrix;
        resultMatrix.resize(Matrix.size());
        normalmatrix* second=static_cast<normalmatrix*>(secondM);
        for (int i=0; i<nAmountPoints;i++){
            resultMatrix[i].resize(Matrix.size());
            for (int j=0; j<nAmountPoints; j++){
                resultMatrix[i][j]=0;
                for (int l=0; l<nAmountPoints; l++){
                resultMatrix[i][j]+=second->Matrix[i][l] * this->Matrix[l][j];
                }
            }
         }
        return new normalmatrix(resultMatrix);
    }


    //функция создания матрицы С=B^(-1/2)*A*B^(-1/2)
    RichardsonSLAU *CreateC(){
        RichardsonSLAU *MatrixB=this->CreateB()->ReB();
        RichardsonSLAU *MatrixC=this->CreateB()->ReB();
        MatrixC->MatrixMatrix(static_cast<RichardsonSLAU*>(this))->MatrixMatrix(MatrixB);
        return MatrixC;
    }

    //функция обращения матрицы(-A)
    RichardsonSLAU * MinesMatrex(){
        for (int i=0; i<Matrix.size(); i++){
            for (int j=0; j<Matrix.size(); j++){
                Matrix[i][j]*=-1;
            }
        }
        return this;
    }


    //функция умножения матрицы на вектор
    vector <double> MultMatrixVector( vector <double> y){
        vector <double> MultVector;
        MultVector.resize(y.size());
        for (int i=0; i<y.size(); i++){
            MultVector[i]=0;
            for (int j=0; j<y.size(); j++){
                MultVector[i]+=Matrix[i][j]*y[j];
            }
        }
        return MultVector;
    }


    //функция вывода матрицы
    void WriteMatrix (){
        for (int i=0; i<Matrix.size(); i++){
            for (int j=0; j<Matrix.size(); j++){
                cout<<Matrix[i][j]<<"  ";
            }
            cout<<endl;
        }
    }
};

#endif // NORMALMATRIX_H
