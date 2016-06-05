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


    //чтение матрицы с клавиатуры
    void ReadMatrix(int N){
        vector <vector <double> > Matrix;
        Matrix.resize(N);
        for (int i=0;i<N; i++){
            for (int j=0;j<N; j++){
                cin>>Matrix[i][j];
            }
        }
    }


    //создание матрицы B из A
    void CreateB(vector<vector<double> > &Matrix, vector<vector<double> > &MatrixB) {
        double max;
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
    }


    //преобразование матрицы B^(-1/2)
    void ReB(vector<vector<double> >&MatrixB){
        for (int i=0; i<MatrixB.size();i++){
            for (int j=0; j<MatrixB.size(); j++){
                if (i==j) {MatrixB[i][j]=pow(MatrixB[i][j],-0.5);}
                }
            }
    }


    //функция умножения матриц (создание матрицы С)
    void MatrixMatrix (vector<vector<double> > &Matrix, vector<vector<double> > & MatrixB, vector<vector<double> > &MatrixC) {
        int N=Matrix.size();
        vector<vector<double> > MatrixAB;
        MatrixAB.resize(N);
        for (int i=0; i<N; i++){
             MatrixAB[i].resize(N);
        }

        for (int i=0; i<N;i++){
            for (int j=0; j<N; j++){
                MatrixAB[i][j]=0;
                for (int l=0; l<N; l++){
               MatrixAB[i][j]+=MatrixB[i][l]*Matrix[l][j];
                }
            }
         }

        for (int i=0; i<N;i++){
            for (int j=0; j<N; j++){
                MatrixC[i][j]=0;
                for (int l=0; l<N; l++){
               MatrixC[i][j]+=MatrixAB[i][l]*MatrixB[l][j];
                }
            }
         }

    }


    //функция обращения матрицы(-A)
    void MinesMatrex(vector<vector<double> > &Matrix){
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

#endif // RICHARDSONSLAU_H
