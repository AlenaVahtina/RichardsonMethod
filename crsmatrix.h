#ifndef CRSMATRIX_H
#define CRSMATRIX_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include "richardsonslau.h"

using namespace std;

class crsmatrix : public RichardsonSLAU
{
public:
    crsmatrix(vector <int> pointer_, vector <int> cols_, vector <double> values_);

    vector <int> pointer;
    vector <int> cols;
    vector <double> values;

    int getType() {
        return CRSMATRIX;
    }

    //чтение матрицы с клавиатуры в Йельском формате
    void ReadMatrix (int length, int nomberNotNullElemet){

        pointer.resize(length);
        cout<<'\n'<<"pointer"<<'\n';
        for (int i=0; i<length;i++){
            cin>>pointer[i];
        }

        values.resize(nomberNotNullElemet);
        cout<<'\n'<<"values"<<'\n';
        for (int i=0; i<nomberNotNullElemet;i++){
            cin>>values[i];
          }

        cols.resize(nomberNotNullElemet);
        cout<<'\n'<<"cols"<<'\n';
        for (int i=0; i<nomberNotNullElemet;i++){
            cin>>cols[i];
        }
    }


    //создание матрицы B из А
    RichardsonSLAU *CreateB()
    {
     vector <double> valuesB;
     vector <int> colsB;
     vector <int>  pointerB;
     int point=0;
     pointerB.push_back(0);
     for (int i=0; i<pointer.size()-1; i++){
        point=0;
        for (int j=pointer[i]; j<pointer[i+1]; j++){
            if (cols[j]==i){
                    valuesB.push_back(values[j]);
                    colsB.push_back(cols[j]);
                    point++;
            }
        }
        if (point!=0){pointerB.push_back(pointerB[i]+point);}
     }
     return new crsmatrix (pointerB, colsB, valuesB);
    }


    //преобразование матрицы(-1/2)
    RichardsonSLAU * ReB()
    {
     for (int i=0; i<values.size(); i++){
        values[i]=1/sqrt(values[i]);
     }
     return this;
    }


    //фунция транспонирования матрицы
    RichardsonSLAU * TCRSMaatrix()
    {
        vector <vector <int> > IntVectors;
        vector <vector <double> > RealVectors;
        IntVectors.resize(pointer.size()-1);
        RealVectors.resize(pointer.size()-1);
        for (int i=1; i<pointer.size(); i++){
            for (int j=pointer[i-1]; j<pointer[i]; j++){
                IntVectors[cols[j]].push_back(i-1);
                RealVectors[cols[j]].push_back(values[j]);
            }
        }

        pointer.push_back(0);
        for (int i=0; i<IntVectors.size(); i++){
            for (int j=0; j<IntVectors[i].size(); j++){

                    cols.push_back(IntVectors[i][j]);
                    values.push_back(RealVectors[i][j]);

            }
           pointer.push_back(pointer[i]+IntVectors[i].size());
        }

        for (int i=0; i<pointer.size();i++){
                cout<<pointer[i]<<"  ";
          }
        cout<<'\n';
        return this;
    }


    //функция умножение матриц
    RichardsonSLAU *MatrixMatrix (RichardsonSLAU * secondM)
    {
        vector <double>resultValues;
        vector <int> resultCols;
        vector <int> resultPointe;

        crsmatrix* second=static_cast<crsmatrix*>(secondM);
        second->TCRSMaatrix();

        double sum;
        int resultPointer;
        resultPointe.push_back(0);
        for (int i=0; i<pointer.size()-1; i++){
            resultPointer=0;
            for(int j=0; j<pointer.size()-1; j++){
                sum=0;
                for (int k=pointer[i]; k<pointer[i+1]; k++){
                    for (int l=second->pointer[j]; l<second->pointer[j+1];l++){
                        if (cols[k]==second->cols[l]){
                            sum+=values[k]*second->values[l];
                            break;
                        }
                 }
                }
                if (sum!=0){
                    resultValues.push_back(sum);
                    resultCols.push_back(j);
                    resultPointer++;
                }
            }
            resultPointe.push_back(resultPointer+resultPointe[i]);
        }
        return new crsmatrix(resultPointe,resultCols,resultValues);
    }


    //функция создания матрицы С=B^(-1/2)*A*B^(-1/2)
    RichardsonSLAU *CreateC(){
        RichardsonSLAU *MatrixB=this->CreateB()->ReB();
        RichardsonSLAU *MatrixC=this->CreateB()->ReB();
        MatrixC->MatrixMatrix(static_cast<RichardsonSLAU*>(this))->MatrixMatrix(MatrixB);
        return MatrixC;
    }

    //функция обращения матрицы (-А)
    RichardsonSLAU * MinesMatrex(){
        for (int i=0; i<values.size(); i++){
           values[i]*=-1;
        }
        return this;
    }


    //функция умножения матрицы на вектор
    vector <double> MultMatrixVector(vector <double> y)
    {
        int N=pointer.size();
        vector <double> vectorC;
        for (int i=1; i<=N; i++){
            for (int j=pointer[i-1]; j<pointer[i]; j++){
                    vectorC[i-1]+=values[j]*y[cols[j]];
            }
        }

        return vectorC;
    }


    //функция вывода матрицы в Йельском формате
    void WriteMatrix (){

        int length=pointer.size();
        cout<<'\n'<<"pointer"<<'\n';
        for (int i=0; i<length;i++){
            cout<<pointer[i]<<"  ";
        }
        cout<<'\n';

        int nomberNotNullElemet=values.size();
        cout<<'\n'<<"values"<<'\n';
        for (int i=0; i<nomberNotNullElemet;i++){
            cout<<values[i]<<"  ";
          }
        cout<<'\n';

        cout<<'\n'<<"cols"<<'\n';
        for (int i=0; i<nomberNotNullElemet;i++){
            cout<<cols[i]<<"  ";
        }
        cout<<'\n';
    }


};

#endif // CRSMATRIX_H
