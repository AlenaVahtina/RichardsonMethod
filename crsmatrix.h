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
    void ReadMatrix (int length, int nomberNotNullElemet, vector <int> &pointer, vector <int> &cols, vector <double> &values){

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
    RichardsonSLAU CreateB()
    {
     vector <double> valuesB;
     vector <int> colsB;
     vector <int>  pointerB;
     int point=0;
     pointerB.push_back(0);
     for (int i=1; i<pointer.size(); i++){
        point=0;
        for (int j=pointer[i-1]; j<pointer[i]; j++){
            if (cols[j]==i){
                    valuesB.push_back(values[j]);
                    colsB.push_back(cols[j]);
                    point++;
            }
        }
        if (point!=0){pointerB.push_back(pointerB[i-1]+point);}
     }
     return crsmatrix (pointerB, colsB, valuesB);
    }


    //преобразование матрицы(-1/2)
    void ReB()
    {
     for (int i=0; i<values.size(); i++){
        values[i]=1/sqrt(values[i]);
     }
    }


    //фунция транспонирования матрицы
    void TCRSMaatrix(vector <double> values, vector <int> cols, vector <int>  pointer,vector <double> &valuesT, vector <int> &colsT, vector <int>  &pointerT)
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

        pointerT.push_back(0);
        for (int i=0; i<IntVectors.size(); i++){
            for (int j=0; j<IntVectors[i].size(); j++){

                    colsT.push_back(IntVectors[i][j]);
                    valuesT.push_back(RealVectors[i][j]);

            }
           pointerT.push_back(pointerT[i]+IntVectors[i].size());
        }

    for (int i=0; i<pointerT.size();i++){
            cout<<pointerT[i]<<"  ";
      }
        cout<<'\n';
    }


    //функция умножение матриц
    void MatrixMatrix (vector <double> values, vector <int> cols, vector <int>  pointer, vector <double> values2, vector <int> cols2, vector <int> pointer2, vector <double> &resultValues, vector <int> &resultCols, vector <int> &resultPointe )
    {

        vector <double> valuesT;
        vector <int> colsT;
        vector <int> pointerT;
        TCRSMaatrix(values2, cols2, pointer2,valuesT, colsT, pointerT);

        double sum;
        int resultPointer;
        resultPointe.push_back(0);
        for (int i=0; i<pointer.size()-1; i++){
            resultPointer=0;
            for(int j=0; j<pointer.size()-1; j++){
                sum=0;
                for (int k=pointer[i]; k<pointer[i+1]; k++){
                    for (int l=pointerT[j]; l<pointerT[j+1];l++){
                        if (cols[k]==colsT[l]){
                            sum+=values[k]*valuesT[l];
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
    }


    //функция создания матрицы С
    void CreateC(vector <double> values, vector <int> cols, vector <int>  pointer,vector <double> valuesB, vector <int> colsB, vector <int>  pointerB, vector <double> &valuesC, vector <int> &colsC, vector <int>  &pointerC){
        ReB();
        vector <double> valuesGap;
        vector <int> colsGap;
        vector <int> pointerGap;
        MatrixMatrix(values, cols, pointer, valuesB, colsB, pointerB, valuesGap, colsGap, pointerGap);
        MatrixMatrix(valuesGap, colsGap, pointerGap, valuesB, colsB, pointerB, valuesC, colsC, pointerC);
    }

    //функция обращения матрицы (-А)
    void MinesMatrex(){
        for (int i=0; i<values.size(); i++){
           values[i]*=-1;
        }
    }


    //функция умножения матрицы на вектор
    void MultMatrixVector(vector <double> vectorB)
    {
        int N=pointer.size();
        vector <double> vectorC;
        for (int i=1; i<=N; i++){
            for (int j=pointer[i-1]; j<pointer[i]; j++){
                    vectorC[i-1]+=values[j]*vectorB[cols[j]];
            }
        }

        return vectorC;
    }


    //функция вывода матрицы в Йельском формате
    void WriteMatrix ( vector <int> pointer, vector <int> cols, vector <double> values){

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
