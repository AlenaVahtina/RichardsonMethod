#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include "plots.h"

using namespace std;

class Richardson
{
public:
    Richardson();
    void setA(double _a){a=_a;}
    void setB(double _b){b=_b;}
    void setS(double _s){s=_s;}
    vector<double> getErrors(){return deltak;}
    vector<double> VectorB;
    //функция расчета итогового значения
    void ItartionR (vector<double> &y, vector<vector<double> > &Matrix,vector<double> f,int kr){

        vector<double> oldy=y;
        //граничные условия, колличество ячеек и шаг
        N=y.size();
        h=(b-a)/N;
        ya=1;
        yb=0;

        //гамма 1 и гамма 2 (26) границы спектра
        gamma1=4*sin(M_PI*h/2*(b-a))*sin(M_PI*h/2*(b-a))/(h*h);
        gamma2=4*cos(M_PI*h/2*(b-a))*cos(M_PI*h/2*(b-a))/(h*h);
        //посчитать p0
        p0=(1-gamma1/gamma2)/(1+gamma1/gamma2);

        flambda(index, lambda);

        //заполнить масив тао
        tao0=2/(gamma1+gamma2);
        tao.resize(s);
        for (int i=0; i<s; i++){
            tao[i]=tao0/(1+p0*lambda[i]);
        }
        MinesMatrex(Matrix);
        f[1]=ya/(h*h);
        f[N-1]=yb/(h*h);
        s=tao.size();
        Plots plot;
        deltak.resize(s);
        MultVector.resize(N);
        for (int i=0; i<s; i++){
            MultMatrixVector(y, Matrix);
            for (int j=1; j<N; j++){
//                std::cout<<"nomer "<<i<<std::endl;
//                std::cout<<tao[i]<<std::endl;
//                std::cout<<f[j]<<std::endl;
//                std::cout<<MultVector[j]/(h*h)<<std::endl;
                y[j]+=tao[i]*(f[j]-MultVector[j]);
//                std::cout<<y[j]<<"  ";
            }
//            cout<<endl;
            //Yplot(y,"output"+atoi(i)+".png")
            if(i%kr==0)//условие сюда ,а там цикл убрать
                plot.IteratPlot(y,std::string("output")+std::to_string(i)+".png",std::string("output")+std::to_string(i)+".dat");
            deltak[i]=IterError(y,oldy);
        }

        y[0]=ya;
        y[N-1]=yb;
    }


    void ItartionRWithGer (vector<double> &y, vector<vector<double> > &Matrix,vector<double> f,int kr){

        vector<double> oldy=y;
        vector<vector<double> >MatrixB;
        vector<vector<double> >MatrixC;
        //граничные условия, колличество ячеек и шаг
        N=y.size();
        h=(b-a)/N;
        ya=1;
        yb=0;

        //-A
        MinesMatrex(Matrix);
        std::cout<<"mitrix\n";
        for (int i=0; i<N;i++){
            for (int j=0; j<N; j++){
                std::cout<<Matrix[i][j]<<" ";
            }
            std::cout<<endl;
        }
        std::cout<<std::endl;

        //работа с матрицей B
        MatrixB.resize(N);
        for (int i=0; i<N; i++){
             MatrixB[i].resize(N);
        }

        CreateB(Matrix, MatrixB);
        //вывод B
        for (int i=0; i<N;i++){
            for (int j=0; j<N; j++){
                std::cout<<MatrixB[i][j]<<" ";
            }
            std::cout<<endl;
        }
        std::cout<<endl;
        ReB(MatrixB);

     //вывод B
        for (int i=0; i<N;i++){
            for (int j=0; j<N; j++){
                std::cout<<MatrixB[i][j]<<"  ";
            }
            std::cout<<endl;
        }

        //работа с матрицей C
            MatrixC.resize(N);
            for (int i=0; i<N; i++){
                MatrixC[i].resize(N);
            }
            
            MatrixMatrix(Matrix,MatrixB, MatrixC);
            //выводC
            for (int i=0; i<N;i++){
                for (int j=0; j<N; j++){
                    std::cout<<MatrixC[i][j]<<"   ";
                }
                std::cout<<endl;
            }

        //гамма 1(26) верхняя граница спектра
        gamma1=4*sin(M_PI*h/2*(b-a))*sin(M_PI*h/2*(b-a))/(h*h);
//      gamma2=4*cos(M_PI*h/2*(b-a))*cos(M_PI*h/2*(b-a))/(h*h);


       //расчет нижней границы с помощью кругов Герщгорина
        vector<double> R;
        R.resize(N);
        for (int i=0; i<N;i++){
            for (int j=0; j<N; j++){
                if (i!=j) {R[i]+=MatrixC[i][j];}
            }
        }

        for (int i=0; i<=N;i++){
            if (gamma2<fabs(MatrixC[i][i]+R[i])) {
                    gamma2=fabs(MatrixC[i][i]+R[i]);
            }
       }
        gamma2=gamma2/(h*h);
        std::cout<<gamma2<<"UI"<<endl;

        //посчитать p0
        p0=(1-gamma1/gamma2)/(1+gamma1/gamma2);

        flambda(index, lambda);

        //заполнить масив тао
        tao0=2/(gamma1+gamma2);
        tao.resize(s);
        for (int i=0; i<s; i++){
            tao[i]=tao0/(1+p0*lambda[i]);
        }

        f[1]=ya/(h*h);
        f[N-1]=yb/(h*h);
        s=tao.size();
        Plots plot;
        deltak.resize(s);
        MultVector.resize(N);
        for (int i=0; i<s; i++){
            MultMatrixVector(y, Matrix);
            for (int j=0; j<N; j++){
//                std::cout<<"nomer "<<i<<std::endl;
//                std::cout<<tao[i]<<std::endl;
//                std::cout<<f[j]<<std::endl;
//                std::cout<<MultVector[j]/(h*h)<<std::endl;
                y[j]+=tao[i]*(f[j]-MultVector[j]/(h*h));
//                std::cout<<y[j]<<"  ";
            }
//            cout<<endl;
            //Yplot(y,"output"+atoi(i)+".png")
            if(i%kr==0)//условие сюда ,а там цикл убрать
                plot.IteratPlot(y,std::string("output")+std::to_string(i)+".png",std::string("output")+std::to_string(i)+".dat");
            deltak[i]=IterError(y,oldy);
        }

        y[0]=ya;
        y[N-1]=yb;
    }



    //создание матрицы B
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
                if ((i==j) and (Matrix[i][j]!=0) and (!Matrix[i].empty())) {
                    MatrixB[i][j]=Matrix[i][j];
                }
                if ((i==j) and (Matrix[i][j]==0) and (!Matrix[i].empty())) {
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


    //функция умножения матриц
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


private:
    double a,b,ya,yb,h;
    int N; //чтение колличества ячеек
    double s;//колличесство итераций
    vector<double> f;//функция, правая часть
    vector<double> y;//искомое расспределение темпиратур, вектор неизвестных
    vector<vector<double> > Matrix;//матрица оператора A, которую заполняет пользователь
    vector<double> MultVector;//результат умножения A' на y
    vector<int> index;
    vector<double> lambda;
    vector<double> tao;
    vector<double> deltak;

    double gamma1, gamma2, p0, tao0;


    void MinesMatrex(vector<vector<double> > &Matrix){
        for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                Matrix[i][j]*=-1;
            }
        }
    }


    //функция сортировки
    void Rsort(vector<int> &index){
        for (int power=1; power<=log(index.size())/log(2); power++){
            for (int k=pow(2,(power-1))-1;k>=0; k--){
                index[2*k]=index[k];
                index[2*k+1]=pow(2,(power+1))-index[2*k];
            }
        }
    }


    //сделать массив лямбда
    void flambda (vector <int> &index, vector <double> &lambda){
        //заполнить массив index
        index.resize(s);
        for (int i=0; i<s; i++) {
            index[i]=i+1;
        }
        Rsort(index);

        lambda.resize(s);
        for (int i=0; i<s; i++){
            lambda[i]=-cos(M_PI*index[i]/(2*s));
        }
    }


    //функция умножения матрицы A' на у
    void MultMatrixVector(vector<double> &y, vector<vector<double> > &Matrix){
        for (int i=0; i<y.size(); i++){
            MultVector[i]=0;
            for (int j=0; j<y.size(); j++){
                MultVector[i]+=Matrix[i][j]*y[j];
            }
        }
    }


    //функция итерационной ошибки
    double IterError(vector<double> &y, vector<double> &oldy) {
        double maxk;
        for (int i=1; i<=y.size(); i++){
            if (fabs(y[i]-oldy[i])>maxk) {maxk=fabs(y[i]-oldy[i]);}
        }
        oldy=y;
        return maxk;
    }

};


#endif // FUNCTIONS_H
