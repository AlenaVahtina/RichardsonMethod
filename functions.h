#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include "plots.h"
#include "richardsonslau.h"

using namespace std;

#ifndef M_PI
    const double M_PI=3,1415926535897932384626;
#endif

class Richardson
{
public:
    const double epselon=0.000000000001;
    const double eps=0.3;
    Richardson();
    void setA(double _a){a=_a;}
    void setB(double _b){b=_b;}
    void setS(double _s){s=_s;}
    vector<double> getErrors(){return deltak;}
    vector<double> VectorB;


    //итерация для y без конкурирующих процессов, когда В единичная матрица
    void ItartionR (vector<double> &y, vector<vector<double> > &Matrix,vector<double> f,int kr){

        RichardsonSLAU SLAU;
        vector<double> oldy=y;
        //граничные условия, колличество ячеек и шаг
        N=y.size();
        h=(b-a)/N;
        ya=1;
        yb=0;

        //гамма 1 и гамма 2 (26) границы спектра
//        gamma1=4*sin(M_PI*h/2*(b-a))*sin(M_PI*h/2*(b-a))/(h*h);
//        gamma2=4*cos(M_PI*h/2*(b-a))*cos(M_PI*h/2*(b-a))/(h*h);

        gamma1=8/(b-a)*(b-a);
        gamma2=4/(h*h);

        //посчитать p0
        p0=(1-gamma1/gamma2)/(1+gamma1/gamma2);

        flambda(index, lambda);

        //заполнить масив тао
        tao0=2/(gamma1+gamma2);
        tao.resize(s);
        for (int i=0; i<s; i++){
            tao[i]=tao0/(1+p0*lambda[i]);
        }
        SLAU.MinesMatrex(Matrix);
        f[0]=ya/(h*h);
        f[N-1]=yb/(h*h);
        s=tao.size();
        Plots plot;
        deltak.resize(s);
        MultVector.resize(N);
        for (int i=0; i<s; i++){
            SLAU.MultMatrixVector(y, Matrix, MultVector);
            for (int j=0; j<N; j++){
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


    void ItartionR2 (vector<double> &y, vector<vector<double> > Matrix,vector<double> f,int kr, vector<double> &deltak,int s,double gamma1,double gamma2,int N,double h,int ya,int yb, vector<vector<double> > MatrixB){

        RichardsonSLAU SLAU;
        vector<double> oldy=y;

        //заполнить масив тао
        double p0=(1-gamma1/gamma2)/(1+gamma1/gamma2);

        flambda(index, lambda);

        double tao0=2/(gamma1+gamma2);
        vector<double> tao;
        tao.resize(s);
        for (int i=0; i<s; i++){
            tao[i]=tao0/(1+p0*lambda[i]);
        }

        SLAU.MinesMatrex(Matrix);
        f[0]=ya/(h*h);
        f[N-1]=yb/(h*h);

        deltak.resize(s);

        vector<double> MultVector;
        MultVector.resize(N);

        for (int i=0; i<s; i++){
            SLAU.MultMatrixVector(y, Matrix, MultVector);
            for (int j=0; j<N; j++){
                y[j]+=(tao[i]/2000)*(f[j]-MultVector[j]);
                cout<<y[j]<<"  ";
            }
             cout<"\n";
            deltak[i]=IterError(y,oldy);
        }

        y[0]=ya;
        y[N-1]=yb;
    }


    //расчет гамма 1*, 1** и 2 для конкурирующих процессов
    void gammacalculation(double& gamma11,double &gamma12, double & gamma2,vector<vector<double> > &Matrix, vector<vector<double> > &MatrixB){

        RichardsonSLAU SLAU;
//vector<vector<double> >MatrixB;
        vector<vector<double> >MatrixC;

        SLAU.MinesMatrex(Matrix);
        SLAU.CreateB(Matrix, MatrixB);
        SLAU.ReB(MatrixB);

        //работа с матрицей C
            MatrixC.resize(N);
            for (int i=0; i<N; i++){
                MatrixC[i].resize(N);
            }

        SLAU.MatrixMatrix(Matrix,MatrixB, MatrixC);

       //расчет нижней границы(гамма2) с помощью кругов Герщгорина
        vector<double> R;
        R.resize(N);
        for (int i=0; i<N;i++){
            for (int j=0; j<N; j++){
                if (i!=j) {
                    R[i]+=fabs(MatrixC[i][j]);
                }
            }
        }
        std::cout.flush();

        gamma2=0;
        for (int i=0; i<N;i++){
            if (gamma2<fabs(MatrixC[i][i]+R[i])) {
                    gamma2=fabs(MatrixC[i][i]+R[i]);
            }
       }
cout<<gamma2<<"   ";
        //q и р
        q=0.1;
        p=q/4;

        //начала S-цикла
        //гамма1*
        gamma11=q*gamma2;
        std::cout<<a<<" "<<b<<"\n";
        cout<<gamma11<<"  гамма1*  ";
        //гамма1**
        gamma12=p*gamma2;
        cout<<gamma12<<"  гамма1**  ";
    }


    void ItartionRFin(vector<double> &y, vector<vector<double> > &Matrix,vector<double> f,int kr){
        N=y.size();
        h=(b-a)/N;
        ya=1;
        yb=0;

        vector<vector<double> >MatrixB;

       MatrixB.resize(N);
            for (int i=0; i<N; i++){
                 MatrixB[i].resize(N);
            }

        std::string plname="";\

        double gamma11;
        double gamma12;
        double gamma2;
        gammacalculation(gamma11, gamma12, gamma2, Matrix, MatrixB);


        vector<double> y1=y;
        vector<double> y2=y;

        bool ko =false;
        while (true){
           ItartionR2(y1,Matrix,f,kr,deltak1,s,gamma11,gamma2,N,h,ya,yb,MatrixB);
           ItartionR2(y2,Matrix,f,kr,deltak2,s,gamma12,gamma2,N,h,ya,yb,MatrixB);

           int istop=s-1;
           if (ko)break;

           if (deltak1[istop]<deltak2[istop])
               {
               s=2*s;
               //continue;
           }else {
               if((deltak1[istop]-deltak2[istop])>eps*deltak1[istop]){
                   s=2*s;
                   gamma11=gamma11;
                   gamma12=gamma12/4;
                //   continue;
               }
               else {
                   s=2*s;
                   gamma11=gamma12;
                   gamma12=gamma12/4;

                   ko=true;
              //     continue;
               //    y1=y2;
               }
           }
           plname+="0";
           Plots p;
//           p.AveragePlotDoble(deltak1,"y1"+plname+".png","y1"+plname);
//           p.AveragePlotDoble(deltak2,"y2"+plname+".png","y2"+plname);
           p.AveragePlotDoble2(deltak1,deltak2,"y1"+plname+".png","y1"+plname);
         //  break;

       }
       y=y1;
    }


    //итерация для y без конкурирующих процессов, когда В не единичная матрица
    void ItartionRWithGer1 (vector<double> &y, vector<vector<double> > &Matrix,vector<double> f,int kr){

           RichardsonSLAU SLAU;
           vector <vector <double> > MatrixB;
           vector <vector <double> > MatrixC;
           vector <double> oldy=y;
           //граничные условия, колличество ячеек и шаг
           N=y.size();
           h=(b-a)/N;
           ya=1;
           yb=0;

           SLAU.MinesMatrex(Matrix);

           //заполнение матрицы B, в данном случае она единичная
            MatrixB.resize(N);
            for (int i=0; i<N; i++){
                MatrixB[i].resize(N);
            }
            for (int i=0; i<N; i++){
                for (int j=0; j<N; j++){
                    MatrixB[j][j]=1;
                }
            }


            //заполнение матрицы С
                MatrixC.resize(N);
                for (int i=0; i<N; i++){
                    MatrixC[i].resize(N);
                }
             SLAU.MatrixMatrix(Matrix, MatrixB, MatrixC);

             for (int i=0; i<N; i++){
                 for (int j=0; j<N; j++){
                     cout<<MatrixC[j][j]<<" ";
                 }
                 cout<<'\n';
             }
             cout<<'\n';

                //расчет нижней границы(гамма2) с помощью кругов Герщгорина
                 vector<double> R;
                 R.resize(N);
                 for (int i=0; i<N;i++){
                     for (int j=0; j<N; j++){
                         if (i!=j) {
                             R[i]+=fabs(MatrixC[i][j]);
                         }
                     }
                 }
                 std::cout.flush();

                 gamma2=0;
                 for (int i=0; i<N;i++){
                     if (gamma2<fabs(MatrixC[i][i]+R[i])) {
                             gamma2=fabs(MatrixC[i][i]+R[i]);
                     }
                }

           gamma1=gamma2/10;

           cout<<"gamma1 = "<<gamma1<<"  gamma2 = "<<gamma2<<"\n";

           //посчитать p0
           p0=(1-gamma1/gamma2)/(1+gamma1/gamma2);

           flambda(index, lambda);

           //заполнить масив тао
           tao0=2/(gamma1+gamma2);
           tao.resize(s);
           for (int i=0; i<s; i++){
               tao[i]=tao0/(1+p0*lambda[i]);
           }


           for (int i=0; i<s; i++){
               cout<<tao[i]<<"  ";
           }
           cout<<endl;

           f[0]=ya/(h*h);
           f[N-1]=yb/(h*h);
           s=tao.size();
           Plots plot;
           deltak.resize(s);
           MultVector.resize(N);
           for (int i=0; i<s; i++){
               SLAU.MultMatrixVector(y, Matrix, MultVector);
               for (int j=0; j<N; j++){
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


    void ItartionRWithGer2 (vector<double> &y, vector<vector<double> > &Matrix,vector<double> f,int kr){

        RichardsonSLAU SLAU;
        Plots plot;
        vector<double> oldy=y;
        vector<double> y1=y;
        vector<double> y2=y;
        vector<vector<double> >MatrixB;
        vector<vector<double> >MatrixC;

        //граничные условия, колличество ячеек и шаг
        N=y.size();
        h=(b-a)/N;
        ya=1;
        yb=0;
        MultVector.resize(N);
        //-A
        SLAU.MinesMatrex(Matrix);

        //работа с матрицей B
        MatrixB.resize(N);
        for (int i=0; i<N; i++){
             MatrixB[i].resize(N);
        }

        SLAU.CreateB(Matrix, MatrixB);
        SLAU.ReB(MatrixB);

        //работа с матрицей C
            MatrixC.resize(N);
            for (int i=0; i<N; i++){
                MatrixC[i].resize(N);
            }

        SLAU.MatrixMatrix(Matrix,MatrixB, MatrixC);

       //расчет нижней границы(гамма2) с помощью кругов Герщгорина
        vector<double> R;
        R.resize(N);
        for (int i=0; i<N;i++){
            for (int j=0; j<N; j++){
                if (i!=j) {
                    R[i]+=fabs(MatrixC[i][j]);
                }
            }
        }
        std::cout.flush();

        gamma2=0;
        for (int i=0; i<N;i++){
            if (gamma2<fabs(MatrixC[i][i]+R[i])) {
                    gamma2=fabs(MatrixC[i][i]+R[i]);
            }
       }
cout<<gamma2<<"   ";
        //q и р
        q=0.1;
        p=q/4;

        //начала S-цикла
        //гамма1*
        gamma11=q*gamma2;
        std::cout<<a<<" "<<b<<"\n";
        cout<<gamma11<<"  гамма1*  ";
        //гамма1**
        gamma12=p*gamma2;
        cout<<gamma12<<"  гамма1**  ";
        bool transition;
        bool ko=false;
        std::string plname="";
        s=128;
        while (true) {
           //заполнение массива индексов
            flambda(index, lambda);
            //вектор ошибок
            deltak1.resize(s);
            deltak2.resize(s);
            MultVector.resize(N);

            //коррекция функции f
            f[0]=ya/(h*h);
            f[N-1]=yb/(h*h);
            //p0 для первого конкурирующего процесса №1
            p01=(1-gamma11/gamma2)/(1+gamma11/gamma2);
            //заполнить масив тао1
            tao01=2/(gamma11+gamma2);
            tao1.resize(s);
            for (int i=0; i<s; i++){
                tao1[i]=tao01/(1+p01*lambda[i]);
            }

            //p0 для первого конкурирующего процесса №2
            p02=(1-gamma12/gamma2)/(1+gamma12/gamma2);
            //заполнить масив тао2
            tao02=2/(gamma12+gamma2);
            tao2.resize(s);
            for (int i=0; i<s; i++){
                tao2[i]=tao02/(1+p02*lambda[i]);
            }

                    for (int i=0; i<s; i++){
                        SLAU.MultMatrixVector(y1, Matrix, MultVector);
                        for (int i=0; i<N; i++){
                            std::cout<<MultVector[i];
                        }
                        for (int j=0; j<N; j++){
                            y1[j]+=tao1[i]*(f[j]-MultVector[j]/(h*h));
                        }
                        deltak1[i]=fabs(y1[i]-oldy[i]);
                    }
                    for (int i=0; i<s; i++){
                        SLAU.MultMatrixVector(y2, Matrix, MultVector);
                        for (int j=0; j<N; j++){
                            y2[j]+=tao2[i]*(f[j]-MultVector[j]/(h*h));
                        }
                        deltak2[i]=IterError(y2,oldy);
                    }

 int           istop=s-1;
            if (ko)break;
            if (deltak2[istop]<deltak1[istop])
                {
                s=2*s;
                //continue;
            }else {
                if((deltak1[istop]-deltak2[istop])>eps*deltak1[istop]){
                    s=2*s;
                    gamma11=gamma11;
                    gamma12=gamma12/4;
                 //   continue;
                }
                else {
                    s=2*s;
                    gamma11=gamma12;
                    gamma12=gamma12/4;

                    ko=true;
               //     continue;
                    y1=y2;
                }
            }
            plname+="0";
            Plots p;
            p.IteratPlot(deltak1,"y2"+plname,"y1"+plname+".png");
            p.IteratPlot(deltak2,"y1"+plname,"y2"+plname+".png");

        }

        y1[0]=ya;
        y1[N-1]=yb;
        y2[0]=ya;
        y2[N-1]=yb;
        //конец S-цикла
        //вывод у1
        cout<<endl;
        for (int i=0; i<N; i++){
            cout<<y1[i]<<"  ";
        }
        cout<<endl;

        //вывод у2
        cout<<endl;
        for (int i=0; i<N; i++){
            cout<<y2[i]<<"   ";
        }
        cout<<endl;
        std::cout<<s;
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
    vector<double> tao1;
    vector<double> tao2;
    vector<double> deltak;
    vector<double> deltak1;
    vector<double> deltak2;

    double gamma1, gamma2, p0, tao0, tao01, tao02, p, q, gamma12, gamma11, p01, p02;


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


    //функция итерационной ошибки
    double IterError(vector<double> &y, vector<double> &oldy) {
        double maxk;
        for (int i=1; i<y.size(); i++){
            if (fabs(y[i]-oldy[i])>maxk) {maxk=fabs(y[i]-oldy[i]);}
        }
        oldy=y;
        return maxk;
    }

};


#endif // FUNCTIONS_H
