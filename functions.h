#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include "plots.h"
#include "richardsonslau.h"
#include "normalmatrix.h"
#include "crsmatrix.h"

using namespace std;

#ifndef M_PI
    const double M_PI=3,1415926535897932384626;
#endif

class Richardson
{
public:
    Richardson();
    const double EPSELON_SOLUTION=0.000000000001;
    const double EPSELON_ERROU=0.3;
    void setA(double _a){a=_a;}
    void setB(double _b){b=_b;}
    void setS(double _iterationNomber){iterationNomber=_iterationNomber;}
    vector<double> getErrors(){return deltak;}
    vector<double> VectorB; //посмотреть что это за хрень

    //функция расчета итогового значения вектора у для еденичной матрицы и без конкурирующих процессов
    void computeResultVector (vector<double> &y, RichardsonSLAU *SLAU,vector<double> f,int fold){

        vector<double> oldy=y;
        //граничные условия, колличество ячеек и шаг
        nAmountPoints=y.size();
        step=(b-a)/nAmountPoints;
        ya=1;
        yb=0;

        //гамма 1 и гамма 2; границы спектра (первый метод более точный)
        //первый метод вычисления гамма
        //gamma1=4*sin(M_PI*step/2*(b-a))*sin(M_PI*step/2*(b-a))/(step*step);
        //gamma2=4*cos(M_PI*step/2*(b-a))*cos(M_PI*step/2*(b-a))/(step*step);

        //второй метод вычисления гамма
        gamma1=8/(b-a)*(b-a);
        gamma2=4/(step*step);

        //посчитать p0
        p0=(1-gamma1/gamma2)/(1+gamma1/gamma2);

        flambda(index, lambda);

        //заполнить масив тао
        tao0=2/(gamma1+gamma2);
        tao.resize(iterationNomber);
        for (int i=0; i<iterationNomber; i++){
            tao[i]=tao0/(1+p0*lambda[i]);
        }
        SLAU->MinesMatrex();
        f[0]=ya/(step*step);
        f[nAmountPoints-1]=yb/(step*step);
        iterationNomber=tao.size();
        Plots plot;
        deltak.resize(iterationNomber);
        MultVector.resize(nAmountPoints);
        for (int i=0; i<iterationNomber; i++){
            MultVector=SLAU->MultMatrixVector(y);
            for (int j=0; j<nAmountPoints; j++){
                y[j]+=tao[i]*(f[j]-MultVector[j]);
            }
            if(i%fold==0)
                plot.IteratPlot(y,std::string("output")+std::to_string(i)+".png",std::string("output")+std::to_string(i)+".dat");
            deltak[i]=IterError(y,oldy);
        }

        y[0]=ya;
        y[nAmountPoints-1]=yb;
    }

    //вспомогательная функция расчета значения у для конкурирующх процессов
    void supportingComputeResultVector (vector<double> &y, vector<vector<double> > Matrix,vector<double> f,int fold, vector<double> &deltak,int iterationNomber,double gamma1,double gamma2,int nAmountPoints,double step,int ya,int yb){

        RichardsonSLAU *SLAU=new normalmatrix(Matrix);
        vector<double> oldy=y;

        double p0=(1-gamma1/gamma2)/(1+gamma1/gamma2);

        flambda(index, lambda);
        //заполнить масив тао
        double tao0=2/(gamma1+gamma2);
        vector<double> tao;
        tao.resize(iterationNomber);
        for (int i=0; i<iterationNomber; i++){
            tao[i]=tao0/(1+p0*lambda[i]);
        }

        SLAU->MinesMatrex();
        f[0]=ya/(step*step);
        f[nAmountPoints-1]=yb/(step*step);
        iterationNomber=tao.size();
        Plots plot;
        deltak.resize(iterationNomber);
        vector<double> MultVector;
        MultVector.resize(nAmountPoints);
        for (int i=0; i<iterationNomber; i++){
            SLAU->MultMatrixVector(y);
            for (int j=0; j<nAmountPoints; j++){
                y[j]+=tao[i]*(f[j]-MultVector[j]);
            }
            deltak[i]=IterError(y,oldy);
        }

        y[0]=ya;
        y[nAmountPoints-1]=yb;
    }

    //функция расчета гамма1* гамма1** и гамма2 
    void gammacalculation(double& gamma11,double &gamma12, double & gamma2,vector<vector<double> > &Matrix){

        RichardsonSLAU *SLAU=new normalmatrix(Matrix);
        vector<vector<double> >MatrixB;
        vector<vector<double> >MatrixC;

        SLAU->MinesMatrex();

        //работа с матрицей B
        MatrixB.resize(nAmountPoints);
        for (int i=0; i<nAmountPoints; i++){
             MatrixB[i].resize(nAmountPoints);
        }

        SLAU->CreateB();
        SLAU->ReB();

        //работа с матрицей C
            MatrixC.resize(nAmountPoints);
            for (int i=0; i<nAmountPoints; i++){
                MatrixC[i].resize(nAmountPoints);
            }

        //SLAU->MatrixMatrix();

       //расчет нижней границы(гамма2) с помощью кругов Герщгорина
        vector<double> R;
        R.resize(nAmountPoints);
        for (int i=0; i<nAmountPoints;i++){
            for (int j=0; j<nAmountPoints; j++){
                if (i!=j) {
                    R[i]+=fabs(MatrixC[i][j]);
                }
            }
        }
        std::cout.flush();

        gamma2=0;
        for (int i=0; i<nAmountPoints;i++){
            if (gamma2<fabs(MatrixC[i][i]+R[i])) {
                    gamma2=fabs(MatrixC[i][i]+R[i]);
            }
       }
cout<<gamma2<<"   ";
        //q и р
        q=0.1;
        p=q/4;

        //гамма1*
        gamma11=q*gamma2;
        std::cout<<a<<" "<<b<<"\n";
        cout<<gamma11<<"  гамма1*  ";
        //гамма1**
        gamma12=p*gamma2;
        cout<<gamma12<<"  гамма1**  ";
    }

    void ItartionRFin(vector<double> &y, vector<vector<double> > &Matrix,vector<double> f,int fold){
        nAmountPoints=y.size();
        step=(b-a)/nAmountPoints;
        ya=1;
        yb=0;
        gamma1=8/(b-a)*(b-a);
        gamma2=4/(step*step);
        std::string plname="";

        gamma11=gamma1;
        gamma12=gamma1;
        std::cout<<"gamma1 "<<gamma1<<"\n";
        std::cout<<"gamma2 "<<gamma2<<"\n";

        //gammacalculation(gamma11,gamma12,gamma2,Matrix);
        std::cout<<"gamma11 "<<gamma11<<"\n";
        std::cout<<"gamma12 "<<gamma12<<"\n";
        std::cout<<"gamma2 "<<gamma2<<"\n";
        gamma2=40000;
        gamma11=4000;
        gamma12=800;
//        gamma2=gamma2*10000;
        vector<double> y1=y;
        vector<double> y2=y;

bool ko =false;
        while (true){
           supportingComputeResultVector(y1,Matrix,f,fold,deltak1,iterationNomber,gamma11,gamma2,nAmountPoints,step,ya,yb);
           supportingComputeResultVector(y2,Matrix,f,fold,deltak2,iterationNomber,gamma12,gamma2,nAmountPoints,step,ya,yb);

           int istop=iterationNomber-1;
           if (ko)break;
//           std::cout<<"\n";
//           for (int i=0;i<iterationNomber;i++){
//               std::cout << deltak1[i]<<" "<<deltak2[i]<<"\n";
//           }

           if (deltak1[istop]<deltak2[istop])
               {
               iterationNomber=2*iterationNomber;
               //continue;
           }else {
               if((deltak1[istop]-deltak2[istop])>EPSELON_ERROU*deltak1[istop]){
                   iterationNomber=2*iterationNomber;
                   gamma11=gamma11;
                   gamma12=gamma12/4;
                //   continue;
               }
               else {
                   iterationNomber=2*iterationNomber;
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

    void ItartionRWithGer1 (vector<double> &y, vector<vector<double> > &Matrix,vector<double> f,int fold){

           RichardsonSLAU *SLAU=new normalmatrix(Matrix);
           vector <vector <double> > MatrixB;
           vector <vector <double> > MatrixC;
           vector <double> oldy=y;
           //граничные условия, колличество ячеек и шаг
           nAmountPoints=y.size();
           step=(b-a)/nAmountPoints;
           ya=1;
           yb=0;

           SLAU->MinesMatrex();

           //заполнение матрицы B, в данном случае она единичная
            MatrixB.resize(nAmountPoints);
            for (int i=0; i<nAmountPoints; i++){
                MatrixB[i].resize(nAmountPoints);
            }
            for (int i=0; i<nAmountPoints; i++){
                for (int j=0; j<nAmountPoints; j++){
                    MatrixB[j][j]=1;
                }
            }


            //заполнение матрицы С
                MatrixC.resize(nAmountPoints);
                for (int i=0; i<nAmountPoints; i++){
                    MatrixC[i].resize(nAmountPoints);
                }
             //SLAU->MatrixMatrix();

             for (int i=0; i<nAmountPoints; i++){
                 for (int j=0; j<nAmountPoints; j++){
                     cout<<MatrixC[j][j]<<" ";
                 }
                 cout<<'\n';
             }
             cout<<'\n';

                //расчет нижней границы(гамма2) с помощью кругов Герщгорина
                 vector<double> R;
                 R.resize(nAmountPoints);
                 for (int i=0; i<nAmountPoints;i++){
                     for (int j=0; j<nAmountPoints; j++){
                         if (i!=j) {
                             R[i]+=fabs(MatrixC[i][j]);
                         }
                     }
                 }
                 std::cout.flush();

                 gamma2=0;
                 for (int i=0; i<nAmountPoints;i++){
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
           tao.resize(iterationNomber);
           for (int i=0; i<iterationNomber; i++){
               tao[i]=tao0/(1+p0*lambda[i]);
           }


           for (int i=0; i<iterationNomber; i++){
               cout<<tao[i]<<"  ";
           }
           cout<<endl;

           f[0]=ya/(step*step);
           f[nAmountPoints-1]=yb/(step*step);
           iterationNomber=tao.size();
           Plots plot;
           deltak.resize(iterationNomber);
           MultVector.resize(nAmountPoints);
           for (int i=0; i<iterationNomber; i++){
               SLAU->MultMatrixVector(y);
               for (int j=0; j<nAmountPoints; j++){
   //                std::cout<<"nomer "<<i<<std::endl;
   //                std::cout<<tao[i]<<std::endl;
   //                std::cout<<f[j]<<std::endl;
   //                std::cout<<MultVector[j]/(step*step)<<std::endl;
                   y[j]+=tao[i]*(f[j]-MultVector[j]);
   //                std::cout<<y[j]<<"  ";
               }
   //            cout<<endl;
               //Yplot(y,"output"+atoi(i)+".png")
               if(i%fold==0)//условие сюда ,а там цикл убрать
               plot.IteratPlot(y,std::string("output")+std::to_string(i)+".png",std::string("output")+std::to_string(i)+".dat");
               deltak[i]=IterError(y,oldy);
           }

           y[0]=ya;
           y[nAmountPoints-1]=yb;
       }


    void ItartionRWithGer2 (vector<double> &y, vector<vector<double> > &Matrix,vector<double> f,int fold){

        RichardsonSLAU *SLAU=new normalmatrix(Matrix);
        vector <vector <double> > MatrixB;
        vector <vector <double> > MatrixC;
        vector <double> oldy=y;
        //граничные условия, колличество ячеек и шаг
        nAmountPoints=y.size();
        step=(b-a)/nAmountPoints;
        ya=1;
        yb=0;

        SLAU->MinesMatrex();

        //заполнение матрицы B, в данном случае она единичная
         MatrixB.resize(nAmountPoints);
         for (int i=0; i<nAmountPoints; i++){
             MatrixB[i].resize(nAmountPoints);
         }
         SLAU->CreateB();
         SLAU->ReB();

         //заполнение матрицы С
             MatrixC.resize(nAmountPoints);
             for (int i=0; i<nAmountPoints; i++){
                 MatrixC[i].resize(nAmountPoints);
             }
          //SLAU->MatrixMatrix();

          for (int i=0; i<nAmountPoints; i++){
              for (int j=0; j<nAmountPoints; j++){
                  cout<<Matrix[i][j]<<" ";
              }
              cout<<'\n';
          }
          cout<<'\n';

          for (int i=0; i<nAmountPoints; i++){
              for (int j=0; j<nAmountPoints; j++){
                  cout<<MatrixB[i][j]<<" ";
              }
              cout<<'\n';
          }
          cout<<'\n';


          for (int i=0; i<nAmountPoints; i++){
              for (int j=0; j<nAmountPoints; j++){
                  cout<<MatrixC[i][j]<<" ";
              }
              cout<<'\n';
          }
          cout<<'\n';

             //расчет нижней границы(гамма2) с помощью кругов Герщгорина
              vector<double> R;
              R.resize(nAmountPoints);
              for (int i=0; i<nAmountPoints;i++){
                  for (int j=0; j<nAmountPoints; j++){
                      if (i!=j) {
                          R[i]+=fabs(MatrixC[i][j]);
                      }
                  }
              }
              std::cout.flush();

              gamma2=0;
              for (int i=0; i<nAmountPoints;i++){
                  if (gamma2<fabs(MatrixC[i][i]+R[i])) {
                          gamma2=fabs(MatrixC[i][i]+R[i]);
                  }
             }

        gamma1=gamma2*(step*step);

        cout<<"gamma1 = "<<gamma1<<"  gamma2 = "<<gamma2<<"\n";

        //посчитать p0
        p0=(1-gamma1/gamma2)/(1+gamma1/gamma2);

        flambda(index, lambda);

        //заполнить масив тао
        tao0=2/(gamma1+gamma2);
        tao.resize(iterationNomber);
        for (int i=0; i<iterationNomber; i++){
            tao[i]=tao0/(1+p0*lambda[i]);
        }


        for (int i=0; i<iterationNomber; i++){
            cout<<tao[i]<<"  ";
        }
        cout<<endl;

        f[0]=ya/(step*step);
        f[nAmountPoints-1]=yb/(step*step);
        iterationNomber=tao.size();
        Plots plot;
        deltak.resize(iterationNomber);
        MultVector.resize(nAmountPoints);
        for (int i=0; i<iterationNomber; i++){
            SLAU->MultMatrixVector(y);
            for (int j=0; j<nAmountPoints; j++){
//                std::cout<<"nomer "<<i<<std::endl;
//                std::cout<<tao[i]<<std::endl;
//                std::cout<<f[j]<<std::endl;
//                std::cout<<MultVector[j]/(step*step)<<std::endl;
                y[j]+=tao[i]*(f[j]-MultVector[j])/Matrix[j][j];
//                std::cout<<y[j]<<"  ";
            }
//            cout<<endl;
            //Yplot(y,"output"+atoi(i)+".png")
            if(i%fold==0)//условие сюда ,а там цикл убрать
            plot.IteratPlot(y,std::string("output")+std::to_string(i)+".png",std::string("output")+std::to_string(i)+".dat");
            deltak[i]=IterError(y,oldy);
        }

        y[0]=ya;
        y[nAmountPoints-1]=yb;
    }


private:
    double a,b,ya,yb,step;
    int nAmountPoints; //чтение колличества ячеек
    double iterationNomber;//колличесство итераций
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
        index.resize(iterationNomber);
        for (int i=0; i<iterationNomber; i++) {
            index[i]=i+1;
        }
        Rsort(index);

        lambda.resize(iterationNomber);
        for (int i=0; i<iterationNomber; i++){
            lambda[i]=-cos(M_PI*index[i]/(2*iterationNomber));
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
