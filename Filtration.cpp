#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

int main()
{    
    //Эти константы настраиваются
    double K=0.0046;
    double betta=22.13;
    double timeConst=1.0;

    int nAmountPoints=51;
    ofstream dataFiele("data6.txt");
    vector < vector < double > > Matrix;
    double step = double (1.0/nAmountPoints);

    //Заполнение матрицы
    Matrix.resize(nAmountPoints);
    for (int i=0; i<nAmountPoints; i++){
        Matrix[i].resize(nAmountPoints);
    }
    for (int i=0; i<nAmountPoints; i++){
        for (int j=0; j<nAmountPoints; j++){
            if(j>0) Matrix[j-1][j]=-K/(step*step);
            Matrix[j][j]=2*K/(step*step)+(1.0/timeConst);
            if(j<nAmountPoints-1) Matrix[j+1][j]=-K/(step*step);
        }
   }
   Matrix[nAmountPoints/2][nAmountPoints/2]+=betta;


   
    //Заполнение вектора f
    for (int i=0;i<nAmountPoints; i++){
        for (int j=0;j<nAmountPoints; j++){
            dataFiele<<Matrix[i][j]<<" ";
        }
        dataFiele<<std::endl;
    }
    dataFiele.close();





    ofstream dataFiele2("data10.txt");
    vector < double > f;
    f.resize(nAmountPoints);
    for (int i=0; i<nAmountPoints; i++){
        f[i]=1;
    }

    for (int i=0;i<nAmountPoints; i++){
            dataFiele2<<f[i]<<" ";
    }
    dataFiele2.close();


    return 0;
}



