#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>

using namespace std;

int main()
{    
    double k1=2;
    double k2=15;
    int nAmountPoints=51;
    ofstream dataFiele("data5.txt");
    vector < vector < double > > Matrix;
    double step = double (1.0/nAmountPoints);

    cout<<step<<"  "<<'\n';

    Matrix.resize(nAmountPoints);
    for (int i=0; i<nAmountPoints; i++){
        Matrix[i].resize(nAmountPoints);
    }
    //for (int i=0; i<nAmountPoints/2; i++){
        for (int j=0; j<nAmountPoints/2; j++){
            if(j>0) Matrix[j-1][j]=k1/(step*step);
            Matrix[j][j]=-2*k1/(step*step);
            if(j<nAmountPoints-1) Matrix[j+1][j]=k1/(step*step);
        }
   // }

    Matrix[nAmountPoints/2-1][nAmountPoints/2]=k1/(step*step);
    Matrix[nAmountPoints/2][nAmountPoints/2]=-(k2+k1)/(step*step);
    Matrix[nAmountPoints/2+1][nAmountPoints/2]=k2/(step*step);

    //for (int i=nAmountPoints/2+1; i<=nAmountPoints; i++){
        for (int j=nAmountPoints/2+1; j<nAmountPoints; j++){
            if(j>0) Matrix[j-1][j]=k2/(step*step);
            Matrix[j][j]=-2*k2/(step*step);
            if(j<nAmountPoints-1) Matrix[j+1][j]=k2/(step*step);
        }
       // cout<<i<<"  ";
    //}

    for (int i=0;i<nAmountPoints; i++){
        for (int j=0;j<nAmountPoints; j++){
            dataFiele<<Matrix[i][j]<<" ";
        }
        dataFiele<<std::endl;
    }
    dataFiele.close();


    return 0;
}



