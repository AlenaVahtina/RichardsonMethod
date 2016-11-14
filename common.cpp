#include "common.h"

//функция расчета гамма1* гамма1** и гамма2
void Common::gammacalculation(double& gamma11,double &gamma12, double & gamma2,BaseMatrix *SLAU, double p, double q, int nAmountPoints){

    NormalMatrix* matrC=static_cast<NormalMatrix*>(new NormalMatrix(SLAU));//перевод в нормальный тип матрицы матрицу);

   //расчет нижней границы(гамма2) с помощью кругов Герщгорина
    vector<double> R;
    R.resize(nAmountPoints);
    for (int i=0; i<nAmountPoints;i++){
        for (int j=0; j<nAmountPoints; j++){
            if (i!=j) {
                R[i]+=fabs(matrC->Matrix[i][j]);
            }
        }
    }
    std::cout.flush();

    gamma2=0;
    for (int i=0; i<nAmountPoints;i++){
        if (gamma2<fabs(matrC->Matrix[i][i]+R[i])) {
                gamma2=fabs(matrC->Matrix[i][i]+R[i]);
        }
   }
    cout<<gamma2<<"   ";
    //q и р
//    q=0.5;
//    p=q/4;

    //гамма1*
    gamma11=q*gamma2;
    cout<<gamma11<<"  гамма1*  ";
    //гамма1**
    gamma12=p*gamma2;
    cout<<gamma12<<"  гамма1**  ";
}


//функция сортировки
void Common::Rsort(vector<int> &index){
    for (int power=1; power<=log(index.size())/log(2); power++){
        for (int k=pow(2,(power-1))-1;k>=0; k--){
            index[2*k]=index[k];
            index[2*k+1]=pow(2,(power+1))-index[2*k];
        }
    }
}


//сделать массив лямбда
void Common::flambda (vector <int> &index, vector <double> &lambda, double iterationNomber){
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
double Common::IterError(vector<double> &y, vector<double> &oldy) {
    double maxk;
    for (int i=1; i<y.size(); i++){
        if (fabs(y[i]-oldy[i])>maxk) {maxk=fabs(y[i]-oldy[i]);}
    }
    oldy=y;
    return maxk;
}
