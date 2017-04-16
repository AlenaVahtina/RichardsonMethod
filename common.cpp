#include "common.h"

//функция расчета гамма1 и гамма2 для модуля без конкурирующих процессов
void Common::gammacalculation1(double& gamma1, double & gamma2,BaseMatrix *SLAU, int nAmountPoints){


    //гамма 1 и гамма 2; границы спектра (первый метод более точный)
    //первый метод вычисления гамма
    //gamma1=4*sin(M_PI*step/2*(b-a))*sin(M_PI*step/2*(b-a))/(step*step);
    //gamma2=4*cos(M_PI*step/2*(b-a))*cos(M_PI*step/2*(b-a))/(step*step);

    //второй метод вычисления гамма
//    gamma1=8/(b-a)*(b-a);
//    gamma2=4/(step*step);



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
    gamma1=gamma2/10;
    cout<<gamma1<<"   "<<gamma2<<"   ";
}


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
    for (uint i=1; i<y.size(); i++){
        if (fabs(y[i]-oldy[i])>maxk) {maxk=fabs(y[i]-oldy[i]);}
    }
    oldy=y;
    return maxk;
}


//является ли s степенью 2
bool Common::fold2(int s) {
    if(s<=0)
      return false;
  while((s%2)==0)
  {
      s=s/2;

   if(s==2) {
      return true;
   }
  }
  return false;
}


bool Common::appliedCriteria(vector<double>y, vector<double>&y1, double relativE, double absolutlyE){
    int nAmountPoints=y.size();
    double max=y[0];
    for (int i=0;i<nAmountPoints; i++)
    {
        if (y[i]>y1[i]) max=y[i];
                else max=y1[i];

        if (fabs(y[i]-y1[i])<relativE*max+absolutlyE)
            return true;
        else
            return false;
    }
    return false;

}


bool Common::error(double gamma1, double gamma2, int iterationNomber)
{
    if ((0>gamma1) || (gamma1>gamma2))
    {
        cout<<"gamma1 or gamma2 does not meet the requirements"<<'\n';
        return false;
    }

    if (!Common::fold2(iterationNomber))
    {
        cout<<"s not multiply 2"<<'\n';
        return false;
    }
//    return true;
}


// фун. для проверки на пересечение
int Common::cross(vector<double> deltak1, vector<double> deltak2)
{
    int iterationNomber=deltak1.size();
    for (int i=3; i<iterationNomber; i+=4) {
        if (deltak1[i]>deltak2[i]) {
            return i;
        }
    }
    return -1;
}


bool Common::criterion27(vector<double> deltak1, vector<double> deltak2, int istop, double eps)
{
    int iterationNomber=deltak1.size();
    for (int i=istop; i<iterationNomber; i+=4) {
        if ((deltak1[i]-deltak2[i])>deltak1[i]*eps) {
            return true;
        }
    }
    return false;
}

bool Common::issolution(vector<double> deltak, double eps) {
    int iterationNomber=deltak.size();
    for (int i=3; i<iterationNomber; i+=4) {
        if (deltak[i]<eps) {
            return true;
        }
    }
    return false;
}

