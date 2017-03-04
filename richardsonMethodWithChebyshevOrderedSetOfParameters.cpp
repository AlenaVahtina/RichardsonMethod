#include "richardsonMethodWithChebyshevOrderedSetOfParameters.h"

//функция расчета итогового значения вектора у для еденичной матрицы и без конкурирующих процессов
void RichardsonMethod::computeResultVectorForE (vector<double> &y, BaseMatrix *SLAU, vector<double> f,int fold, double gamma1, double gamma2){

    //граничные условия, колличество ячеек и шаг
    nAmountPoints=y.size();
    step=(b-a)/nAmountPoints;
    ya=1;
    yb=0;


    SLAU->minesMatrex();

    Common::gammacalculation1(gamma1,gamma2, SLAU, nAmountPoints);

    if (!Common::error(gamma1, gamma2, iterationNomber))
        {calculate(y,SLAU,f,fold,gamma1, gamma2, true);}
    else {
        return;
    }

}


//функция для не единичной без конкурирущих, т.е для gamma вручную
void RichardsonMethod::computeResultVectorForC (vector<double> &y, BaseMatrix *SLAU, vector<double> f,int fold, double gamma1, double gamma2){

//    vector<double> oldy=y;
    //граничные условия, колличество ячеек и шаг
    nAmountPoints=y.size();
    step=(b-a)/nAmountPoints;
    ya=1;
    yb=0;

    gamma2=2;
    gamma1=0.0002;

    if (!Common::error(gamma1, gamma2, iterationNomber))
        {calculate(y,SLAU,f,fold,gamma1, gamma2, false);}
    else {
        return;
    }

}


void  RichardsonMethod::calculate (vector<double> &y, BaseMatrix *SLAU, vector<double> f,int fold, double gamma1, double gamma2, bool matrixType){


    BaseMatrix *SLAUB;
    if (!matrixType){
        SLAUB=SLAU->createB();
        divider=SLAUB->getElement(0,0);
    }
    else
    {
         divider=1;
         SLAU->minesMatrex();
    }

    //посчитать p0
    p0=(1-gamma1/gamma2)/(1+gamma1/gamma2);

    Common::flambda(index, lambda, iterationNomber);

    //заполнить масив тао
    tao0=2/(gamma1+gamma2);
    tao.resize(iterationNomber);
    for (int i=0; i<iterationNomber; i++){
        tao[i]=tao0/(1+p0*lambda[i]);
    }
    SLAU->minesMatrex();

    vector<double> oldy=y;
    f[0]=ya/(step*step);
    f[nAmountPoints-1]=yb/(step*step);
    iterationNomber=tao.size();

    Plots plot;
    deltak.resize(iterationNomber);
    MultVector.resize(nAmountPoints);
    for (int i=0; i<iterationNomber; i++){
        MultVector=SLAU->multMatrixVector(y);
        for (int j=0; j<nAmountPoints; j++){
            if (!matrixType) {
                divider=SLAUB->getElement(j,j);
                divider=-divider;
            }
            y[j]+=tao[i]*(f[j]-MultVector[j])/divider;
        }
        if(i%fold==0)
            plot.iteratPlot(y,std::string("output")+std::to_string(i)+".png",std::string("output")+std::to_string(i)+".dat");
        deltak[i]=Common::IterError(y,oldy);
    }

    y[0]=ya;
    y[nAmountPoints-1]=yb;
}
