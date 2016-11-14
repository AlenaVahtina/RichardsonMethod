#include "richardsonMethodWithChebyshevOrderedSetOfParameters.h"

//функция расчета итогового значения вектора у для еденичной матрицы и без конкурирующих процессов
void RichardsonMethod::computeResultVectorForE (vector<double> &y, BaseMatrix *SLAU, vector<double> f,int fold){

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

    Common::flambda(index, lambda);

    //заполнить масив тао
    tao0=2/(gamma1+gamma2);
    tao.resize(iterationNomber);
    for (int i=0; i<iterationNomber; i++){
        tao[i]=tao0/(1+p0*lambda[i]);
    }
    SLAU->minesMatrex();
    f[0]=ya/(step*step);
    f[nAmountPoints-1]=yb/(step*step);
    iterationNomber=tao.size();
    Plots plot;
    deltak.resize(iterationNomber);
    MultVector.resize(nAmountPoints);
    for (int i=0; i<iterationNomber; i++){
        MultVector=SLAU->multMatrixVector(y);
        for (int j=0; j<nAmountPoints; j++){
            y[j]+=tao[i]*(f[j]-MultVector[j]);
        }
        if(i%fold==0)
            plot.iteratPlot(y,std::string("output")+std::to_string(i)+".png",std::string("output")+std::to_string(i)+".dat");
        deltak[i]=Common::IterError(y,oldy);
    }

    y[0]=ya;
    y[nAmountPoints-1]=yb;
}


//функция для не единичной без конкурирущих, т.е для gamma вручную
void RichardsonMethod::computeResultVectorForC (vector<double> &y, BaseMatrix *SLAU, vector<double> f,int fold){

    vector<double> oldy=y;
    //граничные условия, колличество ячеек и шаг
    nAmountPoints=y.size();
    step=(b-a)/nAmountPoints;
    ya=1;
    yb=0;

    BaseMatrix *SLAUB=SLAU->createB();
    double divider=SLAUB->getElement(0,0);

//    gammacalculation(gamma11, gamma12, gamma2, SLAU2);


    gamma2=2;
    gamma1=0.0002;
    //посчитать p0
    p0=(1-gamma1/gamma2)/(1+gamma1/gamma2);

    Common::flambda(index, lambda);

    //заполнить масив тао
    tao0=2/(gamma1+gamma2);
    tao.resize(iterationNomber);
    for (int i=0; i<iterationNomber; i++){
        tao[i]=tao0/(1+p0*lambda[i]);
    }

    f[0]=ya/(step*step);
    f[nAmountPoints-1]=yb/(step*step);
    iterationNomber=tao.size();

    Plots plot;
    SLAU->minesMatrex();
    deltak.resize(iterationNomber);
    MultVector.resize(nAmountPoints);
    for (int i=0; i<iterationNomber; i++){
        MultVector=SLAU->multMatrixVector(y);
        for (int j=0; j<nAmountPoints; j++){
            y[j]+=tao[i]*(f[j]-MultVector[j])/-divider;
        }
        if(i%fold==0)
            plot.iteratPlot(y,std::string("output")+std::to_string(i)+".png",std::string("output")+std::to_string(i)+".dat");
        deltak[i]=Common::IterError(y,oldy);
    }

    y[0]=ya;
    y[nAmountPoints-1]=yb;
}
