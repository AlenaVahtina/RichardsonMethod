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
        {calculate(y,SLAU,f,fold,gamma1, gamma2, deltak, true, true, 0, iterationNomber);}
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
        {calculate(y,SLAU,f,fold,gamma1, gamma2, deltak, false, false, 0, iterationNomber);}
    else {
        return;
    }

}





//функция расчета итогового значенияя вектора с канкурирующими процессами для единичной матрицы
void RichardsonMethod::computeResultVectorForEWithRivalProcess(vector<double> &y, BaseMatrix *SLAU,vector<double> f,int fold, double gamma11, double gamma12, double gamma2)
{
//    std::string plname="";
    nAmountPoints=y.size();

    step=(b-a)/nAmountPoints;
    ya=1;
    yb=0;
    gamma1=8/(b-a)*(b-a);
    gamma2=40000;
    gamma11=4000;
    gamma12=800;

    SLAU->minesMatrex();

    vector<double> y1=y;
    vector<double> y2=y;
    y1.resize(y.size());
    y2.resize(y.size());

    bool endOfIteration =false;
    while (true){
        calculate(y1,SLAU,f,fold,gamma11, gamma2, deltak1, true, true, 0, iterationNomber);
        calculate(y2,SLAU,f,fold,gamma12, gamma2, deltak2, true, 0, iterationNomber);

        int istop=iterationNomber-1;

        if (endOfIteration)break;
        break;

        if (deltak1[istop]<deltak2[istop])
            {
            iterationNomber=2*iterationNomber;
        }
        else {
            if((deltak1[istop]-deltak2[istop])>EPSELON_ERROU*deltak1[istop]){
                iterationNomber=2*iterationNomber;
                gamma11=gamma11;
                gamma12=gamma12/4;
            }
            else {
                iterationNomber=2*iterationNomber;
                gamma11=gamma12;
                gamma12=gamma12/4;

                endOfIteration=true;
                y1=y2;
            }
        }
//        plname+="0";
//        Plots p;
//        p.AveragePlotDoble(deltak1,"y1"+plname+".png","y1"+plname);
//        p.AveragePlotDoble(deltak2,"y2"+plname+".png","y2"+plname);
//        p.AveragePlotDoble2(deltak1,deltak2,"y1"+plname+".png","y1"+plname);
    }
}





void  RichardsonMethod::calculate (vector<double> &y, BaseMatrix *SLAU, vector<double> f,int fold, double gamma1, double gamma2, vector<double> &deltak, bool matrixType, bool processType, int iterationNomberFrom, int iterationNomberTo){

    iterationNomberFrom=0;
    iterationNomberTo=iterationNomber;

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
        if((i%fold==0) && (processType==false))
            plot.iteratPlot(y,std::string("output")+std::to_string(i)+".png",std::string("output")+std::to_string(i)+".dat");
        deltak[i]=Common::IterError(y,oldy);
    }

    y[0]=ya;
    y[nAmountPoints-1]=yb;
}
