#include "functions.h"

Richardson::Richardson()
{

}


//функция расчета итогового значения вектора у для еденичной матрицы и без конкурирующих процессов
void Richardson::computeResultVectorForE (vector<double> &y, RichardsonSLAU *SLAU, vector<double> f,int fold){

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
void Richardson::supportingComputeResultVector (vector<double> &y, RichardsonSLAU *SLAU,vector<double> f,int fold, vector<double> &deltak,int iterationNomber,double gamma1,double gamma2,int nAmountPoints,double step,int ya,int yb){

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

//        SLAU->MinesMatrex();
    f[0]=ya/(step*step);
    f[nAmountPoints-1]=yb/(step*step);
    iterationNomber=tao.size();
    deltak.resize(iterationNomber);
    vector<double> MultVector;
    MultVector.resize(nAmountPoints);
    for (int i=0; i<iterationNomber; i++){
        MultVector=SLAU->MultMatrixVector(y);
        for (int j=0; j<nAmountPoints; j++){
            y[j]+=tao[i]*(f[j]-MultVector[j]);
        }
        deltak[i]=IterError(y,oldy);
    }

    y[0]=ya;
    y[nAmountPoints-1]=yb;
}


//функция расчета гамма1* гамма1** и гамма2
void Richardson::gammacalculation(double& gamma11,double &gamma12, double & gamma2,RichardsonSLAU *SLAU){

    normalmatrix* matrC=static_cast<normalmatrix*>(new normalmatrix(SLAU));//перевод в нормальный тип матрицы матрицу);

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


//функция расчета итогового значенияя вектора с канкурирующими процессами для единичной матрицы
void Richardson::computeResultVectorForEWithRivalProcess(vector<double> &y, RichardsonSLAU *SLAU,vector<double> f,int fold){
    std::string plname="";
    nAmountPoints=y.size();

    step=(b-a)/nAmountPoints;
    ya=1;
    yb=0;
    gamma1=8/(b-a)*(b-a);
    gamma2=40000;
    gamma11=4000;
    gamma12=800;

    SLAU->MinesMatrex();

    vector<double> y1=y;
    vector<double> y2=y;
    y1.resize(y.size());
    y2.resize(y.size());

    bool endOfIteration =false;
    while (true){
       supportingComputeResultVector(y1,SLAU,f,fold,deltak1,iterationNomber,gamma11,gamma2,nAmountPoints,step,ya,yb);
       supportingComputeResultVector(y2,SLAU,f,fold,deltak2,iterationNomber,gamma12,gamma2,nAmountPoints,step,ya,yb);

       int istop=iterationNomber-1;

       if (endOfIteration)break;

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
       plname+="0";
       Plots p;
//           p.AveragePlotDoble(deltak1,"y1"+plname+".png","y1"+plname);
//           p.AveragePlotDoble(deltak2,"y2"+plname+".png","y2"+plname);
       p.AveragePlotDoble2(deltak1,deltak2,"y1"+plname+".png","y1"+plname);

   }
   y=y1;
}


//функция расчета итогового значенияя вектора с канкурирующими процессами для не единичной матрицы
void Richardson::computeResultVectorForNotEWithRivalProcess(vector<double> &y, RichardsonSLAU *SLAU,vector<double> f,int fold){
    std::string plname="";
    nAmountPoints=y.size();

    step=(b-a)/nAmountPoints;
    ya=1;
    yb=0;
    gamma1=8/(b-a)*(b-a);
    gamma2=4/(step*step);
    gamma11=gamma1;
    gamma12=gamma1;


    vector<double> y1=y;
    vector<double> y2=y;
    y1.resize(y.size());
    y2.resize(y.size());

    SLAU->MinesMatrex();
    RichardsonSLAU *SLAU2=SLAU->CreateC();
    gammacalculation(gamma11, gamma12, gamma2, SLAU2);
//        gamma11=4000;
//        gamma12=800;
//        gamma2=4/(step*step);

    //SLAU2->WriteMatrix();


    bool notExcept2=false;
    bool endOfIteration =false;

    while (true){

       deltak1.resize(iterationNomber);
       deltak2.resize(iterationNomber);

       supportingComputeResultVector(y1,SLAU2,f,fold,deltak1,iterationNomber,gamma11,gamma2,nAmountPoints,step,ya,yb);
       supportingComputeResultVector(y2,SLAU2,f,fold,deltak2,iterationNomber,gamma12,gamma2,nAmountPoints,step,ya,yb);


       int istop=iterationNomber-1;
       if ((endOfIteration) || (deltak1[istop]<EPSELON_SOLUTION) || (deltak2[istop]<EPSELON_SOLUTION))break;

       cout<<deltak1[istop]<<" "<<deltak2[istop]<<"\n";
       if ((deltak1[istop]<deltak2[istop])&&(deltak2[istop]>EPSELON_SOLUTION*10000))
           {
           iterationNomber=2*iterationNomber;
       }
       else {
           cout<<"else "<<deltak1[istop]-deltak2[istop]<<" "<<deltak1[istop]<<" "<< EPSELON_ERROU*deltak1[istop]<<"\n";
           if(((deltak1[istop]-deltak2[istop])>EPSELON_ERROU*deltak1[istop])&&(deltak2[istop]>EPSELON_SOLUTION*10000))
           {
               iterationNomber=2*iterationNomber;
               gamma11=gamma11;
               gamma12=gamma12/4;
               notExcept2=true;
           }
           else {
               iterationNomber=2*iterationNomber;
               gamma11=gamma12;
               gamma12=gamma12/4;

               endOfIteration=true;
               y1=y2;
           }
       }
       plname+="0";
       Plots p;
//           p.AveragePlotDoble(deltak1,"y1"+plname+".png","y1"+plname);
//           p.AveragePlotDoble(deltak2,"y2"+plname+".png","y2"+plname);
           p.AveragePlotDoble2(deltak1,deltak2,"y1"+plname+".png","y1"+plname);



   }
   y=y1;
}


//функция сортировки
void Richardson::Rsort(vector<int> &index){
    for (int power=1; power<=log(index.size())/log(2); power++){
        for (int k=pow(2,(power-1))-1;k>=0; k--){
            index[2*k]=index[k];
            index[2*k+1]=pow(2,(power+1))-index[2*k];
        }
    }
}


//сделать массив лямбда
void Richardson::flambda (vector <int> &index, vector <double> &lambda){
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
double Richardson::IterError(vector<double> &y, vector<double> &oldy) {
    double maxk;
    for (int i=1; i<y.size(); i++){
        if (fabs(y[i]-oldy[i])>maxk) {maxk=fabs(y[i]-oldy[i]);}
    }
    oldy=y;
    return maxk;
}
