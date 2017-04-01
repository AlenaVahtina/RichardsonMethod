#include "normalmatrix.h"
#include "crsmatrix.h"

NormalMatrix::NormalMatrix(vector <vector <double> > Matrix_)
{
    Matrix = Matrix_;
}


//перевод из Йельского формата в обычный
void TranslateCRSNormal(vector<vector<double> > &Matrix, vector <double> values,vector <int> cols, vector <int>  pointer )
{
    int N=pointer.size();
    int colon;
    for (int i=1; i<=N; i++){
        for (int j=pointer[i-1]; j<pointer[i]; j++){
            colon=cols[j];
            Matrix[i-1][colon]=values[j];
        }
    }
}


NormalMatrix::NormalMatrix(BaseMatrix * another){
    if (another->getType()==BaseMatrix::NORMALMATRIX)
    {
        NormalMatrix* normal=static_cast<NormalMatrix*>(another);
        Matrix = normal->Matrix;
    }
    else if (another->getType()==BaseMatrix::CRSMATRIX){
        //перевод из Йельского в нормальный
        CrsMatrix* crs=static_cast<CrsMatrix*>(another);
        vector <int> toPointer=crs->pointer;
        vector <int> toCols=crs->cols;
        vector <double> toValuse=crs->values;

        TranslateCRSNormal(Matrix, toValuse, toCols, toPointer);
    }
}


//чтение матрицы с клавиатуры
void NormalMatrix::readMatrix(int nAmountPoints){
    vector <vector <double> > Matrix;
    Matrix.resize(nAmountPoints);
    for (int i=0;i<nAmountPoints; i++){
        for (int j=0;j<nAmountPoints; j++){
            cin>>Matrix[i][j];
        }
    }
}


//чтение матрицы из файла
void NormalMatrix::readMatrixFile(int nAmountPoints){
    ifstream dataFiele("data.txt");
    if (!dataFiele) {return;}
    Matrix.resize(nAmountPoints);
    for (int i=0;i<nAmountPoints; i++){
        Matrix[i].resize(nAmountPoints);
        for (int j=0;j<nAmountPoints; j++){
            dataFiele>>Matrix[i][j];
        }
    }
    dataFiele.close();
}


//создание матрицы B из A
BaseMatrix *NormalMatrix::createB() {
    double max;
    vector <vector <double> > MatrixB;
    MatrixB.resize(Matrix.size());
    for (int i=0; i<MatrixB.size();i++){
        MatrixB[i].resize(Matrix.size());
        max=Matrix[i][0];
        for (int j=0; j<Matrix.size(); j++){
            if (Matrix[i][j]!=0) {
                if (Matrix[i][j]>max){max=Matrix[i][j];}
            }
        }
        if (Matrix[i].empty()) {std::cout<<"Error. In your matrix you have line with all 0 "<<endl;}
        for (int j=0; j<MatrixB.size(); j++){
            if ((i==j) && (Matrix[i][j]!=0) && (!Matrix[i].empty())) {
                MatrixB[i][j]=Matrix[i][j];
            }
            if ((i==j) && (Matrix[i][j]==0) && (!Matrix[i].empty())) {
                MatrixB[i][j]=max;
            }
        }
    }

    return  new NormalMatrix (MatrixB);
}


//преобразование матрицы^(-1/2)
BaseMatrix * NormalMatrix::reB(){
    for (int i=0; i<Matrix.size();i++){
        for (int j=0; j<Matrix.size(); j++){
            if (i==j) {Matrix[i][j]=pow(Matrix[i][j],-0.5);}
            }
        }
    return this;
}


//функция умножения матриц
BaseMatrix *NormalMatrix::matrixMatrix (BaseMatrix * secondM)
{
    int nAmountPoints=Matrix.size();
    vector <vector <double> > resultMatrix;
    resultMatrix.resize(Matrix.size());
    NormalMatrix* second=static_cast<NormalMatrix*>(secondM);
    for (int i=0; i<nAmountPoints;i++){
        resultMatrix[i].resize(Matrix.size());
        for (int j=0; j<nAmountPoints; j++){
            resultMatrix[i][j]=0;
            for (int l=0; l<nAmountPoints; l++){
            resultMatrix[i][j]+=second->Matrix[i][l] * this->Matrix[l][j];
            }
        }
     }
    return new NormalMatrix(resultMatrix);
}


//функция создания матрицы С=B^(-1/2)*A*B^(-1/2)
BaseMatrix *NormalMatrix::createC(){
    BaseMatrix *MatrixB=this->createB()->reB();
    BaseMatrix *MatrixC=this->createB()->reB();
    return MatrixC->matrixMatrix(static_cast<BaseMatrix*>(this))->matrixMatrix(MatrixB);
}

//функция обращения матрицы(-A)
BaseMatrix * NormalMatrix::minesMatrex(){
    for (int i=0; i<Matrix.size(); i++){
        for (int j=0; j<Matrix.size(); j++){
            Matrix[i][j]*=-1;
        }
    }
    return this;
}


//функция умножения матрицы на вектор
vector <double> NormalMatrix::multMatrixVector( vector <double> y){
    vector <double> MultVector;
    int ysize=y.size();
    double sum=0;
    vector<double> part;
    MultVector.resize(ysize);
    for (int i=0; i<ysize; i++){
        MultVector[i]=0;
        sum=0;
        part=Matrix[i];
        for (int j=0; j<ysize; j++){
            sum+=part[j]*y[j];
        }
        MultVector[i]=sum;
        sum=0;
    }
    return MultVector;
}


//функция вывода матрицы
void NormalMatrix::writeMatrix (){
    for (int i=0; i<Matrix.size(); i++){
        for (int j=0; j<Matrix.size(); j++){
            cout<<Matrix[i][j]<<"  ";
        }
        cout<<endl;
    }
}


//взять элемент
double NormalMatrix::getElement (int row, int column) const
{
    return Matrix[row][column];
}
