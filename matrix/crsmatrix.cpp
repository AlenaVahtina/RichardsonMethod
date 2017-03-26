#include "crsmatrix.h"
#include "normalmatrix.h"

CrsMatrix::CrsMatrix(vector <int> pointer_, vector <int> cols_, vector <double> values_)
{
    pointer=pointer_;
    cols=cols_;
    values=values_;
}

void translateNormalCRS(vector<vector<double> > Matrix, vector <double>  &values,vector <int> &cols, vector <int>  &pointer )
{
    int k=0;
    int numb;
    int N=Matrix[0].size();
    pointer.resize(N+1);
    pointer[0]=0;
    for (int i=0; i<N;i++){
        if (i!=0) {pointer[i]=pointer[i-1]+numb;}
        numb=0;
        for (int j=0; j<N; j++){
            if (Matrix[i][j]!=0){
                k++;
                numb++;
                values.resize(k);
                values[k-1]=Matrix[i][j];
                cols.resize(k);
                cols[k-1]=j;
            }
        }
    }
    pointer[N]=pointer[N-1]+numb;
  }

CrsMatrix::CrsMatrix(BaseMatrix * another){
    if (another->getType()==BaseMatrix::CRSMATRIX){
        CrsMatrix* crs=static_cast<CrsMatrix*>(another);
        pointer=crs->pointer;
        cols=crs->cols;
        values=crs->values;
    }else if (another->getType()==BaseMatrix::NORMALMATRIX){
        //перевод из matrix в три вектора
        NormalMatrix* normal=static_cast<NormalMatrix*>(another);
        vector<vector<double>> to =normal->Matrix;
        translateNormalCRS(to, values,cols, pointer);
    }
}

//чтение матрицы с клавиатуры в Йельском формате
void CrsMatrix::readMatrix (int length, int nomberNotNullElemet){

    pointer.resize(length);
    cout<<'\n'<<"pointer"<<'\n';
    for (int i=0; i<length;i++){
        cin>>pointer[i];
    }

    values.resize(nomberNotNullElemet);
    cout<<'\n'<<"values"<<'\n';
    for (int i=0; i<nomberNotNullElemet;i++){
        cin>>values[i];
      }

    cols.resize(nomberNotNullElemet);
    cout<<'\n'<<"cols"<<'\n';
    for (int i=0; i<nomberNotNullElemet;i++){
        cin>>cols[i];
    }
}

//чтение матрицы с клавиатуры в Йельском формате
void CrsMatrix::readMatrixFile (int length, int nomberNotNullElemet){

    ifstream dataFiele("data.txt");
    if (!dataFiele) {exit (1);}

    pointer.resize(length);

    for (int i=0; i<length;i++){
        dataFiele>>pointer[i];
    }

    values.resize(nomberNotNullElemet);

    for (int i=0; i<nomberNotNullElemet;i++){
        dataFiele>>values[i];
      }

    cols.resize(nomberNotNullElemet);

    for (int i=0; i<nomberNotNullElemet;i++){
        dataFiele>>cols[i];
    }
}


//создание матрицы B из А
BaseMatrix *CrsMatrix::createB()
{
 vector <double> valuesB;
 vector <int> colsB;
 vector <int>  pointerB;
 int point=0;
 pointerB.push_back(0);
 for (int i=0; i<pointer.size()-1; i++){
    point=0;
    for (int j=pointer[i]; j<pointer[i+1]; j++){
        if (cols[j]==i){
                valuesB.push_back(values[j]);
                colsB.push_back(cols[j]);
                point++;
        }
    }
    if (point!=0){pointerB.push_back(pointerB[i]+point);}
 }
 return new CrsMatrix (pointerB, colsB, valuesB);
}


//преобразование матрицы(-1/2)
BaseMatrix * CrsMatrix::reB()
{
 for (int i=0; i<values.size(); i++){
    values[i]=1/sqrt(values[i]);
 }
 return this;
}


//фунция транспонирования матрицы
BaseMatrix * CrsMatrix::TCRSMaatrix()
{
    vector <vector <int> > int_vectors;
    vector <vector <double> > RealVectors;
    int_vectors.resize(pointer.size()-1);
    RealVectors.resize(pointer.size()-1);
    for (int i=1; i<pointer.size(); i++){
        for (int j=pointer[i-1]; j<pointer[i]; j++){
            int_vectors[cols[j]].push_back(i-1);
            RealVectors[cols[j]].push_back(values[j]);
        }
    }

    pointer.push_back(0);
    for (int i=0; i<int_vectors.size(); i++){
        for (int j=0; j<int_vectors[i].size(); j++){

                cols.push_back(int_vectors[i][j]);
                values.push_back(RealVectors[i][j]);

        }
       pointer.push_back(pointer[i]+int_vectors[i].size());
    }

    for (int i=0; i<pointer.size();i++){
            cout<<pointer[i]<<"  ";
      }
    cout<<'\n';
    return this;
}


//функция умножение матриц
BaseMatrix *CrsMatrix::matrixMatrix (BaseMatrix * secondM)
{
    vector <double>result_values;
    vector <int> result_cols;
    vector <int> result_pointe;

    CrsMatrix* second=static_cast<CrsMatrix*>(secondM);
    second->TCRSMaatrix();

    double sum;
    int resultPointer;
    result_pointe.push_back(0);
    for (int i=0; i<pointer.size()-1; i++){
        resultPointer=0;
        for(int j=0; j<pointer.size()-1; j++){
            sum=0;
            for (int k=pointer[i]; k<pointer[i+1]; k++){
                for (int l=second->pointer[j]; l<second->pointer[j+1];l++){
                    if (cols[k]==second->cols[l]){
                        sum+=values[k]*second->values[l];
                        break;
                    }
             }
            }
            if (sum!=0){
                result_values.push_back(sum);
                result_cols.push_back(j);
                resultPointer++;
            }
        }
        result_pointe.push_back(resultPointer+result_pointe[i]);
    }
    return new CrsMatrix(result_pointe,result_cols,result_values);
}


//функция создания матрицы С=B^(-1/2)*A*B^(-1/2)
BaseMatrix *CrsMatrix::createC(){
    BaseMatrix *MatrixB=this->createB()->reB();
    BaseMatrix *MatrixC=this->createB()->reB();
    MatrixC->matrixMatrix(static_cast<BaseMatrix*>(this))->matrixMatrix(MatrixB);
    return MatrixC;
}

//функция обращения матрицы (-А)
BaseMatrix * CrsMatrix::minesMatrex(){
    for (int i=0; i<values.size(); i++){
       values[i]*=-1;
    }
    return this;
}


//функция умножения матрицы на вектор
vector <double> CrsMatrix::multMatrixVector(vector <double> y)
{
    int N=pointer.size();
    vector <double> vectorC;
    for (int i=1; i<=N; i++){
        for (int j=pointer[i-1]; j<pointer[i]; j++){
                vectorC[i-1]+=values[j]*y[cols[j]];
        }
    }

    return vectorC;
}


//функция вывода матрицы в Йельском формате
void CrsMatrix::writeMatrix (){

    int length=pointer.size();
    cout<<'\n'<<"pointer"<<'\n';
    for (int i=0; i<length;i++){
        cout<<pointer[i]<<"  ";
    }
    cout<<'\n';

    int nomberNotNullElemet=values.size();
    cout<<'\n'<<"values"<<'\n';
    for (int i=0; i<nomberNotNullElemet;i++){
        cout<<values[i]<<"  ";
      }
    cout<<'\n';

    cout<<'\n'<<"cols"<<'\n';
    for (int i=0; i<nomberNotNullElemet;i++){
        cout<<cols[i]<<"  ";
    }
    cout<<'\n';
}


//взять элемент
double CrsMatrix::getElement (int row, int column) const
{
    for (int i=pointer[row]; i<pointer[row+1]; i++){
        if (cols[i]==column){
            return values[i];
        }
    }
}
