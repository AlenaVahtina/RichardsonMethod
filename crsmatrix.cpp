#include "crsmatrix.h"

crsmatrix::crsmatrix(vector <int> pointer_, vector <int> cols_, vector <double> values_)
{
    pointer=pointer_;
    cols=cols_;
    values=values_;
}
