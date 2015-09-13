#include <iostream>
#include "armadillo"
#include <time.h>
#include <math.h>
#include <fstream>

using namespace std;
using namespace arma;


int main()
{


//defining size of and iteration

double n = 1000;
int i = 1;
double h = 1/(n+1);

//starting time
double t0  = clock();

//Defining vectors for tridiagonal decomposition
rowvec a = zeros<rowvec>(n);
for(i = 0; i< n; i++ ){
    a(i) = -1/pow(h,double(2)); //filling vectors
}
a(0) = 0.0;

rowvec b = zeros<rowvec>(n); //diagonal vector
for(i = 0; i< n; i++ ){
    b(i) = 2;
}


rowvec c = zeros<rowvec>(n);
for(i = 0; i< n; i++ ){
    c(i) = -1/pow(h,double(2));
}


c(n-1) = 0;

//Making vectors for memory allocation
rowvec u = zeros<rowvec>(n);
rowvec f = zeros<rowvec>(n);
rowvec x = linspace<rowvec>(0,1,n); //Equally spaced x values

for(i = 0; i<n; i++){
    f(i)= pow(h,double(2))*100*exp(-10*x(i)); //getting function value
}
//cout<<"x = "<<x<<endl;
//cout<<"fi = "<<f<<endl;


rowvec temp = zeros<rowvec>(n);
//forward substitution loop
double btemp = b(0);
u(1) = f(0)/btemp;
for(i=1 ; i < n ; i++) {
    temp(i) = c(i-1)/btemp;
    btemp = b(i)-a(i)*temp(i);
    u(i) = (f(i) - a(i)*u(i-1)/btemp);
    }

//cout<<u<<endl;
//backward substitution loop
for(i=n-2; i>=0; i--){
    u(i) -= temp(i+1)*u(i+1);
}

//cout<<"u= "<<u<<endl;

//getting closed form soulution
rowvec ure = zeros<rowvec>(n);
for(i= 1; i<n; i++){
    ure(i) = 1-(1-exp(-10))*x(i)-exp(-10*x(i));
    }
//calculating error
rowvec epsu = zeros<rowvec>(n);
for(i = 0; i<n; i++){
    epsu(i) = log10(abs((u(i)-ure(i))/ure(i)));
}

//cout<<"ure= "<<ure<<endl;

double t1 = clock(); //stop for the first method
double cpu_time = (t1 - t0)/CLOCKS_PER_SEC; //Getting time used

//cout<<"Time used in sec = "<<cpu_time<<endl;



double t3  = clock(); //starting time for LU


//LU Decomposition
//Creationg matrix A

mat A = eye<mat>(n,n);
cout<<"rows in A.n_rows "<<A.n_rows<<endl;
for(i = 0; i < n; i++){
    A(i,i) = 2; //Filling in values in matrix
}
for(i= 1; i< n;i++){
    A(i, i-1) = -1;
}
for(i= 0; i< n-1; i++){
    A(i, i+1) = -1;
}



mat L, U, P; //Decomposing LU

lu(L, U, P, A);
mat v = solve(trimatu(U), solve(trimatl(L), P*trans(f) )); //Solving for v

//cout<<v<<endl;
double t4 = clock(); //end time for LU
double cpu_time2 = (t4 - t3)/CLOCKS_PER_SEC; //time used
cout<<"Size "<< n<<" Time used = "<<cpu_time<<" Time LU = "<<cpu_time2<<endl; //comparing time used
rowvec epsv = zeros<rowvec>(n);
for(i = 0; i<n; i++){
    epsv(i) = log10(abs((v(i)-ure(i))/ure(i)));
}
cout<<"epsu = "<<max(epsu)<<" epsv = "<<max(epsv)<<endl;


ofstream output("ctopl1000.txt"); //saving to file

output <<"x"<<endl;

for (i=1; i<x.size(); i++)
    {
        output << x(i) << " "; // behaves like cout - cout is also a stream
    }
    output << ";" << endl;

output <<"u Tridiagonal "<<endl;

for (i=1; i<u.size(); i++)
        {
            output << u(i) << " "; // behaves like cout - cout is also a stream
        }
        output << ";" << endl;
output <<"u real"<<endl;

for (i=1; i<ure.size(); i++)
      {
        output << ure(i) << " "; // behaves like cout - cout is also a stream
             }
output << ";" << endl;

output<<"LU v"<<endl;

for (i=1; i<v.size(); i++){
    output<<v(i) << " ";
}
output<<";"<<endl;

output<<"epsu error Tri"<<endl;

for (i=1; i<epsu.size(); i++){
    output<<epsu(i) << " ";
}
output<<";"<<endl;
output<<"epsu error LU"<<endl;

for (i=1; i<epsv.size(); i++){
    output<<epsv(i) << " ";
}
output<<";"<<endl;

return 0;

}
/*
Size 10 Time used = 5.5e-05 Time LU = 0.000175
epsu = inf epsv = inf

Size 100 Time used = 0.000151 Time LU = 0.000858
epsu = inf epsv = inf

Size 1000 Time used = 0.00101 Time LU = 0.07432
epsu = inf epsv = inf
*/
