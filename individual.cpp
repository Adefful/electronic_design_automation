
/*
*
* Author: Zadorozhnyi Pavel
* date:   15.04.21
*
*/
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <iomanip> 
using namespace std;
// https://github.com/ochaplashkin/numerical_analysis

double  EPS     =  0.001, 
        phit    = 0.025,    // v
        I0      = 1.2e-10,  // A
        R1      = 10,       // om
        R2      = 5000,     // om
        C       = 1e-12,    // F
        E_max   = 5.;       // v
int n   = 50;
double finish_time = 3e-6;
double h = finish_time / n;

double E(double t) {
    if (t < 1e-6 || t > 2.2e-6) return 0;
    if (t < 1.1e-6) return E_max/0.1e-6*(t - 1e-6);  
    if (t < 2.1e-6) return E_max;
    return -E_max/0.1e-6*(t - 2.2e-6);  
}

vector<double> iteration(vector<vector<double>> a,vector<double> b)
{
    int N = b.size();
    vector<double> x(N);
    int i,j;
    double norma;
    vector<double> xn(N);
    do{
        norma=0.0;
        for(i=0;i < N;i++) {
            xn[i]=-b[i];
            
            for(j=0;j < N;j++) {
                if(i!=j)
                 xn[i]+=a[i][j]*x[j];
            }
            
            xn[i]/=-a[i][i];
        }
        for(i=0;i < N;i++) {
            if(abs(x[i]-xn[i]) > norma)
            norma=abs(x[i]-xn[i]);
            x[i]=xn[i];
        }
    }
    while(norma > EPS);
 return x;
}

vector<double> kramer(vector<vector<double>> a,vector<double> f) {
    vector<double> ans(f.size());
    ans[0] = (-f[0] * a[1][1] + f[1]* a[0][1])/ (a[0][0]*a[1][1] - a[1][0]*a[0][1]);
    ans[1] = (-f[1] * a[0][0] + f[0]*a[1][0]) / (a[0][0]*a[1][1] - a[1][0]*a[0][1]);
    return ans;
}



vector<double> f_implicit(const vector<vector<double>>& u, double t) {
    vector<double> f(u.back().size());
    double u0 = u.back()[0];
    double u1 = u.back()[1];
    double u0_prev = u[u.size() - 2][0];
    double u1_prev = u[u.size() - 2][1];
    f[0] =  -((E(t) - u0)/R1 - C*((u0 - u0_prev) - u1 + u1_prev)/h - I0*(exp(u0/phit) - 1.));
    f[1] =  -(C*((u0 - u0_prev) - (u1 - u1_prev))/h - u1/R2);
    return f;
}

vector<vector<double>> A_implicit(const vector<vector<double>>& u) {
    vector<vector<double>> A(2,vector<double>(2));
    double u0 = u.back()[0];
    A[0][0] = -1./R1 - C/h - I0/phit*exp(u0/phit); 
    A[0][1] = C/h;
    A[1][0] = C/h;
    A[1][1] = -C/h - 1./R2;
    return A;
}

vector<vector<double>> A_explicit() {
    vector<vector<double>> A(2,vector<double>(2));
    A[0][0] = - C/h; 
    A[0][1] = C/h;
    A[1][0] = C/h;
    A[1][1] = -C/h;
    return A;
}

vector<double> f_explicit(const vector<vector<double>>& u, double t){
    vector<double> f(u.back().size());
    double u0 = u.back()[0];
    double u1 = u.back()[1];
    double u0_prev = u[u.size() - 2][0];
    double u1_prev = u[u.size() - 2][1];
    f[0] =  -((E(t - h) - u0_prev)/R1 - C*((u0 - u1)  - u0_prev + u1_prev)/h - I0*(exp(u0_prev/phit) - 1.));
    f[1] = -(C*((u0 - u1)  - u0_prev + u1_prev)/h - u1_prev/R2);
    return f;
}


int main()
{
    freopen("output.csv", "w", stdout);
    vector<vector<double>> u;
    u.push_back(vector<double>(2));
    for(double t = h; t <= finish_time; t += h) {
        u.push_back(vector<double>(2));
        u.back()[0] = u[u.size() - 2][0];
        u.back()[1] = u[u.size() - 2][1];
        vector<double> u_prev(2);
        u_prev[0] = u[u.size() - 2][0];
        u_prev[1] = u[u.size() - 2][1];
        int iter = 0;
        vector<double> x(2);

        do {
            iter++;
            double u0 = u.back()[0];
            double u1 = u.back()[1];
            double u0_prev = u[u.size() - 2][0];
            double u1_prev = u[u.size() - 2][1];
            vector<vector<double>> A = A_implicit(u);
            //cout << '\t' <<  iter << ' ' << u.back()[0] << endl; 
            vector<double> f = f_implicit(u,t);
            double r = 0.01, a1 = 1, a2 = 1;
            // if (abs(x[0]) > r &&  abs(x[1]) > r)  {
            //     a1 = r/ x[0];
            //     a2 = r/ x[1];
            // }
            // else {
            //     a1 = 1;
            //     a2 = 1;
            // }
            x = iteration(A,f);
            //cout << '\t' << iter << '\t' << x[0] << '\t' << x[1] << endl;
            
            u_prev[0] = u.back()[0];
            u_prev[1] = u.back()[1];
            u.back()[0]  += a1*x[0];
            u.back()[1]  += a2*x[1];
           
        }
        while(abs(u.back()[0] - u_prev[0]) > EPS && abs(u.back()[1] - u_prev[1]) > EPS);
        // vector<vector<double>> A = A_implicit(u);
        // do {
        //     iter++;
        //     double u0 = u.back()[0];
        //     double u1 = u.back()[1];
        //     double u0_prev = u[u.size() - 2][0];
        //     double u1_prev = u[u.size() - 2][1];
            
        //     //cout << '\t' <<  iter << ' ' << u.back()[0] << endl; 
        //     vector<double> f = f_implicit(u,t);
            
        //     vector<double> x = iteration(A,f);
        //     //cout << '\t' << iter << '\t' << x[0] << '\t' << x[1] << endl;
        //     u_prev[0] = u.back()[0];
        //     u_prev[1] = u.back()[1];
        //     u.back()[0]  += x[0];
        //     u.back()[1]  += x[1];
        // }
        // while(abs(u.back()[0] - u_prev[0]) > EPS && abs(u.back()[1] - u_prev[1]) > EPS);

        cout << t  << "," <<  u.back()[0]  << "," << u.back()[1]  << endl;
    }
    
    return 0;
}


