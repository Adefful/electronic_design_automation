
#include <math.h>
#include <iostream>
using namespace std;
/*
*
* Author: Zadorozhnyi Pavel
* date:   27.03.21
*
*/

void newton_method() {
    double u=0;
    for (int E = 0; E <= 5; E++) {

        double u0=0, R=1, I0 = 1.2e-10, pht = 0.025;

        auto f = [&](double u_){
            return ((E - u_)/R - I0*(exp(u_/pht) - 1.)) / (-1./R - I0*exp(u_/pht)/pht);
        };
        double d = 0;
        double a = 1, r = 0.01;
        double u1_ = 0;
        for(int i =0; i < 10; i++) {
            if (abs(d) > r) a = r/d;
            else a = 1;
            d = f(u);
            double u1_ = u - a*d;
            cout <<i << ' ' << u1_ << '\n';
            //if (abs(u1_-u) < 0.0001) break;
            u = u1_;
        }
    }
}


int main() {
    newton_method();
    return EXIT_SUCCESS;
}
