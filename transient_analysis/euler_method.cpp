#include "euler_method.hpp"
#include <math.h>

using namespace std;

void analytical_solution(double h = 1./(R*C)) {
    freopen("analytical_solution.txt", "w", stdout);

    cout << "t \t u \n";
    for(int n = 0; n < 30; n++) {
        double t = n*h;
        double u = U0*exp(-t/(R*C));
        cout << t << '\t' << u << endl;
    }
}
void explicit_euler_method(double h = 1./(R*C)) {
    freopen("explicit_euler_method.txt", "w", stdout);
    double u = U0;
    double edge = 4*15./h;
    cout << "t \t u \n";
    for(double t = 0; t < edge; t += h) {
        cout << t << '\t' << u << endl;
        u *= (1 - h/(R*C));
    }
}
void implicit_euler_method(double h = 1./(R*C)) {
    freopen("implicit_euler_method.txt", "w", stdout);
    double u = U0;
    double edge = 4*15./h;
    cout << "t \t u \n";
    for(double t = 0; t < edge; t += h) {
        cout << t << '\t' << u << endl;
        u *= 1./(1. + h/(R*C));
    }
}
void trapezoidal_rule(double h = 1./(R*C)) {
    freopen("trapezoidal_rule.txt", "w", stdout);
    double u = U0;
    double edge = 4*15./h;
    cout << "t \t u \n";
    for(double t = 0; t < edge; t += h) {
        cout << t << '\t' << u << endl;
        u *= (1. - h/(R*C) + 1./(1. + h/(R*C)))/2.;
    }
}

int main() {
    double h = 2/(R*C);
    analytical_solution(h);
    explicit_euler_method(h);
    implicit_euler_method(h);
    trapezoidal_rule(h);
    
    return EXIT_SUCCESS;
}
