
/*
*
* Author: Zadorozhnyi Pavel
* date:   28.04.21
*
*/

#include <iostream>
#include <functional>
#include <algorithm>
#include <math.h>

using namespace std;

const double    E = 5,
                R = 1,
                I0 = 1.2e-10,
                phit = 0.025,
                EPS = 1e-3;

double func(double x) {
    return (E - x)/R - I0*(exp(x/phit)  - 1.);
}

void bin_search () {
    double l = 0, r = 6;
    int i = 0;
    while(abs(func(r) - func(l)) > EPS) {
        double m = (l + r)/2;
        cout << i++ << '\t' <<  m << ' '<< func(m) << endl;
        if (func(m)*func(r) >= 0) {
            r = m;
        } else {
            l = m;
        }
    }
}

void golden_section_search () {
    
    double gr = (sqrt(5.) + 1.) / 2.;
    double a = 0, b = 6;
    int i = 0;

	double x1, x2;  
	while (abs(func(b) - func(a)) > EPS) {
	    x1 = b - (b - a) / gr; 
        x2 = a + (b - a) / gr;
        cout << i++ << '\t' << (a + b )/2 << ' ' << func((a+b) / 2) << endl;
        if (abs(func(x1)) > abs(func(x2))) 
           a = x1; 
       else 
           b = x2;
	}
}

void step_doubling_method() {
    double l = 0, r = 6;
    double h0 = EPS/2;
    double h = h0; 
   
    int i = 0;
    while(abs(r - l) > EPS) {
        cout << i++ << '\t' <<  h << ' ' << l << ' ' << func(l) <<  endl;
        if (i == 200) break;
        if (abs(func(l)) > abs(func(l + h)) && h + l < r) {
            l += h;
            h *= 2;
        } else {
            r = h + l;
            h = h0;
        }
    }
}

int main () {
    
    cout << "\n Bin Search \n";
    bin_search();
     
     cout << "\n golden_section_search \n";
     golden_section_search();
     cout << "\n Step doubling \n";
     step_doubling_method();
     return 0;
}

