#include <iostream>
/*
*
* Author: Zadorozhnyi Pavel
* date:   13.03.21
*
*/

const double  R = 1, // Om
              C = 1, // pF
              U0 = 5; // V

void analytical_solution(double h );
void explicit_euler_method(double h );
void implicit_euler_method(double h );
void trapezoidal_rule(double h );