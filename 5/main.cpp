#include <iostream>
#include <iomanip>
#include <cmath>
#include "Integration_Scheme_Interval.hpp"

using namespace std;

int main() {
    Point A(-0.08, 0.0, 0.0), B(0.39, 0.0, 0.0);
    int N0 = 1;

    auto f = [](const Point &P) { return exp(P.x()); };
    double Iex = exp(B.x()) - exp(A.x());
    double H = (B.x() - A.x()) / N0;

    const int count = 6;
    const char *names[count] = {"Trapezoidal", "Gauss_1", "Simpson", "Gauss_2", "Gauss_3", "Gauss_4"};
    Integration_Scheme::Integration_Scheme_Type types[count] = {
        Integration_Scheme::Trapezoidal,
        Integration_Scheme::Gauss_1,
        Integration_Scheme::Simpson,
        Integration_Scheme::Gauss_2,
        Integration_Scheme::Gauss_3,
        Integration_Scheme::Gauss_4
    };
    int orders[count] = {2, 2, 4, 4, 6, 8};

    cout << scientific << setprecision(4);
    cout << "f = exp(x), [" << A.x() << ", " << B.x() << "], I* = "
         << setprecision(15) << Iex << setprecision(4) << "\n";

    for (int s = 0; s < count; s++) {
        Integration_Scheme_Interval q(types[s]);
        int k = orders[s];

        double I[4];
        for (int m = 0; m < 4; m++)
            I[m] = q.calculate_integral(A, B, N0 << m, f);

        cout << "\n" << names[s] << " (k = " << k << ")\n";
        cout << left
             << setw(14) << "h"
             << setw(14) << "I*-I^h"
             << setw(14) << "ratio"
             << setw(14) << "Runge"
             << setw(14) << "I^R"
             << setw(14) << "I*-I^R" << "\n";

        double c = pow(2.0, k) - 1.0;
        double p[3];
        for (int m = 0; m < 3; m++) {
            double err0  = Iex - I[m];
            double err1  = Iex - I[m + 1];
            double ratio = err0 / err1;                  
            double runge = (I[m + 1] - I[m]) / c;        
            double IR    = I[m + 1] + runge; 
            p[m] = log2(fabs(ratio));

            cout << left
                 << setw(14) << H / (1 << m)
                 << setw(14) << err0
                 << setw(14) << ratio
                 << setw(14) << runge
                 << setw(14) << IR
                 << setw(14) << Iex - IR << "\n";
        }

        cout << fixed << setprecision(3)
             << "order p = " << p[0] << ", " << p[1] << ", " << p[2]
             << "; theory k = " << k << "\n"
             << scientific << setprecision(4);
    }
    return 0;
}