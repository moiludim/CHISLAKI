#include <vector>
#include <cmath>
#include <iostream>

std::vector<double> generation(double a, double b, double r, int n){
    std::vector<double> result;
    result.reserve(n + 1);

    if (std::abs(r - 1.0) < 1e-12) {
        double h = (b - a) / n;
        for (int i = 0; i <= n; ++i)
            result.push_back(a + i * h);
        return result;
    }

    std::vector<double> prog(n + 1, 0.0);
    double power = 1.0;
    for (int i = 1; i <= n; ++i) {
        power *= r;
        prog[i] = prog[i-1] + power;
    }

    double sum = prog[n];
    double h = (b - a) / sum;

    for (int i = 0; i <= n; ++i)
        result.push_back(a + h * prog[i]);

    return result;
}