#pragma once
#include "Integration_Scheme.hpp"
#include <functional>

class Integration_Scheme_Interval : protected Integration_Scheme {
public:
    Integration_Scheme_Interval(Integration_Scheme_Type type);
    double calculate_integral(const Point &A, const Point &B, int Number_Segments,
                              const std::function<double(const Point &P)> &Func) const;
};