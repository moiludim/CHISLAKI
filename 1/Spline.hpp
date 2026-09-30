#pragma once
#include <vector>
#include "Point.hpp"
class Spline{
    public:
    virtual void Update_Spline(const std::vector<Point> &Point, const std::vector<double> &F_value) = 0;
    virtual void Get_Value(const Point &P, double *Res)const = 0;
};