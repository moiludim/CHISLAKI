#pragma once
#include "Spline.hpp"

class Cubic_Spline : public Spline
{
    private:
        std::vector<Point> Points;
        std::vector<double> a,b,c,d;
    public:
        void Update_Spline(const std::vector<Point> &Point, const std::vector<double> &F_value) override;
        void Get_Value(const Point &P, double *Res)const override;
};