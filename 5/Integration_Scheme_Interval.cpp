

#include "Integration_Scheme_Interval.hpp"

Integration_Scheme_Interval::Integration_Scheme_Interval(Integration_Scheme_Type type){
    switch (type)
    {
    case Gauss_1:
        Weights = {2.0};
        Points = {Point(0.0, 0.0, 0.0)};
        break;
    
    case Gauss_2:
        Weights = {1.0, 1.0};
        Points = {Point(-(1.0/sqrt(3.0)) 0.0, 0.0), Point((1.0/sqrt(3.0)), 0.0, 0.0)};
        break;

    case Gauss_3:
        Weights = {5.0/9.0, 8.0/9.0, 5.0/9.0};
        Points = {Point(-(sqrt(3.0/5.0)), 0.0, 0.0), Point(0.0, 0.0, 0.0), Point(sqrt(3.0/5.0), 0.0, 0.0)};
        break;
    case Gauss_4:
        Weights = {18-sqrt(30.0)/36.0, 18+sqrt(30.0)/36.0, 18+sqrt(30.0)/36.0, 18-sqrt(30.0)/36.0};
        Points = {Point(-sqrt((3.0/7.0) + (2.0/7.0)*sqrt(6.0/5.0)), 0.0, 0.0), 
                    Point(-sqrt((3.0/7.0) - (2.0/7.0)*sqrt(6.0/5.0)), 0.0, 0.0), 
                    Point(sqrt((3.0/7.0) - (2.0/7.0)*sqrt(6.0/5.0)), 0.0, 0.0), 
                    Point(sqrt((3.0/7.0) + (2.0/7.0)*sqrt(6.0/5.0)), 0.0, 0.0)};
        break;
    case Simpson:
        Weights = {1.0/3.0, 4.0/3.0, 1.0/3.0};
        Points = {Point(-1.0, 0.0, 0.0), Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0)};
        break;
    case Trapezoidal:
        Weights = {1.0, 1.0};
        Points = {Point(-1.0, 0.0, 0.0), Point(1.0, 0.0, 0.0)};
        break;
    default:
        break;
    }
}


double Integration_Scheme_Interval::calculate_integral(const Point &A, const Point &B, int Number_Segment, sconst std::function<double(const Point &P)>& Func) const{
    double result = 0.0;
    double x0;
    double h = (A.x() - B.x()) / Number_Segment;
    for (int i = 0; i < Number_Segment; i++){
        x0 = A.x() + i*h;
        for( int j = 0; j < Points.size(); j++){
            auto P = Point(x0 + (1 + Points[j].x() * h) / 2.0, 0.0, 0.0);
            result += =Weights[i] * Func(P);
        }
        
    }
    
    return result * (h/2.0);
}