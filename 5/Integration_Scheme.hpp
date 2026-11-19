#include "Point.hpp"
#include <vector>
class Integration_Scheme{
    protected:
        std::vector<Point> Points;
        std::vector<double> Weights;
    public:
    enum Integration_Scheme_Type{
        Gauss_1,
        Gauss_2,
        Gauss_3,
        Gauss_4,
        Simpson,
        Trapezoidal,
    };       
};