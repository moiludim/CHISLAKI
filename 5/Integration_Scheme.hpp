#include "Point.hpp"
#include <vector>
class Integration_Scheme{
    protected:
        std::vector<Point> points;
        std::vector<double> weights;
    public:
    enum Integration_Scheme_Type{
        Gauss_1,
        Gauss_2,
        Gauss_3,
        Gauss_4,
        Simpson,
        Trapezoidal,
    };       
}