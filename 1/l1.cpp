#include <iostream>
#include "Cubic_Spline.hpp"
#include "genh.hpp"
#include <cmath>
#include <iomanip>
const double a = 0.09;
const double b = 0.36;
const double eps = 0.0001;
const double r = 1.0;
const int n1 = 3;
const int n2 = 6;
const int n3 = 12;

std::vector<double> fmain(const std::vector<double>& x){
    int n = x.size();
    std::vector<double> result(n);
    for(int i = 0; i<n; i++){
        result[i] = tan(x[i]);
    }
    return result;
}
std::vector<double> dfmain(const std::vector<double>& x){
    int n = x.size();
    std::vector<double> result(n);
    for(int i = 0; i<n; i++){
        result[i] = 1/pow(cos(x[i]), 2);
    }
    return result;
}
std::vector<double> ddfmain(const std::vector<double>& x){
    int n = x.size();
    std::vector<double> result(n);
    for(int i = 0; i<n; i++){
        result[i] = 2*sin(x[i])/pow(cos(x[i]), 3);
    }
    return result;
}



//генерация по сетке точек с x и y,z = 0
std::vector<Point> genP(const std::vector<double>& x){
    int n = x.size();
    std::vector<Point> result;
    result.reserve(n);
    for(int i = 0; i<n; i++){
        result.emplace_back(x[i], 0.0, 0.0);
    }
    return result;
}


double max(std::vector<double>& vec){
    if (vec.empty()){throw std::invalid_argument("findMax: empty vector");}
    
    double max = vec[0];

    for(int i = 1; i<vec.size(); i++){
        double temp = vec[i];
        if(max<temp){
            max = temp;
        }
    }
    return max;
}



void test(int n){

    std::vector<double> sample = generation(a,b,r,n);

    std::vector<double> f_value = fmain(sample);

    std::vector<Point> P = genP(sample);

    Cubic_Spline spline;
    spline.Update_Spline(P, f_value);
    double Res[3];
    std::vector<double> sample_new;
    //формируем сетку не попадающая в узлы
    for(int i = 0; i + 1 < sample.size(); i++){
        sample_new.push_back(0.5 * (sample[i] + sample[i+1]));
    }

    std::vector<double> deltabig_f;
    std::vector<double> deltabig_df;
    std::vector<double> deltabig_ddf;
    std::vector<Point> P_new = genP(sample_new);
    std::vector<double> f_new   = fmain(sample_new);
    std::vector<double> df_new  = dfmain(sample_new);
    std::vector<double> ddf_new = ddfmain(sample_new);

    std::cout << std::scientific << std::setprecision(6);
    for(int i = 0; i < sample_new.size(); i++){
        spline.Get_Value(P_new[i], Res);

        std::cout << "x = " << sample_new[i]
                << "\tg(x) = "  << Res[0]     << "\tf(x) = "  << f_new[i]
                << "\tg'(x) = " << Res[1]     << "\tf'(x) = " << df_new[i]
                << "\tg''(x) = "<< Res[2]     << "\tf''(x) = "<< ddf_new[i] << "\n";

        deltabig_f.push_back(std::fabs(Res[0] - f_new[i]));
        deltabig_df.push_back(std::fabs(Res[1] - df_new[i]));
        deltabig_ddf.push_back(std::fabs(Res[2] - ddf_new[i]));
    }
    
    std::cout<<"delta(f-g) = "<<max(deltabig_f)<<"\t";
    std::cout<<"delta(f'-g') = "<<max(deltabig_df)<<"\t";
    std::cout<<"delta(f''-g'') = "<<max(deltabig_ddf)<<"\t\n\n";
}


int main(){

    test(n1);
    test(n2);
    test(n3);

}