#include <iostream>
#include "Cubic_Spline.hpp"
#include "genh.hpp"
#include <cmath>
#include <iomanip>

const double a = 0.09;
const double b = 0.36;
const double eps = 0.0001;
const double r = 1.0;
const int n1 = 12;
const int n2 = 24;
const int n3 = 48;

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

void test3(int n){

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


std::vector<double> dddfmain(const std::vector<double>& x){
    int n = x.size();
    std::vector<double> result(n);
    for(int i = 0; i<n; i++){
        result[i] = 2*(1 + 4 * pow(tan(x[i]),2) + 3 * pow(tan(x[i]), 4));
    }
    return result;
}
//вычисляю шаг для выполнения точности епсилон для грубо если разделить минимальное из этих значений пополам, то можно использовать треточечные шаблоны левый правый и центральный
std::vector<double> h_value_MKR3(const std::vector<double>& uzel){
    std::vector<double> result;
    std::vector<double> dddf = dddfmain(uzel);

    for (int i = 0; i < uzel.size()-1; i++){
        result.push_back(sqrt(6 * eps / std::abs(dddf[i+1])));
    }

    return result;
}
//вычисление шага для двухточеченого шаблона 
std::vector<double> h_value_MKR2(const std::vector<double>& uzel){
    std::vector<double> result;
    std::vector<double> ddf = ddfmain(uzel);

    for (int i = 0; i < uzel.size()-1; i++){
        result.push_back(sqrt(2 * eps / std::abs(ddf[i+1])));
    }

    return result;
}

//метод конечных разностей трехточечный левая разность 
std::vector<double> MKR_Left(std::vector<double>& sample, std::vector<double>& f_value){
    std::vector<double> result;
    for(int i = 1; i<sample.size()-1; i++){
        double h = sample[i+1]-sample[i];
        result.push_back((f_value[i+1] - f_value[i]) / h);
    }
    return result;
}

//метод конечных разностей трехточечный центральный без крайних точек
std::vector<double> MKR_Center(std::vector<double>& sample, std::vector<double>& f_value){
    std::vector<double> result;
    for(int i = 1; i<sample.size()-2; i++){
        double hi = sample[i]-sample[i-1];
        double hi1 = sample[i+1]-sample[i];
        double temp =
            -hi1 / (hi * (hi + hi1)) * f_value[i-1]
            + (hi1 - hi) / (hi * hi1) * f_value[i]
            + hi / (hi1 * (hi + hi1)) * f_value[i+1];
        result.push_back(temp);
    }
    return result;
}

//метод конечных разности для второй производной
std::vector<double> ddMKR(std::vector<double>& sample, std::vector<double>& f_value){
    std::vector<double> result;
//всеравно для регулярной сетки 
    for(int i = 1; i<sample.size()-2; i++){
        double h = sample[i+1]-sample[i];
        result.push_back((f_value[i+1] - 2*f_value[i] + f_value[i-1]) / (h*h));
    }
    return result;
}

void test4(int n){

    std::vector<double> sample = generation(a, b, r, n);
    std::vector<double> f_value = fmain(sample);

    std::vector<double> df_value = dfmain(sample);
    std::vector<double> ddf_value = ddfmain(sample);

    std::vector<double> MKRL = MKR_Left(sample, f_value);
    std::vector<double> MKRC = MKR_Center(sample, f_value);
    std::vector<double> ddMKR_value = ddMKR(sample, f_value);

    std::vector<Point> P = genP(sample);

    Cubic_Spline spline;
    spline.Update_Spline(P, f_value);

    double Res[3];

    int size = sample.size();

    std::cout << std::scientific << std::setprecision(6);

    for(int i = 1; i < size - 2; i++){

        spline.Get_Value(P[i], Res);

        double error_Left =
            std::fabs(MKRL[i] - df_value[i]);

        double error_Center =
            std::fabs(MKRC[i - 1] - df_value[i]);

        double error_Spline_df =
            std::fabs(Res[1] - df_value[i]);

        double error_ddMKR =
            std::fabs(ddMKR_value[i - 1] - ddf_value[i]);

        double error_Spline_ddf =
            std::fabs(Res[2] - ddf_value[i]);

        std::cout << "x = " << sample[i]
                  << "\tf'(x) = " << df_value[i]
                  << "\n";

        std::cout << "\tMKR_Left   = " << MKRL[i]
                  << "\tError = " << error_Left
                  << "\n";

        std::cout << "\tMKR_Center = " << MKRC[i - 1]
                  << "\tError = " << error_Center
                  << "\n";

        std::cout << "\tSpline f'  = " << Res[1]
                  << "\tError = " << error_Spline_df
                  << "\n";

        std::cout << "\tf''(x) = " << ddf_value[i]
                  << "\n";

        std::cout << "\tddMKR      = " << ddMKR_value[i - 1]
                  << "\tError = " << error_ddMKR
                  << "\n";

        std::cout << "\tSpline f'' = " << Res[2]
                  << "\tError = " << error_Spline_ddf
                  << "\n\n";
    }
}


int main(){
   
    //test3(n1);
    //test3(n2);
    //test3(n3);


    test4(n2);
}