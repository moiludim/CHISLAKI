#include <iostream>
#include <vector>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <algorithm>

#include "Point.hpp"
#include "genh.hpp"
#include "Cubic_Spline.hpp"
#include "Smoothing_Spline_1D.hpp"

const int N = 1963;
const double M = 0.97;
const double sigma = 5.24;
const int nsmall = 200;

const int PLOT_POINTS_PER_SEGMENT = 10;

void SendSample(FILE* gp, const std::vector<Point>& sample)
{
    const std::size_t count =
        std::min<std::size_t>(nsmall, sample.size());

    for (std::size_t i = 0; i < count; ++i)
    {
        fprintf(gp, "%.10f %.10f\n", sample[i].x(), sample[i].y());
    }

    fprintf(gp, "e\n");
}

template <typename SplineType>
void SendSpline(
    FILE* gp,
    const std::vector<Point>& nodes,
    const SplineType& spline)
{
    const std::size_t count =
        std::min<std::size_t>(nsmall, nodes.size());

    if (count < 2)
    {
        fprintf(gp, "e\n");
        return;
    }

    double result[3];

    for (std::size_t i = 0; i + 1 < count; ++i)
    {
        const double x1 = nodes[i].x();
        const double x2 = nodes[i + 1].x();

        for (int j = 0; j < PLOT_POINTS_PER_SEGMENT; ++j)
        {
            const double t =
                static_cast<double>(j) / PLOT_POINTS_PER_SEGMENT;

            const double x = x1 + (x2 - x1) * t;
            const Point point(x, 0.0, 0.0);

            spline.Get_Value(point, result);

            fprintf(gp, "%.10f %.10f\n", x, result[0]);
        }
    }

    spline.Get_Value(nodes[count - 1], result);

    fprintf(gp, "%.10f %.10f\n",
            nodes[count - 1].x(), result[0]);

    fprintf(gp, "e\n");
}

void PlotSplines(
    const std::vector<Point>& sample,
    const std::vector<Point>& nodes,
    const Cubic_Spline& interpolationSpline,
    const Smoothing_Spline_1D& smoothing04,
    const Smoothing_Spline_1D& smoothing08,
    const Smoothing_Spline_1D& smoothing099,
    const Smoothing_Spline_1D& smoothing0)
{
    FILE* gp = _popen(
        "\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persistent",
        "w"
    );

    if (!gp)
        throw std::runtime_error("Не удалось запустить gnuplot");

    fprintf(gp, "set terminal wxt size 1500,850 font 'Arial,9'\n");

    fprintf(gp,
        "set title 'Интерполяционный и сглаживающие сплайны' "
        "font 'Arial,11'\n");

    fprintf(gp, "set xlabel 'Номер наблюдения' font 'Arial,9'\n");
    fprintf(gp, "set ylabel 'Случайное значение' font 'Arial,9'\n");

    fprintf(gp, "set xtics font 'Arial,8'\n");

    fprintf(gp, "set ytics font 'Arial,8'\n");

    fprintf(gp, "set grid\n");

    fprintf(gp, "set key outside right font 'Arial,8' spacing 1.1\n");

    fprintf(gp, "set lmargin 9\n");
    fprintf(gp, "set rmargin 5\n");
    fprintf(gp, "set tmargin 3\n");
    fprintf(gp, "set bmargin 4\n");

    fprintf(gp, "set xrange [1:50]\n");
    fprintf(gp, "set yrange [-40:40]\n");

    fprintf(gp,
        "plot '-' using 1:2 with points pt 7 ps 1.0 lc rgb 'black' "
        "title 'Исходная выборка', "
        "'-' using 1:2 with lines lw 1.5 lc rgb 'red' "
        "title 'Интерполяционный сплайн', "
        "'-' using 1:2 with lines lw 1.5 lc rgb 'blue' "
        "title 'Сглаживающий сплайн p=0.4', "
        "'-' using 1:2 with lines lw 1.5 lc rgb 'green' "
        "title 'Сглаживающий сплайн p=0.8', "
        "'-' using 1:2 with lines lw 1.5 lc rgb 'orange' "
        "title 'Сглаживающий сплайн p=0.99', "
        "'-' using 1:2 with lines lw 1.5 lc rgb 'black' "
        "title 'Сглаживающий сплайн p=0'\n"
    );

    SendSample(gp, sample);

    SendSpline(gp, nodes, interpolationSpline);
    SendSpline(gp, nodes, smoothing04);
    SendSpline(gp, nodes, smoothing08);
    SendSpline(gp, nodes, smoothing099);
    SendSpline(gp, nodes, smoothing0);
    fflush(gp);
    _pclose(gp);
}

int main()
{
    try
    {
        std::vector<Point> sample =
            GenerateNormalSequence(N, M, sigma);

        if (sample.size() < 2)
            throw std::runtime_error("Недостаточно точек для сплайнов");

        std::vector<Point> nodes;
        std::vector<double> values;

        nodes.reserve(sample.size());
        values.reserve(sample.size());

        for (const Point& p : sample)
        {
            nodes.emplace_back(p.x(), 0.0, 0.0);
            values.push_back(p.y());
        }

        Cubic_Spline interpolationSpline;
        interpolationSpline.Update_Spline(nodes, values);

        Smoothing_Spline_1D smoothing04(0.4);
        smoothing04.Update_Spline(nodes, values);

        Smoothing_Spline_1D smoothing08(0.8);
        smoothing08.Update_Spline(nodes, values);

        Smoothing_Spline_1D smoothing099(0.99);
        smoothing099.Update_Spline(nodes, values);

        Smoothing_Spline_1D smoothing0(0.0);
        smoothing0.Update_Spline(nodes, values);

        PlotSplines(
        sample,
        nodes,
        interpolationSpline,
        smoothing04,
        smoothing08,
        smoothing099,
        smoothing0
        );

        std::cout << "График построен.\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "Ошибка: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
