#include "Cubic_Spline.hpp"
#include <cmath>
#include <stdexcept>

void Cubic_Spline::Update_Spline(const std::vector<Point> &Points,
                                 const std::vector<double> &F_value)
    {
        //Мы в классе кубического сплайна очищаем private Points и присваеваем параметр Points;
        this->Points = Points;
        //Количество сегментов разбиения
        int n_segment = Points.size() - 1;

        //длина отрезков
		double h_current, h_next;

        
		a.resize(n_segment);
		b.resize(n_segment);
		c.resize(n_segment);
		d.resize(n_segment);
        //вектор fi
        std::vector<double> f(n_segment - 1);

        for (int i = 0; i < n_segment - 1; i++)
        {
            h_current = Points[i + 1].x() - Points[i].x();
            h_next    = Points[i + 2].x() - Points[i + 1].x();

            b[i] = 2 * (h_current + h_next);   // диагональ
            a[i + 1] = h_current;              // поддиагональ
            d[i] = h_next;                     // наддиагональ
            f[i] = 3.0 * ((F_value[i + 2] - F_value[i + 1]) / h_next
                        - (F_value[i + 1] - F_value[i]) / h_current);
        }
        //метод прогонки прямой ход
		for (int j = 1; j < n_segment - 1; j++)
		{
			b[j] -= a[j] / b[j - 1] * d[j - 1];
			f[j] -= a[j] / b[j - 1] * f[j - 1];
		}
        //обратный ход
		c[n_segment - 1] = f[n_segment - 2] / b[n_segment - 2];
		for (int j = n_segment - 2; j > 0; j--)
		{
			c[j] = (f[j - 1] - c[j + 1] * d[j - 1]) / b[j - 1];
		}

		c[0] = 0.0;


        for (int i = 0; i < n_segment - 1; i++)
		{
			h_current = Points[i + 1].x() - Points[i].x();
			a[i] = F_value[i];
			b[i] = (F_value[i + 1] - F_value[i]) / h_current - (c[i + 1] + 2.0 * c[i]) * h_current / 3.0;
			d[i] = (c[i + 1] - c[i]) / h_current / 3.0;
		}

		h_current = Points[n_segment].x() - Points[n_segment - 1].x();
		a[n_segment - 1] = F_value[n_segment - 1];
		b[n_segment - 1] = (F_value[n_segment] - F_value[n_segment - 1]) / h_current - 2.0 * c[n_segment - 1] * h_current / 3.0;
		d[n_segment - 1] = -c[n_segment - 1] / h_current / 3.0;
    }

void Cubic_Spline::Get_Value(const Point &P, double * Res)const
	{
		double eps = 1e-7;
		int Num_Segment = Points.size() - 1;
		for (int i = 0; i < Num_Segment; i++)
		{
			if (P.x() > Points[i].x() && P.x() < Points[i + 1].x() ||
				fabs(P.x() - Points[i].x()) < eps ||
				fabs(P.x() - Points[i + 1].x()) < eps)
			{
				double diff = (P.x() - Points[i].x());
				Res[0] = a[i] + b[i] * diff + c[i] * pow(diff, 2) + d[i] * pow(diff, 3);
				Res[1] = b[i] + 2.0 * c[i] * diff + 3.0 * d[i] * pow(diff, 2);
				Res[2] = 2.0 * c[i] + 6.0 * d[i] * diff;
				return;
			}
		}
		throw std::invalid_argument("The point is not found in the segments...");
	}