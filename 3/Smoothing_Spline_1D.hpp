#pragma once
#include "Spline.hpp"


//сглаживающий сплайн
class Smoothing_Spline_1D : public Spline
{
private:
	//параметр сглаживания
	double SMOOTH;
	//точки сетки
	std::vector<Point> Points;
	//кэфы по по базисным функциям
	std::vector<double> alpha;
	//переход на эталонный отрезок [-1, 1]: 
	//Seg_Num - номер сегмента, Х - абсцисса, Ksi - координата на эталонном отрезке
	void Transition_To_Master_Element(int Seg_Num, const double &X, double &Ksi) const;
	//базисные функции на [-1, 1]:
	//Number - номер функции, Ksi - координата на эталонном отрезке
	double Basis_Function(int Number, const double &Ksi) const;
	//производные базисных функций на [-1, 1]:
	//Number - номер функции, Ksi - координата на эталонном отрезке
	double Der_Basis_Function(int Number, const double & Ksi) const;
public:
	//конструктор: SMOOTH - параметр сглаживания
	Smoothing_Spline_1D(const double &SMOOTH);
	//обновление сплайна
	void Update_Spline(const std::vector<Point> &Points, const std::vector<double> &F_Value) override;
	//получение значения сплайна в точке P
	void Get_Value(const Point &P, double * Res)const override;
};


