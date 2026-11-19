#include "genh.hpp"


std::vector<Point> GenerateNormalSequence(int N, double M, double sigma)
{
    std::vector<Point> sequence;
    sequence.reserve(N);

    std::random_device rd;
    std::mt19937 gen(rd());

    std::normal_distribution<double> dist(M, sigma);

    for (int i = 0; i < N; i++)
    {
        double y = dist(gen);

        sequence.emplace_back(i+1, y, 0.0);
    }

    return sequence;
}

