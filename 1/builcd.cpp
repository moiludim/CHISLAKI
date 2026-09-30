#include <iostream>
#include <cstdlib>

int main(){
    int build_result = system(
        "g++ l1.cpp genh.cpp Point.cpp Cubic_Spline.cpp "
        "-o test.exe"
    );

    if (build_result != 0) {
        std::cout << "\nОшибка сборки! Проверьте, что все .cpp файлы лежат рядом.\n";
        return 1;
    }

    std::cout << "\nСборка успешна, запуск...\n\n";
    system(".\\test.exe");
    return 0;
}