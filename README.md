# Численные методы линейной алгебры
Данный репозиторий содержит реализацию библиотеки для работы с линейной алгеброй. Библиотека определяет базовые классы для работы с комплексными числами, векторами и матрицами, а также некоторые вычислительно-устойчивые алгоритмы, в том числе алгоритм по поиску сингулярного разложения.

## Сборка
Проект собирается с помощью системы `Cmake`, минимальная необходимая версия -- `3.21.2`. Пример инструкций для сборки:
```
cmake -B build-cmake .
cmake --build build-cmake
```

В результате будет собрана библиотека `course_project_lib`. Помимо этого, будет собран исполняемый файл `course_project_test` с автоматическими тестами. 

## Работа с библиотекой

Базовые классы -- `Complex`, `Vector<Type>` и `Matrix<Type>` определены в пространстве имен `::svd_computation`, и содержат различные конструкторы. Пример инстанциирования:

```C++
using namespace svd_computation;

Complex a = Complex(2.5, -4.5);
Complex b = exp_form(2, 0);

Vector<long double> v = {5, 6, 7, 8};

Matrix<long double> A = {{-4.5, -3.6, 7.8}, {5.2, -8.9, 1.6}, {-0.3, 4, 5}};
```

Различные вычислительно-устойчивые алгоритмы также лежат в пространстве имен `::svd_computation`. Среди них:

* `void orthonormalize(std::vector<Vector<Type>>& system, const long double eps)` -- ортонормализация Грамма-Шмидта
* `std::pair<Matrix<Type>, Matrix<Type>> get_QR_decomposition(const Matrix<Type>& A, const long double eps)` -- для вычисления QR-разложения
* `Matrix<Type> bidiagonalize(const Matrix<Type>& A, Matrix<Type>* left_basis, Matrix<Type>* right_basis, const long double eps)` -- для вычисления бидиагонализации
* `std::vector<long double> compute_svd(const Matrix<Type>& A, Matrix<Type>* left_basis, Matrix<Type>* right_basis, const long double eps)` -- для вычисления сингулярного разложения

Больше примеров использования базовых классов и алгоритмов можно найти в [тестах](course-project-second-year/tests).
