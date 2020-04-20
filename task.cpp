#include <stdio.h>
#include <ctime>
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include <chrono>

using namespace std::chrono;

// количество строк в исходной квадратной матрице
const int MATRIX_SIZE = 1500;

/// Функция InitMatrix() заполняет переданную в качестве 
/// параметра квадратную матрицу случайными значениями
/// matrix - исходная матрица СЛАУ
void InitMatrix(double** matrix)
{
	for (int i = 0; i < MATRIX_SIZE; ++i)
	{
		matrix[i] = new double[MATRIX_SIZE + 1];
	}

	for (int i = 0; i < MATRIX_SIZE; ++i)
	{
		for (int j = 0; j <= MATRIX_SIZE; ++j)
		{
			matrix[i][j] = rand() % 2500 + 1;
		}
	}
}

/// Функция SerialGaussMethod() решает СЛАУ методом Гаусса 
/// matrix - исходная матрица коэффиициентов уравнений, входящих в СЛАУ,
/// последний столбей матрицы - значения правых частей уравнений
/// rows - количество строк в исходной матрице
/// result - массив ответов СЛАУ
duration<double> SerialGaussMethod(double **matrix, const int rows, double* result)
{
	int k;
	double koef;
	high_resolution_clock::time_point start = high_resolution_clock::now();
	// прямой ход метода Гаусса
	for (k = 0; k < rows; ++k)
	{
		//
		for (int i = k + 1; i < rows; ++i)
		{
			koef = -matrix[i][k] / matrix[k][k];

			for (int j = k; j <= rows; ++j)
			{
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}
	high_resolution_clock::time_point finish = high_resolution_clock::now();
	duration<double> duration = (finish - start);
	// обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];

	for (k = rows - 2; k >= 0; --k)
	{
		result[k] = matrix[k][rows];

		//
		for (int j = k + 1; j < rows; ++j)
		{
			result[k] -= matrix[k][j] * result[j];
		}

		result[k] /= matrix[k][k];
	}

	return duration;
}

duration<double> ParallGaussMethod(double **matrix, const int rows, double* result)
{
	int k;
	high_resolution_clock::time_point start = high_resolution_clock::now();
	// прямой ход метода Гаусса
	for (k = 0; k < rows; ++k)
	{

		cilk_for(int i = k + 1; i < rows; ++i)
		{
			double koef = -matrix[i][k] / matrix[k][k];

			for (int j = k; j <= rows; ++j)
			{
				matrix[i][j] += koef * matrix[k][j];
			}
		}
	}
	high_resolution_clock::time_point finish = high_resolution_clock::now();
	duration<double> duration = (finish - start);
	// обратный ход метода Гаусса
	result[rows - 1] = matrix[rows - 1][rows] / matrix[rows - 1][rows - 1];


	for (k = rows - 2; k >= 0; --k)
	{
		cilk::reducer_opadd<double> result_k(matrix[k][rows]);

		cilk_for(int j = k + 1; j < rows; ++j)
		{
			result_k -= matrix[k][j] * result[j];
		}

		result[k] = result_k->get_value() / matrix[k][k];

	}

	return duration;
}


int main()
{
	srand((unsigned)time(0));

	int i;

	// кол-во строк в матрице, приводимой в качестве примера
	//const int test_matrix_lines = 4;
	//Меняем на:
	const int test_matrix_lines = MATRIX_SIZE;

	double **test_matrix = new double*[test_matrix_lines];
	double **test_matrix_parall = new double*[test_matrix_lines];

	// цикл по строкам
	for (i = 0; i < test_matrix_lines; ++i)
	{
		// (test_matrix_lines + 1)- количество столбцов в тестовой матрице,
		// последний столбец матрицы отведен под правые части уравнений, входящих в СЛАУ
		test_matrix[i] = new double[test_matrix_lines + 1];
		test_matrix_parall[i] = new double[test_matrix_lines + 1];
	}

	// массив решений СЛАУ
	double *result = new double[test_matrix_lines];
	double *result_parall = new double[test_matrix_lines];

	// инициализация тестовой матрицы
	/*test_matrix[0][0] = 2; test_matrix[0][1] = 5;  test_matrix[0][2] = 4;  test_matrix[0][3] = 1;  test_matrix[0][4] = 20;
	test_matrix[1][0] = 1; test_matrix[1][1] = 3;  test_matrix[1][2] = 2;  test_matrix[1][3] = 1;  test_matrix[1][4] = 11;
	test_matrix[2][0] = 2; test_matrix[2][1] = 10; test_matrix[2][2] = 9;  test_matrix[2][3] = 7;  test_matrix[2][4] = 40;
	test_matrix[3][0] = 3; test_matrix[3][1] = 8;  test_matrix[3][2] = 9;  test_matrix[3][3] = 2;  test_matrix[3][4] = 37; */
	//меняем на:
	InitMatrix(test_matrix);

	duration<double>duration_serial = SerialGaussMethod(test_matrix, test_matrix_lines, result);
	
	duration<double>duration_parall = ParallGaussMethod(test_matrix_parall, test_matrix_lines, result_parall);

	for (i = 0; i < test_matrix_lines; ++i)
	{
		delete[]test_matrix[i];
		delete[]test_matrix_parall[i];
	}

	printf("Solution:\n");
	printf("Serial:\n");
	printf("Time serial is:: %lf sec\n\n", duration_serial.count());
	for (i = 0; i < test_matrix_lines; ++i)
	{
		//printf("x(%d) = %lf\n", i, result[i]);
	}
	
	printf("\n\nParall:\n");
	printf("Time parall is:: %lf sec\n\n", duration_parall.count());
	for (i = 0; i < test_matrix_lines; ++i)
	{
		//printf("x(%d) = %lf\n", i, result_parall[i]);
	}

	printf("Boost is:: %lf\n\n", duration_serial.count()/duration_parall.count());
	
	delete[] result;
	delete[] result_parall;

	return 0;
}