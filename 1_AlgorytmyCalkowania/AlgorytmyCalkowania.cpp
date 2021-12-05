// AlgorytmyCalkowania.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <chrono>
#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <windows.h>

using namespace std;

# define M_PI           3.14159265358979323846  /* pi */
double f(double x)
{
    return 3 * sin(2*M_PI * (x - 0.25)) + 3;
}

double multiCorefIntegralRectangleLeft(double lowerBound, double upperBound, int partitionNumber)
{
    double sum = 0;

    //obliczanie szerokości podprzedziału
    double partitionWidth = (upperBound - lowerBound) / partitionNumber;
    //END obliczanie szerokości podprzedziału

	#pragma omp parallel for default(none) firstprivate(partitionNumber, lowerBound, partitionWidth) reduction(+:sum)
    for (int i = 0; i < partitionNumber; i++)
    {
        sum += partitionWidth * f(lowerBound + partitionWidth * i);
    }

    return sum;
}

double fIntegralRectangleRight(double lowerBound, double upperBound, int partitionNumber)
{
    double sum = 0;

    //obliczanie szerokości podprzedziału
    double partitionWidth = (upperBound - lowerBound) / partitionNumber;
    //END obliczanie szerokości podprzedziału

    for (int i = 1; i <= partitionNumber; i++)
    {
        sum += partitionWidth * f(lowerBound + partitionWidth * i);
    }

    return sum;
}

double multiCorefIntegralRectangleRight(double lowerBound, double upperBound, int partitionNumber)
{
    double sum = 0;

    //obliczanie szerokości podprzedziału
    double partitionWidth = (upperBound - lowerBound) / partitionNumber;
    //END obliczanie szerokości podprzedziału

    #pragma omp parallel for default(none) firstprivate(partitionNumber, lowerBound, partitionWidth) reduction(+:sum)
    for (int i = 1; i <= partitionNumber; i++)
    {
        sum += partitionWidth * f(lowerBound + partitionWidth * i);
    }

    return sum;
}

double multiCorefIntegralTrapezoidal(double lowerBound, double upperBound, int partitionNumber)
{
    double sum = 0;

    //obliczanie szerokości podprzedziału
    double partitionWidth = (upperBound - lowerBound) / partitionNumber;
    //END obliczanie szerokości podprzedziału

	#pragma omp parallel for default(none) firstprivate(partitionNumber, lowerBound, partitionWidth) reduction(+:sum)
    for (int i = 1; i < partitionNumber; i++)
    {
        sum += f(lowerBound + partitionWidth * i);
    }
    sum = (sum + (f(lowerBound) + f(upperBound)) / 2) * partitionWidth;

    return sum;
}

double multiCorefIntegralSimpson(double lowerBound, double upperBound, int partitionNumber)
{
    double sum = 0;
    double sumMid = 0;

    //obliczanie szerokości podprzedziału
    double partitionWidth = (upperBound - lowerBound) / partitionNumber;
    //END obliczanie szerokości podprzedziału

    #pragma omp parallel for default(none) firstprivate(partitionNumber, lowerBound, partitionWidth, upperBound) reduction(+:sum) reduction(+:sumMid)
    for (int i = 1; i <= partitionNumber; i++)
    {
        sumMid += f((lowerBound + partitionWidth * i) - partitionWidth / 2);
        if (i < partitionNumber)
        {
            sum += f(lowerBound + partitionWidth * i);
        }
    }
    double result = partitionWidth / 6 * (f(lowerBound)+f(upperBound)+2*sum+4*sumMid);

    return result;
}

typedef double (*integralFunctionT)(double, double, int); //wskaźnik na funkcję całkującą

void multipleTestExecution(int testNumber, double lowerBound, double upperBound, int partitionNumber, integralFunctionT integralFunction )
{
    //badanie czasu obliczania
    vector < long long > executeTimes;
    auto totalTime = 0;
    for (int i = 1; i <= testNumber; i++)
    {
        cout << '\r' << i << '/' << testNumber;
        //pomiar czasu obliczania całki
        chrono::steady_clock::time_point start = chrono::steady_clock::now();
        ////////////////
        integralFunction(lowerBound, upperBound, partitionNumber);
        ////////////////
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        //END pomiar czasu obliczania całki

        executeTimes.push_back(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
        totalTime += chrono::duration_cast<chrono::nanoseconds>(end - start).count();

        //wyświetlanie wyniku i czasu obliczeń
        //cout << i << "--------------" << endl;
        //printf("Wynik: %f\n", result);
        //cout << "Rdzeni: " << omp_get_max_threads() << " Czas obliczeń: "
        //   << chrono::duration_cast<chrono::nanoseconds>(end - start).count()
        //    << "ns.\n";
        //cout << "Średni czas obliczeń: " << totalTime / i << "ns." << endl << endl;
        //END wyświetlanie wyniku, czasu obliczeń i czasu średniego
    }
    cout << "\r                     \r";
    //wyświetlanie posortowanych czasów obliczania
    //cout << "Średni czas obliczeń: " << totalTime / testNumber << "ns." << endl << endl;
    sort(executeTimes.begin(), executeTimes.end());
    //cout << "Posortowane czasy wykonań [ns]: \n";
    cout << omp_get_max_threads() << ' ';
    for (const auto& i : executeTimes)
        cout << i << ' ';
    cout << endl;
    //END wyświetlanie posortowanych czasów obliczania

    //END badanie czasu obliczania
}


int main()
{
    setlocale(LC_CTYPE, "Polish");

    //parametry programu
    omp_set_dynamic(0);
    double lowerBound = 10;
    double upperBound = 20;
    int partitionNumber[] = { 1000000, 2000000, 4000000, 8000000 };
    int testNumber = 100;
    int maxThreads = 12;
    //END parametry programu

    //czy górna granica jest na pewno większa od dolnej
    if (lowerBound > upperBound)
    {
        double tmp = lowerBound;
        lowerBound = upperBound;
        upperBound = tmp;
    }
    //END czy górna granica jest na pewno większa od dolnej

    //zebranie funkcji w tabelę
    integralFunctionT functionPointers[4];
    functionPointers[0] = multiCorefIntegralRectangleLeft;
    functionPointers[1] = multiCorefIntegralRectangleRight;
    functionPointers[2] = multiCorefIntegralTrapezoidal;
    functionPointers[3] = multiCorefIntegralSimpson;

    for (int j = 0; j <= 3; j++)
    {
        cout << "Metoda całkowania: ";
        switch (j)
        {
        case 0: cout << "Prostokątów z niedomiarem"; break;
        case 1: cout << "Prostokątów z nadmiarem"; break;
        case 2: cout << "Trapezów"; break;
        case 3: cout << "Wzór Simpsona"; break;
        }
        cout << endl;
        integralFunctionT currentFunction = functionPointers[j];

        for (int i = 0; i <= 3; i++)
        {
            cout << partitionNumber[i] << endl;
            for (int j = 1; j <= maxThreads; j++)
            {
                omp_set_num_threads(j);
                multipleTestExecution(testNumber, lowerBound, upperBound, partitionNumber[i], currentFunction);
            }
        }
    }

    getchar();
}