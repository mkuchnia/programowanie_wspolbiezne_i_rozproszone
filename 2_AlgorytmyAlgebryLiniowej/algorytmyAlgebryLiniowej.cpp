// algorytmyAlgebryLiniowej.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include <iostream>
#include <omp.h>
#include <time.h>
#include <vector>
#include <chrono>
#include <algorithm>
////////////////////////
#include <windows.h> 


using namespace std;
//#define N 1000  //rozmiar wektorów i macierzy
const int N = 1000;

float scalarProductVectorVector(float VectorA[N], float VectorB[N])
{
    float sum = 0;
    #pragma omp parallel for schedule(static) default(none) firstprivate(VectorA, VectorB) reduction(+:sum)
    for (int i = 0; i < N; i++)
    {
        sum += VectorA[i] * VectorB[i];
    }
    //cout << "scalarProduct: " << sum << endl;
    return sum;
}

void productMatrixVector(float Matrix[N][N], float Vector[N], float resultVector[N])
{
    #pragma omp parallel for schedule(static) default(none) firstprivate(Matrix, Vector)
    for (int i = 0; i < N; i++)
    {
        float sum = 0;
        for (int j = 0; j < N; j++)
        {
            sum += Matrix[i][j] * Vector[j];
        }
        resultVector[i] = sum;
        //cout << "resultVector[" << i << "]: " << sum << endl;
    }
}

void productMatrixMatrix(float MatrixA[N][N], float MatrixB[N][N], float resultMatrix[N][N])
{
    #pragma omp parallel for schedule(static) default(none) firstprivate(MatrixA, MatrixB)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            float sum = 0;
            for (int k = 0; k < N; k++)
            {
                sum += MatrixA[i][k] * MatrixB[k][j];
            }
            resultMatrix[i][j] = sum;
            //cout << "resultMatrix[" << i << "][" << j << "]: " << sum << endl;
        }
    }
}

float fRand()
{
    float fMin = -1000;
    float fMax = 1000;
    float f = (float)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main()
{
    omp_set_dynamic(0);
    srand(time(NULL));
    int testNumber = 100;
    int maxThreads = 12;

    //float Vector1[3] = { 1, 2, 3 };
    //float Vector2[3] = { 4, 5, 6 };
    //float Matrix1[3][3] = { {1,2,3},
    //                        {4,5,6},
    //                        {7,8,9} };
    //float Matrix2[3][3] = { {9,8,7},
    //                        {6,5,4},
    //                        {3,2,1} };
    static float Vector1[N];
    static float Vector2[N];
    static float Matrix1[N][N];
    static float Matrix2[N][N];
    static float resultVector[N];
    static float resultMatrix[N][N];
    //losowanie
    for (int i = 0; i < N; i++)
    {
        Vector1[i] = fRand();
        Vector2[i] = fRand();
        for (int j = 0; j < N; j++)
        {
            Matrix1[i][j] = fRand();
            Matrix2[i][j] = fRand();
        }
    }
    //END losowanie

    //iloczyn skalarny wektorów
    for (int i = 1; i <= maxThreads; i++)
    {
        omp_set_num_threads(i);
        //badanie czasu obliczania
        vector < long long > executeTimes;
        for (int i = 1; i <= testNumber; i++)
        {
            cout << '\r' << i << '/' << testNumber;
            //pomiar czasu obliczania iloczynu skalarnego
            chrono::steady_clock::time_point start = chrono::steady_clock::now();
            ////////////////
            scalarProductVectorVector(Vector1, Vector2);
            ////////////////
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            //END pomiar czasu obliczania iloczynu skalarnego

            executeTimes.push_back(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
        }
        cout << "\r                     \r";
        //wyświetlanie posortowanych czasów obliczania
        sort(executeTimes.begin(), executeTimes.end());
        cout << omp_get_max_threads() << ' ';
        for (const auto& i : executeTimes)
            cout << i << ' ';
        cout << endl;
        //END wyświetlanie posortowanych czasów obliczania

        //END badanie czasu obliczania
    }
    //END iloczyn skalarny wektorów

    for (int i = 1; i <= maxThreads; i++)
    {
        //cout << "********************************" << endl;
        omp_set_num_threads(i);
        //cout << "Wątków: " << omp_get_max_threads() << endl;
        //badanie czasu obliczania
        vector < long long > executeTimes;
        long long totalTime = 0;
        for (int i = 1; i <= testNumber; i++)
        {
            cout << '\r' << i << '/' << testNumber;
            //pomiar czasu obliczania iloczynu macierzy i wektora
            chrono::steady_clock::time_point start = chrono::steady_clock::now();
            ////////////////
            productMatrixVector(Matrix1, Vector1, resultVector);
            ////////////////
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            //END pomiar czasu obliczania iloczynu macierzy i wektora

            executeTimes.push_back(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
            totalTime += chrono::duration_cast<chrono::nanoseconds>(end - start).count();
            //cout << chrono::duration_cast<chrono::nanoseconds>(end - start).count() << endl;
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

    for (int i = 1; i <= maxThreads; i++)
    {
        //cout << "********************************" << endl;
        omp_set_num_threads(i);
        //cout << "Wątków: " << omp_get_max_threads() << endl;
        //badanie czasu obliczania
        vector < long long > executeTimes;
        long long totalTime = 0;
        for (int i = 1; i <= testNumber; i++)
        {
            cout << '\r' << i << '/' << testNumber;
            //pomiar czasu obliczania iloczynu macierzy
            chrono::steady_clock::time_point start = chrono::steady_clock::now();
            ////////////////
            productMatrixMatrix(Matrix1, Matrix2, resultMatrix);
            ////////////////
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            //END pomiar czasu obliczania iloczynu macierzy

            executeTimes.push_back(chrono::duration_cast<chrono::nanoseconds>(end - start).count());
            totalTime += chrono::duration_cast<chrono::nanoseconds>(end - start).count();
            //cout << chrono::duration_cast<chrono::nanoseconds>(end - start).count() << '\t';
        }
        cout << "\r                     \r";
        //wyświetlanie posortowanych czasów obliczania
        //cout << "Średni czas obliczeń: " << totalTime / testNumber << "ns." << endl << endl;
        sort(executeTimes.begin(), executeTimes.end());
        //cout << "Posortowane czasy wykonań [ns]: \n";
        cout << endl << omp_get_max_threads() << ' ';
        for (const auto& i : executeTimes)
            cout << i << ' ';
        cout << endl;
        //END wyświetlanie posortowanych czasów obliczania

        //END badanie czasu obliczania
    }

    getchar();
}