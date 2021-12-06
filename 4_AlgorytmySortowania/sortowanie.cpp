// sortowanie.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <time.h>
#include <Windows.h>
#include <fstream>
#include <vector>
#include <chrono>
#include <omp.h>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <random>

using namespace std;
static fstream fs;
bool TEST = false;
HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);

float fRand()
{
    float fMin = -1000;
    float fMax = 1000;
    float f = (float)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int intRand()
{
    int fMin = -3;
    int fMax = 3;
    int f = rand() % (fMax - fMin) + fMin;
    return f;
}

void generateRandomArray(int n, float* Array)
{
    for (int i = 0; i < n; i++)
    {
        Array[i] = fRand();
    }
}

void generatePermutationArray(int n, float* Array)
{
    for (int i = 0; i < n; i++)
    {
        Array[i] = i;
    }
    shuffle(Array, Array+n, default_random_engine(time(NULL)));
}

float* copyArray(int n, float* Array)
{
    float* newArray = new float[n];
    for (int i = 0; i < n; i++)
    {
        newArray[i] = Array[i];
    }
    return newArray;
}

void printArray(int n, float* Array)
{
    for (int i = 0; i < n; i++)
    {
        cout << Array[i] << endl;;
    }
}

void printArrays(int n, float* Array, float* sortedArray)
{
    for (int i = 0; i < n; i++)
    {
        cout << setw(20) << left << Array[i] << setw(20) << left << sortedArray[i] << endl;;
    }
}
//FUNKCJE SORTUJ¥CE/////////////////////////////
void bubbleSort(int n, float* Array)
{
    for (int i = 0; i < n -1; i++)
    {
        for (int j = 0; j < n - i - 1; j++)
        {
            if (Array[j] > Array[j + 1])
            {
                float temp = Array[j];
                Array[j] = Array[j+1];
                Array[j + 1] = temp;
            }
        }
    }
}

void bubbleSort(float* Array, int left, int right)
{
    bubbleSort(right - left, &Array[left]);
}

void pipeBubbleSort(int n, float* Array)
{
    //tablica postêpu
    int* doneTable;
    doneTable = new int[n];
    for (int i = 0; i < n; i++) doneTable[i] = 0;

    #pragma omp parallel for schedule(dynamic,1) firstprivate(n) shared(doneTable, Array)
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < n - i - 1; j++)
        {
            while (doneTable[j] != i || doneTable[j] != doneTable[j + 1]) { Sleep(0); }
            if (Array[j] > Array[j + 1])
            {
                float temp = Array[j];
                Array[j] = Array[j + 1];
                Array[j + 1] = temp;
            }
            doneTable[j] = i + 1;
        }
    }
    delete[] doneTable;
}

//prawie identyczna jak pipeBubbleSort
void oddEvenSort(int n, float* Array)
{
    //tablica postêpu
    int* doneTable;
    doneTable = new int[n];
    for (int i = 0; i < n; i++) doneTable[i] = 0;

    //kolejne fazy
    //#pragma omp parallel for schedule(dynamic, 1) firstprivate(n) shared(doneTable, Array)
    for (int i = 0; i < n; i++)
    {
        if (i % 2 == 0)
        //jeœli faza parzysta
        {
            #pragma omp parallel for firstprivate(n) shared(doneTable, Array)
            for (int j = 0; j < n;j += 2)
            {
                while (doneTable[j] != i || doneTable[j] != doneTable[j + 1]) { Sleep(0); }
                if (Array[j] > Array[j + 1])
                {
                    float temp = Array[j];
                    Array[j] = Array[j + 1];
                    Array[j + 1] = temp;
                }
                doneTable[j] = i + 1;
                doneTable[j + 1] = i + 1;
            }
            //jeœli nieparzysta liczba elementów ustaw ostatni jako obs³u¿ony
            if (n % 2 == 1) {
                while (doneTable[n - 1] != i) { Sleep(0); }
                doneTable[n - 1] = i + 1;
            }
        }
        else
        //jeœli faza nieparzysta
        {
            //ustaw pierwszy element jako obs³u¿ony
            while (doneTable[0] != i) { Sleep(0); }
            doneTable[0] = i + 1;
            //#pragma omp parallel for schedule(static) firstprivate(n) shared(doneTable, Array)
            for (int j = 1; j < n-1;j += 2)
            {
                while (doneTable[j] != i || doneTable[j] != doneTable[j + 1]) { Sleep(0); }
                if (Array[j] > Array[j + 1])
                {
                    float temp = Array[j];
                    Array[j] = Array[j + 1];
                    Array[j + 1] = temp;
                }
                doneTable[j] = i + 1;
                doneTable[j + 1] = i + 1;
            }
            //jeœli parzysta liczba elementów ustaw ostatni jako obs³u¿ony
            if (n % 2 == 0) {
                while (doneTable[n - 1] != i) { Sleep(0); }
                doneTable[n - 1] = i + 1;
            }
        }
    }
    delete[] doneTable;
}

void quickSort(float* Array, int left, int right)
{
    float v = Array[(left + right) / 2];
    int i, j;
    float x;
    i = left;
    j = right;
    do
    {
        while (Array[i] < v) i++;
        while (Array[j] > v) j--;
        if (i <= j)
        {
            x = Array[i];
            Array[i] = Array[j];
            Array[j] = x;
            i++; j--;
        }
    } while (i <= j);
    if (j > left) quickSort(Array, left, j);
    if (i < right) quickSort(Array, i, right);
}

void quickSort(int n, float* Array)
{
    quickSort(Array, 0, n-1);
}

void quickSortNaive(float* Array, int left, int right)
{
    float v = Array[(left + right) / 2];
    int i, j;
    float x;
    i = left;
    j = right;
    do
    {
        while (Array[i] < v) i++;
        while (Array[j] > v) j--;
        if (i <= j)
        {
            x = Array[i];
            Array[i] = Array[j];
            Array[j] = x;
            i++; j--;
        }
    } while (i <= j);

    #pragma omp parallel shared(Array) firstprivate (left, j)
    {
        #pragma omp sections nowait
        {
            #pragma omp section
            {
                if (j > left) quickSortNaive(Array, left, j);
            }
            #pragma omp section
            {
                if (i < right) quickSortNaive(Array, i, right);
            }
        }
    }
    //#pragma omp parallel shared(Array) firstprivate (left, j)
    //{
    //    #pragma omp single nowait
    //    {
    //        quickSortNaive(Array, left, j);
    //    }
    //    #pragma omp single nowait
    //    {
    //        quickSortNaive(Array, i, right);
    //    }
    //}
        
    
}

void quickSortNaive(int n, float* Array)
{
    quickSortNaive(Array, 0, n - 1);
}

void rankSort(int n, float* Array)
{
    float* localArray = copyArray(n, Array);
    for (int i = 0;i < n;i++)
    {
        int x = 0;
        for (int j = 0;j < n;j++)
        {
            if (localArray[i] > localArray[j]) x++;
        }
        Array[x] = localArray[i];
    }
    delete[] localArray;
}

void rankSort(float* Array, int left, int right)
{
    rankSort(right - left, &Array[left]);
}

void rankSortMulti(int n, float* Array)
{
    float* localArray = copyArray(n, Array);
    #pragma omp parallel for firstprivate(n) shared(localArray, Array)
    for (int i = 0;i < n;i++)
    {
        int x = 0;
        for (int j = 0;j < n;j++)
        {
            if (localArray[i] > localArray[j]) x++;
        }
        Array[x] = localArray[i];
    }
    delete[] localArray;
}

void countingSort(int n, float* Array)
{
    float* localArray = copyArray(n, Array);

    int* histogram = new int[n];
    for (int i = 0; i < n; i++) histogram[i] = 0;
    for (int i = 0; i < n; i++) histogram[(int)localArray[i]]++;
    //suma prefiksowa
    for (int i = 1; i < n; i++) histogram[i] = histogram[i] + histogram[i - 1];
    //umieszczanie elementów
    for (int i = n - 1;i >= 0;i--)
    {
        Array[histogram[(int)localArray[i]] - 1] = localArray[i];
        histogram[(int)localArray[i]]--;
    }
    delete[] localArray;
    delete[] histogram;
}

void countingSort(float* Array, int left, int right)
{
    countingSort(right - left, &Array[left]);
}

void countingSortMulti(int n, float* Array)
{
    float* localArray = copyArray(n, Array);

    int* histogram = new int[n];
    #pragma omp parallel for firstprivate(n) shared(histogram)
    for (int i = 0; i < n; i++) histogram[i] = 0;
    #pragma omp parallel for firstprivate(n) shared(histogram, localArray)
    for (int i = 0; i < n; i++) histogram[(int)localArray[i]]++;
    //suma prefiksowa
    for (int i = 1; i < n; i++) histogram[i] = histogram[i] + histogram[i - 1];
    //umieszczanie elementów
    #pragma omp parallel for firstprivate(n) shared(histogram, localArray, Array)
    for (int i = n - 1;i >= 0;i--)
    {
        Array[histogram[(int)localArray[i]] - 1] = localArray[i];
        histogram[(int)localArray[i]]--;
    }
    delete[] localArray;
    delete[] histogram;
}

void divSingleRank(int n, float* Array)
{
    //przydzielenie liczby zadañ procesom
    int p = omp_get_max_threads();
    int* nPerThread = new int[p];
    for (int i = 0; i < p;i++) nPerThread[i] = n / p;
    if (n % p != 0) for (int i = 0; i < n % p; i++) nPerThread[i]++;
    //dolne granice procesów
    int* leftThread = new int[p];
    leftThread[0] = 0;
    for (int i = 1; i < p; i++) leftThread[i] = leftThread[i - 1] + nPerThread[i - 1];
    //sortowanie podci¹gów
    for (int i = 0; i < p; i++)
    {
        //posortuj
        rankSort(Array, leftThread[i], leftThread[i]+ nPerThread[i]);
    }
    //agregacja podci¹gów
    //scalanie drzewiaste, dopóki nie uzyska siê jednego ci¹gu
    int aggregations = p;
    for (int i = 2; aggregations != 1; i *= 2)
    {
        float* localArray = copyArray(n, Array);
        for (int j = 0; j < p - i / 2; j += i)
        {
            //scal
            int ArrayLeftIndex = leftThread[j]; int ArrayLeftMax = leftThread[j] + nPerThread[j];
            int ArrayRightIndex = leftThread[j + i / 2]; int ArrayRightMax = leftThread[j + i / 2] + nPerThread[j + i / 2];
            int MainArrayIndex = ArrayLeftIndex;
            while ((ArrayLeftIndex < ArrayLeftMax || ArrayRightIndex < ArrayRightMax) && MainArrayIndex < n)
            {
                if (ArrayLeftIndex >= ArrayLeftMax)
                {
                    Array[MainArrayIndex] = localArray[ArrayRightIndex];
                    MainArrayIndex++;
                    ArrayRightIndex++;
                }
                else if (ArrayRightIndex >= ArrayRightMax)
                {
                    Array[MainArrayIndex] = localArray[ArrayLeftIndex];
                    MainArrayIndex++;
                    ArrayLeftIndex++;
                }
                else
                {
                    if (localArray[ArrayLeftIndex] <= localArray[ArrayRightIndex])
                    {
                        Array[MainArrayIndex] = localArray[ArrayLeftIndex];
                        MainArrayIndex++;
                        ArrayLeftIndex++;
                    }
                    else
                    {
                        Array[MainArrayIndex] = localArray[ArrayRightIndex];
                        MainArrayIndex++;
                        ArrayRightIndex++;
                    }
                }
            }
            nPerThread[j] += nPerThread[j + i / 2];
            nPerThread[j + i / 2] = 0;
            aggregations--;
        }
        delete[] localArray;
    }

    delete[] nPerThread;
    delete[] leftThread;
}

void divSingleBubble(int n, float* Array)
{
    //przydzielenie liczby zadañ procesom
    int p = omp_get_max_threads();
    int* nPerThread = new int[p];
    for (int i = 0; i < p;i++) nPerThread[i] = n / p;
    if (n % p != 0) for (int i = 0; i < n % p; i++) nPerThread[i]++;
    //dolne granice procesów
    int* leftThread = new int[p];
    leftThread[0] = 0;
    for (int i = 1; i < p; i++) leftThread[i] = leftThread[i - 1] + nPerThread[i - 1];
    //sortowanie podci¹gów
    for (int i = 0; i < p; i++)
    {
        //posortuj
        bubbleSort(Array, leftThread[i], leftThread[i] + nPerThread[i]);
    }
    //agregacja podci¹gów
    //scalanie drzewiaste, dopóki nie uzyska siê jednego ci¹gu
    int aggregations = p;
    for (int i = 2; aggregations != 1; i *= 2)
    {
        float* localArray = copyArray(n, Array);
        for (int j = 0; j < p - i / 2; j += i)
        {
            //scal
            int ArrayLeftIndex = leftThread[j]; int ArrayLeftMax = leftThread[j] + nPerThread[j];
            int ArrayRightIndex = leftThread[j + i / 2]; int ArrayRightMax = leftThread[j + i / 2] + nPerThread[j + i / 2];
            int MainArrayIndex = ArrayLeftIndex;
            while ((ArrayLeftIndex < ArrayLeftMax || ArrayRightIndex < ArrayRightMax) && MainArrayIndex < n)
            {
                if (ArrayLeftIndex >= ArrayLeftMax)
                {
                    Array[MainArrayIndex] = localArray[ArrayRightIndex];
                    MainArrayIndex++;
                    ArrayRightIndex++;
                }
                else if (ArrayRightIndex >= ArrayRightMax)
                {
                    Array[MainArrayIndex] = localArray[ArrayLeftIndex];
                    MainArrayIndex++;
                    ArrayLeftIndex++;
                }
                else
                {
                    if (localArray[ArrayLeftIndex] <= localArray[ArrayRightIndex])
                    {
                        Array[MainArrayIndex] = localArray[ArrayLeftIndex];
                        MainArrayIndex++;
                        ArrayLeftIndex++;
                    }
                    else
                    {
                        Array[MainArrayIndex] = localArray[ArrayRightIndex];
                        MainArrayIndex++;
                        ArrayRightIndex++;
                    }
                }
            }
            nPerThread[j] += nPerThread[j + i / 2];
            nPerThread[j + i / 2] = 0;
            aggregations--;
        }
        delete[] localArray;
    }

    delete[] nPerThread;
    delete[] leftThread;
}

void divSingleQuick(int n, float* Array)
{
    //przydzielenie liczby zadañ procesom
    int p = omp_get_max_threads();
    int* nPerThread = new int[p];
    for (int i = 0; i < p;i++) nPerThread[i] = n / p;
    if (n % p != 0) for (int i = 0; i < n % p; i++) nPerThread[i]++;
    //dolne granice procesów
    int* leftThread = new int[p];
    leftThread[0] = 0;
    for (int i = 1; i < p; i++) leftThread[i] = leftThread[i - 1] + nPerThread[i - 1];
    //sortowanie podci¹gów
    for (int i = 0; i < p; i++)
    {
        //posortuj
        quickSort(Array, leftThread[i], leftThread[i] + nPerThread[i] -1);
    }
    //agregacja podci¹gów
    //scalanie drzewiaste, dopóki nie uzyska siê jednego ci¹gu
    int aggregations = p;
    for (int i = 2; aggregations != 1; i *= 2)
    {
        float* localArray = copyArray(n, Array);
        for (int j = 0; j < p - i / 2; j += i)
        {
            //scal
            int ArrayLeftIndex = leftThread[j]; int ArrayLeftMax = leftThread[j] + nPerThread[j];
            int ArrayRightIndex = leftThread[j + i / 2]; int ArrayRightMax = leftThread[j + i / 2] + nPerThread[j + i / 2];
            int MainArrayIndex = ArrayLeftIndex;
            while ((ArrayLeftIndex < ArrayLeftMax || ArrayRightIndex < ArrayRightMax) && MainArrayIndex < n)
            {
                if (ArrayLeftIndex >= ArrayLeftMax)
                {
                    Array[MainArrayIndex] = localArray[ArrayRightIndex];
                    MainArrayIndex++;
                    ArrayRightIndex++;
                }
                else if (ArrayRightIndex >= ArrayRightMax)
                {
                    Array[MainArrayIndex] = localArray[ArrayLeftIndex];
                    MainArrayIndex++;
                    ArrayLeftIndex++;
                }
                else
                {
                    if (localArray[ArrayLeftIndex] <= localArray[ArrayRightIndex])
                    {
                        Array[MainArrayIndex] = localArray[ArrayLeftIndex];
                        MainArrayIndex++;
                        ArrayLeftIndex++;
                    }
                    else
                    {
                        Array[MainArrayIndex] = localArray[ArrayRightIndex];
                        MainArrayIndex++;
                        ArrayRightIndex++;
                    }
                }
            }
            nPerThread[j] += nPerThread[j + i / 2];
            nPerThread[j + i / 2] = 0;
            aggregations--;
        }
        delete[] localArray;
    }

    delete[] nPerThread;
    delete[] leftThread;
}

void divMultiRank(int n, float* Array)
{
    //przydzielenie liczby zadañ procesom
    int p = omp_get_max_threads();
    int* nPerThread = new int[p];
    for (int i = 0; i < p;i++) nPerThread[i] = n / p;
    if (n % p != 0) for (int i = 0; i < n % p; i++) nPerThread[i]++;
    //dolne granice procesów
    int* leftThread = new int[p];
    leftThread[0] = 0;
    for (int i = 1; i < p; i++) leftThread[i] = leftThread[i - 1] + nPerThread[i - 1];
    //sortowanie podci¹gów
#pragma omp parallel for schedule(static, 1) firstprivate(p) shared(Array, leftThread, nPerThread)
    for (int i = 0; i < p; i++)
    {
        //posortuj
        rankSort(Array, leftThread[i], leftThread[i] + nPerThread[i]);
    }
    //agregacja podci¹gów
    //scalanie drzewiaste, dopóki nie uzyska siê jednego ci¹gu
    int aggregations = p;
    for (int i = 2; aggregations != 1; i *= 2)
    {
        float* localArray = copyArray(n, Array);
#pragma omp parallel for firstprivate(i, n) shared(localArray, Array, aggregations, nPerThread, leftThread) schedule(static, 1)
        for (int j = 0; j < p - i / 2; j += i)
        {
            //scal
            int ArrayLeftIndex = leftThread[j]; int ArrayLeftMax = leftThread[j] + nPerThread[j];
            int ArrayRightIndex = leftThread[j + i / 2]; int ArrayRightMax = leftThread[j + i / 2] + nPerThread[j + i / 2];
            int MainArrayIndex = ArrayLeftIndex;
            while ((ArrayLeftIndex < ArrayLeftMax || ArrayRightIndex < ArrayRightMax) && MainArrayIndex < n)
            {
                if (ArrayLeftIndex >= ArrayLeftMax)
                {
                    Array[MainArrayIndex] = localArray[ArrayRightIndex];
                    MainArrayIndex++;
                    ArrayRightIndex++;
                }
                else if (ArrayRightIndex >= ArrayRightMax)
                {
                    Array[MainArrayIndex] = localArray[ArrayLeftIndex];
                    MainArrayIndex++;
                    ArrayLeftIndex++;
                }
                else
                {
                    if (localArray[ArrayLeftIndex] <= localArray[ArrayRightIndex])
                    {
                        Array[MainArrayIndex] = localArray[ArrayLeftIndex];
                        MainArrayIndex++;
                        ArrayLeftIndex++;
                    }
                    else
                    {
                        Array[MainArrayIndex] = localArray[ArrayRightIndex];
                        MainArrayIndex++;
                        ArrayRightIndex++;
                    }
                }
            }
            nPerThread[j] += nPerThread[j + i / 2];
            nPerThread[j + i / 2] = 0;
#pragma omp atomic
            aggregations--;
        }
        delete[] localArray;
    }

    delete[] nPerThread;
    delete[] leftThread;
}

void divMultiBubble(int n, float* Array)
{
    //przydzielenie liczby zadañ procesom
    int p = omp_get_max_threads();
    int* nPerThread = new int[p];
    for (int i = 0; i < p;i++) nPerThread[i] = n / p;
    if (n % p != 0) for (int i = 0; i < n % p; i++) nPerThread[i]++;
    //dolne granice procesów
    int* leftThread = new int[p];
    leftThread[0] = 0;
    for (int i = 1; i < p; i++) leftThread[i] = leftThread[i - 1] + nPerThread[i - 1];
    //sortowanie podci¹gów
#pragma omp parallel for schedule(static, 1) firstprivate(p) shared(Array, leftThread, nPerThread)
    for (int i = 0; i < p; i++)
    {
        //posortuj
        bubbleSort(Array, leftThread[i], leftThread[i] + nPerThread[i]);
    }
    //agregacja podci¹gów
    //scalanie drzewiaste, dopóki nie uzyska siê jednego ci¹gu
    int aggregations = p;
    for (int i = 2; aggregations != 1; i *= 2)
    {
        float* localArray = copyArray(n, Array);
#pragma omp parallel for firstprivate(i, n) shared(localArray, Array, aggregations, nPerThread, leftThread) schedule(static, 1)
        for (int j = 0; j < p - i / 2; j += i)
        {
            //scal
            int ArrayLeftIndex = leftThread[j]; int ArrayLeftMax = leftThread[j] + nPerThread[j];
            int ArrayRightIndex = leftThread[j + i / 2]; int ArrayRightMax = leftThread[j + i / 2] + nPerThread[j + i / 2];
            int MainArrayIndex = ArrayLeftIndex;
            while ((ArrayLeftIndex < ArrayLeftMax || ArrayRightIndex < ArrayRightMax) && MainArrayIndex < n)
            {
                if (ArrayLeftIndex >= ArrayLeftMax)
                {
                    Array[MainArrayIndex] = localArray[ArrayRightIndex];
                    MainArrayIndex++;
                    ArrayRightIndex++;
                }
                else if (ArrayRightIndex >= ArrayRightMax)
                {
                    Array[MainArrayIndex] = localArray[ArrayLeftIndex];
                    MainArrayIndex++;
                    ArrayLeftIndex++;
                }
                else
                {
                    if (localArray[ArrayLeftIndex] <= localArray[ArrayRightIndex])
                    {
                        Array[MainArrayIndex] = localArray[ArrayLeftIndex];
                        MainArrayIndex++;
                        ArrayLeftIndex++;
                    }
                    else
                    {
                        Array[MainArrayIndex] = localArray[ArrayRightIndex];
                        MainArrayIndex++;
                        ArrayRightIndex++;
                    }
                }
            }
            nPerThread[j] += nPerThread[j + i / 2];
            nPerThread[j + i / 2] = 0;
#pragma omp atomic
            aggregations--;
        }
        delete[] localArray;
    }

    delete[] nPerThread;
    delete[] leftThread;
}

void divMultiQuick(int n, float* Array)
{
    //przydzielenie liczby zadañ procesom
    int p = omp_get_max_threads();
    int* nPerThread = new int[p];
    for (int i = 0; i < p;i++) nPerThread[i] = n / p;
    if (n % p != 0) for (int i = 0; i < n % p; i++) nPerThread[i]++;
    //dolne granice procesów
    int* leftThread = new int[p];
    leftThread[0] = 0;
    for (int i = 1; i < p; i++) leftThread[i] = leftThread[i - 1] + nPerThread[i - 1];
    //sortowanie podci¹gów
#pragma omp parallel for schedule(static, 1) firstprivate(p) shared(Array, leftThread, nPerThread)
    for (int i = 0; i < p; i++)
    {
        //posortuj
        quickSort(Array, leftThread[i], leftThread[i] + nPerThread[i] -1);
    }
    //agregacja podci¹gów
    //scalanie drzewiaste, dopóki nie uzyska siê jednego ci¹gu
    int aggregations = p;
    for (int i = 2; aggregations != 1; i *= 2)
    {
        float* localArray = copyArray(n, Array);
#pragma omp parallel for firstprivate(i, n) shared(localArray, Array, aggregations, nPerThread, leftThread) schedule(static, 1)
        for (int j = 0; j < p - i / 2; j += i)
        {
            //scal
            int ArrayLeftIndex = leftThread[j]; int ArrayLeftMax = leftThread[j] + nPerThread[j];
            int ArrayRightIndex = leftThread[j + i / 2]; int ArrayRightMax = leftThread[j + i / 2] + nPerThread[j + i / 2];
            int MainArrayIndex = ArrayLeftIndex;
            while ((ArrayLeftIndex < ArrayLeftMax || ArrayRightIndex < ArrayRightMax) && MainArrayIndex < n)
            {
                if (ArrayLeftIndex >= ArrayLeftMax)
                {
                    Array[MainArrayIndex] = localArray[ArrayRightIndex];
                    MainArrayIndex++;
                    ArrayRightIndex++;
                }
                else if (ArrayRightIndex >= ArrayRightMax)
                {
                    Array[MainArrayIndex] = localArray[ArrayLeftIndex];
                    MainArrayIndex++;
                    ArrayLeftIndex++;
                }
                else
                {
                    if (localArray[ArrayLeftIndex] <= localArray[ArrayRightIndex])
                    {
                        Array[MainArrayIndex] = localArray[ArrayLeftIndex];
                        MainArrayIndex++;
                        ArrayLeftIndex++;
                    }
                    else
                    {
                        Array[MainArrayIndex] = localArray[ArrayRightIndex];
                        MainArrayIndex++;
                        ArrayRightIndex++;
                    }
                }
            }
            nPerThread[j] += nPerThread[j + i / 2];
            nPerThread[j + i / 2] = 0;
#pragma omp atomic
            aggregations--;
        }
        delete[] localArray;
    }

    delete[] nPerThread;
    delete[] leftThread;
}
//END FUNKCJE SORTUJ¥CE/////////////////////////////
typedef void (*sortFunctionT)(int, float*); //wskaŸnik na funkcjê sortuj¹c¹

void multipleTestExecution(int n, float* Array, sortFunctionT sortFunction, int maxTests)
{
    //badanie czasu obliczania
    vector < long long > executeTimes;
    for (int i = 1; i <= maxTests; i++)
    {
        if (!TEST) cout << '\r' << i << '/' << maxTests;
        //kopia lokalna wygenerowanej tablicy
        float* localArray = copyArray(n, Array);
        //pomiar czasu obliczania
        chrono::steady_clock::time_point start = chrono::steady_clock::now();
        sortFunction(n, localArray);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        //END pomiar czasu obliczania
        executeTimes.push_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
        //poka¿ wyniki
        if (TEST) printArrays(n, Array, localArray);
        delete[] localArray;
    }
    if (!TEST) cout << "\r                     \r";
    //wyœwietlanie posortowanych czasów obliczania
    sort(executeTimes.begin(), executeTimes.end());
    if (!TEST) cout << setw(10) << left << omp_get_max_threads(); fs << setw(10) << left << omp_get_max_threads();
    if (!TEST) for (const auto& i : executeTimes)
    {
        cout << setw(10) << left << i; fs << setw(10) << left << i;
    }
    cout << endl; fs << endl;
    //END wyœwietlanie posortowanych czasów obliczania

    //END badanie czasu obliczania
}

int main()
{
    chrono::steady_clock::time_point start = chrono::steady_clock::now();
    setlocale(LC_CTYPE, "Polish");
    omp_set_dynamic(0);
    srand(time(NULL));
    ostringstream filenameStream;
    filenameStream << "wyniki\\out" << time(NULL) << ".txt";
    fs.open(filenameStream.str(), fstream::in | fstream::out | fstream::trunc);
    //ZMIENNE OGÓLNE
    //wektor w¹tków
    vector < int > threads;
    //liczba testów
    int maxTests = 10;
    if (TEST) maxTests = 1;
    //wektory funkcji
    vector < sortFunctionT > sortFunctions;
    vector < string > sortFunctionsNames;
    //wektor wielkoœci
    vector < int > n;
    //END ZMIENNE OGÓLNE
    //SORTOWANIE PERMUTACJI, WSZYSTKIE ALGORYTMY
    threads.clear();
    for (int i = 1; i <= 12; i++) threads.push_back(i);     //i <= 50 dla divFunkcji
    //wektory funkcji
    sortFunctions.clear();
    sortFunctionsNames.clear();
    { sortFunctions.push_back(rankSortMulti); sortFunctionsNames.push_back("rankSortMulti"); }
    { sortFunctions.push_back(countingSortMulti); sortFunctionsNames.push_back("countingSortMulti"); }
    { sortFunctions.push_back(bubbleSort); sortFunctionsNames.push_back("bubbleSort"); }
    { sortFunctions.push_back(pipeBubbleSort); sortFunctionsNames.push_back("pipeBubbleSort"); }
    { sortFunctions.push_back(oddEvenSort); sortFunctionsNames.push_back("oddEvenSort"); }
    { sortFunctions.push_back(quickSort); sortFunctionsNames.push_back("quickSort"); }
    { sortFunctions.push_back(divSingleRank); sortFunctionsNames.push_back("divSingleRank"); }
    { sortFunctions.push_back(divSingleBubble); sortFunctionsNames.push_back("divSingleBubble"); }
    { sortFunctions.push_back(divSingleQuick); sortFunctionsNames.push_back("divSingleQuick"); }
    { sortFunctions.push_back(divMultiRank); sortFunctionsNames.push_back("divMultiRank"); }
    { sortFunctions.push_back(divMultiBubble); sortFunctionsNames.push_back("divMultiBubble"); }
    { sortFunctions.push_back(divMultiQuick); sortFunctionsNames.push_back("divMultiQuick"); }
    //wektor wielkoœci
    n.clear();
    n = { 1250, 2500, 5000, 10000 };    //ma³e wartoœci, dla najs³abszych algorytmów
    if (TEST) n = { 1 };
    //dla ka¿dej wielkoœci danych
    for (int i = 0; i < n.size(); i++)
    {
        cout << n[i] << endl; fs << n[i] << endl;
        float* Array = new float[n[i]];
        generatePermutationArray(n[i], Array);
        //dla ka¿dego algorytmu
        for (int j = 0; j < sortFunctions.size(); j++)
        {
            cout << sortFunctionsNames[j] << endl; fs << sortFunctionsNames[j] << endl;
            //dla ka¿dej liczby w¹tków
            for (int k = 0; k < threads.size(); k++)
            {
                omp_set_num_threads(threads[k]);
                multipleTestExecution(n[i], Array, sortFunctions[j], maxTests);
                //funkcje jednow¹tkowe wykonaj tylko na jednym w¹tku
                if (sortFunctionsNames[j] == "bubbleSort") break;
                if (sortFunctionsNames[j] == "quickSort") break;
            }
        }
        delete[] Array;
    }
    //END SORTOWANIE PERMUTACJI, WSZYSTKIE ALGORYTMY
    ////////////////////////////////////////////////////////////////////////////////////////

    fs.close();
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    int min = chrono::duration_cast<chrono::minutes>(end - start).count();
    int sec = chrono::duration_cast<chrono::seconds>(end - start).count()%60;
    int milisec = chrono::duration_cast<chrono::milliseconds>(end - start).count()%1000;
    cout << "Zakoñczono: " << min << "m " << sec << "s " << milisec << "ms " << endl;
    getchar();
}


