
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <Windows.h>

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
    int f = rand()%(fMax-fMin)+fMin;
    return f;
}

void generateRandomMatrix(int n, float ** Matrix)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n + 1; j++)
        {
            Matrix[i][j] = fRand();
            if (TEST) Matrix[i][j] = intRand();
        }
    }
}

void printMatrix(int n, float** Matrix)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n + 1; j++)
        {
            cout << setw(10) << left << Matrix[i][j];
        }
        cout << endl;
    }
}

void printMatrix(int n, float** Matrix, int row, int col)
{
    SetConsoleCursorPosition(hConsole, { 0, 0 });
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n + 1; j++)
        {
            if (i==row && j==col) SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 12);
            cout << setw(10) << left << Matrix[i][j];
            if (i == row && j == col) SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7);
        }
        cout << endl;
    }
    getchar();
}

void printMatrix(int n, int** Matrix, int row, int col)
{
    SetConsoleCursorPosition(hConsole, { 0, 0 });
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n + 1; j++)
        {
            if (i == row && j == col) SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 12);
            cout << setw(10) << left << Matrix[i][j];
            if (i == row && j == col) SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7);
        }
        cout << endl;
    }
    getchar();
}

void GaussElimination(int n, float ** AB)
{
    for (int k = 0; k < n; k++)
    {
        //jeśli pierwszy rozpatrywany element na przekątnej jest równy 0, to zamieniaj z kolejnymi wierszami (o ile istnieją) do skutku lub do końca wierszy
        //int nextRow = k + 1;
        //while (AB[k][k] == 0 && nextRow < n)
        //{
        //    for (int j = 0; j < n; j++)
        //    {
        //        swap(AB[k][j], AB[nextRow][j]);
        //    }
        //    nextRow += 1;
        //}
        if (AB[k][k] == 0)
        {
            break;
        }
        //OBLICZENIA
        for (int j = k + 1; j < n + 1; j++)
        {
            AB[k][j] = AB[k][j] / AB[k][k];
        }
        AB[k][k] = 1;
        //ELIMINACJA
        for (int i = k + 1; i < n; i++)
        {
            for (int j = k + 1; j < n + 1; j++)
                AB[i][j] = AB[i][j] - AB[i][k] * AB[k][j];
            AB[i][k] = 0;
        }
    }
}

void rowBlockGaussElimination(int n, float** AB)
{
    for (int k = 0; k < n; k++)
    {
        #pragma omp parallel for firstprivate(k, n) shared(AB)
        for (int j = k + 1; j < n + 1; j++)
        {
            AB[k][j] = AB[k][j] / AB[k][k];
        }
        AB[k][k] = 1;

        #pragma omp parallel for firstprivate(k, n) shared(AB)
        for (int i = k + 1; i < n; i++)
        {
            for (int j = k + 1; j < n + 1; j++)
                AB[i][j] = AB[i][j] - AB[i][k] * AB[k][j];
            AB[i][k] = 0;
        }
    }
}

void rowCyclicGaussElimination(int n, float** AB)
{
    for (int k = 0; k < n; k++)
    {
        for (int i = k + 1; i < n + 1; i++)
        {
            AB[k][i] = AB[k][i] / AB[k][k];
        }
        AB[k][k] = 1;

        #pragma omp parallel for schedule(dynamic, 1) firstprivate(k, n) shared(AB)
        for (int i = k + 1; i < n; i++)
        {
            for (int j = k + 1; j < n + 1; j++)
                AB[i][j] = AB[i][j] - AB[i][k] * AB[k][j];
            AB[i][k] = 0;
        }
    }
}

void pipeGaussElimination(int n, float** AB)
{
    //tablica przechowująca postęp działania
    int** doneTable; 
    doneTable = new int* [n];
    for (int j = 0; j < n; j++) doneTable[j] = new int[n + 1];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n + 1; j++) doneTable[i][j] = 0;
    }

    //tablica z kolejnością wykonywania
    int** executeOrderTable;
    executeOrderTable = new int* [(n * (1 + n)) / 2];
    for (int j = 0; j < (n * (1 + n)) / 2; j++) executeOrderTable[j] = new int[2];

    int tmp = 0;
    for (int k = 0; k < n; k++)
    {
        int var = k;
        while (var < n)
        {
            executeOrderTable[tmp][0] = var;    //który wiersz przypada wątkowi
            executeOrderTable[tmp][1] = k;      //przez który wiersz będzie dzielenie
            //cout << var << ": " << k << endl;;
            var++;
            tmp++;
        }
    }
    //tablice zostały przygotowane
    int k = 0;
    #pragma omp parallel for schedule(dynamic, 1) firstprivate(n) shared(AB, doneTable, executeOrderTable)
    for (k = 0; k < (n*(1+n))/2; k++)
    {
        int row = executeOrderTable[k][0]; //wiersz przydzielony wątkowi
        int iteration = executeOrderTable[k][1]; //wiersz przez który wątek będzie dzielił - nr iteracji



        //cout << row << iteration << endl;
        if (row == iteration) //jeśli wiersz jest pierwszy w kolejnej iteracji
        {
            //obliczanie
            //jesli element nie został jeszcze obliczony przez poprzednie iteracje to czekaj
            //bool flag = true;
            //do
            //{
            //    //TU TKWI PROBLEM
            //    if (doneTable[row][row] == iteration)
            //    {
            //        flag = false;
            //    }
            //    else
            //    {
            //        //!!!!!!!!!!!!!!!!!cout << 'z' << doneTable[row][row] << endl; - z tą komendą działa
            //        //cout << 'z'; //- z tą też
            //        //iteration = iteration; //- z tą nie
            //        //if (doneTable[row][row] == iteration)         //Z TYM
            //        //{                                             //TEŻ
            //        //    flag = false;                             //NIE
            //        //}                                             //DZIAŁA
            //        //int i;    // też nie
            //        //int i = doneTable[row][row] / iteration;  //też nie
            //        //cout; //też nie
            //        cout << ""; // działa
            //        flag = true; 
            //    }
            //} while (flag);
            while (doneTable[row][row] != iteration) { Sleep(0); }
            float ratio = AB[row][row];
            for (int i = row; i < n + 1; i++)
            {
                //if (i == row) cout << 'y' << doneTable[row][i] << endl;
                //jesli element nie został jeszcze obliczony przez poprzednie iteracje to czekaj
                //bool flag;
                //do
                //{
                //    flag = true;
                //    if (doneTable[row][i] == iteration)
                //    {
                //        flag = false;
                //        //cout << 'y';
                //    }
                //    else
                //    {
                //        flag = true;
                //        //cout << 'x';
                //        //cout << iteration << endl;
                //    }
                //} while (flag);
                while (doneTable[row][i] != iteration) { Sleep(0); }
                AB[row][i] = AB[row][i] / ratio;
                doneTable[row][i]=iteration+1;    //obliczony element otrzymuje status o jeden większy niż poprzednio
                //if (i==row) cout << 'x' << doneTable[row][i] << endl;
                //#pragma omp critical
                //{
                //    printMatrix(n, doneTable, row, i);
                //};
            }
        }
        else
        {
            //eliminacja
            //dopóki elementy nie zostały obliczone przez poprzednie iteracje to czekaj
            //bool flag = true;
            //do
            //{
            //    if (doneTable[row][iteration] == iteration && doneTable[iteration][iteration] == iteration+1)
            //    {
            //        flag = false;
            //        //cout << 'y';
            //    }
            //    else
            //    {
            //        //TU UTYKAJĄ WĄTKI
            //        flag = true;
            //        //cout << 'x';
            //        //cout << iteration+1 << doneTable[iteration][iteration] << endl;   //widoczna różnica, przyczyna oczekiwania   //TEŻ WYMAGANE DO DZIAŁANIA
            //        cout << "";
            //    }
            //} while (flag);
            while (doneTable[row][iteration] != iteration || doneTable[iteration][iteration] != iteration + 1) { Sleep(0); }
            float ratio = AB[row][iteration] / AB[iteration][iteration];
            //dla wszystkich pozostałych elementów w wierszu
            for (int j = iteration; j < n + 1; j++)
            {
                ////dopóki elementy nie zostały obliczone przez poprzednie iteracje to czekaj
                //bool flag;
                //do
                //{
                //    if (doneTable[row][j] == iteration && doneTable[iteration][j] == iteration + 1)
                //    {
                //        flag = false;
                //        //cout << 'y';
                //    }
                //    else
                //    {
                //        flag = true;
                //        //cout << 'x';
                //        //cout << iteration << endl;
                //    }
                //} while (flag);
                while (doneTable[row][j] != iteration || doneTable[iteration][j] != iteration + 1) { Sleep(0); }
                AB[row][j] = AB[row][j] - ratio * AB[iteration][j];
                doneTable[row][j]=iteration+1;
                //#pragma omp critical
                //{
                //    printMatrix(n, doneTable, row, j);
                //}
            }
        }
        
    }

    for (int j = 0; j < n; j++) delete[] doneTable[j];
    delete[] doneTable;
    for (int j = 0; j < (n * (1 + n)) / 2; j++) delete[] executeOrderTable[j];
    delete[] executeOrderTable;
}

void multiGaussElimination(int n, float** AB)
{
    for (int k = 0; k < n; k++)
    {
        for (int i = k + 1; i < n + 1; i++)
        {
            AB[k][i] = AB[k][i] / AB[k][k];
        }
        AB[k][k] = 1;

        #pragma omp parallel for schedule(dynamic, 1) firstprivate(k, n) shared(AB)
        for (int i = k + 1; i < n; i++)
        {
            for (int j = k + 1; j < n + 1; j++)
                AB[i][j] = AB[i][j] - AB[i][k] * AB[k][j];
            AB[i][k] = 0;
        }
    }
}

float diagonalMultiplication(int n, float ** AB)
{
    float result = 1;
    for (int i = 0; i < n; i++)
        result *= AB[i][i];
    return result;
}

void generateRandomCorrectMatrix(int n, float ** Matrix)
{
    float det = 0;
    while (det == 0)
    {
        generateRandomMatrix(n, Matrix);
        //kopia lokalna wygenerowanej macierzy
        float ** localMatrix;
        localMatrix = new float* [n];
        for (int i = 0; i < n; i++) localMatrix[i] = new float[n + 1];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n + 1; j++)
            {
                localMatrix[i][j] = Matrix[i][j];
            }
        }
        GaussElimination(n, localMatrix);
        det = diagonalMultiplication(n, localMatrix);
        for (int j = 0; j < n; j++) delete[] localMatrix[j];
        delete[] localMatrix;
    }
}

void GaussSolve(int n, float** AB, float* X)
{
    pipeGaussElimination(n, AB);
    if (TEST) printMatrix(n, AB);
    for (int i = n - 1; i >= 0; i--)
    {
        float s = AB[i][n];
        for (int j = n - 1; j > i; j--)
        {
            s -= AB[i][j] * X[j];
        }
        X[i] = s / AB[i][i];
    }
}

void multipleTestExecution(int n, float** AB, float* X, int maxTests)
{
    //badanie czasu obliczania
    vector < long long > executeTimes;
    for (int i = 1; i <= maxTests; i++)
    {
        if (!TEST) cout << '\r' << i << '/' << maxTests;
        //kopia lokalna wygenerowanej macierzy
        float** localMatrix;
        localMatrix = new float* [n];
        for (int j = 0; j < n; j++) localMatrix[j] = new float[n + 1];
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < n + 1; k++)
            {
                localMatrix[j][k] = AB[j][k];
            }
        }
        //pomiar czasu obliczania
        chrono::steady_clock::time_point start = chrono::steady_clock::now();
        GaussSolve(n, localMatrix, X);
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        //END pomiar czasu obliczania
        for (int j = 0; j < n; j++) delete[] localMatrix[j];
        delete[] localMatrix;
        executeTimes.push_back(chrono::duration_cast<chrono::milliseconds>(end - start).count());
        //pokaż wyniki
        if (TEST) for (int j = 0; j < n; j++) cout << "X" << j << ": " << X[j] << endl;
    }
    if (!TEST) cout << "\r                     \r";
    //wyświetlanie posortowanych czasów obliczania
    sort(executeTimes.begin(), executeTimes.end());
    if (!TEST) cout << setw(10) << left << omp_get_max_threads(); fs << setw(10) << left << omp_get_max_threads();
    if (!TEST) for (const auto& i : executeTimes)
    {
        cout << setw(10) << left << i; fs << setw(10) << left << i;
    }
    cout << endl; fs << endl;
    //END wyświetlanie posortowanych czasów obliczania

    //END badanie czasu obliczania
}

int main()
{
    chrono::steady_clock::time_point start = chrono::steady_clock::now();
    setlocale(LC_CTYPE, "Polish");
    srand(time(NULL));
    int matrixSizes[] = { 160, 240, 320, 400 };
    if (TEST) { matrixSizes[0] = 3; matrixSizes[1] = 4; matrixSizes[2] = 5; matrixSizes[3] = 6; }
    int maxTests = 100;
    if (TEST) maxTests = 1;
    int maxThreads = 12;
    if (TEST) maxThreads = 1;
    fs.open("out.txt", fstream::in | fstream::out | fstream::trunc);

    //pętla dla wszystkich rozmiarów macierzy
    for (int i = 0; i < 4 ; i++)
    {
        int n = matrixSizes[i];
        cout << n << endl; fs << n << endl;
        //dynamiczna macierz rozszerzona
        float ** AB;
        AB = new float* [n];
        for (int j = 0; j < n; j++) AB[j] = new float[n + 1];
        //dynamiczna macierz niewiadomych
        float * X;
        X = new float[n];

        generateRandomCorrectMatrix(n, AB);
        if (TEST) printMatrix(n, AB);
        //pętla dla wszyskich liczb wątków
        for (int j = 1; j <= maxThreads;j++)
        {
            omp_set_num_threads(j);
            if (TEST) omp_set_num_threads(4);
            multipleTestExecution(n, AB, X, maxTests);
        }
        //usuwanie macierzy dynamicznych
        for (int j = 0; j < n; j++) delete[] AB[j];
        delete[] AB;
        delete[] X;
    }
    
    fs.close();
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    int sec = chrono::duration_cast<chrono::seconds>(end - start).count();
    int min = chrono::duration_cast<chrono::minutes>(end - start).count();
    int milisec = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "Zakończono: " << min << "m " << sec << "s " << milisec << "ms " << endl;
    getchar();
}
