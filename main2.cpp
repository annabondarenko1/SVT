#include <iostream>
#include <cmath>
#include "inmost.h"

using namespace INMOST;
using namespace std;

void
create_matrix_and_vector(int N, Sparse::Matrix& A, Sparse::Vector& b, Sparse::Vector& real_sol)
{
    for(int i = 0; i < (N - 1) * (N - 1); i++) {

        double x = (double)(i % (N - 1) + 1) / N;
        double y = (double)(i / (N - 1) + 1) / N;

        b[i] = 26 * std::sin(5 * x) * std::cos(5 * y) / N / N;

        real_sol[i] = std::sin(5 * x) * std::cos(5 * y);

        A[i][i] = 4;

        if(i % (N - 1) + 1 < (N - 1)) {
            A[i][i + 1] = -1;
        } else {
            b[i] += std::sin(5) * std::cos(5 * y);
        }

        if(i % (N - 1) - 1 >= 0) {
            A[i][i - 1] = -1;
        }

        if(i + N - 1 < (N - 1) * (N - 1)) {
            A[i][i + N - 1] = -1;
        } else {
            b[i] += std::sin(5 * x) * std::cos(5);
        }

        if(i - N + 1 >= 0) {
            A[i][i - N + 1] = -1;
        }
    }
}

void
test_function(int N, Sparse::Vector sol, Sparse::Vector real_sol)
{
    double val, l2_norm = 0, c_norm = 0;

    for(int i = 0; i < (N - 1) * (N - 1); i++) {
        val = abs(sol[i] - real_sol[i]);
        l2_norm += val * val;

        if (val > c_norm) {
            c_norm = val;
        }
    }

    std::cout << "l2_norm: " << sqrt(l2_norm) / N << std::endl;
    std::cout << "c_norm: " << c_norm << std::endl;
}

int
main()
{
    int N;
    std::cin >> N;

    Sparse::Matrix A;
    Sparse::Vector b;
    Sparse::Vector sol;
    Sparse::Vector real_sol;

    A.SetInterval(0, (N - 1) * (N - 1));
    b.SetInterval(0, (N - 1) * (N - 1));
    sol.SetInterval(0, (N - 1) * (N - 1));
    real_sol.SetInterval(0, (N - 1) * (N - 1));

    create_matrix_and_vector(N, A, b, real_sol);

    Solver S(Solver::INNER_MPTILU2);

    S.SetParameter("absolute_tolerance", "1e-10");
    S.SetParameter("relative_tolerance", "1e-6");
    
    S.SetMatrix(A);
    
    bool solved = S.Solve(b, sol);

    cout << "num.iters: " << S.Iterations() << endl;
    cout << "prec.time: " << S.PreconditionerTime() << endl;
    cout << "iter.time: " << S.IterationsTime() << endl;
    if(!solved){
        cout << "Linear solver failure!" << endl;
        cout << "Reason: " << S.ReturnReason() << endl;
    }

    test_function(N - 1, sol, real_sol);
    return 0;
}
