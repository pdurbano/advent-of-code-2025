#ifndef _SNF_H_
#define _SNF_H_

#include <algorithm>
#include <iostream>
#include <optional>
#include <vector>

//
// Support Classes
//

struct GCD {
    GCD(int a, int b)
    {
        int x0 = a < 0 ? -1 : 1;
        int y0 = 0;

        int x1 = 0;
        int y1 = b < 0 ? -1 : 1;

        int c0 = a * x0;
        int c1 = b * y1;

        while (c1 != 0) {
            if (c1 > c0) {
                c1 = c1 - c0;
                x1 = x1 - x0;
                y1 = y1 - y0;
            }
            else {
                int temp;

                temp = c1;
                c1 = c0 - c1;
                c0 = temp;

                temp = x1;
                x1 = x0 - x1;
                x0 = temp;

                temp = y1;
                y1 = y0 - y1;
                y0 = temp;
            }
        }

        bet = c0;
        sig = x0;
        tau = y0;
    }

    int bet;
    int sig;
    int tau;
};

struct Mat {
    struct Row {
        Row(int *p) : p(p) {}
        int &operator[](int col) { return p[col]; }

      private:
        int *p;
    };

    struct ConstRow {
        ConstRow(const int *p) : p(p) {}
        const int &operator[](int col) const { return p[col]; }

      private:
        const int *p;
    };

    Mat() : M(0), N(0) {}

    Mat(size_t rows, size_t cols) : M(rows), N(cols), mem(rows * cols, 0) {}

    Mat(const std::vector<std::vector<int>> &plan, size_t rows) : M(rows), N(plan.size()), mem(rows * plan.size(), 0)
    {
        for (size_t n = 0; n < N; ++n)
            for (auto m : plan[n])
                mem[m * N + n] = 1;
    }

    Mat(const std::vector<int> *plan, size_t plan_size, size_t rows) : M(rows), N(plan_size), mem(rows * plan_size, 0)
    {
        for (size_t n = 0; n < N; ++n)
            for (auto m : plan[n])
                mem[m * N + n] = 1;
    }

    ConstRow operator[](int k) const { return ConstRow(&mem[k * N]); }
    Row operator[](int k) { return Row(&mem[k * N]); }

    void print() const
    {
        const int *p = mem.data();
        for (size_t m = 0; m < M; ++m) {
            std::cout << "|";
            for (size_t n = 0; n < N - 1; ++n)
                std::cout << *p++ << ", ";
            std::cout << *p++ << "|\n";
        }
    }

    Mat operator*(const Mat &o) const
    {
        Mat ret(M, o.N);
        auto K = N;
        if (K != o.M)
            return {};

        for (size_t m = 0; m < ret.M; ++m)
            for (size_t n = 0; n < ret.N; ++n)
                for (size_t k = 0; k < K; ++k)
                    ret[m][n] += (*this)[m][k] * o[k][n];

        return ret;
    }

    Mat T() const
    {
        Mat ret(N, M);
        for (size_t m = 0; m < ret.M; ++m)
            for (size_t n = 0; n < ret.N; ++n)
                ret[m][n] = (*this)[n][m];
        return ret;
    }

    static Mat eye(size_t N)
    {
        Mat ret(N, N);
        for (size_t n = 0; n < N; ++n)
            ret[n][n] = 1;
        return ret;
    }

    static Mat swap_two(int a, int b, size_t N)
    {
        auto ret = eye(N);
        if (a == b)
            return ret;

        ret[a][a] = 0;
        ret[a][b] = 1;
        ret[b][a] = 1;
        ret[b][b] = 0;

        return ret;
    }

    static Mat two_row_op(int a, int b, const Mat &op, size_t N)
    {
        auto ret = eye(N);

        if (op.M < 2 || op.N < 2)
            return ret;

        if (a == b) {
            ret[a][b] = op[0][0];
            return ret;
        }

        ret[a][a] = op[0][0];
        ret[a][b] = op[0][1];
        ret[b][a] = op[1][0];
        ret[b][b] = op[1][1];

        return ret;
    }

    static Mat add_col(int a, int b, int scale_a, size_t N)
    {
        auto ret = eye(N);
        if (a == b) {
            ret[a][b] = 1 + scale_a;
            return ret;
        }

        ret[a][b] = scale_a;

        return ret;
    }

    static Mat add_row(int a, int b, int scale_a, size_t N)
    {
        auto ret = eye(N);
        if (a == b) {
            ret[b][a] = 1 + scale_a;
            return ret;
        }

        ret[b][a] = scale_a;

        return ret;
    }

    size_t M, N;

  private:
    std::vector<int> mem;
};

//
// Smith Normal Form
//
struct SNF {
    SNF(Mat A)
    {
        auto M = A.M;
        auto N = A.N;

        S = Mat::eye(M);
        T = Mat::eye(N);

        int j = -1;
        for (int t = 0; t < M; ++t) {
            int j_prev = j;

            // Step I: Choosing a pivot
            {
                bool found = false;
                for (j = j_prev + 1; j < N; ++j) {
                    for (int k = t; k < M; ++k) {
                        if (A[k][j] != 0) {
                            found = true;
                            if (k != t) {
                                // Move pivot to row t.
                                auto L = Mat::swap_two(k, t, M);
                                A = L * A;
                                S = L * S;
                            }
                            break;
                        }
                    }

                    if (found)
                        break;
                }

                if (!found)
                    // Only zeros remain, nothing left to pivot on
                    break;
            }

            int a_tj = A[t][j];

            // Loop over steps II and III for rows and columns until pivot for row t is complete
            bool done = false;
            while (!done) {
                // Columns - Step II: Improving the pivot
                for (int k = j_prev + 1; k < N; ++k) {
                    if (k == j)
                        continue;

                    int a_tk = A[t][k];
                    if (a_tk == 0)
                        continue;
                    if (a_tk % a_tj == 0)
                        continue;

                    auto [bet, sig, tau] = GCD(a_tj, a_tk);
                    int alp = a_tj / bet;
                    int gam = a_tk / bet;

                    Mat R0(2, 2);
                    R0[0][0] = sig;
                    R0[0][1] = -gam;
                    R0[1][0] = tau;
                    R0[1][1] = alp;
                    auto R = Mat::two_row_op(j, k, R0, N);
                    A = A * R;
                    T = T * R;
                    a_tj = A[t][j];
                }

                // Columns - Step III: Eliminating entries
                for (int k = j_prev + 1; k < N; ++k) {
                    if (k == j)
                        continue;

                    int a_tk = A[t][k];
                    if (a_tk == 0)
                        continue;

                    int f = a_tk / a_tj;
                    auto R = Mat::add_col(j, k, -f, N);
                    A = A * R;
                    T = T * R;
                    a_tj = A[t][j];
                }

                // Rows - Step II: Improving the pivot
                for (int k = t + 1; k < M; ++k) {
                    int a_kj = A[k][j];
                    if (a_kj == 0)
                        continue;
                    if (a_kj % a_tj == 0)
                        continue;

                    auto [bet, sig, tau] = GCD(a_tj, a_kj);
                    int alp = a_tj / bet;
                    int gam = a_kj / bet;

                    Mat L0(2, 2);
                    L0[0][0] = sig;
                    L0[0][1] = tau;
                    L0[1][0] = -gam;
                    L0[1][1] = alp;
                    auto L = Mat::two_row_op(t, k, L0, M);
                    A = L * A;
                    S = L * S;
                    a_tj = A[t][j];
                }

                // Rows - Step III: Eliminating entries
                for (int k = t + 1; k < M; ++k) {
                    int a_kj = A[k][j];
                    if (a_kj == 0)
                        continue;

                    int f = a_kj / a_tj;
                    auto L = Mat::add_row(t, k, -f, M);
                    A = L * A;
                    S = L * S;
                    a_tj = A[t][j];
                }

                // Confirm pivot is complete
                done = true;
                for (int k = j_prev + 1; k < N; ++k) {
                    if (k == j)
                        continue;

                    if (A[t][k] != 0) {
                        done = false;
                        break;
                    }
                }

                if (!done)
                    continue;

                for (int k = t + 1; k < M; ++k) {
                    if (A[k][j] != 0) {
                        done = false;
                        break;
                    }
                }
            }
        }

        // Shift null columns right
        for (int n = 0; n < N - 1; ++n) {
            bool all_zero = true;
            for (int m = 0; m < M; ++m) {
                if (A[m][n] != 0) {
                    all_zero = false;
                    break;
                }
            }
            if (!all_zero)
                continue;

            all_zero = true;
            for (int k = n + 1; k < N; ++k) {
                for (int m = 0; m < M; ++m) {
                    if (A[m][k] != 0) {
                        all_zero = false;
                        auto R = Mat::swap_two(n, k, N);
                        A = A * R;
                        T = T * R;
                        break;
                    }
                }
                if (!all_zero)
                    break;
            }
            if (all_zero)
                break;
        }

        B = A;
    }

    Mat B;
    Mat S;
    Mat T;
};

namespace Solver {

struct Result {
    Mat X;
    Mat Y;
    int free_start;
    int free_end;
};

std::optional<Result> solve(const Mat &B, const Mat &S, const Mat &T, const std::vector<int> &c)
{
    Mat C(c.size(), 1);
    for (int n = 0; n < c.size(); ++n)
        C[n][0] = c[n];

    auto D = S * C;

    auto max_rank = std::min(B.M, B.N);
    for (int i = max_rank; i < D.M; ++i)
        if (D[i][0] != 0)
            return {};

    int free_start = max_rank;
    int free_end = B.N;

    // Construction initializes all locations to zero so we don't necessarily have to write them all below.
    Mat Y(free_end, 1);

    for (int i = 0; i < max_rank; ++i) {
        if (B[i][i] != 0) {
            Y[i][0] = D[i][0] / B[i][i];
            if (D[i][0] != Y[i][0] * B[i][i])
                return {};
        }
        else {
            if (D[i][0] != 0)
                return {};
            if (i < free_start)
                free_start = i;
        }
    }

    auto X = T * Y;

    return Result{X, Y, free_start, free_end};
}

}; // namespace Solver

#endif