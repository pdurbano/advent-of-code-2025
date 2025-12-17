/*
 * clang++ -O3 -std=c++20 day10.cpp -o day10
 */

#include "snf.h"

#include <charconv>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

std::pair<std::vector<std::vector<int>>, std::vector<int>> parse_line(const std::string &line)
{
    constexpr bool verbose = false;

    if (verbose)
        std::cout << line << "\n";

    const char *begin = line.data();
    const char *end = begin + line.size();

    std::vector<std::vector<int>> buttons;
    std::vector<int> joltage;
    char last_delim = '\0';
    bool curly_braces = false;
    while (begin < end) {
        char next1_delim;
        char next2_delim;
        char term1_delim;
        char term2_delim;
        switch (last_delim) {
        case '{':
        case '(':
        case ',':
            next1_delim = ',';
            next2_delim = '\0';
            term1_delim = ')';
            term2_delim = '}';
            break;
        case ')':
        case '}':
        default:
            next1_delim = '(';
            next2_delim = '{';
            term1_delim = '\0';
            term2_delim = '\0';
            break;
        }

        while (begin < end &&
               (*begin != next1_delim && *begin != next2_delim && *begin != term1_delim && *begin != term2_delim))
            ++begin;

        if (begin >= end - 1)
            break;
        last_delim = *begin;
        ++begin;

        if (last_delim == term1_delim)
            continue;

        if (last_delim == term2_delim)
            break;

        if (last_delim == '(')
            buttons.emplace_back();
        else if (last_delim == '{')
            curly_braces = true;

        int value;
        auto [ptr, ec] = std::from_chars(begin, end, value);
        if (ec == std::errc()) {
            if (!curly_braces)
                buttons.back().emplace_back(value);
            else
                joltage.emplace_back(value);
            begin = ptr;
        }
        else {
            std::cerr << "Parse error\n";
            break;
        }
    }

    if (verbose) {
        for (const auto &button : buttons) {
            std::cout << "(";
            for (auto value : button)
                std::cout << value << ",";
            std::cout << ")\n";
        }
        std::cout << "{";
        for (auto value : joltage)
            std::cout << value << ",";
        std::cout << "}\n";
    }

    return {buttons, joltage};
}

bool check_match(const std::vector<int> *buttons, size_t num_buttons, const std::vector<int> &joltage, int pushes)
{
    for (auto value : joltage)
        if (value > pushes)
            return false;

    if (pushes == 0)
        return true;

    if (num_buttons == 0)
        return false;

    // Use Smith normal form of data matrix to test feasibility before recursing further.
    Mat A(buttons, num_buttons, joltage.size());
    auto [B, S, T] = SNF(A);

    auto soln = Solver::solve(B, S, T, joltage);
    if (!soln)
        return false;

    auto [X, _, free_start, free_end] = *soln;
    if (free_start >= free_end) {
        // no free variables
        int sum = 0;
        for (int m = 0; m < X.M; ++m) {
            if (X[m][0] < 0)
                return false;
            sum += X[m][0];
        }

        if (sum == pushes)
            return true;

        return false;
    }

    const auto &button = buttons[num_buttons - 1];

    int start = std::numeric_limits<int>::max();
    for (auto wire : button) {
        auto jolts = joltage[wire];
        if (jolts < start)
            start = jolts;
    }

    for (int n = start; n >= 0; --n) {
        std::vector<int> new_joltage(joltage);

        for (auto wire : button)
            new_joltage[wire] -= n;

        if (check_match(buttons, num_buttons - 1, new_joltage, pushes - n))
            return true;

        // If this is the last button we've already failed. No reason to back off after optimal # of pushes.
        if (num_buttons == 1)
            return false;
    }

    return false;
}

int main()
{
    std::ifstream file("day_10_input.txt");
    if (!file) {
        std::cerr << "Failed to open file.\n";
        return 1;
    }

    int total = 0;
    int line_num = 0;
    std::string line;
    while (std::getline(file, line)) {
        std::cout << "[" << line_num << "]  " << line << "\n";
        line_num++;

        auto [buttons, joltage] = parse_line(line);

        int best_case = -1;
        int worst_case = 0;
        for (auto value : joltage) {
            worst_case += value;
            if (value > best_case)
                best_case = value;
        }

        bool found = false;

        // Search every feasible number of pushes in increasing order until we find one that works
        for (int n = best_case; n < worst_case; ++n) {
            if (check_match(buttons.data(), buttons.size(), joltage, n)) {
                found = true;
                std::cout << n << " pushes\n";
                total += n;
                break;
            }
        }

        if (!found) {
            std::cerr << "Failed to solve the problem on this line.\n";
            return 1;
        }
    }

    std::cout << "Part 2 Result: " << total << "\n";

    return 0;
}
