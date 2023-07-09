#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <complex>
#include <algorithm>
using namespace std;

int N;
vector<vector<long long>> table;
vector<vector<long long>> dp;

int main() {
    cin >> N;
    table = vector<vector<long long>>(981, vector<long long>(N));
    for (auto& a : table) {
        for (auto& x : a) {
            cin >> x;
        }
    }
    long long result = -1;
    dp = vector<vector<long long>>(981, vector<long long>(1<<N, -1));
    const int M = (1<<N) - 1;
    dp[0][0] = 0;
    for (int x = 0; x < 981; ++ x) {
        for (int p = 0; p < M-1; ++ p) {
            if (x > 0 && dp[x-1][p] > dp[x][p]) {
                dp[x][p] = dp[x-1][p];
            }
            if (dp[x][p] == -1) continue;
            for (int i = 0; i < N; ++ i) {
                if (p & (1<<i)) continue;
                int pp = p|(1<<i);
                auto next = dp[x][p] + table[x][i];
                if (pp == M) {
                    if (next > result) {
                        result = next;
                    }
                } else if (x + 10 < 981) {
                    if (next > dp[x+10][pp]) {
                        dp[x+10][pp] = next;
                    }
                }
            }
        }
    }

    cout << result << endl;
}
