#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <complex>
#include <algorithm>
using namespace std;

struct Attendee {
    double x;
    double y;
    vector<double> tastes;
};

struct Problem {
    double roomWidth;
    double roomHeight;
    double stageWidth;
    double stageHeight;
    double stageBottom;
    double stageLeft;
    vector<int> musicians;
    vector<Attendee> attendees;
};

long long calcIdeal(const Problem& problem, int taste) {
    double bestScore = -1;
    pair<int, int> bestPos;
    for (int x = (int)problem.stageLeft + 10; x <= (int)problem.stageLeft + (int)problem.stageWidth - 10; ++ x) {
        for (int y = (int)problem.stageBottom + 10; y <= (int)problem.stageBottom + (int)problem.stageHeight - 10; ++ y) {
            double score = 0;
            for (auto& a : problem.attendees) {
                if (a.tastes[taste] <= 0) continue;
                auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
                score += a.tastes[taste] / d2;
            }
            if (score > bestScore) {
                bestScore = score;
                bestPos = {x, y};
            }
        }
    }
    auto [x, y] = bestPos;
    long long score = 0;
    for (auto& a : problem.attendees) {
        if (a.tastes[taste] <= 0) continue;
        auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
        score += (long long)ceil(1000000 * a.tastes[taste] / d2); // NOLINT(cppcoreguidelines-narrowing-conversions)
    }
    return score;
}

int main() {
    Problem problem;
    cin >> problem.roomWidth >> problem.roomHeight;
    cin >> problem.stageWidth >> problem.stageHeight;
    cin >> problem.stageLeft >> problem.stageBottom;
    int musicianN, tasteN, attendeeN;
    cin >> musicianN >> tasteN;
    problem.musicians.resize(musicianN);
    for (auto& m : problem.musicians) {
        cin >> m;
    }
    cin >> attendeeN;
    problem.attendees.resize(attendeeN);
    for (auto& a : problem.attendees) {
        cin >> a.x >> a.y;
        a.tastes.resize(tasteN);
        for (auto& t : a.tastes) {
            cin >> t;
        }
    }

    map<int, int> tasteCount;
    for (auto taste : problem.musicians) {
        tasteCount[taste]++;
    }

    long long ideal = 0;
    for (auto x : tasteCount) {
        ideal += calcIdeal(problem, x.first) * x.second;
    }
    cout << ideal << endl;
}
