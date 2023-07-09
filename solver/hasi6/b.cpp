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

struct Pillar {
    double x;
    double y;
    double r;
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
    vector<Pillar> pillars;
};

bool isBlocked(double x0, double y0, double x1, double y1, double x2, double y2, double radius) {
    complex<double> p0(x0, y0), p1(x1, y1), p2(x2, y2);
    if (real(conj(p1 - p0) * (p2 - p0)) < 0) return false;
    if (real(conj(p0 - p1) * (p2 - p1)) < 0) return false;
    double t = real(conj(p2 - p0) * (p0 - p1)) / norm(p0 - p1);
    return abs(p2 - (p0 + (p0 - p1) * t)) < radius;
}

map<int, int> tasteCount;
map<int, long long> tasteValue;

void calcIdeal(const Problem& problem, int taste) {
    double x = problem.stageLeft + 10 + rand() % ((int)problem.stageWidth - 19);
    double y = problem.stageBottom + 10 + rand() % ((int)problem.stageHeight - 19);
    long long score = 0;
    for (auto& a : problem.attendees) {
        if (a.tastes[taste] <= 0) continue;
        for (auto pr : problem.pillars) {
            if (isBlocked(a.x, a.y, x, y, pr.x, pr.y, pr.r)) {
                goto next2;
            }
        }
        {
            auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
            score += (long long)ceil(1000000 * a.tastes[taste] / d2); // NOLINT(cppcoreguidelines-narrowing-conversions)
        }
        next2:;
    }
    if (!tasteValue.count(taste) || score > tasteValue[taste]) {
        tasteValue[taste] = score;
        if (tasteValue.size() == tasteCount.size()) {
            long long res = 0;
            for (auto [t, tasteCnt] : tasteCount) {
                res += tasteValue[t] * tasteCnt;
            }
            cout << res << endl;
        }
    }
}

int main() {
    Problem problem;
    cin >> problem.roomWidth >> problem.roomHeight;
    cin >> problem.stageWidth >> problem.stageHeight;
    cin >> problem.stageLeft >> problem.stageBottom;
    int musicianN, tasteN, attendeeN, pillarN;
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
    cin >> pillarN;
    problem.pillars.resize(pillarN);
    for (auto& p : problem.pillars) {
        cin >> p.x >> p.y >> p.r;
    }

    for (auto taste : problem.musicians) {
        tasteCount[taste]++;
    }

    for (auto x : tasteCount) {
        cout << x.second << " ";
    }
    cout << endl;

    for (;;) {
        for (auto x : tasteCount) {
            calcIdeal(problem, x.first);
        }
    }
}
