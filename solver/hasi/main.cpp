#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
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

bool isBlocked(double x0, double y0, double x1, double y1, double x2, double y2) {
    complex<double> p0(x0, y0), p1(x1, y1), p2(x2, y2);
    if (imag(conj(p1 - p0) * (p2 - p0)) < 0) return false;
    if (imag(conj(p0 - p1) * (p2 - p1)) < 0) return false;
    double t = real(conj(p2 - p0) * (p0 - p1)) / norm(p0 - p1);
    return abs(p2 - (p0 + (p0 - p1) * t)) <= 5;
}

pair<bool, long long> calcScore(const Problem& problem, vector<pair<double, double>> placements) {
    if (placements.size() != problem.musicians.size()) {
        return {false, 0};
    }
    for (unsigned i = 0; i < placements.size(); i++) {
        auto [x, y] = placements[i];
        if (!(
            problem.stageBottom + 10 <= x && x <= problem.stageBottom + problem.stageHeight - 10 &&
            problem.stageLeft + 10 <= y && y <= problem.stageLeft + problem.stageWidth - 10
        )) {
            return {false, 0};
        }
        for (unsigned j = 0; j < i; j++) {
            auto [x2, y2] = placements[j];
            if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
                return {false, 0};
            }
        }
    }
    double score = 0;
    for (auto& a : problem.attendees) {
        for (unsigned i = 0; i < placements.size(); i++) {
            auto [x, y] = placements[i];
            for (unsigned j = 0; j < placements.size(); j++) if (i != j) {
                auto [x2, y2] = placements[j];
                if (isBlocked(a.x, a.y, x, y, x2, y2)) {
                    goto next;
                }
            }
            {
                auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
                score += (long long)ceil(1000000 * a.tastes[problem.musicians[i]] / d2); // NOLINT(cppcoreguidelines-narrowing-conversions)
            }
            next:;
        }
    }
    return {true,score};
}

vector<pair<double, double>> solve(const Problem& problem) {
    vector<pair<double, double>> res(problem.musicians.size());
    double x = problem.stageBottom + 10;
    double y = problem.stageLeft + 10;
    for (auto & re : res) {
        re = {x, y};
        x += 10;
        if (x > problem.stageBottom + problem.stageHeight - 10) {
            x = problem.stageBottom + 10;
            y += 10;
        }
    }
    return res;
}

int main() {
    Problem problem;
    cin >> problem.roomWidth >> problem.roomHeight;
    cin >> problem.stageWidth >> problem.stageHeight;
    cin >> problem.stageBottom >> problem.stageLeft;
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

    auto placement = solve(problem);
    auto res = calcScore(problem, placement);
    if (!res.first) throw runtime_error("invalid placement");
    cerr << "score = " << res.second << endl;
}
