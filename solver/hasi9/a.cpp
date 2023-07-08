#include <iostream>
#include <fstream>
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

pair<bool, long long> calcScore(const Problem& problem, vector<pair<double, double>> placements) {
    if (placements.size() != problem.musicians.size()) {
        return {false, 0};
    }
    for (unsigned i = 0; i < placements.size(); i++) {
        auto [x, y] = placements[i];
        if (!(
                problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
                problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10
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
    vector<double> factor(placements.size(), 1);
    if (!problem.pillars.empty()) {
        for (unsigned i = 0; i < placements.size(); i++) {
            auto [x1, y1] = placements[i];
            for (unsigned j = 0; j < placements.size(); j++) {
                if (j != i && problem.musicians[i] == problem.musicians[j]) {
                    auto [x2, y2] = placements[j];
                    auto d2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
                    factor[i] += 1 / sqrt(d2);
                }
            }
        }
    }
    double score = 0;
    for (auto& a : problem.attendees) {
        for (unsigned i = 0; i < placements.size(); i++) {
            auto [x, y] = placements[i];
            for (auto pr : problem.pillars) {
                if (isBlocked(a.x, a.y, x, y, pr.x, pr.y, pr.r)) {
                    goto next;
                }
            }
            for (unsigned j = 0; j < placements.size(); j++) {
                if (i != j) {
                    auto [x2, y2] = placements[j];
                    if (isBlocked(a.x, a.y, x, y, x2, y2, 5)) {
                        goto next;
                    }
                }
            }
            {
                auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
                auto s = (long long)ceil(1000000 * a.tastes[problem.musicians[i]] / d2);
                if (factor[i] > 1) {
                    s = (long long)ceil(factor[i] * s); // NOLINT(cppcoreguidelines-narrowing-conversions)
                }
                score += s; // NOLINT(cppcoreguidelines-narrowing-conversions)
            }
            next:;
        }
    }
    return {true,score};
}


long long calcScore9(const Problem& problem, vector<pair<double, double>> placements, unsigned idx, int taste) {
    long long score = 0;
    auto [x, y] = placements[idx];
    for (auto& a : problem.attendees) {
        for (unsigned j = 0; j < placements.size(); j++) {
            if (j == idx) continue;
            auto [x2, y2] = placements[j];
            if (isBlocked(a.x, a.y, x, y, x2, y2, 5)) {
                goto next;
            }
        }
        {
            auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
            score += (long long)ceil(1000000 * a.tastes[taste] / d2);
        }
        next:;
    }
    return score;
}

bool checkPlacements(double x, double y, const vector<pair<double, double>>& placements) {
    for (auto [x2, y2] : placements) {
        if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
            return false;
        }
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc < 1) {
        cerr << "Usage: " << argv[0] << " <input file>" << endl;
        return 1;
    }

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

    vector<pair<double, double>> placements;
    {
        auto y = problem.stageBottom + problem.stageHeight - 10;
        for (int x = (int)problem.stageLeft + 10; x <= problem.stageLeft + problem.stageWidth - 10; x += 10) {
            placements.emplace_back(x, y);
        }
    }
    {
        auto x = problem.stageLeft + problem.stageWidth - 10;
        for (int y = (int)problem.stageBottom + 10; y <= problem.stageBottom + problem.stageHeight - 10; y += 10) {
            if (checkPlacements(x, y, placements)) {
                placements.emplace_back(x, y);
            }
        }
    }
    for (int x = (int)problem.stageLeft + 10; x <= problem.stageLeft + problem.stageWidth - 10; x += 1) {
        for (int y = (int)problem.stageBottom + 10; y <= problem.stageBottom + problem.stageHeight - 10; y += 1) {
            if (!checkPlacements(x, y, placements)) continue;
            for (auto& a : problem.attendees) {
                for (auto [x2, y2] : placements) {
                    if (isBlocked(a.x, a.y, x, y, x2, y2, 5)) {
                        goto next;
                    }
                }
            }
            cerr << "leak " << x << " " << y << endl;
            next:;
        }
    }
    map<int, int> tasteCount;
    for (auto taste : problem.musicians) {
        tasteCount[taste]++;
    }
    cout << placements.size() << endl;
    for (auto [x, y] : placements) {
        cout << x << " " << y << endl;
    }
    for (auto [taste, cnt] : tasteCount) {
        cout << taste << " " << cnt << endl;
        for (unsigned idx = 0; idx < placements.size(); idx++) {
            cout << calcScore9(problem, placements, idx, taste) << " ";
        }
        cout << endl;
    }
}
