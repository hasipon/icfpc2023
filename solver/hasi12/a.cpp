#include <iostream>
#include <iomanip>
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

pair<bool, vector<long long>> calcScore(const Problem& problem, const vector<pair<double, double>>& placements) {
    if (placements.size() != problem.musicians.size()) {
        return {false, {}};
    }
    for (unsigned i = 0; i < placements.size(); i++) {
        auto [x, y] = placements[i];
        if (!(
                problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
                problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10
        )) {
            return {false, {}};
        }
        for (unsigned j = 0; j < i; j++) {
            auto [x2, y2] = placements[j];
            if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
                return {false, {}};
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
    vector<long long> score(placements.size());
    for (auto& a : problem.attendees) {
        for (unsigned i = 0; i < placements.size(); i++) {
            auto [x, y] = placements[i];
            for (unsigned j = 0; j < placements.size(); j++) {
                if (i != j) {
                    auto [x2, y2] = placements[j];
                    if (isBlocked(a.x, a.y, x, y, x2, y2, 5)) {
                        goto next;
                    }
                }
            }
            for (auto pr : problem.pillars) {
                if (isBlocked(a.x, a.y, x, y, pr.x, pr.y, pr.r)) {
                    goto next;
                }
            }
            {
                auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
                auto iik = (long long)ceil(1000000 * a.tastes[problem.musicians[i]] / d2);
                if (factor[i] > 1) {
                    score[i] += (long long)ceil(10 * factor[i] * iik); // NOLINT(cppcoreguidelines-narrowing-conversions)
                } else {
                    score[i] += 10 * iik;
                }
            }
            next:;
        }
    }
    return {true,score};
}

long long sumScore(const vector<long long>& score) {
    long long sum = 0;
    for (auto s : score) {
        if (s > 0) {
            sum += s;
        }
    }
    return sum;
}

bool check(const Problem& problem, const vector<pair<double, double>>& placements, unsigned idx, double x, double y) {
    if (!(problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 && problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10)) return false;
    for (unsigned j = 0; j < placements.size(); j++) {
        if (j == idx) continue;
        auto [x2, y2] = placements[j];
        if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100 + 1e-6) {
            return false;
        }
    }
    return true;
}

pair<double, double> move(const Problem& problem, const vector<pair<double, double>>& placements, int idx) {
    double theta = rand() * 2 * M_PI / RAND_MAX;
    double L = 0, R = 1 + (rand() * 9.0) / RAND_MAX;
    while (R - L >= 1e-4) {
        auto M = (L + R) / 2;
        auto x = placements[idx].first + M * cos(theta);
        auto y = placements[idx].second + M * sin(theta);
        if (check(problem, placements, idx, x, y)) {
            L = M;
        } else {
            R = M;
        }
    }
    if (L < 2e-3) return {-1, -1};
    return {placements[idx].first + L * cos(theta), placements[idx].second + L * sin(theta)};
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

    ifstream fin(argv[1]);
    int n;
    fin >> n;
    vector<pair<double, double>> placements(n);
    for (auto& p : placements) {
        fin >> p.first >> p.second;
    }

    auto [ok, aScore] = calcScore(problem, placements);
    if (!ok) {
        cerr << "Invalid placements" << endl;
        return 1;
    }
    auto score = sumScore(aScore);
    cerr << "score = " << score << endl;

    for (;;) {
        int idx = rand() % placements.size();
        auto prev = placements[idx];
        placements[idx] = move(problem, placements, idx);
        if (placements[idx].first < 0) {
            placements[idx] = prev;
            continue;
        }
        auto [ok2, aScore2] = calcScore(problem, placements);
        if (!ok2) {
            cerr << "Invalid placements" << endl;
            return 1;
        }
        auto score2 = sumScore(aScore2);
        if (score2 > score) {
            score = score2;
            cerr << "score = " << score << endl;
            cout << "{\"placements\":[";
            for (unsigned i = 0; i < placements.size(); i++) {
                if (i > 0) cout << ",";
                cout << "{\"x\":" << std::setprecision(16) << placements[i].first << ",\"y\":" << std::setprecision(16) << placements[i].second << "}";
            }
            cout << "],\"volumes\":[";
            for (unsigned i = 0; i < aScore.size(); i++) {
                if (i > 0) cout << ",";
                cout << (aScore[i] > 0 ? "10" : "0");
            }
            cout << "]}" << endl;
        } else {
            placements[idx] = prev;
        }
    }
}
