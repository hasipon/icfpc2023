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

int main(int argc, char** argv) {
    if (argc < 2) {
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
    int nPlacements;
    fin >> nPlacements;
    vector<pair<double, double>> placements(nPlacements);
    for (auto& p : placements) {
        fin >> p.first >> p.second;
    }
    vector<pair<int, int>> tasteCnt;
    vector<vector<long long>> scoreTable;
    {
        int taste, nTaste;
        while (fin >> taste >> nTaste) {
            tasteCnt.emplace_back(taste, nTaste);
            scoreTable.emplace_back(nPlacements);
            for (auto& s : scoreTable.back()) {
                fin >> s;
            }
        }
    }

    vector<int> tasteIdxes;
    for (unsigned i = 0; i < tasteCnt.size(); i++) {
        for (int j = 0; j < tasteCnt[i].second; j++) {
            tasteIdxes.push_back((int)i);
        }
    }

    long long bestScore = 0;
    for (int iter = 0; iter < 1000000000; ++ iter) {
        random_shuffle(tasteIdxes.begin(), tasteIdxes.end());

        long long score = 0;
        for (int i = 0; i < nPlacements; i++) {
            score += scoreTable[tasteIdxes[i]][i];
        }

        for (int tt = 0; tt < 10000000; ++ tt) {
            int i1 = rand() % nPlacements;
            int i2 = rand() % tasteIdxes.size();
            if (i1 == i2) continue;
            long long score2 = score;
            score2 -= scoreTable[tasteIdxes[i1]][i1];
            score2 += scoreTable[tasteIdxes[i2]][i1];
            if (i2 < nPlacements) {
                score2 -= scoreTable[tasteIdxes[i2]][i2];
                score2 += scoreTable[tasteIdxes[i1]][i2];
            }
            if (score2 > score) {
                swap(tasteIdxes[i1], tasteIdxes[i2]);
                score = score2;
            }
        }
        if (score > bestScore) {
            bestScore = score;

            long long scoreCheck = 0;
            for (int i = 0; i < nPlacements; i++) {
                scoreCheck += scoreTable[tasteIdxes[i]][i];
            }
            if (score != scoreCheck) {
                cerr << "score mismatch: " << score << " != " << scoreCheck << endl;
                return 1;
            }

            cout << score << " ";
            for (int i = 0; i < nPlacements; i++) {
                cout << tasteCnt[tasteIdxes[i]].first << " ";
            }
            cout << endl;
        }
    }
}
