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

bool checkPlacements(double x, double y, const vector<pair<double, double>>& placements) {
    for (auto [x2, y2] : placements) {
        if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
            return false;
        }
    }
    return true;
}

pair<int, int> makeRandPos(const Problem& problem, const vector<pair<double, double>>& placements) {
    for (;;) {
        double x0 = problem.stageLeft + 10 + rand() % ((int)problem.stageWidth - 19);
        double y0 = problem.stageBottom + 10 + rand() % ((int)problem.stageHeight - 19);
        if (checkPlacements(x0, y0, placements)) return {x0, y0};
    }
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
    vector<pair<double, double>> planPos(nPlacements);
    for (auto& p : planPos) {
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

    vector<int> tastes {62,8,8,85,85,85,85,85,85,85,85,85,85,85,85,8,8,8,85,49,8,8,8,8,8,8,8,18,18,18,18,18,18,18,18,18,18,49,49,52,49,49,49,49,49,49,49,49,49,49,55,55,55,55,106,106,106,106,106,106,106,106,64,106,106,106,106,106,106,59,59,59,59,59,59,59,59,59,59,59,19,59,59,64,64,64,64,64,19,19,19,19,19,19,19,19,19,64,64,64,64,64,64,19,62,62,62,62,62,62,62,62,62,62,104,104,104,104,104,104,104,104,43,63,63,63,74,63,63,5,5,35,5,5,5,74,74,74,11,94,94,2,94,2,94,94,2,2,94,70,70,70,70,70,70,70,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,32,27,27,27,27,27,5,5,5,5,5,5,5,5,5,32,5,5,5,5,5,5,5,5,13,29,49,86,57,52,52,52,52,52,52,93,93,93,93,93,93,93,93,93,99,99,99,99,99,99,99,99,99,99,99,99,99,99,86,86,13,13,39,2,20,20,2,2,12,20,20,20,20,52,52,52,52,52,52,52,52,52,52,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,36,36,36,36,36,36,36,36,36,36,36,36,36,36,36,1,1,1,1,1,95,95,95,36,67,67,104,104,76,76,104,104,67,67,67,56,56,6,104,17,17,17,17,29,29,29,17,29,29,29,29,29,29,29,29,29,29,74,74,74,74,74,74,74,46,46,46,46,46,46,74,1,1,1,1,1,1,1,1,1,17,17,17,17,17,17,53,53,53,53,53,53,53,53,90,90,90,90,90,90,90,90,90,53,90,53,53,79,79,90,90,90,90,90,90,90,90,53,53,53,53,16,16,16,16,16,16,90,16,16,16,61,61,61,61,61,61,61,61,61,61,61,61,47,47,47,47,72,72,72,72,83,83,57,57,57,57,57,83,83,83,83,83,83,83,83,83,47,102,62,47,47,47,106,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,10,10,10,10,10,10,10,10,10,10,56,56,56,30,0,0,0,0,0,0,0,0,57,57,44,44,44,44,44,44,44,44,44,44,47,47,47,47,47,39,30,30,30,30,30,30,98,98,30,30,77,77,77,77,41,41,41,41,41,41,41,41,64,64,64,41,41,41,78,78,78,78,78,78,84,84,84,84,84,84,84,84,84,84,84,84,84,84,84,79,79,79,71,71,71,71,79,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,30,30,30,30,30,30,30,30,30,28,28,28,104,104,104,104,104,104,60,60,60,60,60,60,60,60,60,28,28,28,28,28,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,102,102,102,102,102,28,28,28,28,39,39,39,39,39,73,38,45,77,11,87,71,23,69,65,42,6,7,98,34,91,65,35,16,35,69,31,91,8,27,38,9,82,45,75,95,98,96,33,71,40,80,89,80,87,65,75,23,50,97,37,3,42,69,50,94,96,91,103,25,25,7,103,101,75,98,14,23,55,54,50,33,14,77,91,96,98,65,77,22,4,24,38,103,4,82,31,103,98,26,9,26,50,75,105,42,78,65,11,82,11,22,96,73,14,72,38,43,31,69,88,80,12,87,34,48,7,91,67,26,9,9,27,105,96,75,51,51,51,51,51,51,51,79,79,79,79,79,79,79,79,79,79,51,54,54,81,54,54,54,70,70,70,70,70,70,70,70,70,21,21,21,21,21,21,21,21,21,21,21,21,70,70,70,70,70,70,75,75,57,57,57,57,57,54,54,54,54,54,86,103,6,27,88,39,23,11,42,14,46,42,101,23,98,24,77,40,80,7,48,40,26,65,42,43,12,65,39,58,6,46,101,105,96,38,98,23,89,72,98,72,68,58,87,98,97,65,38,91,65,69,26,11,103,6,3,89};
    map<int, vector<unsigned>> plan;
    for (unsigned i = 0; i < tastes.size(); i++) {
        plan[tastes[i]].push_back(i);
    }

    vector<pair<double, double>> placements(problem.musicians.size());
    vector<pair<double, double>> used;
    vector<unsigned> noPlan;
    for (unsigned i = 0; i < problem.musicians.size(); ++ i) {
        if (plan.count(problem.musicians[i]) && !plan[problem.musicians[i]].empty()) {
            auto idx = plan[problem.musicians[i]].back();
            plan[problem.musicians[i]].pop_back();
            placements[i] = planPos[idx];
            used.push_back(placements[i]);
        } else {
            noPlan.push_back(i);
        }
    }
    for (auto i : noPlan) {
        placements[i] = makeRandPos(problem, used);
        used.push_back(placements[i]);
    }

    cout << "{\"placements\":[";
    for (unsigned i = 0; i < placements.size(); i++) {
        if (i > 0) cout << ",";
        cout << "{\"x\":" << placements[i].first << ",\"y\":" << placements[i].second << "}";
    }
    cout << "]}" << endl;

    auto [ok, score] = calcScore(problem, placements);
    if (!ok) {
        cerr << "Invalid placements" << endl;
        return 1;
    }
    cerr << "score = " << score << endl;
}
