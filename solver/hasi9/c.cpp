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

    vector<int> tastes {61,232,170,349,41,41,333,192,159,245,68,81,346,266,312,184,184,184,409,418,338,348,332,207,253,295,160,233,388,326,429,233,218,265,265,206,173,205,41,76,410,384,302,146,252,88,140,140,413,173,413,168,91,319,375,263,32,62,13,354,402,402,208,152,64,110,432,39,437,199,411,46,219,188,83,144,364,366,321,386,386,394,408,363,256,299,256,131,247,247,225,234,103,191,260,244,249,244,109,136,54,244,397,244,32,211,279,14,374,25,158,345,393,135,258,235,166,99,213,55,55,273,143,261,214,307,58,289,187,358,172,81,81,92,92,433,44,72,72,343,343,174,148,419,87,200,200,425,425,387,138,398,398,420,370,47};
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
