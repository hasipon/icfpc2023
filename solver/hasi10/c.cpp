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

    vector<int> tastes {110,172,172,110,131,131,265,48,254,254,254,313,318,60,194,313,84,318,160,160,352,42,42,46,46,110,79,79,316,139,200,112,112,112,139,112,112,56,316,163,163,163,299,256,256,111,111,308,111,111,107,107,10,224,224,7,7,7,278,39,39,224,357,107,223,223,230,92,92,92,92,338,338,113,68,68,68,34,34,166,208,63,63,34,308,308,208,83,83,83,83,83,308,206,206,206,198,198,198,198,198,373,373,373,192,283,257,29,317,29,310,310,310,29,169,169,237,272,350,350,350,322,115,115,322,322,282,234,234,75,282,27,297,297,297,4,309,116,28,28,27,28,28,28,10,202,192,227,227,192,239,239,281,35,70,35,35,233,227,17,313,313,291,291,291,271,106,106,106,106,328,225,225,149,288,165,149,149,330,279,147,260,147,12,125,125,199,125,193,193,171,108,108,108,171,314,29,237,133,29,314,51,30,181,181,51,305,51,51,355,355,355,71,146,36,370,19,19,195,195,123,123,195,268,33,116,347,136,22,22,32,32,32,98,98,305,305,305,335,164,93,314,261,261,261,298,298,298,298,344,344,246,335,335,335,124,372,302,102,302,124,124,307,62,62,62,62,186,335,158,158,158,21,158,337,337,21,103,103,179,179,179,179,179,179,21,21,132,132,236,236,247,247,337,337,352,343,343,343,180,180,180,356,76,204,151,151,151,151,23,23,218,218,204,154,289,289,289,154,248,248,76,232,232,154,66,124,249,249,66,66,105,105,368,266,266,5,5,57,57,95,95,243,243,178,266,126,126,95,340,340,340,340,287,246,246,238,287,263,263,263,135,135,246,246,135,135,77,135,77,77,340,340,284,129,129,129,129,284,217,217,45,196,196,196,148,201,201,148};
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
