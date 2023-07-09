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

    vector<int> tastes {103,466,199,479,291,110,347,89,326,333,88,251,548,300,156,420,165,340,42,408,541,100,246,269,552,497,281,317,18,521,173,368,195,410,243,416,531,73,405,241,353,396,186,426,90,420,93,349,451,201,417,119,91,224,7,63,448,216,17,367,480,373,222,464,148,77,512,303,31,273,299,81,389,359,152,11,311,239,85,113,304,56,378,390,127,50,245,215,486,240,146,205,476,264,252,271,316,32,292,279,235,260,181,297,407,198,406,357,343,471,139,133,176,136,414,82,546,358,218,255,435,159,337,135,440,1,70,444,556,459,58,431,44,432,45,314,60,204,342,232,76,425,319,94,59,155,308,310,47,192,277,509,57,470,262,140,182,469,305,399,522,12,345,499,516,415,551,543,441,30,75,151,131,402,15,26,481,329,398,83,187,147,234,461,321,473,508,61,233,456,130,259,539,214,433,542,203,352,527,533,104,9,79,363,379,434,419,328,350,554,49,294,534,189,385,474,154,2,518,109,194,295,266,453,8,180,261,377,538,290,24,362,157,37,62,244,324,92,361,365,383,141,320,540,33,196,500,507,66,344,411,395,475,334,446,153,442,87,123,126,537,412,166,338,467,501,219,376,238,510,97,270,423,478,558,28,125,355,163,52,366,529,226,191,514,391,5,183,437,306,96,375,247,462,504,65,485,175,447,179,10,101,112,525,220,212,275,185,284,465,107,172,488,108,138,178,322,487,121,211,170,388,236,207,296,134,221,424,524,351,150,490,536,137,64,0,452,553,400,228,202,339,346,515,248,454,454,439,114,449,495,394,413,280,356,404,116,458,313,168,256,23,332,225,393,430,129,429,38,382,206,327,197,161,526,84,48,315,210,268,257,3,177,335,520,387,184,421,249,428,80,164,164,74,14,397,25,494,143,285,188,86,549,436,511,418,242,544};
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
