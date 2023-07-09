#include <iostream>
#include <vector>
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

/*
bool isBlocked(double x0, double y0, double x1, double y1, double x2, double y2) {
    complex<double> p0(x0, y0), p1(x1, y1), p2(x2, y2);
    if (real(conj(p1 - p0) * (p2 - p0)) < 0) return false;
    if (real(conj(p0 - p1) * (p2 - p1)) < 0) return false;
    double t = real(conj(p2 - p0) * (p0 - p1)) / norm(p0 - p1);
    return abs(p2 - (p0 + (p0 - p1) * t)) < 5;
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

double calcScore2(const Problem& problem, int taste, double x, double y) {
    double score = 0;
    for (auto& a : problem.attendees) {
        auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
        score += a.tastes[taste] / d2;
    }
    return score;
}

double yamaScore(const Problem& problem, int taste, double x0, double y0) {
    const double D = 0.5;
    const double dx[4] = {+D, -D, 0, 0};
    const double dy[4] = {0, 0, +D, -D};
    double xx = x0;
    double yy = y0;
    double score = calcScore2(problem, taste, xx, yy);
    for (;;) {
        int k = -1;
        for (int i = 0; i < 4; ++ i) {
            double x = xx + dx[i];
            double y = yy + dy[i];
            if (
                    problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
                    problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10
                    ) {
                double s = calcScore2(problem, taste, x, y);
                if (s > score) {
                    score = s;
                    k = i;
                }
            }
        }
        if (k == -1) break;
        xx += dx[k];
        yy += dy[k];
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

pair<int, int> makeStart(const Problem& problem, const vector<pair<double, double>>& placements) {
    for (;;) {
        double x0 = problem.stageLeft + 10 + rand() % ((int)problem.stageWidth - 19);
        double y0 = problem.stageBottom + 10 + rand() % ((int)problem.stageHeight - 19);
        if (checkPlacements(x0, y0, placements)) return {x0, y0};
    }
}

double calcScore3(const Problem& problem, const vector<unsigned>& perm, const vector<pair<double, double>>& placements, const set<pair<unsigned, unsigned>>& blocking, int taste, double x, double y) {
    double score = 0;
    for (unsigned k = 0; k < problem.attendees.size(); ++ k) {
        auto& a = problem.attendees[k];
        bool blocked = false;
        for (unsigned jj = 0; jj < placements.size(); ++ jj) {
            auto j = perm[jj];
            auto [x1, y1] = placements[j];
            if (!blocked && isBlocked(a.x, a.y, x, y, x1, y1)) {
                blocked = true;
            }
            if (!blocking.count({j, k}) && isBlocked(a.x, a.y, x1, y1, x, y)) {
                auto d2 = (a.x - x1) * (a.x - x1) + (a.y - y1) * (a.y - y1);
                score -= a.tastes[problem.musicians[j]] / d2;
            }
        }
        if (!blocked) {
            auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
            score += a.tastes[taste] / d2;
        }
    }
    return score;
}

pair<double, double> yama(const Problem& problem, int taste, pair<int, int> start, const vector<pair<double, double>>& placements) {
    const double D = 1;
    const double dx[4] = {+D, -D, 0, 0};
    const double dy[4] = {0, 0, +D, -D};
    double xx = start.first;
    double yy = start.second;
    double score = calcScore2(problem, taste, xx, yy);
    for (;;) {
        int k = -1;
        for (int i = 0; i < 4; ++ i) {
            double x = xx + dx[i];
            double y = yy + dy[i];
            if (
                    problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
                    problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10
                    ) {
                if (!checkPlacements(x, y, placements)) continue;
                double s = calcScore2(problem, taste, x, y);
                if (s > score) {
                    score = s;
                    k = i;
                }
            }
        }
        if (k == -1) break;
        xx += dx[k];
        yy += dy[k];
    }
    return {xx, yy};
}

pair<double, double> yama3(const Problem& problem, int taste, double x0, double y0, double score, const vector<pair<double, double>>& placements, const vector<unsigned>& perm, const set<pair<unsigned, unsigned>>& blocking) {
    const double D = 1;
    const double dx[4] = {+D, -D, 0, 0};
    const double dy[4] = {0, 0, +D, -D};
    double xx = x0;
    double yy = y0;
    for (;;) {
        int k = -1;
        for (int i = 0; i < 4; ++ i) {
            double x = xx + dx[i];
            double y = yy + dy[i];
            if (problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 && problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10) {
                if (!checkPlacements(x, y, placements)) continue;
                double s = calcScore3(problem, perm, placements, blocking, taste, x, y);
                if (s > score) {
                    score = s;
                    k = i;
                }
            }
        }
        if (k == -1) break;
        xx += dx[k];
        yy += dy[k];
    }
    return {xx, yy};
}

vector<double> calcWeightTaste(const Problem& problem){
    int tasteN = problem.attendees[0].tastes.size();
    vector<double> weight(problem.attendees[0].tastes.size(), 0);
    for(auto a : problem.attendees){
        double diffX = min(abs(a.x - problem.stageLeft), abs(a.x - problem.stageLeft + problem.stageWidth));
        double diffY = min(abs(a.y - problem.stageBottom), abs(a.y - problem.stageBottom + problem.stageHeight));
        for(int i = 0; i<a.tastes.size(); i++){
            weight[i] += (a.tastes[i] * 1000000.) / (diffX * diffX + diffY * diffY);
        }
    }
    return weight;
}

vector<unsigned > calcPerm(const Problem &p){
    auto weightTaste = calcWeightTaste(p);
    vector<pair<pair<double, double>, int> > sortV;
    for(int i = 0; i<p.musicians.size(); i++){
        auto m = p.musicians[i];
        auto s = yamaScore(p, p.musicians[i], p.stageLeft + p.stageWidth / 2, p.stageBottom + p.stageHeight / 2);
        sortV.emplace_back(make_pair(make_pair(weightTaste[m], s), i));
    }
    sort(sortV.begin(), sortV.end());
    vector<unsigned > perm;
    for(auto v : sortV){
        perm.emplace_back(v.second);
    }
    return perm;
}

vector<pair<double, double>> solve(const Problem& problem) {
    unsigned N = problem.musicians.size();
    vector<pair<double, double>> res(N);
    vector<unsigned > perm = calcPerm(problem);

    vector<pair<double, double>> placements;
    set<pair<unsigned, unsigned>> blocking;
    for (auto i : perm) {
        int taste = problem.musicians[i];
        auto best = yama(problem, problem.musicians[i], makeStart(problem, placements), placements);
        double bestScore = calcScore3(problem, perm, placements, blocking, taste, best.first, best.second);
        for (int tt = 0; tt < 5; ++ tt) {
            auto pos = yama(problem, problem.musicians[i], makeStart(problem, placements), placements);
            auto s = calcScore3(problem, perm, placements, blocking, taste, pos.first, pos.second);
            if (s > bestScore) {
                bestScore = s;
                best = pos;
            }
        }
        auto [x2, y2] = yama3(problem, taste, best.first, best.second, bestScore, placements, perm, blocking);
        res[i] = {x2, y2};
        placements.push_back(res[i]);
        for (unsigned k = 0; k < problem.attendees.size(); ++ k) {
            auto& a = problem.attendees[k];
            for (unsigned jj = 0; jj < placements.size()-1; jj++) {
                auto j = perm[jj];
                if (blocking.count({j, k})) continue;
                auto [x, y] = placements[jj];
                if (isBlocked(a.x, a.y, x, y, x2, y2)) {
                    blocking.insert({j, k});
                }
            }
        }
    }
    return res;
}
*/

using grid = pair<int, int>;
using grid_score = pair<double, grid>;

vector<grid_score> gridRanking(const Problem &problem) {
    constexpr int stageGrid = 10;

    vector<int> tasteCount(problem.attendees[0].tastes.size());
    for (auto taste: problem.musicians) {
        tasteCount[taste]++;
    }

    vector<double> attendeesContributions(problem.attendees.size());
    for (int aud = 0; aud < problem.attendees.size(); aud++) {
        const auto &audience = problem.attendees[aud];
        for (int i = 0; i < audience.tastes.size(); i++) {
            attendeesContributions[aud] = audience.tastes[i] * tasteCount[i];
        }
    }

    const auto x_grids = int(ceil((problem.stageWidth - 20) / stageGrid));
    const auto y_grids = int(ceil((problem.stageHeight - 20) / stageGrid));
    vector<grid_score> grids;
    grids.reserve(x_grids * y_grids);

    for (int x_grid = 0; x_grid < x_grids; x_grid++) {
        for (int y_grid = 0; y_grid < y_grids; y_grid++) {
            const int x = problem.stageLeft + 10 + x_grid * stageGrid;
            const int y = problem.stageBottom + 10 + y_grid * stageGrid;
            double score = 0;
            for (int aud = 0; aud < problem.attendees.size(); aud++) {
                const auto &audience = problem.attendees[aud];
                const auto distance = (x - audience.x) * (x - audience.x) + (y - audience.y) * (y - audience.y);
                score += attendeesContributions[aud] / distance;
            }
            grids.emplace_back(score, make_pair(x, y));
        }
    }

    sort(grids.begin(), grids.end(), greater<>{});

    return grids;
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

    const auto grids = gridRanking(problem);
    for (const auto &grid: grids) {
        //cout << "score: " << grid.first << "(" << grid.second.first << "," << grid.second.second << ")" << endl;
        cout << grid.first << "," << grid.second.first << "," << grid.second.second << endl;
    }

    /*

    auto placement = solve(problem);
    cout << "{\"placements\":[";
    for (unsigned i = 0; i < placement.size(); i++) {
        if (i > 0) cout << ",";
        cout << "{\"x\":" << placement[i].first << ",\"y\":" << placement[i].second << "}";
    }
    cout << "]}" << endl;

    auto res = calcScore(problem, placement);
    if (!res.first) throw runtime_error("invalid placement");
    cerr << "score = " << res.second << endl;
     */
}
