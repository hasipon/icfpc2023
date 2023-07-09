#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <set>
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

using grid = pair<int, int>;
using grid_score = pair<double, grid>;

vector<vector<grid_score>> gridRanking(const Problem &problem) {
    constexpr int stageGrid = 10;

    auto x_grids = int(ceil((problem.stageWidth - 20) / stageGrid));
    auto y_grids = int(ceil((problem.stageHeight - 20) / stageGrid));
    if (x_grids == 0) {
        x_grids = 1;
    }
    if (y_grids == 0) {
        y_grids = 1;
    }

    vector<vector<grid_score>> gridss(problem.attendees[0].tastes.size());
    for (auto &grids: gridss) {
        grids.reserve(x_grids * y_grids);
    }

    for (int x_grid = 0; x_grid < x_grids; x_grid++) {
        const int x = problem.stageLeft + 10 + x_grid * stageGrid;
        for (int y_grid = 0; y_grid < y_grids; y_grid++) {
            const int y = problem.stageBottom + 10 + y_grid * stageGrid;

            vector<double> score(problem.attendees[0].tastes.size());
            for (int aud = 0; aud < problem.attendees.size(); aud++) {
                const auto &audience = problem.attendees[aud];

                bool blocked = false;
                for (auto pr : problem.pillars) {
                    if (isBlocked(audience.x, audience.y, x, y, pr.x, pr.y, pr.r)) {
                        blocked = true;
                        break;
                    }
                }
                if (!blocked) {
                    const auto distance = (x - audience.x) * (x - audience.x) + (y - audience.y) * (y - audience.y);
                    for (int taste = 0; taste < audience.tastes.size(); taste++) {
                        score[taste] += audience.tastes[taste] / distance;
                    }
                }
            }
            for (int taste = 0; taste < problem.attendees[0].tastes.size(); taste++) {
                gridss[taste].emplace_back(score[taste], make_pair(x, y));
            }
        }
    }

    for (auto &grids:gridss) {
        sort(grids.begin(), grids.end(), greater<>{});
    }

    return gridss;
}

vector<pair<double, double>> solve(const Problem& problem, vector<vector<grid_score>> gridss) {
    // ビームの方針。
    // musician を全部配置しないといけない。taste ごとに考慮するべき順番はわかったがどれを置くべきかどうかはわかっていないのでそれを探索する
    // taste ごとに greedy におけるところを preference に基づいてきめる
    // ただ、既に置いてしまっている場合にはもう置けないことと、近くにおいている場合に死角になるために置かない方がいいという話がある
    // 置く時に上位3件ぐらい試してそれでスコアが高いやつ
    // ブロックするか確認することにする。

    vector<vector<int>> musicianByTaste(problem.attendees[0].tastes.size());
    for (int i = 0; i < problem.musicians.size(); i++) {
        musicianByTaste[problem.musicians[i]].push_back(i);
    }

    vector<pair<double, double>> res(problem.musicians.size());
    vector<int> progress(gridss.size());
    vector<int> used(gridss.size());

    set<pair<int, int>> placed;

    for (int i = 0; i < problem.musicians.size(); i++) {
        double score_max = -1e10;
        int score_max_taste = -1;
        pair<int, int> score_max_position;
        for (int j = 0; j < gridss.size(); j++) {
            auto use = used[j];
            if (musicianByTaste[j].size() <= use) {
                // no more musician for the taste.
                continue;
            }

            for (;; progress[j]++) {
                const auto [score, point] = gridss[j][progress[j]];
                if (placed.find(point) == placed.end()) {
                    if (score_max < score) {
                        score_max = score;
                        score_max_taste = j;
                        score_max_position = point;
                    }
                    break;
                }
            }
        }

        // Place taste(score_max_taste) to score_max_position.
        res[musicianByTaste[score_max_taste][used[score_max_taste]]] = make_pair(score_max_position.first, score_max_position.second);

        placed.insert(score_max_position);
        used[score_max_taste]++;
        progress[score_max_taste]++;
    }

    return res;
}

int main() {
    Problem problem;
    cin >> problem.roomWidth >> problem.roomHeight;
    cin >> problem.stageWidth >> problem.stageHeight;
    cin >> problem.stageLeft >> problem.stageBottom;
    int musicianN, tasteN, attendeeN, pillarN;
    cin >> musicianN >> tasteN;
    problem.musicians.resize(musicianN);
    for (auto &m: problem.musicians) {
        cin >> m;
    }
    cin >> attendeeN;
    problem.attendees.resize(attendeeN);
    for (auto &a: problem.attendees) {
        cin >> a.x >> a.y;
        a.tastes.resize(tasteN);
        for (auto &t: a.tastes) {
            cin >> t;
        }
    }
    cin >> pillarN;
    problem.pillars.resize(pillarN);
    for (auto &p: problem.pillars) {
        cin >> p.x >> p.y >> p.r;
    }

    auto gridss = gridRanking(problem);

    auto placement = solve(problem, gridss);

    cout << "{\"placements\":[";
    for (unsigned i = 0; i < placement.size(); i++) {
        if (i > 0) cout << ",";
        cout << "{\"x\":" << placement[i].first << ",\"y\":" << placement[i].second << "}";
    }
    cout << "]}" << endl;

    auto res = calcScore(problem, placement);
    if (!res.first) throw runtime_error("invalid placement");
    cerr << "score = " << res.second << endl;
}
