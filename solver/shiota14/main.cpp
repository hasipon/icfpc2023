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

pair<bool, long long> calcScore(const Problem& problem, const vector<pair<int, int>> &placements) {
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
                    auto d2 = double(x1 - x2) * double(x1 - x2) + double(y1 - y2) * double(y1 - y2);
                    factor[i] += 1 / sqrt(d2);
                }
            }
        }
    }
    long long score = 0;
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

constexpr auto GRID = 50;
using grid_index_t = vector<vector<vector<int>>>;

bool checkBlocked(const Problem &problem, const grid_index_t &grid_index, const Attendee &a, int x1, int y1, int i, const vector<pair<int, int>> &placements) {
    auto x_grids = int(problem.roomWidth / GRID);
    auto y_grids = int(problem.roomHeight / GRID);

    complex<double> aud(a.x, a.y);
    //cout << "aud" << aud << endl;
    complex<double> k(x1, y1);
    //cout << "k" << k << endl;
    auto ak = k - aud;
    complex<double> rotation(0, -1);
    auto normal = ak * rotation;
    normal /= abs(normal);
    //cout << "normal(" << normal.real() << "," << normal.imag() << ")" << endl;

    auto left_top = k - 5.*normal;
    //cout << "left_top" << left_top << endl;
    auto right_top = k + 5.*normal;
    //cout << "right_top" << right_top << endl;
    auto left_bottom = aud - 5.*normal;
    //cout << "left_bottom" << left_bottom << endl;
    auto right_bottom = aud + 5.*normal;
    //cout << "right_bottom" << right_bottom << endl;

    auto incl_left_bottom = complex<double>(
            min({left_top.real(), right_top.real(), left_bottom.real(), right_bottom.real()}),
            min({left_top.imag(), right_top.imag(), left_bottom.imag(), right_bottom.imag()})
    );
    if (incl_left_bottom.real() < problem.stageLeft) {
        incl_left_bottom.real(problem.stageLeft);
    }
    if (incl_left_bottom.imag() < problem.stageBottom) {
        incl_left_bottom.imag(problem.stageBottom);
    }
    //cout << "incl_left_bottom" << incl_left_bottom << endl;

    auto incl_right_top = complex<double>(
            max({left_top.real(), right_top.real(), left_bottom.real(), right_bottom.real()}),
            max({left_top.imag(), right_top.imag(), left_bottom.imag(), right_bottom.imag()})
    );
    if (problem.stageLeft + problem.stageWidth < incl_right_top.real()) {
        incl_right_top.real(problem.stageLeft + problem.stageWidth);
    }
    if (problem.stageBottom + problem.stageHeight < incl_right_top.imag()) {
        incl_right_top.imag(problem.stageBottom + problem.stageHeight);
    }
    //cout << "incl_right_top" << incl_right_top << endl;

    // TODO: Reduce the number of grids to check by consider rotation.
    for (int x_grid = incl_left_bottom.real(); x_grid <= incl_right_top.real() + GRID; x_grid += GRID) {
        for (int y_grid = incl_left_bottom.imag(); y_grid <= incl_right_top.imag() + GRID; y_grid += GRID) {
            auto x_grid_id = x_grid / GRID;
            auto y_grid_id = y_grid / GRID;
            if (0 <= x_grid_id && x_grid_id < x_grids && 0 <= y_grid_id && y_grid_id < y_grids) {
                //cout << "grid_id : " << x_grid_id << " : " << y_grid_id << endl;

                for (auto j : grid_index[x_grid_id][y_grid_id]) {
                    if (i != j) {
                        auto [x2, y2] = placements[j];
                        if (isBlocked(a.x, a.y, x1, y1, x2, y2, 5)) {
                            return true;
                        }
                    }
                }
            }
        }
    }

    return false;
}

long long calcScoreForOneTaste(
        const Problem& problem,
        const vector<vector<int>> &musicianByTaste,
        const vector<pair<int, int>> &placements,
        int taste,
        const grid_index_t &grid_index
) {
    vector<double> factor(placements.size(), 1);
    if (!problem.pillars.empty()) {
        for (auto i : musicianByTaste[taste]) {
            const auto &[x1, y1] = placements[i];
            if (x1 == -1 && y1 == -1) {
                // Since musician is determined from first so not determined means break immediately.
                break;
            }

            for (auto j : musicianByTaste[taste]) {
                if (i == j) {
                    continue;
                }

                const auto &[x2, y2] = placements[j];
                if (x2 == -1 && y2 == -1) {
                    // Since musician is determined from first so not determined means break immediately.
                    break;
                }

                auto d2 = double(x1 - x2) * double(x1 - x2) + double(y1 - y2) * double(y1 - y2);
                factor[i] += 1 / sqrt(d2);
            }
        }
    }

    long long score = 0;
    for (auto& a : problem.attendees) {
        for (auto i: musicianByTaste[taste]) {
            const auto &[x1, y1] = placements[i];
            if (x1 == -1 && y1 == -1) {
                // Since musician is determined from first so not determined means break immediately.
                break;
            }

            if (checkBlocked(problem, grid_index, a, x1, y1, i, placements)) {
                goto next;
            }

            for (auto pr: problem.pillars) {
                if (isBlocked(a.x, a.y, x1, y1, pr.x, pr.y, pr.r)) {
                    goto next;
                }
            }
            {
                auto d2 = (a.x - x1) * (a.x - x1) + (a.y - y1) * (a.y - y1);
                auto s = (long long) ceil(1000000 * a.tastes[taste] / d2);
                if (factor[i] > 1) {
                    s = (long long) ceil(factor[i] * s); // NOLINT(cppcoreguidelines-narrowing-conversions)
                }
                score += s; // NOLINT(cppcoreguidelines-narrowing-conversions)
            }
            next:;
        }
    }
    return score;
}

using grid = pair<int, int>;
using grid_score = pair<double, grid>;

vector<vector<grid_score>> gridRanking(const Problem &problem, const vector<grid> &placements, const grid_index_t &grid_index) {
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
        int x = problem.stageLeft + 10 + x_grid * stageGrid;
        if(x > problem.stageLeft + problem.stageWidth / 2){
            x += int(problem.stageWidth ) % 10;
        }
        for (int y_grid = 0; y_grid < y_grids; y_grid++) {
             int y = problem.stageBottom + 10 + y_grid * stageGrid;
            if(y > problem.stageBottom + problem.stageHeight / 2){
                y += int(problem.stageHeight ) % 10;
            }

            vector<double> score(problem.attendees[0].tastes.size());
            for (int aud = 0; aud < problem.attendees.size(); aud++) {
                const auto &audience = problem.attendees[aud];

                if (checkBlocked(problem, grid_index, audience, x, y, -1, placements)) {
                    goto next;
                }

                for (auto pr : problem.pillars) {
                    if (isBlocked(audience.x, audience.y, x, y, pr.x, pr.y, pr.r)) {
                        goto next;
                    }
                }
                {
                    const auto distance = (x - audience.x) * (x - audience.x) + (y - audience.y) * (y - audience.y);
                    for (int taste = 0; taste < audience.tastes.size(); taste++) {
                        score[taste] += audience.tastes[taste] / distance;
                    }
                }
                next:;
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


vector<pair<int, int>> solve(const Problem& problem) {
    vector<vector<int>> musicianByTaste(problem.attendees[0].tastes.size());
    for (int i = 0; i < problem.musicians.size(); i++) {
        musicianByTaste[problem.musicians[i]].push_back(i);
    }

    auto x_grids = int(problem.roomWidth / GRID);
    auto y_grids = int(problem.roomHeight / GRID);

    vector<vector<vector<int>>> grid_index(x_grids, vector<vector<int>>(y_grids));
    vector<pair<int, int>> res(problem.musicians.size(), make_pair(-1, -1));
    auto gridss = gridRanking(problem, res, grid_index);

    vector<int> progress(gridss.size());
    vector<int> used(gridss.size());
    vector<long long> last_score(gridss.size());

    set<pair<int, int>> placed;

    for (int i = 0; i < problem.musicians.size(); i++) {
        /*
        if (i != 0 && i % 10 == 0) {
            // update heat map.
            gridss = gridRanking(problem, res, grid_index);
            progress = used;
        }
         */

        // TODO: 未実装。これらの盤面を評価して、B件だけとっておく。
        // constexpr int B = 3;

        vector<pair<double, int>> taste_candidates(gridss.size());
        for (int j = 0; j < gridss.size(); j++) {
            if (musicianByTaste[j].size() <= used[j]) {
                // no more musician for the taste.
                taste_candidates[j] = make_pair(-1e10, j);
                continue;
            }

            for (;; progress[j]++) {
                const auto [score, point] = gridss[j][progress[j]];
                if (placed.find(point) == placed.end()) {
                    taste_candidates[j] = make_pair(score, j);
                    goto next;
                }
            }
            next:;
        }

        sort(taste_candidates.begin(), taste_candidates.end(), greater<>());

        // taste は heatmap からくる良さそうなものを test_trials 個 試す
        //constexpr int taste_trials = 50;
        using score_taste_id_progress_id = tuple<double, double, int, int>;
        vector<score_taste_id_progress_id> candidates;

        for (int taste_trial = 0; taste_trial < taste_candidates.size(); taste_trial++) {
            const auto [taste_score, j] = taste_candidates[taste_trial];
            if (taste_score == -1e10) {
                continue;
            }

            // calcScoreForOneTaste は ある taste に関する貢献しか計算しないので、差分で ranking して top1 以外の taste が選ばれやすくする
            auto previous_score = last_score[j];

            // ある taste について 上位M件の配置を試す
            constexpr int M = 1;

            int found_count = 0;
            for (int jj = progress[j]; jj < gridss[j].size() && found_count < M; jj++) {
                const auto &[s, point] = gridss[j][jj];
                // cout << "FOUND " << j << " to " << point.first << "," << point.second << " score " << s << endl;
                if (placed.find(point) == placed.end()) {
                    res[musicianByTaste[j][used[j]]] = point;

                    {
                        // update grid_index
                        auto x_in_grid = int(point.first / GRID);
                        auto y_in_grid = int(point.second / GRID);

                        grid_index[x_in_grid][y_in_grid].push_back(musicianByTaste[j][used[j]]);
                    }

                    auto score = calcScoreForOneTaste(problem, musicianByTaste, res, j, grid_index);

                    {
                        // restore grid_index
                        auto x_in_grid = int(point.first / GRID);
                        auto y_in_grid = int(point.second / GRID);

                        grid_index[x_in_grid][y_in_grid].pop_back();
                    }

                    candidates.emplace_back(score - previous_score, score, j, jj);
                    res[musicianByTaste[j][used[j]]] = make_pair(-1, -1);
                    found_count++;
                }
            }
        }

        const auto& [score_diff, score, taste_id, progress_id]  = *max_element(candidates.begin(), candidates.end());
        const auto& [_, point] = gridss[taste_id][progress_id];
        res[musicianByTaste[taste_id][used[taste_id]]] = point;
        placed.insert(point);
        used[taste_id]++;
        progress[taste_id] = progress_id + 1;
        last_score[taste_id] = score;

        {
            auto x_in_grid = int(point.first / GRID);
            auto y_in_grid = int(point.second / GRID);

            grid_index[x_in_grid][y_in_grid].push_back(i);
        }
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

    auto placement = solve(problem);

    cout << "{\"placements\":[";
    for (unsigned i = 0; i < placement.size(); i++) {
        if (i > 0) cout << ",";
        cout << "{\"x\":" << placement[i].first << ",\"y\":" << placement[i].second << "}";
    }
    cout << "],";
    cout << "\"volumes\":[";
    for (unsigned i = 0; i < placement.size(); i++) {
        if (i > 0) cout << ",";
        cout << "10";
    }
    cout << "]}" << endl;

    auto res = calcScore(problem, placement);
    if (!res.first) throw runtime_error("invalid placement");
    cerr << "score = " << res.second << endl;
}
