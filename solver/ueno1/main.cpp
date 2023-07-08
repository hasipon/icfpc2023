#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <random>
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

bool isBlocked(double x0, double y0, double x1, double y1, double x2, double y2) {
    complex<double> p0(x0, y0), p1(x1, y1), p2(x2, y2);
    if (real(conj(p1 - p0) * (p2 - p0)) < 0) return false;
    if (real(conj(p0 - p1) * (p2 - p1)) < 0) return false;
    double t = real(conj(p2 - p0) * (p0 - p1)) / norm(p0 - p1);
    return abs(p2 - (p0 + (p0 - p1) * t)) <= 5;
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

    constexpr auto grid = 50;
    auto x_grids = int(problem.roomWidth / grid);
    auto y_grids = int(problem.roomHeight / grid);

    // cout << "x_grids : " << x_grids << endl;
    // cout << "y_grids : " << y_grids << endl;

    vector<vector<vector<int>>> grid_index(x_grids, vector<vector<int>>(y_grids));
    for (auto i = 0; i < placements.size(); i++) {
        auto [x, y] = placements[i];

        auto x_in_grid = int(x / grid);
        auto y_in_grid = int(y / grid);

        grid_index[x_in_grid][y_in_grid].push_back(i);
    }

    double score = 0;
    for (auto& a : problem.attendees) {
        for (unsigned i = 0; i < placements.size(); i++) {
            complex<double> aud(a.x, a.y);
            //cout << "aud" << aud << endl;
            auto [x, y] = placements[i];
            complex<double> k(x, y);
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
            for (int x_grid = incl_left_bottom.real(); x_grid <= incl_right_top.real() + grid; x_grid += grid) {
                for (int y_grid = incl_left_bottom.imag(); y_grid <= incl_right_top.imag() + grid; y_grid += grid) {
                    auto x_grid_id = x_grid / grid;
                    auto y_grid_id = y_grid / grid;
                    if (0 <= x_grid_id && x_grid_id < x_grids && 0 <= y_grid_id && y_grid_id < y_grids) {
                        //cout << "grid_id : " << x_grid_id << " : " << y_grid_id << endl;

                        for (auto j : grid_index[x_grid_id][y_grid_id]) {
                            if (i != j) {
                                auto [x2, y2] = placements[j];
                                if (isBlocked(a.x, a.y, x, y, x2, y2)) {
                                    goto next;
                                }
                            }
                        }
                    }
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

bool checkPlacements(double x, double y, const vector<pair<double, double>>& placements) {
    for (auto [x2, y2] : placements) {
        if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
            return false;
        }
    }
    return true;
}

pair<double, double> yama(const Problem& problem, int taste, double x0, double y0, const vector<pair<double, double>>& placements) {
    const int dx[4] = {+1, -1, 0, 0};
    const int dy[4] = {0, 0, +1, -1};
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

pair<int, int> makeStart(const Problem& problem, const vector<pair<double, double>>& placements) {
    for (;;) {
        double x0 = problem.stageLeft + 10 + rand() % ((int)problem.stageWidth - 19);
        double y0 = problem.stageBottom + 10 + rand() % ((int)problem.stageHeight - 19);
        if (checkPlacements(x0, y0, placements)) return {x0, y0};
    }
}

vector<pair<double, double>> solve(const Problem& problem) {
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    vector<pair<double, double>> res(problem.musicians.size());
    vector<unsigned> perm(problem.musicians.size());
    for (unsigned i = 0; i < perm.size(); ++ i) perm[i] = i;
    shuffle(perm.begin(), perm.end(), engine);
    vector<pair<double, double>> placements;
    for (auto i : perm) {
        auto [x0, y0] = makeStart(problem, placements);
        res[i] = yama(problem, problem.musicians[i], x0, y0, placements);
        placements.push_back(res[i]);
    }
    return res;
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
}
