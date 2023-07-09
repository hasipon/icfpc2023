#include <iostream>
#include <vector>
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

                bool blocked = false;
                for (auto pr : problem.pillars) {
                    if (isBlocked(audience.x, audience.y, x, y, pr.x, pr.y, pr.r)) {
                        blocked = true;
                        break;
                    }
                }
                if (!blocked) {
                    const auto distance = (x - audience.x) * (x - audience.x) + (y - audience.y) * (y - audience.y);
                    score += attendeesContributions[aud] / distance;
                }
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

    const auto grids = gridRanking(problem);
    for (const auto &grid: grids) {
        cout << grid.first << "," << grid.second.first << "," << grid.second.second << endl;
    }
}
