#define LOCAL_DEBUG 1
#define ENABLE_GV 1

#include "inada1.h"
#define GV_JS
#include "gv.hpp"
#include <limits.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <complex>
#include <set>
#include <map>
#include <random>
#include <filesystem>
#include <fstream>
#include <cassert>
#include <algorithm>

using namespace std;
namespace fs = std::filesystem;

#if LOCAL_DEBUG
std::mt19937 g_rand(12345);
#else
std::random_device rd;
std::mt19937 g_rand(rd());
#endif

#include <iostream>
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

// 問題仕様に沿ったスコア計算 (激重)
pair<bool, long long> calcScore(const Problem& problem, vector<pair<double, double>> placements) {
    if (placements.size() != problem.musicians.size()) {
        return { false, 0 };
    }
    for (unsigned i = 0; i < placements.size(); i++) {
        auto [x, y] = placements[i];
        if (!(
            problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
            problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10
            )) {
            return { false, 0 };
        }
        for (unsigned j = 0; j < i; j++) {
            auto [x2, y2] = placements[j];
            if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
                return { false, 0 };
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
    return { true,score };
}

void gvPlacements(const vector<pair<double, double> >& placements) {
    for (const auto& p : placements) {
        gvCircle(p.first, p.second, 5, gvRGB(255, 0, 0));
    }
}

// 外周に沿った placement を返却する
vector<pair<double, double>> makeEdgePlacement(const Problem& problem) {
    const auto x = problem.stageLeft + 10.0;
    const auto y = problem.stageBottom + 10.0;
    const auto w = problem.stageWidth - 20.0;
    const auto h = problem.stageHeight - 20.0;

    const auto speed = 10.0;
    auto dx = 0.0;
    auto dy = speed;
    auto left = x;
    auto right = x + w;
    auto top = y;
    auto bottom = y + h;
    auto cx = right;
    auto cy = y;

    vector<pair<double, double>> placements;
    for (int i = 0; i < problem.musicians.size(); i++) {
        placements.emplace_back(cx, cy);
        cx += dx;
        cy += dy;


        if (cy > bottom) {
            right -= 10.0;
            cy = bottom;
            cx = right;

            dx = -speed;
            dy = 0.0;
        }

        if (cx < left) {
            bottom -= 10.0;
            cx = left;
            cy = bottom;

            dy = -speed;
            dx = 0.0;
        }

        if (cy < top) {
            left += 10.0;
            cy = top;
            cx = left;

            dx = speed;
            dy = 0.0;
        }

        if (cx > right) {
            top += 10.0;
            cx = right;
            cy = top;

            dx = 0.0;
            dy = speed;
        }
    }

    return placements;
}

// 1つのTastesとAttendeesのみを考慮した簡易スコア計算
double calcScore2(const Problem& problem, int taste, double x, double y) {
    double score = 0;
    for (auto& a : problem.attendees) {
        auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
        score += a.tastes[taste] / d2;
    }
    return score;
}

// calcScore2 を使った山登りを行ってスコアだけ返す
double yamaScore(const Problem& problem, int taste, double x0, double y0) {
    const double D = 0.5;
    const double dx[4] = { +D, -D, 0, 0 };
    const double dy[4] = { 0, 0, +D, -D };
    double xx = x0;
    double yy = y0;
    double score = calcScore2(problem, taste, xx, yy);
    for (;;) {
        // --
        gvNewTime();
        gvRect(problem.stageLeft, problem.stageBottom, problem.stageWidth, problem.stageHeight, gvRGB(64, 64, 128));
        gvCircle(xx, yy, 5, gvRGB(255, 0, 0));
        gvOutput("Score: %.6lf", score);
        // --
        int k = -1;
        for (int i = 0; i < 4; ++i) {
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

// (x,y) にplacementを追加できるか (Stageの縁は考慮しない)
bool checkPlacements(double x, double y, const vector<pair<double, double>>& placements) {
    for (auto [x2, y2] : placements) {
        if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
            return false;
        }
    }
    return true;
}

pair<bool, long long> calcScoreWithCache(
    const Problem& problem,
    const vector<pair<double, double>>& placements,
	map<pair<double, double>, vector<double>>& cache,
    bool ignore_factor
    ) {

    if (cache.empty()) {
        if (placements.size() != problem.musicians.size()) {
            return { false, 0 };
        }
        for (unsigned i = 0; i < placements.size(); i++) {
            auto [x, y] = placements[i];
            if (!(
                problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
                problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10
                )) {
                return { false, 0 };
            }
            for (unsigned j = 0; j < i; j++) {
                auto [x2, y2] = placements[j];
                if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
                    return { false, 0 };
                }
            }
        }
    }

	vector<double> factor(placements.size(), 1);
	if (!ignore_factor && !problem.pillars.empty()) {
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

    if (cache.empty()) {
        for (unsigned i = 0; i < placements.size(); i++) {
            if (cache[placements[i]].empty()) {
                cache[placements[i]].resize(problem.attendees[0].tastes.size());
            }
        }
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

				for (unsigned t = 0; t < a.tastes.size(); t++)
				{
					const auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
					auto s = (long long)ceil(1000000 * a.tastes[t] / d2);
					cache[placements[i]][t] += s;
				}
			next:;
			}
		}
    }

    double score = 0;
	for (unsigned i = 0; i < placements.size(); i++) {
        if (cache.find(placements[i]) == cache.end()) {
            cerr << "Broken cache" << endl;
            throw 1;
        }
        auto s = cache[placements[i]][problem.musicians[i]];
		if (factor[i] > 1) {
			s = (long long)ceil(factor[i] * s); // NOLINT(cppcoreguidelines-narrowing-conversions)
		}
        score += s;
	}

    return { true,score };
}

long long swapDeltaScore(
    const Problem& problem,
    const vector<pair<double, double>>& placements,
    map<pair<double, double>, vector<double>>& cache,
    bool ignore_factor,
    int swap_i, int swap_j
) {
    assert(!cache.empty());
    assert(cache.size() == placements.size());
    if (!ignore_factor && !problem.pillars.empty()) {
        throw 1;
    }

    double score = 0;
	score -= cache[placements[swap_i]][problem.musicians[swap_i]];
	score -= cache[placements[swap_j]][problem.musicians[swap_j]];

	score += cache[placements[swap_i]][problem.musicians[swap_j]];
	score += cache[placements[swap_j]][problem.musicians[swap_i]];
    return score;
}

vector<pair<double, double>> solve(const Problem& problem) {
    map<pair<double, double>, vector<double>> cache;

    auto placements = makeEdgePlacement(problem);
    auto res = calcScore(problem, placements);
    cerr << "calcScore: " << res.second << endl;
    res = calcScoreWithCache(problem, placements, cache, false);
    cerr << "calcScoreWithCache1: " << res.second << endl;

    shuffle(placements.begin(), placements.end(), g_rand);

    res = calcScore(problem, placements);
    cerr << "calcScore: " << res.second << endl;
    res = calcScoreWithCache(problem, placements, cache, false);
    cerr << "calcScoreWithCache2: " << res.second << endl;

    auto score = res.second;
    for (int ite = 0; ite < 10; ite++) {
		for (int i = 0; i < placements.size(); i++) {
			for (int j = i + 1; j < placements.size(); j++) {
                auto ds = swapDeltaScore(problem, placements, cache, false, i, j);
                if (0 < ds) {
                    score += ds;
					std::swap(placements[i], placements[j]);
                }
			}
		}
        cerr << "score:" << score << endl;
    }

    return placements;
}

void readProblem(std::istream& is, Problem& problem) {
    is >> problem.roomWidth >> problem.roomHeight;
    is >> problem.stageWidth >> problem.stageHeight;
    is >> problem.stageLeft >> problem.stageBottom;
    int musicianN, tasteN, attendeeN, pillarN;
    is >> musicianN >> tasteN;
    problem.musicians.resize(musicianN);
    for (auto& m : problem.musicians) {
        is >> m;
    }
    is >> attendeeN;
    problem.attendees.resize(attendeeN);
    for (auto& a : problem.attendees) {
        is >> a.x >> a.y;
        a.tastes.resize(tasteN);
        for (auto& t : a.tastes) {
            is >> t;
        }
    }
    is >> pillarN;
    problem.pillars.resize(pillarN);
    for (auto& p : problem.pillars) {
        is >> p.x >> p.y >> p.r;
    }
}

void writePlacementsJSON(std::ostream& os, const std::vector<pair<double, double> >& placements) {
    cout << "{\"placements\":[";
    for (unsigned i = 0; i < placements.size(); i++) {
        if (i > 0) cout << ",";
        cout << "{\"x\":" << placements[i].first << ",\"y\":" << placements[i].second << "}";
    }
    cout << "]}" << endl;
}

int main(int argc, char* argv[]) {
#if LOCAL_DEBUG
    fs::current_path(R"(c:\projects\hasipon\icfpc2023\solver\inada)");

    int problem_id = 1;
    Problem problem;
    {
        ifstream ifs("../../problems.kyopro/" + to_string(problem_id) + ".kyopro");
		readProblem(ifs, problem);
    }
#else
    Problem problem;
    readProblem(cin, problem);
#endif

    cerr << "pillar:" << problem.pillars.size() << endl;

    auto placement = solve(problem);
    gvPlacements(placement);

    auto res = calcScore(problem, placement);
    writePlacementsJSON(cout, placement);
#if LOCAL_DEBUG
    {
		ofstream ofs(to_string(problem_id) + "-inada1-" + to_string(res.second) + ".json");
        writePlacementsJSON(ofs, placement);
    }
#endif
    if (!res.first) throw runtime_error("invalid placement");
    cerr << "score = " << res.second << endl;
}
