#define LOCAL_DEBUG 1
#define ENABLE_GV 1

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

pair<bool, vector<long long>> calcScore(const Problem& problem, const vector<pair<double, double>>& placements) {
    if (placements.size() != problem.musicians.size()) {
        return { false, {} };
    }
    for (unsigned i = 0; i < placements.size(); i++) {
        auto [x, y] = placements[i];
        if (!(
            problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
            problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10
            )) {
            return { false, {} };
        }
        for (unsigned j = 0; j < i; j++) {
            auto [x2, y2] = placements[j];
            if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
                return { false, {} };
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
    vector<long long> score(placements.size());
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
                auto iik = (long long)ceil(1000000 * a.tastes[problem.musicians[i]] / d2);
                if (factor[i] > 1) {
                    score[i] += (long long)ceil(10 * factor[i] * iik); // NOLINT(cppcoreguidelines-narrowing-conversions)
                }
                else {
                    score[i] += 10 * iik;
                }
            }
        next:;
        }
    }
    return { true,score };
}

long long sumScore(const vector<long long>& score) {
    long long sum = 0;
    for (auto s : score) {
        if (s > 0) {
            sum += s;
        }
    }
    return sum;
}

void gvStage(const Problem& problem) {
    gvRect(problem.stageLeft , problem.stageBottom, problem.stageWidth, problem.stageHeight, gvRGB(128,128,0));
}

void gvPlacements(const vector<pair<double, double> >& placements) {
    for (const auto& p : placements) {
        gvCircle(p.first, p.second, 5, gvRGB(255, 0, 0));
    }
}

void gvCandidates(const multimap<double, pair<double, double> >& candidates) {
    for (const auto& p : candidates) {
        gvCircle(p.second.first, p.second.second, 5, gvRGB(255, 0, 0));
        gvText(p.second.first, p.second.second, 3, "%.lf", p.first);
    }
}

void gvHeatmap(const map<pair<double, double>, double>& heatmap) {
    for (const auto& p : heatmap) {
        gvCircle(p.first.first, p.first.second, 1, gvRGB(0, 255, 0));
    }
}

pair<bool, long long> calcScoreWithCache(
    const Problem& problem,
    const vector<pair<double, double>>& placements,
    map<pair<double, double>, vector<long long>>& cache,
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

    long long score = 0;
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
    map<pair<double, double>, vector<long long>>& cache,
    bool ignore_factor,
    int swap_i, int swap_j
) {
    assert(!cache.empty());
    assert(cache.size() == placements.size());
    if (!ignore_factor && !problem.pillars.empty()) {
        throw 1;
    }

    long long score = 0;
    score -= cache[placements[swap_i]][problem.musicians[swap_i]];
    score -= cache[placements[swap_j]][problem.musicians[swap_j]];

    score += cache[placements[swap_i]][problem.musicians[swap_j]];
    score += cache[placements[swap_j]][problem.musicians[swap_i]];
    return score;
}

vector<pair<double, double>> solve(const Problem& problem, const map<pair<double, double>, double>& heatmap) {
    // step1. 外周に沿って配置, スコアに影響のないものは削除
	map<pair<double, double>, int> p_tastes;
	map<pair<double, double>, long long> p_score;
	map<int, int> free_tastes;
    {
        multimap<double, pair<double, double>> candidates;
        auto add_candidate = [&](double x, double y) {
            if (!(problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
                problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10
                )) {
                cerr << "Reject" << x << " " << y << endl;
                return;
            }
            auto rx = std::round(x / 10) * 10;
            auto ry = std::round(y / 10) * 10;
            double ddx[] = { 0, 10, -10, 0, 0, 10, 10, -10, -10 };
            double ddy[] = { 0, 0, 0, 10, -10, 10, -10, 10, -10 };
            for (int i = 0; i < 9; i++) {
                auto it = heatmap.find(make_pair(rx + ddx[i], ry + ddy[i]));
                if (it != heatmap.end()) {
                    candidates.emplace(it->second, make_pair(x, y));
                    return;
                }
            }
            cerr << " not found " << x << " " << y << endl;
        };

        const auto x = problem.stageLeft + 10.0;
        const auto y = problem.stageBottom + 10.0;
        const auto w = problem.stageWidth - 20.0;
        const auto h = problem.stageHeight - 20.0;

        const auto speed = 10.0;
        const auto left = x;
        const auto right = x + w;
        const auto top = y;
        const auto bottom = y + h;

        auto dx = 0.0;
        auto dy = speed;
        auto cx = right;
        auto cy = y;
        double ex = long long(w) % 10;
        double ey = long long(h) % 10;

        for (;;) {
            add_candidate(cx, cy);
            if (0 < dy) add_candidate(cx, cy + ey);
            if (0 > dy) add_candidate(cx, cy - ey);
            if (0 < dx) add_candidate(cx + ex, cy);
            if (0 > dx) add_candidate(cx - ex, cy);

            cx += dx;
            cy += dy;

            if (cy > bottom) {
                cy = bottom;
                cx = right;

                dx = -speed;
                dy = 0.0;
            }

            if (cx < left) {
                cx = left;
                cy = bottom;

                dy = -speed;
                dx = 0.0;
            }

            if (cy < top) {
                cy = top;
                cx = left;

                dx = speed;
                dy = 0.0;
            }

            if (cx > right) {
                cx = right;
                cy = top;

                dx = 0.0;
                dy = speed;
                break;
            }
        }

        vector<pair<double, double> > placements;
        for (auto it = candidates.rbegin(); it != candidates.rend(); ++it) {
            bool ok = true;
            for (unsigned i = 0; i < placements.size(); i++) {
                auto [x, y] = placements[i];
                auto [x2, y2] = it->second;
                if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                placements.push_back(it->second);
                if (placements.size() == problem.musicians.size()) break;
            }
        }

        gvNewTime();
        gvStage(problem);
        gvPlacements(placements);

        auto checkPlacements = [&](double x, double y) -> bool {
            for (unsigned i = 0; i < placements.size(); ++i) {
                auto [x2, y2] = placements[i];
                if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
                    return false;
                }
            }
            return true;
        };

        set<pair<double, double>> edge_points(placements.begin(), placements.end());

        while (placements.size() < problem.musicians.size()) {
            double x0 = problem.stageLeft + 10 + g_rand() % ((int)problem.stageWidth - 19);
            double y0 = problem.stageBottom + 10 + g_rand() % ((int)problem.stageHeight - 19);
            if (checkPlacements(x0, y0)) {
                placements.emplace_back(x0, y0);
            }
        }

        gvNewTime();
        gvStage(problem);
        gvPlacements(placements);

        {
            map<pair<double, double>, vector<long long>> cache;
            auto res = calcScoreWithCache(problem, placements, cache, true);
            auto score = res.second;
            cerr << "score: " << score << endl;

            bool updated = true;
            while (updated) {
                updated = false;
                for (int i = 0; i < placements.size(); i++) {
                    for (int j = i + 1; j < placements.size(); j++) {
                        if (problem.musicians[i] == problem.musicians[j]) continue;
                        auto ds = swapDeltaScore(problem, placements, cache, true, i, j);
                        if (0 < ds) {
                            score += ds;
                            std::swap(placements[i], placements[j]);
                            updated = true;
                        }
                    }
                }
            }
            cerr << "score: " << score << endl;
        }

		gvNewTime();
		auto res = calcScore(problem, placements);
        for (int i = 0; i < placements.size(); i++) {
            const auto& p = placements[i];
            if (0 < res.second[i] && edge_points.find(p) != edge_points.end()) {
                p_tastes.emplace(p, problem.musicians[i]);
                p_score.emplace(p, res.second[i]);
                gvCircle(p.first, p.second, 5, gvRGB(0, 0, 255));
                gvText(p.first, p.second, 3, "%.lf", double(res.second[i]));
            } else {
                free_tastes[problem.musicians[i]]++;
            }
        }

		gvStage(problem);
		gvPlacements(placements);
    }

    // step2. 1で得た人の周りに残りの人を配置
	{
        multimap<double, pair<double, double>> candidates;
        auto add_candidate = [&](double x, double y) {
            if (!(problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
                problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10
                )) {
                cerr << "Reject" << x << " " << y << endl;
                return;
            }
            candidates.emplace(1, make_pair(x, y));
        };
        for (auto& elem: p_score) {
            auto [x, y] = elem.first;
            auto dx = 1.0L / 2;
            auto dy = sqrtl(3) / 2;
        }
	}

	gvNewTime();
	gvClose();
    exit(0);

	vector <pair<double, double>> dummy;
    return dummy;
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

void writePlacementsJSON(std::ostream& os, const std::vector<pair<double, double> >& placements, const vector<long long>& aScore) {
    os << "{\"placements\":[";
    for (unsigned i = 0; i < placements.size(); i++) {
        if (i > 0) os << ",";
        os << "{\"x\":" << placements[i].first << ",\"y\":" << placements[i].second << "}";
    }
    os << "],\"volumes\":[";
    for (unsigned i = 0; i < aScore.size(); i++) {
        if (i > 0) os << ",";
        os << (aScore[i] > 0 ? "10" : "0");
    }
    os << "]}" << endl;
}

void readHeatmap(std::istream& is, std::map<pair<double, double>, double>& heatmap) {
    string str;
    while (std::getline(is, str)) {
        double score, x, y;
        if (sscanf_s(str.c_str(), "%lf,%lf,%lf", &score, &x, &y) == 3) {
            heatmap[make_pair(x, y)] = score;
        }
    }
}

int main(int argc, char* argv[]) {
#if LOCAL_DEBUG
    fs::current_path(R"(c:\projects\hasipon\icfpc2023\solver\inada)");
#endif
    int problem_id = atoi(getenv("PROBLEM_ID"));

    Problem problem;
    {
        ifstream ifs(string(getenv("REPO_ROOT")) + "/problems.kyopro/" + to_string(problem_id) + ".kyopro");
        readProblem(ifs, problem);
    }

    std::map<pair<double, double>, double> heatmap;
    {
        ifstream ifs(string(getenv("REPO_ROOT")) + "/heatmap/" + to_string(problem_id) + ".csv");
        readHeatmap(ifs, heatmap);
    }

    auto placement = solve(problem, heatmap);
    gvStage(problem);
    gvPlacements(placement);

    auto res = calcScore(problem, placement);
    writePlacementsJSON(cout, placement, res.second);
#if LOCAL_DEBUG
    {
        ofstream ofs(to_string(problem_id) + "-inada1-" + to_string(sumScore(res.second)) + ".json");
        writePlacementsJSON(ofs, placement, res.second);
    }
#endif
    if (!res.first) throw runtime_error("invalid placement");
    cerr << "score = " << sumScore(res.second) << endl;
}
