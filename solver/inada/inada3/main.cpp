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
#include <queue>
#include <algorithm>

using namespace std;
namespace fs = std::filesystem;

#if LOCAL_DEBUG
std::mt19937 g_rand(12345);
#else
std::random_device rd;
std::mt19937 g_rand(rd());
#endif

// from: https://github.com/ei1333/library
// Unlicense license
template< typename flow_t, typename cost_t >
struct PrimalDual {
    struct edge {
        int to;
        flow_t cap;
        cost_t cost;
        int rev;
        bool isrev;
    };

    vector< vector< edge > > graph;
    vector< cost_t > potential, min_cost;
    vector< int > prevv, preve;
    const cost_t INF;

    PrimalDual(int V) : graph(V), INF(numeric_limits< cost_t >::max()) {}

    void add_edge(int from, int to, flow_t cap, cost_t cost) {
        graph[from].emplace_back(edge { to, cap, cost, (int)graph[to].size(), false });
        graph[to].emplace_back(edge { from, 0, -cost, (int)graph[from].size() - 1, true });
    }

    cost_t min_cost_flow(int s, int t, flow_t f) {
        int V = (int)graph.size();
        cost_t ret = 0;
        using Pi = pair< cost_t, int >;
        priority_queue< Pi, vector< Pi >, greater< Pi > > que;
        potential.assign(V, 0);
        preve.assign(V, -1);
        prevv.assign(V, -1);

        while (f > 0) {
            min_cost.assign(V, INF);
            que.emplace(0, s);
            min_cost[s] = 0;
            while (!que.empty()) {
                Pi p = que.top();
                que.pop();
                if (min_cost[p.second] < p.first) continue;
                for (int i = 0; i < (int)graph[p.second].size(); i++) {
                    edge& e = graph[p.second][i];
                    cost_t nextCost = min_cost[p.second] + e.cost + potential[p.second] - potential[e.to];
                    if (e.cap > 0 && min_cost[e.to] > nextCost) {
                        min_cost[e.to] = nextCost;
                        prevv[e.to] = p.second, preve[e.to] = i;
                        que.emplace(min_cost[e.to], e.to);
                    }
                }
            }
            if (min_cost[t] == INF) return -1;
            for (int v = 0; v < V; v++) potential[v] += min_cost[v];
            flow_t addflow = f;
            for (int v = t; v != s; v = prevv[v]) {
                addflow = min(addflow, graph[prevv[v]][preve[v]].cap);
            }
            f -= addflow;
            ret += addflow * potential[t];
            for (int v = t; v != s; v = prevv[v]) {
                edge& e = graph[prevv[v]][preve[v]];
                e.cap -= addflow;
                graph[v][e.rev].cap += addflow;
            }
        }
        return ret;
    }

    void output() {
        for (int i = 0; i < graph.size(); i++) {
            for (auto& e : graph[i]) {
                if (e.isrev) continue;
                auto& rev_e = graph[e.to][e.rev];
                if (e.cap == 0)
					cerr << i << "->" << e.to << " (flow: " << rev_e.cap << "/" << rev_e.cap + e.cap << ")" << endl;
            }
        }
    }
};

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
    map<pair<double, double>, vector<long long>>& cache) {
    if (!problem.pillars.empty()) throw 1;
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
        score += cache[placements[i]][problem.musicians[i]];
    }

    return { true,score };
}

long long swapDeltaScore(
    const Problem& problem,
    const vector<pair<double, double>>& placements,
    map<pair<double, double>, vector<long long>>& cache,
    int swap_i, int swap_j
) {
    assert(!cache.empty());
    assert(cache.size() == placements.size());
    if (!problem.pillars.empty()) throw 1;

    long long score = 0;
    score -= cache[placements[swap_i]][problem.musicians[swap_i]];
    score -= cache[placements[swap_j]][problem.musicians[swap_j]];

    score += cache[placements[swap_i]][problem.musicians[swap_j]];
    score += cache[placements[swap_j]][problem.musicians[swap_i]];
    return score;
}

long long getPointTasteScore(pair<double, double> point, int taste, map<pair<double, double>, vector<long long>>& cache) {
    return cache[point][taste];
}

vector<pair<double, double>> solve(const Problem& problem, const map<pair<double, double>, double>& heatmap) {
    // step1. 外周に沿って配置, スコアに影響のないものは削除
	set<pair<double, double>> fixed_points;
    {
        multimap<double, pair<double, double>> candidates;
        auto add_candidate = [&](double x, double y) {
            if (!(problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
                problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10
                )) {
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
        const double nx = long long(w) / 10;
        const double ny = long long(h) / 10;
        const double ex = long long(w) % 10;
        const double ey = long long(h) % 10;

        const auto speed_x = 10.0 + ex / nx;
        const auto speed_y = 10.0 + ey / ny;
        const auto left = x;
        const auto right = x + w;
        const auto top = y;
        const auto bottom = y + h;

        auto dx = 0.0;
        auto dy = speed_y;
        auto cx = right;
        auto cy = y;

        for (;;) {
            add_candidate(cx, cy);

            cx += dx;
            cy += dy;

            if (cy > bottom) {
                cy = bottom;
                cx = right;

                dx = -speed_x;
                dy = 0.0;
            }

            if (cx < left) {
                cx = left;
                cy = bottom;

                dy = -speed_y;
                dx = 0.0;
            }

            if (cy < top) {
                cy = top;
                cx = left;

                dx = speed_x;
                dy = 0.0;
            }

            if (cx > right) {
                cx = right;
                cy = top;

                dx = 0.0;
                dy = speed_y;
                break;
            }
        }

        vector<pair<double, double> > placements;
        auto checkPlacements = [&placements](double x, double y) -> bool {
            for (unsigned i = 0; i < placements.size(); ++i) {
                auto [x2, y2] = placements[i];
                if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
                    return false;
                }
            }
            return true;
        };
        for (auto it = candidates.rbegin(); it != candidates.rend(); ++it) {
            auto [x, y] = it->second;
            if (checkPlacements(x, y)) {
                placements.push_back(it->second);
                if (placements.size() == problem.musicians.size()) {
                    break;
                }
            }
        }

        gvNewTime();
        gvStage(problem);
        gvPlacements(placements);


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
            auto res = calcScoreWithCache(problem, placements, cache);
            auto score = res.second;
            cerr << "score: " << score << endl;

            bool updated = true;
            while (updated) {
                updated = false;
                for (int i = 0; i < placements.size(); i++) {
                    for (int j = i + 1; j < placements.size(); j++) {
                        if (problem.musicians[i] == problem.musicians[j]) continue;
                        auto ds = swapDeltaScore(problem, placements, cache, i, j);
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
                fixed_points.emplace(p);
                gvCircle(p.first, p.second, 5, gvRGB(0, 0, 255));
                gvText(p.first, p.second, 3, "%.lf", double(res.second[i]));
            }
        }

		gvNewTime();
		gvStage(problem);
		gvPlacements(placements);
    }

    // step2. 1で得た人の周りに残りの人を配置
	{
        const auto w = problem.stageWidth - 20.0;
        const auto h = problem.stageHeight - 20.0;
        const double nx = long long(w) / 10;
        const double ny = long long(h) / 10;
        const double ex = long long(w) % 10;
        const double ey = long long(h) % 10;
        const auto speed_x = 10.0 + ex / nx;
        const auto speed_y = 10.0 + ey / ny;

        set<pair<double, double> > candidates_set;
        auto add_candidate = [&](double x, double y) {
			if (!(problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
				problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10)) {
                return;
            }
            candidates_set.emplace(x, y);
        };

        for (auto& fixed_p: fixed_points) {
            auto [x, y] = fixed_p;

            for (int t = 0; t < 2; t++) {
                auto sx = t == 0 ? speed_x : -speed_x;
                auto sy = t == 0 ? speed_y : -speed_y;
				if (fixed_points.find(make_pair(x + sx, y)) != fixed_points.end()) {
					auto cx = x + sx / 2;
					auto dy = sqrt(100.0 - sx * sx / 4) + 0.01;
					add_candidate(cx, y + dy);
					add_candidate(cx, y - dy);
				}

				if (fixed_points.find(make_pair(x, y + sy)) != fixed_points.end()) {
					auto cy = y + sy / 2;
					auto dx = sqrt(100.0 - sx * sx / 4) + 0.01;
					add_candidate(x + dx, cy);
					add_candidate(x - dx, cy);
				}
            }
        }

		for (const auto& p : candidates_set) {
			gvCircle(p.first, p.second, 5, gvRGB(0, 0, 255));
		}

		gvNewTime();
		gvStage(problem);
        // gvPlacements(placements);

        vector<pair<double, double>> placements(fixed_points.begin(), fixed_points.end());
        auto checkPlacements = [&placements](double x, double y) -> bool {
            for (unsigned i = 0; i < placements.size(); ++i) {
                auto [x2, y2] = placements[i];
                if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
                    return false;
                }
            }
            return true;
        };
        for (const auto& p : candidates_set) {
            if (checkPlacements(p.first, p.second)) {
                placements.push_back(p);
                if (placements.size() == problem.musicians.size()) break;
            }
        }
        set<pair<double, double>> edge_points(placements.begin(), placements.end());
        while (placements.size() < problem.musicians.size()) {
            double x0 = problem.stageLeft + 10 + g_rand() % ((int)problem.stageWidth - 19);
            double y0 = problem.stageBottom + 10 + g_rand() % ((int)problem.stageHeight - 19);
            if (checkPlacements(x0, y0)) {
                placements.emplace_back(x0, y0);
            }
        }

        {
            map<pair<double, double>, vector<long long>> cache;
            auto res = calcScoreWithCache(problem, placements, cache);
            auto score = res.second;
            cerr << "score: " << score << endl;

			bool updated = true;
			while (updated) {
				updated = false;
				for (int i = 0; i < placements.size(); i++) {
					for (int j = i + 1; j < placements.size(); j++) {
						if (problem.musicians[i] == problem.musicians[j]) continue;
						auto ds = swapDeltaScore(problem, placements, cache, i, j);
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

        {
            map<pair<double, double>, vector<long long>> cache;
            auto res = calcScoreWithCache(problem, placements, cache);
            auto score = res.second;
            cerr << "score: " << score << endl;

            // Flow
            vector<int> taste_count(problem.attendees[0].tastes.size());
            for (int taste : problem.musicians) {
                taste_count[taste]++;
            }
            auto points(placements);
            sort(points.begin(), points.end());
			PrimalDual<int, double> g(problem.attendees[0].tastes.size() + points.size() + 2);
            const int e_p0 = problem.attendees[0].tastes.size();
            const int e_s = problem.attendees[0].tastes.size() + placements.size();
            const int e_t = e_s + 1;

            int sum_taste_count = 0;
            // s -> taste
            for (int i = 0; i < problem.attendees[0].tastes.size(); i++) {
                g.add_edge(e_s, i, taste_count[i], 0);
                sum_taste_count += taste_count[i];
            }
            cerr << "sum_taste_count" << sum_taste_count << endl;

            // taste -> point
            for (int i = 0; i < problem.attendees[0].tastes.size() ; i++) {
                for (int j = 0; j < points.size(); ++j) {
                    g.add_edge(i, e_p0 + j, 1, -getPointTasteScore(points[j], i, cache));
                }
            }

            // point -> t
            for (int i = 0; i < points.size(); i++) {
                g.add_edge(e_p0 + i, e_t, 1, 0);
            }

            auto min_cost = g.min_cost_flow(e_s, e_t, points.size());
            cerr << "min_cost:" << min_cost << endl;
            cerr << "-min_cost:" << -long long(min_cost) << endl;

            cerr << "Points size" << points.size() << endl;

            multimap<int, pair<double, double>> taste_to_point;
            set<pair<double, double> > used;
            for (int i = 0; i < problem.attendees[0].tastes.size(); i++) {
                for (auto& e : g.graph[i]) {
                    if (e.isrev) continue;
					// auto& rev_e = g.graph[e.to][e.rev];
					if (e.cap == 0) {
						if (e_p0 <= e.to && e.to < e_p0 + points.size()) {
							const auto taste = i;
							const auto& point = points[e_p0 + e.to];
							if (used.find(point) == used.end()) {
								used.emplace(point);
								taste_to_point.emplace(taste, point);
							} else {
								throw 1;
							}
                        }
                        else {
                            cerr << "e_p0" << " = " << e_p0 << endl;
                            cerr << "e_s" << " = " << e_s << endl;
                            cerr << "e_t" << " = " << e_t << endl;
                            cerr << i << " -> " << e.to << endl;
                            throw 1;
                        }
					}
                }
            }

            if (taste_to_point.size() != points.size()) {
                cerr << taste_to_point.size() << endl;
                throw 1;
            }

            vector<pair<double, double> > best_placements(problem.musicians.size());
            for (int i = 0; i < problem.musicians.size(); ++i) {
                auto it = taste_to_point.find(problem.musicians[i]);
                if (it == taste_to_point.end()) {
                    throw 1;
                }
                best_placements[i] = it->second;
                taste_to_point.erase(it);
            }

            return best_placements;
        }

        {
			gvNewTime();
			auto res = calcScore(problem, placements);
            fixed_points.clear();
            gvNewTime();
            gvStage(problem);
            gvOutput("FixedPoints");
			for (int i = 0; i < placements.size(); i++) {
				const auto& p = placements[i];
				if (0 < res.second[i] && edge_points.find(p) != edge_points.end()) {
					fixed_points.emplace(p);
					gvCircle(p.first, p.second, 5, gvRGB(0, 0, 255));
                    gvText(p.first, p.second, 3, "%d", problem.musicians[i]);
				}
			}
        }

		gvNewTime();
		gvStage(problem);
        gvPlacements(placements);
		return placements;
	}

	gvNewTime();
	gvClose();
    exit(0);

	// vector <pair<double, double>> dummy;
    // return dummy;
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
        ofstream ofs(to_string(problem_id) + "-inada3-" + to_string(sumScore(res.second)) + ".json");
        writePlacementsJSON(ofs, placement, res.second);
    }
#endif
    if (!res.first) throw runtime_error("invalid placement");
    cerr << "score = " << sumScore(res.second) << endl;
}
