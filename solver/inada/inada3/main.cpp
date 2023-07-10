#define LOCAL_DEBUG 0
#define ENABLE_GV 0

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
#include "json.hpp"
using json = nlohmann::json;
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

vector<pair<double, double>> solve(const Problem& problem, vector<pair<double, double>> placements) {
    // Prev score
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
    for (int i = 0; i < problem.attendees[0].tastes.size(); i++) {
        for (int j = 0; j < points.size(); ++j) {
            g.add_edge(i, e_p0 + j, 1, -getPointTasteScore(points[j], i, cache));
        }
    }

    // point -> t
    for (int i = 0; i < points.size(); i++) {
        g.add_edge(e_p0 + i, e_t, 1, 0);
    }

    auto min_cost = g.min_cost_flow(e_s, e_t, points.size());
    cerr << "-min_cost:" << -long long(min_cost) << endl;

    multimap<int, pair<double, double>> taste_to_point;
    for (int i = 0; i < problem.attendees[0].tastes.size(); i++) {
        for (auto& e : g.graph[i]) {
            if (e.isrev) continue;
            if (e.cap == 0 && e_p0 <= e.to && e.to < e_p0 + points.size()) {
                const auto taste = i;
                const auto& point = points[e.to - e_p0];
                taste_to_point.emplace(taste, point);
            }
        }
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

    if (getenv("SOLUTION") != nullptr) {
        std::ifstream f(getenv("SOLUTION"));
        json data = json::parse(f);
        vector<pair<double, double>> placements;
        for (auto& e : data["placements"]) {
            double x, y;
            e["x"].get_to(x);
            e["y"].get_to(y);
            placements.emplace_back(x, y);
        }

        auto best_placements = solve(problem, placements);
        if (!best_placements.empty()) {
            auto res = calcScore(problem, best_placements);
            writePlacementsJSON(cout, best_placements, res.second);
            if (!res.first) throw runtime_error("invalid placement");
            cerr << "score = " << sumScore(res.second) << endl;
        }
    }

    for (int i = 1; i < argc; i++) {
        std::ifstream f(argv[i]);
        json data = json::parse(f);
        vector<pair<double, double>> placements;
        for (auto& e : data["placements"]) {
            double x, y;
            e["x"].get_to(x);
            e["y"].get_to(y);
            placements.emplace_back(x, y);
        }

        auto best_placements = solve(problem, placements);
        if (!best_placements.empty()) {
            auto res = calcScore(problem, best_placements);
            writePlacementsJSON(cout, best_placements, res.second);
            if (!res.first) throw runtime_error("invalid placement");
            cerr << "score = " << sumScore(res.second) << endl;
        }
    }

#if LOCAL_DEBUG
    {
        ofstream ofs(to_string(problem_id) + "-inada3-" + to_string(sumScore(res.second)) + ".json");
        writePlacementsJSON(ofs, placement, res.second);
    }
#endif
}
