#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <complex>
#include <algorithm>

// github.com/Johniel/contests
// solver/johniel8/a.cpp

#include <bits/stdc++.h>

#define each(i, c) for (auto& i : c)
#define unless(cond) if (!(cond))
// #define endl "\n"

using namespace std;

template<typename P, typename Q> ostream& operator << (ostream& os, pair<P, Q> p) { os << "(" << p.first << "," << p.second << ")"; return os; }
template<typename P, typename Q> istream& operator >> (istream& is, pair<P, Q>& p) { is >> p.first >> p.second; return is; }
template<typename T> ostream& operator << (ostream& os, vector<T> v) { os << "("; for (auto& i: v) os << i << ","; os << ")"; return os; }
template<typename T> istream& operator >> (istream& is, vector<T>& v) { for (auto& i: v) is >> i; return is; }
template<typename T> ostream& operator << (ostream& os, set<T> s) { os << "#{"; for (auto& i: s) os << i << ","; os << "}"; return os; }
template<typename K, typename V> ostream& operator << (ostream& os, map<K, V> m) { os << "{"; for (auto& i: m) os << i << ","; os << "}"; return os; }
template<typename E, size_t N> istream& operator >> (istream& is, array<E, N>& a) { for (auto& i: a) is >> i; return is; }
template<typename E, size_t N> ostream& operator << (ostream& os, array<E, N>& a) { os << "[" << N << "]{"; for (auto& i: a) os << i << ","; os << "}"; return os; }

template<typename T> inline T setmax(T& a, T b) { return a = std::max(a, b); }
template<typename T> inline T setmin(T& a, T b) { return a = std::min(a, b); }

__attribute__((constructor)) static void ___initio(void) { ios_base::sync_with_stdio(false); cin.tie(nullptr); cout.setf(ios_base::fixed); cout.precision(15); return ; }

using lli = long long int;
using ull = unsigned long long;
using str = string;
template<typename T> using vec = vector<T>;

constexpr array<int, 8> di({0, 1, -1, 0, 1, -1, 1, -1});
constexpr array<int, 8> dj({1, 0, 0, -1, 1, -1, -1, 1});
constexpr lli mod = 1e9 + 7;
// constexpr lli mod = 998244353;

namespace geo {
  typedef complex<double> point;
  struct box { point mn, mx; };
  const double EPS = 1e-11;
  bool eq(double a, double b)
  {
    return fabs(a - b) < EPS;
  }
  point normal(point v)
  {
    return v * point(0, -1);
  }
  double dot(point a, point b)
  {
    return (a.real() * b.real() + a.imag() * b.imag());
  }
  double cross(point a, point b)
  {
    return (a.real() * b.imag() - a.imag() * b.real());
  }
  double distance_lp(point a1, point a2, point b)
  {
    if(dot(a2-a1, b-a1) < EPS)return abs(b-a1);
    if(dot(a1-a2, b-a2) < EPS)return abs(b-a2);
    return abs(cross(a2-a1, b-a1)) / abs(a2-a1);
  }
  double distance_bp(box b, point p)
  {
    point a1(b.mn.real(), b.mx.imag());
    point a2(b.mx.real(), b.mn.imag());
    return min(
      {
        distance_lp(b.mn, a1, p),
        distance_lp(b.mx, a1, p),
        distance_lp(b.mn, a2, p),
        distance_lp(b.mx, a2, p),
      });
  }
  double distance_pp(double x1, double y1, double x2, double y2)
  {
    double x = x1 - x2;
    double y = y1 - y2;
    return sqrt( x * x + y * y );
  }
  double distance_pp(pair<double, double> a, pair<double, double> b)
  {
    return distance_pp(a.first, a.second, b.first, b.second);
  }
};

uint32_t xorshift(void)
{
  // https://shindannin.hatenadiary.com/entry/2021/03/06/115415
  static uint32_t x = 123456789;
  static uint32_t y = 362436069;
  static uint32_t z = 521288629;
  static uint32_t w = 88675123;
  uint32_t t;

  t = x ^ (x << 11);
  x = y; y = z; z = w;
  return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
}

struct Attendee {
    double x;
    double y;
    vector<double> tastes;
  geo::point position(void) const { return geo::point(x, y); }
};
ostream& operator << (ostream& os, Attendee a) {
  os << "Attendee" << make_pair(make_pair(a.x, a.y), a.tastes);
  return os;
}


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
    // lli score = 0;
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

class HenkakuState{
public:
    double dir;
    double dist;
    bool isStart;
    bool isCenter;
    int id;
    double memoDir;

    bool operator<(const HenkakuState &s) const{
        return dir < s.dir;
    }
    bool operator>(const HenkakuState &s) const{
        return dir > s.dir;
    }
    double endDir() const {
        if(isStart){
            return dir + memoDir;
        }else {
            return dir;
        }
    }
};

double evalPlacement(const Problem &p, const vector<pair<double, double>>& placements, double x, double y, int tasteId){
    double nearestDist = 10000000000000000000000.0;
    double nearestDir = 0.0;
    auto center = make_pair(x, y);
    auto nodes = vector<HenkakuState>();

    // ミュージシャンのふちを追加
    for(int pid = 0; pid < placements.size(); pid++){
        auto p = placements[pid];
        double dx = p.first - center.first;
        double dy = p.second - center.second;
        double dist = sqrt(dx * dx + dy * dy);
        double asinRes;
        if(dist <= 5.0){
            asinRes = M_PI / 2.0;
        } else {
            asinRes = asin(5.0 / dist);
        }
        double dir = atan2(dx, dy) - asinRes;
        nodes.emplace_back(HenkakuState{dir, dist, true, false, pid, asinRes*2.0});
        if(dist < nearestDist){
            nearestDist = dist;
            nearestDir = dir;
        }
    }
    // 客の中心を追加
    for(int i = 0; i<p.attendees.size(); i++){
        auto a = p.attendees[i];
        double dx = a.x - center.first;
        double dy = a.y - center.second;
        double dist = sqrt(dx * dx +dy*dy);
        double dir = atan2(dx, dy);
        nodes.emplace_back(HenkakuState{dir,  dist, false, true, i});
    }
    for(int i =0; i<nodes.size(); i++){
        nodes[i].dir -= nearestDir;
        if(nodes[i].dir < 0.0){
            nodes[i].dir += M_PI * 2.0;
        }
    }
    sort(nodes.begin(), nodes.end());
    priority_queue<HenkakuState, vector<HenkakuState>, greater<> > heap;
    double result = 0;
    for(auto node : nodes){
        double nextDir = node.dir;
        while(!heap.empty() && heap.top().endDir() < nextDir){
            heap.pop();
        }
        if(node.isStart){
            heap.push(node);
        }else{
            auto a = p.attendees[node.id];
            if(heap.empty() || node.dist < heap.top().dist){
                result += 1000000.0 * a.tastes[tasteId] / node.dist / node.dist;
            }
        }
    }
    return result;
}

double calcScore3(const Problem& problem, const vector<unsigned>& perm, const vector<pair<double, double>>& placements, const set<pair<unsigned, unsigned>>& blocking, int taste, double x, double y) {
    return evalPlacement(problem, placements, x, y, taste);

    /*

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
                score -= 1000000* a.tastes[problem.musicians[j]] / d2;
            }
        }
        if (!blocked) {
            auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
            score += 1000000*a.tastes[taste] / d2;
        }
    }
    return score;
     */
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
    map<double, vector<unsigned>> yy;
    for (unsigned i = 0; i < N; ++ i) {
        auto s = yamaScore(problem, problem.musicians[i], problem.stageLeft + problem.stageWidth / 2, problem.stageBottom + problem.stageHeight / 2);
        yy[-s].push_back(i);
    }
    vector<unsigned> perm;
    for (auto& p : yy) {
        random_shuffle(p.second.begin(), p.second.end());
        for (auto i : p.second) perm.push_back(i);
    }
    vector<pair<double, double>> placements;
    set<pair<unsigned, unsigned>> blocking;
    for (auto i : perm) {
        int taste = problem.musicians[i];
        auto best = yama(problem, problem.musicians[i], makeStart(problem, placements), placements);
        double bestScore = calcScore3(problem, perm, placements, blocking, taste, best.first, best.second);
        for (int tt = 0; tt < 100; ++ tt) {
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

vector<pair<double, double>> solve2(const Problem& problem) {
    unsigned N = problem.musicians.size();
    vector<pair<double, double>> res(N);
    vector<unsigned > perm = calcPerm(problem);

    vector<pair<double, double>> placements;
    set<pair<unsigned, unsigned>> blocking;
    for (auto i : perm) {
        int taste = problem.musicians[i];
        auto best = yama(problem, problem.musicians[i], makeStart(problem, placements), placements);
        double bestScore = calcScore3(problem, perm, placements, blocking, taste, best.first, best.second);
        for (int tt = 0; tt < 100; ++ tt) {
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

namespace min_cost_flow {
  using Cap = long long int;
  using Cost = long long int;
  const Cost inf = (1LL << 60);

  struct E {
    int src, dst;
    Cap cap, flow;
    Cost cost;
    int rev;
    E(int s, int d, Cap cap_, Cost cost_, int r) : src(s), dst(d), cap(cap_), cost(cost_), rev(r), flow(0) {}
    Cap residue(void) const { return cap - flow; }
  };
  typedef vector<E> Es;
  typedef vector<Es> G;

  // max node conut
  const int N = 1000 * 5000 + 3;

  pair<int, int> path[N]; // path[dst]:=(src,index of e)
  Cost dist[N];
  Cost potential[N];

  void add_edge(G& g, const int src, const int dst, const Cap cap, const Cost cost)
  {
    assert(src < g.size());
    assert(dst < g.size());
    const int s = g[src].size();
    const int d = g[dst].size();
    g[src].push_back(E(src, dst, cap, +cost, d));
    g[dst].push_back(E(dst, src, 0,   -cost, s));
    return ;
  }

  bool sssp(const G &g, const int src, const int snk)
  {
    const int size = g.size();
    fill(dist, dist + size, inf);
    dist[src] = 0;
    path[src] = {src, -1};
    using S = pair<Cost, int>;
    priority_queue<S, vector<S>, greater<S>> q;
    for (q.push({0, src}); q.size();) {
      const S next = q.top();
      q.pop();
      if (dist[next.second] != next.first) continue;
      if (next.second == snk) break;
      for (int i = 0; i < g[next.second].size(); ++i) {
        const E& e = g[next.second][i];
        if (e.residue() <= 0) continue;
        const Cost rcost = e.cost + (potential[e.src] - potential[e.dst]);
        if (dist[e.dst] > rcost + dist[e.src]) {
          dist[e.dst] = rcost + dist[e.src];
          q.push({dist[e.dst], e.dst});
          path[e.dst] = make_pair(e.src, i);
        }
      }
    }
    return dist[snk] != inf;
  }

  pair<Cost, Cap> run(G& g, const int src, const int snk, Cap req)
  {
    assert(src < g.size());
    assert(snk < g.size());

    const int size = g.size();
    fill(potential, potential + size, 0);
    pair<Cost, Cap> result = {0, 0};
    while (0 < req && sssp(g, src, snk)) {
      for (int i = 0; i < size; ++i) {
        potential[i] += dist[i];
      }
      Cap mn = req;
      for (int i = snk; i != path[i].first; i = path[i].first) {
        const int v = path[i].first;
        const int e = path[i].second;
        mn = min(mn, g[v][e].residue());
      }
      for (int i = snk; i != path[i].first; i = path[i].first) {
        const int v = path[i].first;
        const int e = path[i].second;
        result.first += mn * g[v][e].cost;
        g[v][e].flow += mn;
        g[g[v][e].dst][g[v][e].rev].flow -= mn;
      }
      req -= mn;
      result.second += mn;
    }
    return result;
  }
};
namespace mcf = min_cost_flow;

lli calcInstBaseScore(const Problem& problem, const pair<double, double> pos, const int inst) {
  const double x = pos.first;
  const double y = pos.second;
  double score = 0;
  for (auto& a : problem.attendees) {
    auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
    score += (long long)ceil(1000000 * a.tastes[inst] / d2); // NOLINT(cppcoreguidelines-narrowing-conversions)
  }
  return score;
}

bool is_inside(const Problem& problem, double x, double y)
{
  return problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
    problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10;
}

bool is_inside(const Problem& problem, pair<double, double> p)
{
  return is_inside(problem, p.first, p.second);
}


pair<vector<pair<double, double>>, vec<double>> solveF(const Problem& problem)
{
  map<int, vec<int>> musicianIndexByTaste;
  for (int i = 0; i < problem.musicians.size(); ++i) {
    musicianIndexByTaste[problem.musicians[i]].push_back(i);
  }

  const int src = mcf::N - 1;
  const int snk = mcf::N - 2;

  const int M = problem.musicians.size();

  const pair<double, double> upperRight = make_pair(problem.stageLeft + problem.stageWidth, problem.stageBottom + problem.stageHeight);
  const pair<double, double> bottomLeft = make_pair(problem.stageLeft, problem.stageBottom);

  vec<double> volumes(problem.musicians.size(), 10.0);

  vec<pair<double, double>> candidates;
  for (int i = 2; ; ++i) {
    pair<double, double> p = upperRight;
    p.first -= 10.0 * i;
    p.second -= 10.0;
    if (!is_inside(problem, p)) break;
    candidates.push_back(p);
  }
  for (int j = 1; ; ++j) {
    pair<double, double> p = upperRight;
    p.first -= 10.0;
    p.second -= 10.0 * j;
    if (!is_inside(problem, p)) break;
    candidates.push_back(p);
  }

  for (int i = 2; ; ++i) {
    pair<double, double> p = bottomLeft;
    p.first += 10.0 * i;
    p.second += 10.0;
    if (!is_inside(problem, p)) break;
    candidates.push_back(p);
  }
  candidates.pop_back();
  for (int j = 1; ; ++j) {
    pair<double, double> p = bottomLeft;
    p.first += 10.0;
    p.second += 10.0 * j;
    if (!is_inside(problem, p)) break;
    candidates.push_back(p);
  }
  candidates.pop_back();

  const int enabled = candidates.size();

  for (int i = -10; i <= 10; ++i) {
    for (int j = -10; j <= 10; ++j) {
      pair<double, double> p = bottomLeft;
      p.first += problem.stageWidth / 2.0;
      p.second += problem.stageHeight / 2.0;
      p.first += 10.0 * i;
      p.second += 10.0 * j;
      if (is_inside(problem, p)) candidates.push_back(p);
    }
  }

  const int T = problem.attendees[0].tastes.size();
  lli mx = 0;
  for (int t = 0; t < T; ++t) {
    for (int i = 0; i < candidates.size(); ++i) {
      int j = t + candidates.size();
      mx = max(mx, calcInstBaseScore(problem, candidates[i], t));
    }
  }
  mcf::G g(mcf::N);
  for (int t = 0; t < T; ++t) {
    mcf::add_edge(g, src, t, musicianIndexByTaste[t].size(), mx);
  }
  for (int i = 0; i < candidates.size(); ++i) {
    mcf::add_edge(g, T + i, snk, 1, mx);
  }
  for (int t = 0; t < T; ++t) {
    for (int i = 0; i < candidates.size(); ++i) {
      mcf::add_edge(g, t, T + i, 1, mx - calcInstBaseScore(problem, candidates[i], t));
    }
  }
  cerr << mx << ' ' << candidates.size() << ' ' << T << endl;
  auto result = mcf::run(g, src, snk, 1LL << 60);
  cerr << result << endl;
  assert(0 < result.second);

  vector<pair<double, double>> placements(M, make_pair(-100, -100));
  for (int i = 0; i < mcf::N; ++i) {
    each (e, g[i]) {
      if (0 < e.flow) {
        if (e.src == src || e.dst == snk) continue;
        // cerr << make_pair(e.src, e.dst) << make_pair(e.flow, e.cap) << endl;
        if (0 <= e.src && e.src < T && T <= e.dst && e.dst < T + candidates.size()) {
          int a = e.src;
          int b = e.dst - T;
          // cerr << make_pair(a, b) << endl;
          assert(musicianIndexByTaste[a].size());
          placements[musicianIndexByTaste[a].back()] = candidates[b];
          if (enabled <= b) {
            volumes[musicianIndexByTaste[a].back()] = 0.0;
          }
          musicianIndexByTaste[a].pop_back();
        }
      }
    }
  }
  each (p, placements) {
    if (p.first < 0) p = makeStart(problem, placements);
  }
  return {placements, volumes};
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

    geo::box box;
    box.mn = geo::point(problem.stageLeft, problem.stageBottom);
    box.mx = geo::point(problem.stageLeft + problem.stageWidth, problem.stageBottom + problem.stageHeight);
    sort(problem.attendees.begin(), problem.attendees.end(), [&] (auto x, auto y) {
      return geo::distance_bp(box, x.position()) < geo::distance_bp(box, y.position());
    });
    while (true) {
      auto a = problem.attendees.back();
      if (geo::distance_bp(box, a.position()) < 300) break;
      problem.attendees.pop_back();
    }

    auto [placement, volumes] = solveF(problem);
    auto res = calcScore(problem, placement);
    cerr << "score = " << res.second << endl;

    cout << "{\"placements\":[";
    for (unsigned i = 0; i < placement.size(); i++) {
        if (i > 0) cout << ",";
        cout << "{\"x\":" << placement[i].first << ",\"y\":" << placement[i].second << "}";
    }
    cout << "],\"volumes\":[";
    for (unsigned i = 0; i < volumes.size(); i++) {
      if (i > 0) cout << ",";
      cout << volumes[i];
    }
    cout << "]}" << endl;

}
