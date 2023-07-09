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
  os << make_pair(make_pair(a.x, a.y), a.tastes);
  return os;
}

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

using Placement = pair<double, double>;
using Placements = vec<Placement>;

bool is_inside(const Problem& problem, double x, double y)
{
  return problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
    problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10;
}

bool is_inside(const Problem& problem, pair<double, double> p)
{
  return is_inside(problem, p.first, p.second);
}

bool isBlocked(double x0, double y0, double x1, double y1, double x2, double y2, double radius) {
  complex<double> p0(x0, y0), p1(x1, y1), p2(x2, y2);
  if (real(conj(p1 - p0) * (p2 - p0)) < 0) return false;
  if (real(conj(p0 - p1) * (p2 - p1)) < 0) return false;
  double t = real(conj(p2 - p0) * (p0 - p1)) / norm(p0 - p1);
  return abs(p2 - (p0 + (p0 - p1) * t)) < radius;
}

pair<bool, long long> calcScore(const Problem& problem, const vector<pair<double, double>>& placements) {
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
        cerr << make_pair(i, placements[i]) << "and" << make_pair(j, placements[j]) << "are conflicted"<< endl;
        return {false, 0};
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


void show_result(const Problem& problem, const vec<pair<double, double>> placements)
{
  auto score = calcScore(problem, placements);
  cerr << "score = " << score << endl;
  cout << "{\"placements\":[";
  for (unsigned i = 0; i < placements.size(); i++) {
    if (i > 0) cout << ",";
    cout << "{\"x\":" << placements[i].first << ",\"y\":" << placements[i].second << "}";
  }
  cout << "]}" << endl;
  return ;
}


lli swapMusicians(const Problem& problem, vector<pair<double, double>> placements, const int i, const int j) {
  const int instI = problem.musicians[i];
  const int instJ = problem.musicians[j];
  if (instI == instJ) return 0;
  const vec<int> insts({instI, instJ});
  double score = 0;
  for (int _ = 0; _ < 2; ++_) {
    map<int, vec<int>> musicianIndexesByInstrument;
    for (int k = 0; k < problem.musicians.size(); ++k) {
      musicianIndexesByInstrument[problem.musicians[k]].push_back(k);
    }

    vector<double> factor(placements.size(), 1);
    if (!problem.pillars.empty()) {
      each (inst, insts) {
        each (a, musicianIndexesByInstrument[inst]) {
          each (b, musicianIndexesByInstrument[inst]) {
            if (a == b) continue;
            auto [x1, y1] = placements[a];
            auto [x2, y2] = placements[b];
            auto d2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
            factor[i] += 1.0 / sqrt(d2);
          }
        }
      }
    }
    static map<pair<int, pair<double, double>>, bool> memo;
    for (int ai = 0; ai < problem.attendees.size(); ++ai) {
      auto& a = problem.attendees[ai];
      each (inst, insts) {
        each (I, musicianIndexesByInstrument[inst]) {
          const auto key = make_pair(ai, placements[I]);
          auto [x, y] = placements[I];
          if (memo.count(key)) {
            if (memo[key]) {goto next;}
            else           {goto valid;}
          }
          for (unsigned J = 0; J < placements.size(); J++) {
            if (I != J) {
              auto [x2, y2] = placements[J];
              if (isBlocked(a.x, a.y, x, y, x2, y2, 5)) {
                memo[key] = true;
                goto next;
              }
            }
          }
          for (auto pr : problem.pillars) {
            if (isBlocked(a.x, a.y, x, y, pr.x, pr.y, pr.r)) {
              memo[key] = true;
              goto next;
            }
          }
          {
            memo[key] = false;
            valid:;
            auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
            auto s = (long long)ceil(1000000 * a.tastes[inst] / d2);
            if (factor[inst] > 1) {
              s = (long long)ceil(factor[inst] * s); // NOLINT(cppcoreguidelines-narrowing-conversions)
            }
            if (_) score += s; // NOLINT(cppcoreguidelines-narrowing-conversions)
            else   score -= s; // NOLINT(cppcoreguidelines-narrowing-conversions)
          }
          next:;
        }
      }
    }
    if (_) break;
    swap(placements[i], placements[j]);
  }
  return score;
}

bool checkPlacements(double x, double y, const vector<pair<double, double>>& placements) {
  for (unsigned i = 0; i < placements.size(); ++ i) {
    auto [x2, y2] = placements[i];
    if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
      return false;
    }
  }
  return true;
}

pair<double, double> makeStart(const Problem& problem, const vector<pair<double, double>>& placements) {
  for (;;) {
    double x0 = problem.stageLeft + 10 + rand() % ((int)problem.stageWidth - 19);
    double y0 = problem.stageBottom + 10 + rand() % ((int)problem.stageHeight - 19);
    if (checkPlacements(x0, y0, placements)) return {x0, y0};
  }
}

Placements solve10(const Problem& problem)
{
  map<int, vec<int>> musicianIndexesByInstrument;
  for (int k = 0; k < problem.musicians.size(); ++k) {
    musicianIndexesByInstrument[problem.musicians[k]].push_back(k);
  }
  each (i, musicianIndexesByInstrument) cerr << i.first << ' ' << i.second.size() << endl;

  const pair<double, double> stageBottomLeft = make_pair(problem.stageLeft, problem.stageBottom);
  vec<pair<double, double>> placements(problem.musicians.size());
  vec<int> fixed3;

  pair<double, double> curr = stageBottomLeft;
  curr.first += 20.0;
  curr.second += 10.0;
  assert(is_inside(problem, curr));
  pair<double, double> bbb = curr;
  int index3 = 0;
  {
    {
      cerr << "bbb" << bbb << endl;
      bbb.first -= sqrt(10.0 * 10.0 - 5.0 * 5.0) + geo::EPS;
      bbb.second += 10.0/2.0;
      placements[musicianIndexesByInstrument[3][index3++]] = bbb;
    }

    for (int k = 0; k < 8; ++k) {
      placements[musicianIndexesByInstrument[3][index3++]] = curr;
      curr.first += 10.0;
    }
    curr.second += 10.0;
    pair<double, double> corner = curr;
    corner.first -= 10.0;
    curr = corner;
    curr.second -= 5.0;
    curr.first += sqrt(10.0 * 10.0 - 5.0 * 5.0) + geo::EPS;
    for (int k = 0; k < 10; ++k) {
      if (k != 7 && k != 6) {
        placements[musicianIndexesByInstrument[3][index3++]] = curr;
      }
      curr.second += 10.0;
    }
    // assert(musicianIndexesByInstrument[3].size() == index3);
    vec<pair<double, double>> v;
    for (int i = 0; i < 9; ++i) {
      for (int j = 0; j < 8; ++j) {
        pair<double, double> b = corner;
        b.first -= 10 * j;
        b.second += 10 * i;
        v.push_back(b);
        assert(is_inside(problem, b));
      }
    }
    // sort(v.begin(), v.end(), [&] (auto a, auto b) { return geo::distance_pp(a, corner) < geo::distance_pp(b, corner); });
    for (int inst = 0, k = 0; inst < 3; ++inst) {
      each (j, musicianIndexesByInstrument[inst]) {
        placements[j] = v[k++];
      }
    }

    while (index3 < musicianIndexesByInstrument[3].size()) {
      fixed3.push_back(musicianIndexesByInstrument[3][index3++]);
    }
    // assert(fixed3.size() == 4);
    // placements[fixed3[0]] = make_pair(320,993);
    // while (!is_inside(problem, placements[fixed3[0]])) {
    //   // cerr << placements[fixed3[0]] << endl;
    //   placements[fixed3[0]].second -= 1.0;
    // }

    // placements[fixed3[1]] = placements[fixed3[0]];
    // placements[fixed3[1]].first -= 10.0;

    const pair<double, double> stageUpperLeft = make_pair(problem.stageLeft, problem.stageBottom + problem.stageHeight);
    placements[fixed3[0]] = bbb;
    placements[fixed3[0]].second += 10.0;

    placements[fixed3[1]] = bbb;
    placements[fixed3[1]].second += 20.0;

    placements[fixed3[2]] = bbb;
    placements[fixed3[2]].second += 30.0;

    placements[fixed3[3]] = bbb;
    placements[fixed3[3]].second += 40.0;

    placements[fixed3[4]] = bbb;
    placements[fixed3[4]].second += 50.0;

    placements[fixed3[5]] = bbb;
    placements[fixed3[5]].second += 60.0;

    cerr << v.size() << endl;
    // cerr << v << endl;
    // cerr << placements << endl;
    cerr << "corner:" << corner << endl;
    assert(27 + 24 + 20 <= v.size());
    assert(is_inside(problem, corner));
    assert(musicianIndexesByInstrument[3].size() == index3);
    each (i, placements) {
      unless (is_inside(problem, i)) {
        cerr << i << endl;
        assert(is_inside(problem, i));
      }
    }
  }

  while (true) {
    break;

    Placements curr = placements;
    for (int i = 0; i < curr.size(); ++i) {
      unless (count(fixed3.begin(), fixed3.end(), i)) curr[i].second += 0.5;
    }
    bool f = true;

    for (int i = 0; i < curr.size(); ++i) {
      f = f && is_inside(problem, curr[i]);
      each (j, fixed3) {
        if (i == j) continue;
        auto [x1, y1] = curr[i];
        auto [x2, y2] = curr[j];
        f = f && (100 <= (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
      }
    }
    unless (f) break;
    placements = curr;
  }
  return placements;
}

int main(int argc, char** argv) {
  Problem problem;
  {
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
  }

  geo::box box;
  box.mn = geo::point(problem.stageLeft, problem.stageBottom);
  box.mx = geo::point(problem.stageLeft + problem.stageWidth, problem.stageBottom + problem.stageHeight);
  sort(problem.attendees.begin(), problem.attendees.end(), [&] (auto x, auto y) {
    return geo::distance_bp(box, x.position()) < geo::distance_bp(box, y.position());
  });
  for (int i = 0; i < 30; ++i) {
    if (0 < *max_element(problem.attendees[i].tastes.begin(), problem.attendees[i].tastes.end())) {
      cerr << i << ":" << problem.attendees[i] << "," << geo::distance_bp(box, problem.attendees[i].position()) << endl;
    }
  }

  auto placements = solve10(problem);
  show_result(problem, placements);

  return 0;
}



// void ______()
// {
//   vector<pair<double, double>> placements;
//   {
//     pair<double, double> curr = make_pair(problem.stageLeft, problem.stageBottom);
//     curr.first += 10.0;
//     curr.second += 10.0;
//     int idx = 0;
//     for (int dir = 0; dir < 4; ++dir) {
//       const int dx[] = {0, +1, 0, -1};
//       const int dy[] = {+1, 0, -1, 0};
//       while (placements.size() < problem.musicians.size()) {
//         if (placements.empty() || placements.back() != curr) placements.push_back(curr);
//         pair<double, double> next = curr;
//         next.first  += dx[dir] * 10.0;
//         next.second += dy[dir] * 10.0;
//         unless (is_inside(problem, next)) break;
//         curr = next;
//       }
//     }
//     auto [x1, y1] = placements.front();
//     auto [x2, y2] = placements.back();
//     auto d2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
//     if (d2 < 100) placements.pop_back();
//     while (placements.size() < problem.musicians.size()) {
//       placements.push_back(makeStart(problem, placements));
//     }
//   }
//   auto score = calcScore(problem, placements);
//   assert(score.first);
//   for (int _ = 0; _ < 5000; ++_) {
//     int i = xorshift() % problem.musicians.size();
//     int j = xorshift() % problem.musicians.size();
//     if (problem.musicians[i] == problem.musicians[j]) continue;
//     lli diff = swapMusicians(problem, placements, i, j);
//     cerr << _ << ' ' << make_pair(i, j) << ' ' << diff << ", " << score << endl;
//     if (0 < diff) {
//       swap(placements[i], placements[j]);
//       score.second += diff;
//       if (_ % 20 == 0) score = calcScore(problem, placements);
//     }
//   }
//   cerr << score << endl;
// }
