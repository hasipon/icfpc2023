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

using Placement = pair<double, double>;
using Placements = vec<Placement>;
using Volume = double;

Placement operator - (Placement a, Placement b) { return {a.first - b.first, a.second - b.second}; }
Placement operator + (Placement a, Placement b) { return {a.first + b.first, a.second + b.second}; }

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
  geo::point position(void) const { return geo::point(x, y); }};
ostream& operator << (ostream& os, Attendee a) {
  os << "Attendee" << make_pair(make_pair(a.x, a.y), a.tastes);
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

bool is_inside(const Problem& problem, double x, double y)
{
  return problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
    problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10;
}

bool is_inside(const Problem& problem, Placement p)
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

pair<bool, vector<long long>> calcScore(const Problem& problem, const vector<pair<double, double>>& placements) {
    if (placements.size() != problem.musicians.size()) {
        return {false, {}};
    }
    for (unsigned i = 0; i < placements.size(); i++) {
        auto [x, y] = placements[i];
        if (!(
                problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
                problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10
        )) {
          cerr << make_pair(i, placements[i]) << " is out of the stage" << endl;
            return {false, {}};
        }
        for (unsigned j = 0; j < i; j++) {
            auto [x2, y2] = placements[j];
            if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
              cerr << make_pair(i, placements[i]) << " and " << make_pair(j, placements[j]) << " are conflicted" << endl;
                return {false, {}};
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
    int a_idx = -1;
    for (auto& a : problem.attendees) {
      ++a_idx;
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
                // if (a_idx == 8 && problem.musicians[i] == 13) cerr << make_pair(a_idx, i) << a.tastes[problem.musicians[i]] << endl;
                // if (i == 13) cerr << make_pair(a_idx, i) << a << endl;
                if (factor[i] > 1) {
                    score[i] += (long long)ceil(10 * factor[i] * iik); // NOLINT(cppcoreguidelines-narrowing-conversions)
                } else {
                    score[i] += 10 * iik;
                }
            }
            next:;
        }
    }
    return {true,score};
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

bool checkPlacements2(double x, double y, const vector<pair<double, double>>& placements, unsigned idx) {
    for (unsigned i = 0; i < placements.size(); ++ i) {
        if (i == idx) continue;
        auto [x2, y2] = placements[i];
        if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
            return false;
        }
    }
    return true;
}

pair<int, int> makeStart2(const Problem& problem, const vector<pair<double, double>>& placements, int idx) {
    for (;;) {
        double x0 = problem.stageLeft + 10 + rand() % ((int)problem.stageWidth - 19);
        double y0 = problem.stageBottom + 10 + rand() % ((int)problem.stageHeight - 19);
        if (checkPlacements2(x0, y0, placements, idx)) return {x0, y0};
    }
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

pair<int, int> makeStart(const Problem& problem, const vector<pair<double, double>>& placements) {
    for (;;) {
        double x0 = problem.stageLeft + 10 + rand() % ((int)problem.stageWidth - 19);
        double y0 = problem.stageBottom + 10 + rand() % ((int)problem.stageHeight - 19);
        if (checkPlacements(x0, y0, placements)) return {x0, y0};
    }
}


pair<Placements, vec<Volume>> solve14(const Problem& problem) {
  constexpr double eps = 1e-11;
  map<int, vec<int>> musicianIndexesByInstrument;
  for (int k = 0; k < problem.musicians.size(); ++k) {
    musicianIndexesByInstrument[problem.musicians[k]].push_back(k);
  }
  for (int k = 0; k < problem.musicians.size(); ++k) cerr << make_pair(k, problem.musicians[k]); cerr << endl;
  // (0,1) (1,3) (2,1) (3,3) (4,3) (5,2) (6,4) (7,2) (8,1) (9,2) (10,4) (11,3) (12,1) (13,2) (14,2)
  const int M = problem.musicians.size();
  Placements placements(M, {-100.0, -100.0});
  vec<Volume> volumes(M, 0.0);

// 5:Attendee((633,420),(100,-1237,-9375,-2009,-549,-6490,25,-5368,-7747,-3660,-7313,-6043,-765,-9502,-8033,)),7
// 8:Attendee((725,389),(100,-8956,-8144,-2330,-8553,-5303,-6105,-8065,-8400,-6954,-2930,-3143,-7431,66,-3861,)),14 <<<
// 17:Attendee((617,474),(100,18,-2022,-5457,-2152,-7627,-8498,-8065,-4641,-6700,-6120,-321,-9587,-9669,-3018,)),23
// 19:Attendee((615,486),(100,-5369,-987,-5144,-7748,-9674,-651,-9189,59,-4614,-7112,-647,-5960,-5138,-769,)),25
// 27:Attendee((642,372),(100,24,-2097,-3545,59,-5614,-546,-9095,-4324,-9452,-9202,-2707,-1347,-4026,-4045,)),31
// 37:Attendee((597,574),(100,-4555,-649,-6526,-1940,-9924,83,-3871,-5887,-3289,-8768,-9859,32,-8555,-8328,)),43.1045

  int z = musicianIndexesByInstrument[0][0];
  placements[z] = make_pair(650.0, 413.0);
  volumes[z] = 10.0;

  const int X = 725.0;
  const int Y = 413.0;
  placements[musicianIndexesByInstrument[13][0]] = make_pair(X, Y + 10.0);
  volumes[musicianIndexesByInstrument[13][0]] = 10.0;

  int x1 = musicianIndexesByInstrument[1][0];
  int x2 = musicianIndexesByInstrument[1][1];
  placements[x1] = make_pair(X - 5.0 - eps, Y);
  placements[x2] = make_pair(X + 5.0 + eps, Y);

  placements[musicianIndexesByInstrument[3][0]] = make_pair(X - 10.0, Y + 10.0);
  placements[musicianIndexesByInstrument[3][1]] = make_pair(X + 10.0, Y + 10.0);
  // placements[musicianIndexesByInstrument[3][2]] = make_pair(X + 00.0, Y + 20.0);
  placements[musicianIndexesByInstrument[13][1]] = make_pair(X + 00.0, Y + 20.0);

  placements[musicianIndexesByInstrument[7][0]] = make_pair(X - 10.0, Y + 20.0);
  placements[musicianIndexesByInstrument[7][1]] = make_pair(X + 10.0, Y + 20.0);

  int y1 = musicianIndexesByInstrument[9][0];
  int y2 = musicianIndexesByInstrument[9][1];
  placements[y1] = placements[x1] + Placement(-10.0, 0);
  placements[y2] = placements[x2] + Placement(+10.0, 0);

  Placement curr;
  Placements restart(
    {
      placements[y2],
      placements[musicianIndexesByInstrument[3][1]],
      placements[musicianIndexesByInstrument[7][1]],
      placements[musicianIndexesByInstrument[7][1]] + Placement(10.0, 10.0),
      placements[musicianIndexesByInstrument[7][1]] + Placement(20.0, 20.0),
      placements[musicianIndexesByInstrument[7][1]] + Placement(30.0, 30.0),
    });
  int cnt = 0;
  each (p, placements) {
    if (p.first < 0) {
      unless (is_inside(problem, curr)) {
        curr = restart[cnt++] + Placement(+10.0, 0);
      }
      p = curr;
      curr = curr + Placement(+10.0, 0);
    }
  }

  return {placements, volumes};
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
  for (int i = 0; i < 600; ++i) {
    unless (650 <= problem.attendees[i].x && problem.attendees[i].x < 786 || 420 <= problem.attendees[i].y && problem.attendees[i].y < 520) continue;
    if (0 < *max_element(problem.attendees[i].tastes.begin()+1, problem.attendees[i].tastes.end())) {
      cerr << i << ":" << problem.attendees[i] << "," << geo::distance_bp(box, problem.attendees[i].position()) << endl;
    }
  }
  map<int, vec<int>> musicianIndexesByInstrument;
  for (int k = 0; k < problem.musicians.size(); ++k) {
    musicianIndexesByInstrument[problem.musicians[k]].push_back(k);
  }
  each (i, musicianIndexesByInstrument) cerr << make_pair(i.first, i.second.size()) << ' '; cerr << endl;


    auto [placements, volumes] = solve14(problem);

    auto score = calcScore(problem, placements);
    cerr << "score = " << score << endl;
    vec<lli> volumeScores;
    for (int i = 0; i < placements.size(); ++i) {
      if (0 < volumes[i]) {
        volumeScores.push_back(score.second[i]);
      } else {
        volumeScores.push_back(0);
      }
    }
    cerr << volumeScores << accumulate(volumeScores.begin(), volumeScores.end(), 0LL) << endl;
    cerr << volumes << endl;

    cout << "{\"placements\":[";
    for (unsigned i = 0; i < placements.size(); i++) {
        if (i > 0) cout << ",";
        cout << "{\"x\":" << placements[i].first << ",\"y\":" << placements[i].second << "}";
    }
    cout << "],\"volumes\":[";
    for (unsigned i = 0; i < volumes.size(); i++) {
        if (i > 0) cout << ",";
        cout << (volumes[i] > 0 ? "10" : "0");
    }
    cout << "]}" << endl;
}
