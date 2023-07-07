// github.com/Johniel/contests
// solution/johniel1/main.cpp

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
using point = pair<double, double>;

constexpr array<int, 8> di({0, 1, -1, 0, 1, -1, 1, -1});
constexpr array<int, 8> dj({1, 0, 0, -1, 1, -1, -1, 1});

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
  double x, y;
  vec<double> tastes;
  Attendee() {}
};

struct Problem {
  double roomWidth, roomHeight;
  double stageWidth, stageHeight;
  point stageBottomLeft;
  int instruments;
  vec<int> musicians;
  vec<Attendee> attendees;
};

istream& operator >> (istream& is, Problem& p)
{
  is >> p.roomWidth >> p.roomHeight;
  is >> p.stageWidth >> p.stageHeight;
  is >> p.stageBottomLeft;
  int ms;
  is >> ms >> p.instruments;
  p.musicians.resize(ms);
  is >> p.musicians;
  int as;
  is >> as;
  p.attendees.resize(as);
  each (a, p.attendees) {
    is >> a.x >> a.y;
    a.tastes.resize(p.instruments);
    is >> a.tastes;
  }
  return is;
}

void show_placements(vec<point> ps)
{
  cout << "{" << endl;
  cout << "\"placements\": [" << endl;
  for (int i = 0; i < ps.size(); ++i) {
    cout << "{\"x\": " << ps[i].first << ", \"y\": " << ps[i].second << "}";
    if (i + 1 < ps.size()) cout << ",";
    cout << endl;
  }
  cout << "]" << endl;
  cout << "}";
  return ;
}

bool is_inside(pair<double, double> a, pair<double, double> b, pair<double, double> p)
{
  return a.first <= p.first && p.first <= b.first && a.second <= p.second && p.second <= b.second;
}

int main(int argc, char *argv[])
{
  Problem p;
  cin >> p;

  vec<point> candidates;

  point stagemn = p.stageBottomLeft;
  point stagemx = point(stagemn.first + p.stageWidth, stagemn.second + p.stageHeight);
  point b = p.stageBottomLeft;
  b.first += 10.0;
  b.second += 10.0;

  candidates.push_back(b);
  for (int i = 0; i < p.musicians.size(); ++i) {
    for (int j = 0; j < p.musicians.size(); ++j) {
      point a;
      a.first = b.first + (i * 10.0);
      a.second = b.second + (j * 10.0);
      if (is_inside(stagemn, stagemx, a)) {
        candidates.push_back(a);
      } else {
        if (j == 0) i = (1 << 29);
        j = (1 << 29);
      }
    }
  }

  return 0;
}
