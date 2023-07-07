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
  double score = 0;
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

pair<bool, long long> calcScoreNth(const Problem& problem, const vector<pair<double, double>>& placements, int nth) {
  if (placements.size() != problem.musicians.size()) {
    return {false, 0};
  }
  const unsigned i = nth;
  {
    auto [x, y] = placements[i];
    if (!(
          problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
          problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10
          )) {
      return {false, 0};
    }
    for (unsigned j = 0; j < placements.size(); j++) {
      if (i == j) continue;
      auto [x2, y2] = placements[j];
      if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
        return {false, 0};
      }
    }
  }
  double score = 0;
  for (auto& a : problem.attendees) {
    const unsigned i = nth;
    {
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
  return
    a.first + 10.0 <= p.first &&
    p.first + 10.0 <= b.first &&
    a.second + 10.0 <= p.second &&
    p.second + 10.0 <= b.second;
}

lli get_time(void) { return time(NULL); }

int main(int argc, char *argv[])
{
  const lli m_startTime = get_time();
  Problem problem;
  {
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
  }

  vec<point> candidates;

  point stagemn = make_pair(problem.stageLeft, problem.stageBottom);
  point stagemx = point(stagemn.first + problem.stageWidth, stagemn.second + problem.stageHeight);
  point b = stagemn;
  b.first += 10.0;
  b.second += 10.0;

  for (int i = 0; i < problem.musicians.size(); ++i) {
    for (int j = 0; j < problem.musicians.size(); ++j) {
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

  random_device seed_gen;
  mt19937 engine(seed_gen());
  // shuffle(candidates.begin(), candidates.end(), engine);

  vec<point> curr;
  while (curr.size() < problem.musicians.size()) {
    curr.push_back(candidates.back());
    candidates.pop_back();
  }
  pair<bool, lli> currScore = calcScore(problem, curr);

  vec<point> best = curr;
  pair<bool, lli> bestScore = calcScore(problem, best);
  assert(bestScore.first);

  // clog << x << endl;
  const lli msec = 30;
  const lli saTimeStart    = get_time();            // 焼きなまし開始時刻。get_timeは、時間を返す関数( {( {
  const lli saTimeEnd      = m_startTime+msec;     // 焼きなまし終了時刻。m_startTimeはこのプログラム自体の開始時間
  lli saTimeCurrent        = saTimeStart;          // 現在の時刻

  while ((saTimeCurrent = get_time()) < saTimeEnd) {
    vec<point> next = curr;
    pair<bool, lli> nextScore = currScore;
    const lli diffTime = saTimeEnd - saTimeCurrent;
    const lli diffTime2 = diffTime * diffTime;

    for (int _ = (int)max({3.0, sqrt(diffTime2), log2(diffTime2)}); --_; ) {
      const int replaced = xorshift() % curr.size();
      const int deployed = xorshift() % candidates.size();

      pair<bool, lli> r = calcScoreNth(problem, next, replaced);
      unless (r.first) {
        // --_;
        continue;
      }
      swap(candidates[deployed], next[replaced]);
      pair<bool, lli> d = calcScoreNth(problem, next, replaced);
      unless (d.first) {
        swap(candidates[deployed], next[replaced]);
        // --_;
        continue;
      }
      nextScore.second -= r.second;
      nextScore.second += d.second;
    }

    unless (nextScore.first) continue;

    // https://shindannin.hatenadiary.com/entry/20121224/1356364040
    const lli T = saTimeEnd-saTimeStart;
    const lli t = saTimeCurrent-saTimeStart;
    const lli R = 10000;
    const bool FORCE_NEXT = R*(T-t)>T*(xorshift()%R);

    cerr << currScore << "," << nextScore << "," << FORCE_NEXT << endl;
    if (currScore < nextScore || FORCE_NEXT) {
      curr = next;
      currScore = nextScore;
    }
    if (bestScore < nextScore) {
      best = next;
      bestScore = nextScore;
    }
  }

  show_placements(best);
  cerr << bestScore << endl;

  return 0;
}
