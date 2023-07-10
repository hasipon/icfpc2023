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

void show_result(Placements placements, vec<Volume> volumes)
{
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
  return ;
}

Placement operator += (Placement& a, Placement b) {
  a.first += b.first;
  a.second += b.second;
  return a;
}
Placement operator -= (Placement& a, Placement b) {
  a.first -= b.first;
  a.second -= b.second;
  return a;
}
Placement operator + (Placement a, Placement b) { return a += b; }
Placement operator - (Placement a, Placement b) { return a -= b; }

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

  pair<double, double> rot(double x, double y, double th)
  {
    double a = x * cos(th) - y * sin(th);
    double b = x * sin(th) + y * cos(th);
    return make_pair(a, b);
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

struct Pillar {
  double x;
  double y;
  double r;
  geo::point position(void) const { return geo::point(x, y); }
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

  map<int, vec<int>> musician_indexes_by_instrument(void) const {
    map<int, vec<int>> musicianIndexesByInstrument;
    for (int k = 0; k < musicians.size(); ++k) {
      musicianIndexesByInstrument[musicians[k]].push_back(k);
    }
    return musicianIndexesByInstrument;
  }
};

istream& operator >> (istream& is, Problem& problem)
{
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
  return is;
}

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

pair<bool, vector<lli>> calcScore(const Problem& problem, const Placements& placements) {
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
  vector<lli> score(placements.size());
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
        auto iik = (lli)ceil(1000000 * a.tastes[problem.musicians[i]] / d2);
        if (factor[i] > 1) {
          score[i] += (lli)ceil(10 * factor[i] * iik); // NOLINT(cppcoreguidelines-narrowing-conversions)
        } else {
          score[i] += 10 * iik;
        }
      }
      next:;
    }
  }
  return {true,score};
}

map<int, vec<lli>> swapCache;
map<pair<pair<int, int>, int>, vec<lli>> swappedCache;
map<pair<int, int>, bool> blockedCache;
pair<bool, vec<lli>> calcSwapScore(const Problem& problem, Placements placements, const int target1, const int target2)
{
  assert(problem.musicians[target1] != problem.musicians[target2]);
  if (placements.size() != problem.musicians.size()) return {false, {}};
  auto fn = [&placements, &problem] (
    const int target,
    const bool D,
    map<int, vec<lli>>& swapCache,
    map<pair<pair<int, int>, int>, vec<lli>>& swappedCache,
    pair<int, int> key,
    map<pair<int, int>, bool>& blockedCache)
  {
    if (D && swapCache.count(target)) {
      return make_pair(true, swapCache[target]);
    }
    auto K = make_pair(key, target);
    if (!D && swappedCache.count(K)) {
      return make_pair(true, swappedCache[K]);
    }
    const unsigned i = target; {
      auto [x, y] = placements[i];
      if (!(problem.stageBottom + 10 <= y && y <= problem.stageBottom + problem.stageHeight - 10 &&
            problem.stageLeft + 10 <= x && x <= problem.stageLeft + problem.stageWidth - 10)) {
        cerr << make_pair(i, placements[i]) << " is out of the stage" << endl;
        return make_pair(false, vec<lli>());
      }
      for (unsigned j = 0; j < i; j++) {
        auto [x2, y2] = placements[j];
        if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
          cerr << make_pair(i, placements[i]) << " and " << make_pair(j, placements[j]) << " are conflicted" << endl;
          return make_pair(false, vec<lli>());
        }
      }
    }
    vector<double> factor(placements.size(), 1);
    if (!problem.pillars.empty()) {
      const unsigned i = target; {
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
    vector<lli> score(placements.size(), 0);
    int a_idx = -1;
    for (auto& a : problem.attendees) {
      ++a_idx;
      const unsigned i = target; {
        auto [x, y] = placements[i];
        pair<int, int> KEY = make_pair(a_idx, i);
        if (blockedCache.count(KEY)) {
          if (blockedCache[KEY]) goto next;
          else                   goto notblocked;
        }
        for (unsigned j = 0; j < placements.size(); j++) {
          if (i != j) {
            auto [x2, y2] = placements[j];
            if (isBlocked(a.x, a.y, x, y, x2, y2, 5)) {
              blockedCache[key] = true;
              goto next;
            }
          }
        }
        for (auto pr : problem.pillars) {
          if (isBlocked(a.x, a.y, x, y, pr.x, pr.y, pr.r)) {
            blockedCache[KEY] = true;
            goto next;
          }
        }
        {
          blockedCache[KEY] = false;
          notblocked:;
          auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
          auto iik = (lli)ceil(1000000 * a.tastes[problem.musicians[i]] / d2);
          if (factor[i] > 1) {
            score[i] += (lli)ceil(10 * factor[i] * iik); // NOLINT(cppcoreguidelines-narrowing-conversions)
          } else {
            score[i] += 10 * iik;
          }
        }
        next:;
      }
    }
    if (D) swapCache[target] = score;
    else swappedCache[K] = score;
    return make_pair(true, score);
  };
  pair<bool, vec<lli>> a1 = fn(target1, true, swapCache, swappedCache, make_pair(target1, target2), blockedCache);
  pair<bool, vec<lli>> a2 = fn(target2, true, swapCache, swappedCache, make_pair(target1, target2), blockedCache);
  unless (a1.first && a2.first) return {false, {}};
  unless (a1.second[target1] || a2.second[target2]) {
    return {true, vec<lli>(problem.musicians.size(), 0)};
  }
  swap(placements[target1], placements[target2]);
  pair<bool, vec<lli>> b1 = fn(target1, false, swapCache, swappedCache, make_pair(target1, target2), blockedCache);
  pair<bool, vec<lli>> b2 = fn(target2, false, swapCache, swappedCache, make_pair(target2, target2), blockedCache);
  unless (b1.first && b2.first) return {false, {}};
  vec<lli> scores(problem.musicians.size());
  for (int i = 0; i < problem.musicians.size(); ++i) {
    scores[i] = b1.second[i] + b2.second[i] - a1.second[i] - a2.second[i];
  }
  return {true, scores};
}

lli sumScore(const vector<lli>& score) {
  lli sum = 0;
  for (auto s : score) {
    if (0 < s) sum += s;
  }
  return sum;
}

bool checkPlacements2(double x, double y, const Placements& placements, unsigned idx) {
  for (unsigned i = 0; i < placements.size(); ++ i) {
    if (i == idx) continue;
    auto [x2, y2] = placements[i];
    if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
      return false;
    }
  }
  return true;
}

pair<int, int> makeStart2(const Problem& problem, const Placements& placements, int idx) {
  for (;;) {
    double x0 = problem.stageLeft + 10 + rand() % ((int)problem.stageWidth - 19);
    double y0 = problem.stageBottom + 10 + rand() % ((int)problem.stageHeight - 19);
    if (checkPlacements2(x0, y0, placements, idx)) return {x0, y0};
  }
}

bool checkPlacements(double x, double y, const Placements& placements)
{
  for (unsigned i = 0; i < placements.size(); ++ i) {
    auto [x2, y2] = placements[i];
    if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
      return false;
    }
  }
  return true;
}

pair<int, int> makeStart(const Problem& problem, const Placements& placements)
{
  for (;;) {
    double x0 = problem.stageLeft + 10 + rand() % ((int)problem.stageWidth - 19);
    double y0 = problem.stageBottom + 10 + rand() % ((int)problem.stageHeight - 19);
    if (checkPlacements(x0, y0, placements)) return {x0, y0};
  }
}

pair<Placements, vec<Volume>> read_output(const Problem& problem, const str filepath)
{
  const int M = problem.musicians.size();
  Placements placements(M);
  vec<Volume> volumes(M);
  assert(M);
  {
    ifstream fin(filepath);
    str s;
    getline(fin, s);
    each (c, s) unless (isdigit(c) || c == '.') c = ' ';
    istringstream sin(s);
    each (p, placements) assert(sin >> p);
    each (volume, volumes) assert(sin >> volume);
  }
  return {placements, volumes};
}

bool vibrate(const Problem& problem, Placements& placements, double rate)
{
  bool moved = false;
  for (int i = 0; i < placements.size(); ++i) {
    auto& p = placements[i];
    const double x = p.first + ((xorshift() % 20) - 10.0) * rate;
    const double y = p.second + ((xorshift() % 20) - 10.0) * rate;
    if (is_inside(problem, x, y) && checkPlacements2(x, y, placements, i)) {
      p.first = x;
      p.second = y;
      moved = true;
    }
  }
  return moved;
}

void pull_placements(const Problem& problem, Placements& placements, const vec<Volume>& volumes, double rate)
{
  for (int i = 0; i < placements.size(); ++i) {
    auto& p = placements[i];
    double vx = 0;
    double vy = 0;
    double vz = 0;
    each (a, problem.attendees) {
      auto dx = p.first - a.x;
      auto dy = p.second - a.y;
      auto d2 = dx * dx + dy * dy;
      auto speed = a.tastes[i] / d2;
      vx -= dx * speed;
      vy -= dy * speed;
    }
    double x = p.first + sqrt(vx * vx) * rate;
    double y = p.second + sqrt(vy * vy) * rate;
    if (is_inside(problem, x, y) && checkPlacements2(x, y, placements, i)) {
      p.first = x;
      p.second = y;
    }
  }
  return ;
}

pair<Placements, vec<Volume>> improve(const Problem& problem)
{
  const int M = problem.musicians.size();
  const int INST = set<int>(problem.musicians.begin(), problem.musicians.end()).size();
  // 90-shohei15-5-6702-9-hasi12.2673-hasi12.3454-hasi12.806.json
  auto [placements, volumes] = read_output(problem, "../../solutions/90-shohei15-2-9760-90-ueno5-hasi13.json.json-hasi12.309-hasi12.2673-hasi12.3454.json");
  auto musicianIndexesByInstrument = problem.musician_indexes_by_instrument();
  const auto cs = calcScore(problem, placements);
  assert(cs.first);
  vec<lli> bestScores = cs.second;
  lli currScore = sumScore(bestScores);
  lli bestScore = currScore;
  Placements currPlacements = placements;
  Placements bestPlacements = currPlacements;
  for (int loop = 0; loop < 1; ) {
    bool updated = false;
    for (int weight = 1; weight < 4; ++weight) {
      cerr << make_pair(loop, weight) << ": " << endl;
      while (!vibrate(problem, placements, 1.0 / weight)) ;
      for (int inst1 = 0; inst1 < INST; ++inst1) {
        for (int inst2 = inst1 + 1; inst2 < INST; ++inst2) {
          cerr << "# " << make_pair(inst1, inst2) << "" << endl;
          each (i, musicianIndexesByInstrument[inst1]) {
            // cerr << ">> " << make_pair(i, '_') << "" << endl;
            each (j, musicianIndexesByInstrument[inst2]) {
              pair<bool, vec<lli>> swapScore = calcSwapScore(problem, placements, i, j);
              cerr << ">> " << make_pair(i, j) << "" << endl;
              // cerr << swapScore << endl;
              assert(swapScore.first);
              lli diff = sumScore(swapScore.second);
              if (swapScore.first && 0 < diff) {
                // cerr << "updated:" << make_pair(i, j) << "+" << diff << ","<< currScore << "->" << currScore+diff << endl;
                currScore += diff;
                swapCache.clear();
                swappedCache.clear();
                blockedCache.clear();
                swap(placements[i], placements[j]);
              }
              if (bestScore < currScore) {
                bestScore = currScore;
                bestPlacements = currPlacements;
                updated = true;
              }
            }
          }
        }
      }
    }
    if (updated) {
      for (int i = 0; i < placements.size(); ++i) {
        if (0 < bestScores[i]) volumes[i] = 10.0;
        else                   volumes[i] =  0.0;
      }
      show_result(bestPlacements, volumes);
    }
  }
  return make_pair(bestPlacements, volumes);
}

int main(int argc, char** argv)
{
  Problem problem;
  cin >> problem;

  auto [placements, volumes] = improve(problem);
  // show_result(placements, volumes);
  cerr << calcScore(problem, placements) << endl;

  return 0;
}
