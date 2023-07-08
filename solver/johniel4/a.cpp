#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <complex>
#include <algorithm>
#include <numeric>

using namespace std;

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
    return abs(p2 - (p0 + (p0 - p1) * t)) < 5;
}

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace geometry = boost::geometry;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// https://kazuki-nagasawa.hatenablog.com/entry/memo_20140604_boost_rtree
typedef double              distance_t;
typedef std::vector<size_t> indexes_t;
typedef std::vector<double>          coordinates_t;
typedef std::vector< coordinates_t > coordinates_set_t;
namespace geometry = boost::geometry;
typedef geometry::model::point<distance_t, 2, geometry::cs::cartesian> point_t;
typedef geometry::model::box<point_t> box_t;
typedef geometry::model::segment<point_t> segment_t;
typedef std::pair<box_t, size_t> box_value_t;
typedef std::pair<point_t, size_t> point_value_t;
typedef std::pair<segment_t, pair<size_t, size_t>> segment_value_t;

typedef geometry::index::rtree<segment_value_t, geometry::index::quadratic<16>> rtree_segment_t;
typedef geometry::index::rtree<point_value_t, geometry::index::quadratic<16>> rtree_point_t;

template<typename T> using vec = vector<T>;

#define each(i, c) for (auto& i : c)
#define unless(cond) if (!(cond))
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
rtree_point_t rtp;
rtree_segment_t rts;
complex<double> rot(const complex<double> &p, double theta) { return p * polar(1.0, theta); }
bool isBlocked_RT(const Attendee& attendee, double x, double y, const vector<pair<double, double>>& placements)
{
  point_t m(x, y);
  point_t a(attendee.x, attendee.y);
  vector<pair<point_t, int>> nearest;
  rtp.query(bgi::nearest(segment_t(m, a), 1), back_inserter(nearest));
  if (nearest.empty()) return false;
  const int k = nearest[0].second;
  return isBlocked(attendee.x, attendee.y, x, y, placements[k].first, placements[k].second);
}

vec<pair<size_t, size_t>> listNewBlockings(
  const Problem& problem,
  double x, double y,
  const vector<pair<double, double>>& placements)
{
  vec<segment_value_t> nearest;
  const int M = problem.musicians.size();
  const int A = problem.attendees.size();

  point_t mn = point_t(x - 9, y - 9);
  point_t mx = point_t(x + 9, y + 9);
  rts.query(bgi::intersects(box_t(mn, mx)), back_inserter(nearest));
  vec<pair<size_t, size_t>> v;
  cerr << "nearest.size():" << nearest.size() << endl;
  each (i, nearest) {
    int a_index = i.second.first;
    int m_index = i.second.second;
    auto a = problem.attendees[a_index];
    auto p = placements[m_index];
    if (!isBlocked(a.x, a.y, x, y, p.first, p.second)) {
      v.push_back(i.second);
    }
  }
  return v;
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

bool checkPlacements(double x, double y, const vector<pair<double, double>>& placements) {
  vector<pair<point_t, int>> nearest;
  rtp.query(bgi::nearest(point_t(x, y), 1), back_inserter(nearest));
  if (nearest.empty()) return true;
  auto [x2, y2] = placements[nearest[0].second];
  if ((x - x2) * (x - x2) + (y - y2) * (y - y2) < 100) {
    return false;
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

double calcScoreDiff(
  const Problem& problem,
  const vector<pair<double, double>>& placements,
  const set<pair<unsigned, unsigned>>& blocking,
  int taste, double x, double y)
{
    double score = 0;

    // {{a index, m index}}
    vec<pair<size_t, size_t>> newblockings = listNewBlockings(problem, x, y, placements);
    each (b, newblockings) {
      int a_index = b.first;
      int m_index = b.second;
      auto& a = problem.attendees[a_index];
      if (!blocking.count({m_index, a_index})) {
        auto x1 = placements.at(m_index).first;
        auto y1 = placements.at(m_index).second;
        auto d2 = (a.x - x1) * (a.x - x1) + (a.y - y1) * (a.y - y1);
        score -= a.tastes[problem.musicians[m_index]] / d2;
      }
    }

    for (unsigned k = 0; k < problem.attendees.size(); ++ k) {
        auto& a = problem.attendees[k];
        if (!isBlocked_RT(a, x, y, placements)) {
            auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
            score += a.tastes[taste] / d2;
        }
    }
    return score;
}

// {score diff, new position}
pair<double, pair<double, double>> yama(
  const Problem& problem,
  int taste,
  pair<int, int> start,
  const vector<pair<double, double>>& placements,
  const set<pair<unsigned, unsigned>>& blocking
  ) {
    const double D = 1;
    const double dx[4] = {+D, -D, 0, 0};
    const double dy[4] = {0, 0, +D, -D};
    double xx = start.first;
    double yy = start.second;
    double scoreD = calcScoreDiff(problem, placements, blocking, taste, xx, yy);
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
                double s = calcScoreDiff(problem, placements, blocking, taste, x, y);
                if (s > scoreD) {
                    scoreD = s;
                    k = i;
                }
            }
        }
        if (k == -1) break;
        xx += dx[k];
        yy += dy[k];
    }
    return {scoreD, {xx, yy}};
}

vector<pair<double, double>> solve(const Problem& problem) {
  const unsigned N = problem.musicians.size();
    vector<pair<double, double>> res(N);
    vector<unsigned> perm(N);
    iota(perm.begin(), perm.end(), 0);
    random_shuffle(perm.begin(), perm.end());
    vector<pair<double, double>> placements;
    set<pair<unsigned, unsigned>> blocking;
    double bestScore = 0;
    int _ = 0;
    for (auto i : perm) {
      cerr << _++ << ' ' << perm.size() << make_pair(rtp.size(), rts.size()) << endl;
      int taste = problem.musicians[i];
        auto curr = yama(problem, problem.musicians[i], makeStart(problem, placements), placements, blocking);
        for (int tt = 0; tt < 5; ++ tt) {
          auto next = yama(problem, problem.musicians[i], makeStart(problem, placements), placements, blocking);
          setmax(curr, next);
        }
        auto [x2, y2] = curr.second;
        res[i] = {x2, y2};
        bestScore += curr.first;

        // cerr << _++ << ' ' << perm.size() << ' ' << curr << endl;

        for (int j = 0; j < problem.attendees.size(); ++j) {
          point_t m(x2, y2);
          point_t a(problem.attendees[j].x, problem.attendees[j].y);
          // cerr << j << ' ' << problem.attendees.size() << endl;
          if (!isBlocked_RT(problem.attendees[j], x2, y2, placements)) {
            pair<int, int> k = make_pair(j, placements.size());
            rts.insert(make_pair(segment_t(a, m), k));
          }
        }
        rtp.insert(make_pair(point_t(x2, y2), placements.size()));
        placements.push_back(res[i]);

        // {{a index, m index}}
        vec<pair<size_t, size_t>> newblockings = listNewBlockings(problem, x2, y2, placements);
        each (b, newblockings) blocking.insert({b.first, b.second});
    }
    return res;
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

    auto placement = solve(problem);
    cout << "{\"placements\":[";
    for (unsigned i = 0; i < placement.size(); i++) {
        if (i > 0) cout << ",";
        cout << "{\"x\":" << placement[i].first << ",\"y\":" << placement[i].second << "}";
    }
    cout << "]}" << endl;

    auto res = calcScore(problem, placement);
    if (!res.first) throw runtime_error("invalid placement");
    cerr << "score = " << res.second << endl;
}
