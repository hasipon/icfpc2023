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

vec<pair<size_t, size_t>> listBlockings(const Problem& problem, double x, double y, const vector<pair<double, double>>& placements)
{
  vec<segment_value_t> nearest;
  rts.query(bgi::nearest(point_t(x, y), 1 << 29), back_inserter(nearest));
  vec<pair<size_t, size_t>> v;
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

pair<double, double> yama(const Problem& problem, int taste, double x0, double y0, const vector<pair<double, double>>& placements) {
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

vector<pair<double, double>> solve(const Problem& problem) {
  vector<pair<double, double>> res(problem.musicians.size());
  vector<unsigned> perm(problem.musicians.size());
  for (unsigned i = 0; i < perm.size(); ++ i) perm[i] = i;
  random_shuffle(perm.begin(), perm.end());
  vector<pair<double, double>> placements;
  for (auto i : perm) {
    auto [x0, y0] = makeStart(problem, placements);
    res[i] = yama(problem, problem.musicians[i], x0, y0, placements);

    rtp.insert(make_pair(point_t(res[i].first, res[i].second), placements.size()));
    for (int j = 0; j < problem.attendees.size(); ++j) {
      point_t m(res[i].first, res[i].second);
      point_t a(problem.attendees[j].x, problem.attendees[j].y);
      pair<int, int> k = make_pair(j, placements.size());
      rts.insert(make_pair(segment_t(a, m), k));
    }
    placements.push_back(res[i]);
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
