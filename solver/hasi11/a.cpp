#include <iostream>
#include <fstream>
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

long long calcScore11(const Problem& problem, double x, double y, int taste) {
    long long score = 0;
    for (auto& a : problem.attendees) {
        auto d2 = (a.x - x) * (a.x - x) + (a.y - y) * (a.y - y);
        score += (long long)ceil(1000000 * a.tastes[taste] / d2); // NOLINT(cppcoreguidelines-narrowing-conversions)
    }
    return score;
}

int main(int argc, char** argv) {
    if (argc < 1) {
        cerr << "Usage: " << argv[0] << " <input file>" << endl;
        return 1;
    }

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

    cout << problem.musicians.size() << endl;
    int y = (int)problem.stageBottom + 10;
    for (int x = (int)problem.stageLeft + 10; x <= problem.stageLeft + problem.stageWidth - 10; ++ x) {
        for (auto taste : problem.musicians) {
            cout << calcScore11(problem, x, y, taste) << " ";
        }
        cout << endl;
    }
}
