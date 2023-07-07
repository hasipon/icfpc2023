#include <iostream>

#define REP(i,b,n) for(int i=b;i<(int)n;i++)
#define rep(i,n)   REP(i,0,n)
#define dbg(x) cout << __LINE__ << ' ' << #x << " = " << (x) << endl

using namespace std;

class Pos {
public:
    double x, y;
};

typedef long long ll;
typedef vector<double> vd;
typedef vector<int> vi;


const double MUSICIAN_OFFSET = 10;

class Problem {
public:
    double roomW, roomH, stageW, stageH, stageX, stageY;
    int musicianN, tasteN;
    vi musicianTaste;
    int attendeeN;
    vector<vd> attendeeTasteWeight;
    vector<Pos> attendeePos;
};


Problem input() {

    Problem p;

    cin >> p.roomW >> p.roomH;
    cin >> p.stageW >> p.stageH;
    cin >> p.stageX >> p.stageY;

    cin >> p.musicianN >> p.tasteN;
    p.musicianTaste = vi(p.musicianN);

    rep(i, p.musicianN){
        cin >> p.musicianTaste[i];
    }

    cin >> p.attendeeN;
    p.attendeeTasteWeight = vector<vd>(p.attendeeN);
    p.attendeePos = vector<Pos>(p.attendeeN);

    rep(i, p.attendeeN){
        cin >> p.attendeePos[i].x >> p.attendeePos[i].y;
        p.attendeeTasteWeight[i] = vd(p.tasteN);
        rep(j, p.tasteN){
            cin >> p.attendeeTasteWeight[i][j];
            cout << p.attendeeTasteWeight[i][j] << ' ';
        }
    }
    return p;
}

int main() {

    Problem p = input();
    // offset
    p.stageX += MUSICIAN_OFFSET;
    p.stageY += MUSICIAN_OFFSET;
    p.stageH -= MUSICIAN_OFFSET;
    p.stageW -= MUSICIAN_OFFSET;

    return 0;
}
