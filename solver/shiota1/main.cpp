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
typedef pair<int, int> pii
typedef vector<double> vd;


const double MUSICIAN_OFFSET = 10;

double roomW, roomH, stageW, stageH, stageX, stageY
int musicianN, tasteN;
vd musicianTaste(musicianN);
int attendeeN;
vi attendeeTasteWeight(attendeeN);
vector<Pos> attendeePos(attendeeN);

void input(){
    cin >> roomW >> roomH;
    cin >> stageW >> stageH;
    cin >> stageX >> stageY;

    cin >> musicianN >> tasteN;

    rep(i, musicianN){
        cin >> musicianTaste[i];
    }

    cin >> attendeeN;

    rep(i, attendeeN){
        cin >> attendeePos[i].x >> attendeePos[i].y;
        rep(j, tasteN){
            cin >> attendeeTasteWeight[j];
        }
    }
}

int main() {

    // offset
    stageX += MUSICIAN_OFFSET;
    stageY += MUSICIAN_OFFSET;
    stageH -= MUSICIAN_OFFSET;
    stageW -= MUSICIAN_OFFSET;

    //


    std::cout << "Hello, World!" << std::endl;
    return 0;
}
