#ifndef __DEFS__
#define __DEFS__

#include <bits/stdc++.h>

using namespace std;

#define INF 1e32
typedef float cost_t;

string sol_name(string in_name, cost_t f_cost) {
    string out_name = in_name;
    out_name += ".f"s + to_string(f_cost) + ".sol";
    return out_name;
}

bool hasEnding (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

inline double toRadians(const double degree) { 
    return M_PI / 180 * degree; 
} 
  
double GeoDistance(double lat1, double long1, double lat2, double long2) { // in km
    lat1 = toRadians(lat1); 
    long1 = toRadians(long1); 
    lat2 = toRadians(lat2); 
    long2 = toRadians(long2); 
      
    // Haversine Formula 
    auto dlong = long2 - long1; 
    auto dlat = lat2 - lat1; 
  
    auto ans = pow(sin(dlat / 2), 2) + cos(lat1) * cos(lat2) * pow(sin(dlong / 2), 2); 
    ans = 2 * asin(sqrt(ans)); 
    ans = ans * 6371.009; 
    return ans; 
} 

struct Timer {
    chrono::time_point<chrono::steady_clock> ts;
    Timer(){
        Reset();
    }
    void Reset(){
        ts = chrono::steady_clock::now();
    }
    long long ElapsedMs(){
        return chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now()-ts).count();
    }
    void PrintTime(string tag){
        cerr<<"[T] "<<tag<<": "<<ElapsedMs()<<"ms"<<endl;
    }
};

struct Facility {
    int id;
    cost_t cost;
    vector<cost_t> dist;
};

class FLInstance {
    public:
    vector<Facility> facilities;
    // Facility/client count
    int f = 0, c = 0;

    virtual cost_t FCDistance(int fidx, int cidx) const = 0;
    // {
    //     return facilities[fidx].dist[cidx]; 
    // }

    // // Client to client distance apx using an intermediate facility
    // cost_t CCDistanceApx(int cidx_a, int cidx_b) const {
    //     if (cidx_a == cidx_b) return 0;
    //     cost_t dist = INF;
    //     for (int fid = 0; fid < f; fid++)
    //         dist = min(dist, FCDistance(fid, cidx_a)+FCDistance(fid, cidx_b));
    //     return dist;
    // }
    // The _literature_ format for the instances only contains the facility to client distances. To match that with our instances, we lose the client-client distances too. As long as all the clients are also facilities, this is equivalent (just slow).
    virtual cost_t CCDistance(int cidx_a, int cidx_b) const = 0;
    // {
    //     return CCDistanceApx(cidx_a,cidx_b);
    // }

    virtual pair<int, cost_t> ClosestFacility(int cid) const {
        pair<int,cost_t> closest = {-1, INF};
        for (int fid = 0; fid < f; fid++) {
            auto d = FCDistance(fid, cid);
            if (d < closest.second)
                closest = {fid, d}; 
        }
        return closest;
    }

    template <typename Iterable>
    pair<int, cost_t> ClosestFacility(int cid, const Iterable& fs) const {
        pair<int,cost_t> closest = {-1, INF};
        for (const int& fid : fs) {
            auto d = FCDistance(fid, cid);
            if (d < closest.second)
                closest = {fid, d}; 
        }
        return closest;
    }
    template <typename Iterable, typename Iterable2>
    cost_t Potential(const Iterable& cs, const Iterable2& fs) const {
        cost_t p = 0;
        for (const auto& c : cs) {
            p += ClosestFacility(c, fs).second;
        }
        return p;
    }
    template <typename Iterable>
    vector<int> Ball(int center_cid, cost_t radius, const Iterable& cs) const {
        vector<int> b;
        for (const auto& c : cs) {
            auto d = CCDistance(center_cid, c);
            if (d > radius) continue;
            b.push_back(c);
        }
        return b;
    }

    bool IsUniformFL() const {
        if (facilities.empty()) return false;
        auto c = facilities.begin()->cost;
        for (const auto& f : facilities)
            if (c != f.cost)
                return false;
        return true;
    }
    void SetFacilityCost(cost_t f_cost) {
        for (auto& f : facilities)
            f.cost = f_cost;
    }  
};

struct Point {
    cost_t x = 0, y = 0;

    static cost_t DistLP2(const Point& a, const Point& b) {
        return sqrt(pow(a.x-b.x,2)+pow(a.y-b.y,2));
    }
    static cost_t DistGeo(const Point& a, const Point& b) {
        return GeoDistance(a.x,a.y,b.x,b.y);
    }

    static function<cost_t(const Point&,const Point&)> Dist(const string& metric) {
        if (metric == "LP_2")
            return Point::DistLP2;
        else if (metric == "GEO")
            return Point::DistGeo;
        else {
            cout<<"Unknown metric "<<metric<<endl;
            exit(1);
        }
    }
};

struct PointGen {
    vector<cost_t> xs;

    static cost_t DistLP2(const PointGen& a, const PointGen& b) {
        assert(a.xs.size() == b.xs.size());
        cost_t r = 0;
        for (int i = 0; i < a.xs.size(); i++)
            r += pow(a.xs[i]-b.xs[i],2);
        return sqrt(r);
    }
    static function<cost_t(const PointGen&,const PointGen&)> Dist(const string& metric) {
        if (metric == "LP_2")
            return PointGen::DistLP2;
        else {
            cout<<"Unknown metric "<<metric<<endl;
            exit(1);
        }
    }
};

template<typename Point_t>
class FLFPts : public FLInstance {
    public:
    vector<Point_t> facilities_pts, clients_pts;
    function<cost_t(Point_t,Point_t)> distance_fn;

    FLFPts(const string& fname);

    int AddFacility(const Point_t& pt, cost_t cost = 1) {
        int id = facilities.size();
        facilities.emplace_back();
        facilities.back().id = id;
        facilities.back().cost = cost;
        facilities_pts.push_back(pt);
        f++;
        return id;
    }

    bool populated_facs = false;
    void PopulateFacs(bool force = false) {
        if (populated_facs && !force) return;
        cerr<<"Populating matrix:"<<endl;
        for (int i = 0; i < f; i++) {
            if (i % 100 == 0) {
                cerr<<"\r\tProgress: "<<i<<"/"<<f;
                fflush(stderr);
            }
            facilities[i].dist.resize(c);
            for (int j = 0; j < c; j++)
                facilities[i].dist[j] = distance_fn(facilities_pts[i], clients_pts[j]);
        }
        cerr<<"Populated"<<endl;
        populated_facs = true;
    }
    void ClearFacs() {
        populated_facs = false;
        for (int i = 0; i < f; i++)
            facilities[i].dist.clear();
    }

    virtual cost_t FCDistance(int fidx, int cidx) const {
        if (populated_facs) return facilities[fidx].dist[cidx];
        return distance_fn(facilities_pts[fidx], clients_pts[cidx]);
    }

    virtual cost_t CCDistance(int cidx_a, int cidx_b) const override {
        return distance_fn(clients_pts[cidx_a], clients_pts[cidx_b]);
    }

    pair<int, cost_t> ClosestFacility(const Point_t& pt) const {
        pair<int,cost_t> closest = {-1, INF};
        for (int fid = 0; fid < f; fid++) {
            auto d = distance_fn(facilities_pts[fid], pt);
            if (d < closest.second)
                closest = {fid, d}; 
        }
        return closest;
    }
};

template<>
FLFPts<Point>::FLFPts(const string& fname) {
    cerr<<"Using Point"<<endl;
    cerr<<"Loading "<<fname<<endl;
    ifstream fin(fname);
    if (!fin.good()) {
        throw "Unable to open file"s;
    }
    string line, metric;
    for (int i = 0; i < 1 && getline(fin, line); i++) {
        if (line.empty() || line[0] == '#') { i--; continue; }
        stringstream ss(line);
        ss>>f>>c;
    }
    cerr<<"\tF="<<f<<" | C="<<c<<endl;
    facilities.resize(f);
    facilities_pts.resize(f);
    clients_pts.resize(c);
    for (int i = 0; i < 1 && getline(fin, line); i++) {
        if (line.empty() || line[0] == '#') { i--; continue; }
        stringstream ss(line);
        ss>>metric;
    }
    for (int i = 0; i < f && getline(fin, line); i++) {
        if (line.empty() || line[0] == '#') { i--; continue; }
        stringstream ss(line);
        facilities[i].id = i;
        ss>>facilities[i].cost;
        ss>>facilities_pts[i].x>>facilities_pts[i].y;
    }
    for (int i = 0; i < c && getline(fin, line); i++) {
        if (line.empty() || line[0] == '#') { i--; continue; }
        stringstream ss(line);
        ss>>clients_pts[i].x>>clients_pts[i].y;
    }
    distance_fn = Point::Dist(metric);
}

template<>
FLFPts<PointGen>::FLFPts(const string& fname) {
    cerr<<"Using PointGen"<<endl;
    cerr<<"Loading "<<fname<<endl;
    ifstream fin(fname);
    if (!fin.good()) {
        cerr<<"Unable to open file"<<endl;
        exit(1);
    }
    string line, metric;
    for (int i = 0; i < 1 && getline(fin, line); i++) {
        if (line.empty() || line[0] == '#') { i--; continue; }
        stringstream ss(line);
        ss>>f>>c;
    }
    facilities.resize(f);
    facilities_pts.resize(f);
    clients_pts.resize(c);
    int dim;
    for (int i = 0; i < 1 && getline(fin, line); i++) {
        if (line.empty() || line[0] == '#') { i--; continue; }
        stringstream ss(line);
        ss>>metric>>dim;
    }
    for (int i = 0; i < f && getline(fin, line); i++) {
        if (line.empty() || line[0] == '#') { i--; continue; }
        stringstream ss(line);
        facilities[i].id = i;
        ss>>facilities[i].cost;
        facilities_pts[i].xs.resize(dim);
        for (auto& x : facilities_pts[i].xs)
            ss>>x;
    }
    for (int i = 0; i < c && getline(fin, line); i++) {
        if (line.empty() || line[0] == '#') { i--; continue; }
        stringstream ss(line);
        clients_pts[i].xs.resize(dim);
        for (auto& x : clients_pts[i].xs)
            ss>>x;
    }
    distance_fn = PointGen::Dist(metric);
}

// Input file for the _simple_ FL format. For online algos, the clients must be processed in 0..c-1 order
class FLFSimple : public FLInstance {
    public:

    FLFSimple(const string& fname) {
        cerr<<"Loading "<<fname<<endl;
        ifstream fin(fname);
        if (!fin.good()) {
            cerr<<"Unable to open file"<<endl;
            exit(1);
        }
        string line;
        int input_line = 0;
        while (getline(fin, line)) {
            if (line.empty()) continue;
            if (line[0] == '#') continue;
            stringstream ss(line);
            if (input_line == 0) {
                ss>>f>>c;
                facilities.resize(f);
            } else {
                int fid = input_line-1;
                ss>>facilities[fid].id>>facilities[fid].cost;
                facilities[fid].dist.resize(c);
                for (auto& d : facilities[fid].dist) ss>>d;
            }
            input_line++;
        };
        fin.close();
    }
};

template<typename T>
FLFPts<T>* RequireFLPts(const FLInstance* instance) {
    FLFPts<T>* instance_pts = dynamic_cast<FLFPts<T>*>(const_cast<FLInstance*>(instance));
    assert(instance_pts != nullptr);
    return instance_pts;
}

class Solution {
    public:
    string name;
    vector<pair<int,cost_t>> connected;
    vector<pair<int,cost_t>> opened;
    cost_t service_cost = 0, facility_cost = 0, cost = 0;
    bool fake = false;

    long long elapsed_ms = -1;

    Solution(){}
    Solution(cost_t service_cost, cost_t facility_cost) : service_cost(service_cost), facility_cost(facility_cost), cost(service_cost+facility_cost), fake(true) {}
    Solution(string fname){
        ifstream fin(fname);
        if (!fin.good()) {
            cerr<<"Unable to open file"<<endl;
            exit(1);
        }
        string line;
        int input_line = 0;
        while (getline(fin, line)) {
            if (line.empty()) continue;
            stringstream ss(line);
            if (line[0] == '#' && name.empty()) {
                ss>>name>>name>>cost>>facility_cost>>service_cost>>elapsed_ms>>elapsed_ms;
                continue;
            }
            if (input_line == 0) {
                int opened_size;
                ss>>opened_size;
                opened.resize(opened_size);
            } else {
                int fid = input_line-1;
                ss>>opened[fid].first>>opened[fid].second;
            }
            input_line++;
        };
    }
    void Connect(FLInstance* instance) {
        connected.clear();
        const auto opened_vec = GetOpened();
        cost_t service_cost_t = 0;
        for (int i = 0; i < instance->c; i++) {
            connected.push_back(instance->ClosestFacility(i, opened_vec));
            service_cost_t += connected.back().second;
        }
        if (service_cost > 0 && service_cost_t > 0)
            assert(abs(service_cost-service_cost_t)/min(service_cost,service_cost_t) < 1e-4);
    }
    Solution(const vector<int>& connected_, const FLInstance* instance, const string& name) : name(name) {
        set<int> opened_set;
        for (int c = 0; c < connected_.size(); c++) {
            const auto& f = connected_[c];

            auto fc_dist = instance->FCDistance(f,c);
            connected.push_back({f, fc_dist});
            service_cost += fc_dist;
            
            if (opened_set.count(f)) continue;
            auto f_cost = instance->facilities[f].cost;
            opened_set.insert(f);
            opened.push_back({f, f_cost});
            facility_cost += f_cost;
        }
        cost = facility_cost + service_cost;
    }
    vector<int> GetOpened() const {
        vector<int> opened_vec;
        for (const auto& x : opened)
            opened_vec.push_back(x.first);
        return opened_vec;
    }
    void Print(bool print_opened = false, bool print_connection = false) const {
        cerr<<"[+] "<<name<<" (rt: "<<elapsed_ms<<"ms)"<<endl;
        cerr<<"\tCost: "<<cost<<" (service: "<<service_cost<<", facility_c: "<<facility_cost<<")"<<endl;
        cerr<<"\tOpened facilities: "<<opened.size()<<endl;
        if (print_opened)
            for (const auto& x : opened)
                cerr<<"\t\t"<<x.first<<" (cost="<<x.second<<")"<<endl;
        if (print_connection) {
            cerr<<"\tConnections:"<<endl;
            int i = 0;
            for (int i = 0; i < connected.size(); i++)
                cerr<<"\t\t"<<i<<"->"<<connected[i].first<<"("<<connected[i].second<<")"<<endl;
        }
    }
    void Dump(const string& path) {
        ofstream fout(path);
        fout<<"#\t"<<name<<"\t"<<cost<<"\t"<<facility_cost<<"\t"<<service_cost<<"\t"<<opened.size()<<"\t"<<elapsed_ms<<endl;
        fout<<opened.size()<<endl;
        for (const auto& op : opened)
            fout<<op.first<<" "<<op.second<<"\n";
        fout.close();
    }
};

// Interface for FL algos, .run should be called only once (because rng status)
class FLAlgo {
    protected:
    std::default_random_engine rng;
    std::uniform_real_distribution<double> uar_unit{0.0, 1.0};
    double GetRandom() {
        return uar_unit(rng);
    }
    virtual Solution Solve(const FLInstance* instance) = 0;
    
    public:
    string name, name_tmp;
    Solution Run(const FLInstance* instance) {
        Timer timer;
        auto sol = Solve(instance);
        sol.elapsed_ms = timer.ElapsedMs();
        return sol;
    }
    void Reseed(int seed){
        rng.seed(seed);
    }
    void SaveName() {
        name_tmp = name;
    }
    void RestoreName() {
        name = name_tmp;
    }
};

class FLRandom : public FLAlgo {
    public:
    int n;
    FLRandom(int n) : n(n) {
        name = "random(" + to_string(n) + ")";
    }
    protected:
    virtual Solution Solve(const FLInstance* instance) override {
        set<int> lucky;
        while (lucky.size() < min(n,instance->f))
            lucky.insert(rand()%instance->f);
        vector<int> connected;
        vector<int> opened(lucky.begin(), lucky.end());
        for (int c = 0; c < instance->c; c++) {
            connected.push_back(instance->ClosestFacility(c, opened).first);
        }
        return Solution(connected, instance, name);
    }
};

#endif