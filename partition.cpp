#ifndef __PARTITION__
#define __PARTITION__

#include <bits/stdc++.h>
#include "defs.cpp"

using namespace std;

struct FLPartitionRoot {
    
    struct Quadrant {
        Point center;
        cost_t side;

        cost_t cost = 0;
        Point facility;

        Quadrant *parent = nullptr, *subq[4] = {nullptr};

        bool IsOpen() {
            return subq[0] != nullptr;
        }
        void Open(const Point& p) {
            facility = p;
            subq[0] = new Quadrant{{center.x-side/4, center.y-side/4}, side/2}; // UL
            subq[1] = new Quadrant{{center.x+side/4, center.y-side/4}, side/2}; // UR
            subq[2] = new Quadrant{{center.x-side/4, center.y+side/4}, side/2}; // DL
            subq[3] = new Quadrant{{center.x+side/4, center.y+side/4}, side/2}; // DR
            for (int i = 0; i < 4; i++) subq[i]->parent = this;
        } 
        Quadrant* Locate(const Point& p) {
            if (!IsOpen()) return this;
            if (p.x <= center.x && p.y <= center.y) return subq[0]->Locate(p);
            if (p.x >= center.x && p.y <= center.y) return subq[1]->Locate(p);
            if (p.x <= center.x && p.y >= center.y) return subq[2]->Locate(p);
            if (p.x >= center.x && p.y >= center.y) return subq[3]->Locate(p);
            assert(false);
        }
    } *root;

    cost_t f = -1, a;
    function<cost_t(Point,Point)> distance_fn;
    void Init(const Point& center, cost_t f_, cost_t a_, const function<cost_t(Point,Point)>& distance_fn_) {
        f = f_;
        a = a_;
        distance_fn = distance_fn_;
        root = new Quadrant{center, f / static_cast<cost_t>(sqrt(2))};
    }
    
    bool Serve(const Point& client) {
        if (!root->IsOpen()) {
            root->Open(client);
            return true;
        }
        auto q = root->Locate(client);
        cost_t dist = distance_fn(client, q->parent->facility);
        for (Quadrant* parent_it = q->parent; parent_it != nullptr; parent_it = parent_it->parent)
            dist = min(dist, distance_fn(client, parent_it->facility));
        q->cost += dist;
        if (q->cost >= a*f) {
            q->Open(client);
            return true;
        }
        return false;
    }
};

struct FLPartitionRootCollection {
    Point center;
    cost_t f, a;

    cost_t side;
    map<pair<int,int>, FLPartitionRoot> collection;

    void Init() {
        side = f / sqrt(2);
    }

    pair<int,int> GetRoot(const Point& client) {
        return {round((client.x-center.x)/side), round((client.y-center.y)/side)};
    }
    Point GetCenter(const Point& client) {
        auto root = GetRoot(client);
        return {center.x + root.first*side, center.y + root.second*side};
    }

    bool Serve(const Point& client) { // Return true if opened there
        auto& root = collection[GetRoot(client)];
        if (root.f < 0) root.Init(GetCenter(client), f, a, Point::DistLP2);
        return root.Serve(client);
    }

};

// Treat everything as in R^2
class FLPartition : public FLAlgo {
    public:
    cost_t a;
    FLPartition(cost_t a) : a(a){
        name = "partition("s + to_string(a) + ")";
    }

    protected:
    virtual Solution Solve(const FLInstance* instance) override {
        assert(instance->IsUniformFL());
        cost_t f_cost = instance->facilities.begin()->cost;

        auto instance_pts = RequireFLPts<Point>(instance);

        FLPartitionRootCollection roots{instance_pts->clients_pts[0], f_cost, a};
        roots.Init();

        vector<int> connected;
        vector<int> opened;
        for (int c = 0; c < instance_pts->c; c++) {
            if (roots.Serve(instance_pts->clients_pts[c])) {
                auto f = instance->ClosestFacility(c).first;
                opened.push_back(f);
            }
            connected.push_back(instance->ClosestFacility(c, opened).first);
        }
        return Solution(connected, instance, name);
    }
};

#endif