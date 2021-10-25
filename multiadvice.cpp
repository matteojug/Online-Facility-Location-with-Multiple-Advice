#ifndef __MULTIADV__
#define __MULTIADV__

#include <bits/stdc++.h>
#include "defs.cpp"

using namespace std;

struct AdviceSolution {
    vector<int> F;
};
class FLAdvice : public FLAlgo {
    protected:
    vector<AdviceSolution> advices;

    public:
    FLAdvice() {}
    FLAdvice(int seed){
        Reseed(seed);
    }
    void AddAdvice(const AdviceSolution& advice) {
        advices.push_back(advice);
    }
    set<int> Suggested() {
        set<int> all;
        for (auto& advice : advices)
            for (auto& f : advice.F)
                all.insert(f);
        return all;
    }
    set<int> SuggestedDeduped(const FLInstance* instance) {
        auto suggested = Suggested();
        auto instance_pts = RequireFLPts<PTS_T>(instance);

        auto f_dist = [&](int a, int b) {
            return instance_pts->distance_fn(instance_pts->facilities_pts[a],instance_pts->facilities_pts[b]);
        };
        vector<int> suggested_v = {suggested.begin(), suggested.end()};
        set<int> dups;
        for (int i = 0; i < suggested_v.size(); i++)
            for (int j = i+1; j < suggested_v.size(); j++)
                if (f_dist(suggested_v[i], suggested_v[j]) == 0)
                    dups.insert(suggested_v[j]);
        for (auto& x : dups)
            suggested.erase(x);
        return suggested;
    }
    int Advices() {
        return advices.size();
    }
};

class FLAdviceTrust : public FLAdvice {
    public:
    FLAdviceTrust() {
        name = "trust";
    }
    protected:
    virtual Solution Solve(const FLInstance* instance) override {
        assert(instance->IsUniformFL());
        auto suggested = SuggestedDeduped(instance);
        assert(!suggested.empty());

        vector<int> connected;
        for (int c = 0; c < instance->c; c++) {
            auto best_f = instance->ClosestFacility(c, suggested);
            connected.push_back(best_f.first);
        }
        return Solution(connected, instance, name);
    }
};

class HST {
    public:
    struct Node {
        int id, i, leaves, leaf = -1;
        float w_edges;
        vector<shared_ptr<Node>> children;
        shared_ptr<Node> p = nullptr;
    };

    private:
    vector<int> v;
    float min_dist_orig = INF, max_dist_orig = 0, scaler;
    function<float(int,int)> dist_fn_orig;
    std::uniform_real_distribution<float> uar_unit{0.0, 1.0};

    float SampleBeta(default_random_engine& rng, float eps = 1e-5) {
        auto p = uar_unit(rng);
        float a = 1, b = 2, c;
        while (b-a > eps) {
            c = (a+b)/2;
            if (log(c)/log(2) > p)
                b = c;
            else
                a = c;
        }
        return (a+b)/2;
    }

    shared_ptr<Node> root = nullptr;
    int node_id = 0;
    unordered_map<int,vector<int>> pivots;
    unordered_map<int,shared_ptr<Node>> nodes;

    shared_ptr<Node> Partition(const vector<int> S, int i) {
        auto node = make_shared<Node>();
        node->id = node_id++;
        node->i = i+1;
        node->leaves = S.size();
        node->w_edges = pow(2,node->i);
        if (S.size() == 1) {
            nodes[S.front()] = node;
            node->leaf = S.front();
            return node;
        }

        unordered_map<int,vector<int>> sub;
        for (auto& x : S) 
            sub[pivots[x][i]].push_back(x);
        assert(sub.count(-1) == 0);
        for (auto& kv : sub) {
            auto child = Partition(kv.second, i-1);
            child->p = node;
            node->children.push_back(child);
        }
        
        return node;
    }

    public:
    HST(const vector<int>& v_, const function<float(int,int)>& dist_fn_orig, default_random_engine& rng, bool verbose=false) : v(v_), dist_fn_orig(dist_fn_orig) {
        assert(v.size() > 1);
        shuffle(v.begin(), v.end(), rng);

        float beta = SampleBeta(rng);
        for (int i = 0; i < v.size(); i++)
            for (int j = i+1; j < v.size(); j++) {
                auto d = dist_fn_orig(v[i],v[j]);
                assert(d != 0);
                max_dist_orig = max(max_dist_orig, d);
                min_dist_orig = min(min_dist_orig, d);
            }
        // the smallest distance is strictly more than 1
        scaler = (1+1e-6) / min_dist_orig;
        if (verbose)
        cerr<<"hst: v:"<<v.size()<<", dists:"<<min_dist_orig<<","<<max_dist_orig<<" -> "<<min_dist_orig*scaler<<","<<max_dist_orig*scaler<<endl;

        int delta = ceil(log(max_dist_orig*scaler)/log(2));
        if (verbose) cerr<<"hst delta:"<<delta<<endl;

        for (auto& x : v) {
            auto& p = pivots[x];
            p.resize(delta);
            int curr = 0;
            for (int i = delta-1; i >= 0; i--) {
                float beta_i = pow(2,i-1)*beta;
                while (curr < v.size() && dist_fn_orig(x, v[curr])*scaler > beta_i)
                    curr++;
                if (curr == v.size()) {
                    p[i] = -1;
                    break;
                }
                p[i] = curr;
            }
        }
        root = Partition(v, delta-1);
    }
    
    float Dist(int a, int b) {
        if (v.size() == 1) return 0;

        // auto real_d = dist_fn_orig(a,b);
        float d = 0;
        shared_ptr<Node> pa = nodes[a], pb = nodes[b];
        while (pa->i < pb->i) {
            pa = pa->p;
            d += pa->w_edges;
        }
        while (pa->i > pb->i) {
            pb = pb->p;
            d += pb->w_edges;
        }
        while (pa->id != pb->id) {
            pa = pa->p;
            d += pa->w_edges;

            pb = pb->p;
            d += pb->w_edges;
        }
        d /= scaler;
        // cout<<a<<","<<b<<" => rd:"<<real_d<<" vs d:"<<d<<endl;
        return d;
    }
    pair<float,float> Distortion() {
        float max_d = 0, avg_d = 0;
        int avg_dc = 0;
        for (int i = 0; i < v.size(); i++)
            for (int j = i+1; j < v.size(); j++) {
                auto d = Dist(v[i],v[j]) / dist_fn_orig(v[i],v[j]);
                max_d = max(max_d, d);
                avg_d += d;
                avg_dc++;
            }
        avg_d /= avg_dc;
        // cerr<<"HST: max_d="<<max_d<<", avg_d="<<avg_d<<endl;
        return {max_d, avg_d};
    }
    shared_ptr<Node> Root() {
        return root;
    }
    shared_ptr<Node> Leaf(int x) {
        return nodes[x];
    }
    pair<int,float> ClosestTo(shared_ptr<Node> subtree, int target) {
        // cout<<subtree->id<<"("<<subtree->i<<") vs "<<target<<endl;
        if (subtree->leaf != -1)
            return {subtree->leaf, Dist(subtree->leaf, target)};
        pair<int,float> closest = {-1,INF};
        for (auto& c : subtree->children) {
            // cout<<c->id<<"("<<c->i<<")"<<endl;
            auto d = ClosestTo(c, target);
            if (d.second < closest.second)
                closest = d;
        }
        return closest;
    }
};

class FLAdviceHST : public FLAdvice {
    public:
    int hst_count;
    FLAdviceHST(int seed, int hst_count = 1) : FLAdvice(seed), hst_count(hst_count) {
        name = "advice_hst";
        if (hst_count != 1) {
            if (name.back() == ')') name.back() = ',';
            else name += '(';
            name += to_string(hst_count) + ")";
        }
    }
    protected:
    virtual Solution Solve(const FLInstance* instance) override {
        assert(instance->IsUniformFL());
        cost_t f_cost = instance->facilities.begin()->cost;
        auto suggested = SuggestedDeduped(instance);
        assert(!suggested.empty());
        if (suggested.size() == 1) {
            // single node hst is a pain, just revert to trusty
            vector<int> connected;
            for (int c = 0; c < instance->c; c++) {
                auto best_f = instance->ClosestFacility(c, suggested);
                connected.push_back(best_f.first);
            }
            return Solution(connected, instance, name);
        }

        auto instance_pts = RequireFLPts<PTS_T>(instance);
        auto f_dist = [&](int a, int b) {
            return instance_pts->distance_fn(instance_pts->facilities_pts[a],instance_pts->facilities_pts[b]);
        };
        
        cerr<<"hst size:"<<suggested.size()<<" | instance size:"<<instance->c<<endl;
        Timer timer;
        HST hst({suggested.begin(),suggested.end()}, f_dist, rng);
        auto hst_distortion = hst.Distortion();
        for (int i = 1; i < hst_count; i++) {
            HST hst_t({suggested.begin(),suggested.end()}, f_dist, rng);
            auto hst_distortion_t = hst_t.Distortion();
            if (hst_distortion_t < hst_distortion) {
                hst = std::move(hst_t);
                hst_distortion = hst_distortion_t;
            }
        }
        cerr<<"hst_distortion:"<<hst_distortion.first<<","<<hst_distortion.second<<endl;
        timer.PrintTime("Got HST");

        unordered_map<int, int> node_opened;
        unordered_map<int, float> node_potential;

        vector<int> connected;
        set<int> opened;
        for (int c = 0; c < instance->c; c++) {
            // cout<<"fc:"<<c<<endl;
            auto q_p = instance->ClosestFacility(c, suggested).first;
            auto node = hst.Leaf(q_p);
            set<int> root_path;
            while (node->p != nullptr) {
                root_path.insert(node->id);
                node = node->p;
            }
            assert(node == hst.Root());
            // cout<<"\trootpath size:"<<root_path.size()<<endl;
            while (node->leaf == -1) {
                // cout<<"\tnode:"<<node->id<<"("<<node->i<<")"<<endl;
                // no open nodes
                if (!node_opened[node->id]) {
                    shared_ptr<HST::Node> next;
                    for (auto& c : node->children) {
                        // cout<<c->id<<endl;
                        auto closest = hst.ClosestTo(c, q_p);
                        if (closest.second > f_cost/3.0) continue;
                        if (next == nullptr || c->leaves > next->leaves)
                            next = c;
                    }
                    assert(next != nullptr);
                    node = next;
                    continue;
                }
                // opened node in q_p subtree (x)
                shared_ptr<HST::Node> x;
                for (auto& c : node->children) {
                    if (!root_path.count(c->id)) continue;
                    x = c;
                    break;
                }
                if (x != nullptr && node_opened[x->id]) {
                    node = x;
                    continue;
                }
                // find y (closest child != x)
                cost_t y_dist = INF;
                for (auto& q : opened)
                    y_dist = min(y_dist, hst.Dist(q_p,q));
                node_potential[x->id] += y_dist;
                if (node_potential[x->id] >= f_cost) {
                    node = x;
                    continue;
                }
                break;
            }
            // cout<<"outloop:"<<node->leaf<<endl;

            // ended in a leaf
            if (node->leaf != -1) {
                opened.insert(node->leaf);
                do {
                    node_opened[node->id]++;
                    node = node->p;
                } while (node != nullptr);
            }
            
            connected.push_back(instance->ClosestFacility(c, opened).first);
        }
        return Solution(connected, instance, name);
    }
};

class FLMahdian : public FLAlgo {
    vector<FLAlgo*> algos;
    public:
    FLMahdian(const vector<FLAlgo*>& algos, const string& algo_name = "") : algos(algos) {
        name = "mahdian(";
        if (algo_name != "")
            name += algo_name + ")";
        else {
            for (auto& a : algos)
                name += a->name + ",";
            name[name.size()-1] = ')';
        }
    }
    protected:
    struct SubAlgo {
        Solution sol;
        cost_t cost;
        set<int> opened;
    };
    virtual Solution Solve(const FLInstance* instance) override {
        vector<SubAlgo> sols;
        for (auto& a : algos)
            sols.push_back({a->Run(instance), 0, {}});

        vector<int> connected;
        for (int c = 0; c < instance->c; c++) {
            for (auto& sol : sols) {
                auto f = sol.sol.connected[c];
                sol.cost += f.second;
                if (sol.opened.count(f.first)) continue;
                sol.opened.insert(f.first);
                sol.cost += instance->facilities[f.first].cost;
            }
            int best = 0;
            for (int i = 1; i < sols.size(); i++)
                if (sols[i].cost < sols[best].cost)
                    best = i;
            connected.push_back(sols[best].sol.connected[c].first);
        }
        
        // for (int i = 0; i < algos.size(); i++)
        //     cout<<algos[i]->name<<":"<<sols[i].cost<<" ("<<sols[i].opened.size()<<endl;
        
        return Solution(connected, instance, name);
    }
};

class FLMahdian2 : public FLAlgo {
    FLAlgo *algo_a, *algo_b;
    float gamma;

    public:
    FLMahdian2(FLAlgo *a, FLAlgo *b, float gamma = 2) : algo_a(a), algo_b(b), gamma(gamma) {
        name = "mahdian["s + to_string(gamma) + "](" + a->name + "," + b->name + ")";
    }
    protected:
    struct SubAlgo {
        Solution sol;
        cost_t cost;
        set<int> opened;
    };
    virtual Solution Solve(const FLInstance* instance) override {
        vector<SubAlgo> sols = {
            {algo_a->Run(instance), 0, {}},
            {algo_b->Run(instance), 0, {}}
        };

        vector<int> connected;
        for (int c = 0; c < instance->c; c++) {
            for (auto& sol : sols) {
                auto f = sol.sol.connected[c];
                sol.cost += f.second;
                if (sol.opened.count(f.first)) continue;
                sol.opened.insert(f.first);
                sol.cost += instance->facilities[f.first].cost;
            }
            int best;
            if (sols[0].cost <= (gamma-1)*sols[1].cost)
                best = 0;
            else
                best = 1;
            connected.push_back(sols[best].sol.connected[c].first);
        }
        
        // for (int i = 0; i < algos.size(); i++)
        //     cout<<algos[i]->name<<":"<<sols[i].cost<<" ("<<sols[i].opened.size()<<endl;
        
        return Solution(connected, instance, name);
    }
};

#endif