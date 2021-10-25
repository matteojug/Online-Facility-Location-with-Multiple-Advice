#include <bits/stdc++.h>

#ifdef USE_PTS_GEN
#define PTS_T PointGen
#else
#define PTS_T Point
#endif

#include "defs.cpp"
#include "mpi.cpp"
#include "meyerson.cpp"
#include "fotakis.cpp"
#include "partition.cpp"
#include "multiadvice.cpp"

using namespace std;

map<string,FLAlgo*> algos;
vector<string> algos_as_oracle;

#define ALPHA_STEPS 10
#define MAHDIAN_GAMMAS {2., 1.75, 1.5, 1.25, 1.}
// #define MAHDIAN_GAMMAS {2.}
#define HST_COUNT {1, 5, 10, 50, 100}
// #define HST_COUNT {1}

void InitAlgos(int seed){
    auto algos_offline = vector<FLAlgo*>({
        new FLRandom(1),
        new FLRandom(10),
        new FLRandom(100),
        new FLRandom(1000),
    });
    auto algos_online = vector<FLAlgo*>({
        // new FLMeyerson(seed),
        new FLFotakis(0.1),
        new FLFotakis(0.5),
        new FLFotakis(1),
        new FLFotakis(5),
        new FLFotakis(10),
        
        #ifndef USE_PTS_GEN
        new FLPartition(0.001),
        new FLPartition(0.005),
        new FLPartition(0.01),
        new FLPartition(0.05),
        new FLPartition(0.1),
        new FLPartition(0.5),
        new FLPartition(1),
        new FLPartition(5),
        #endif
    });

    int steps = ALPHA_STEPS;
    for (int i = 0; i < steps; i++) {
        float fi = (i+1)/static_cast<float>(steps);
        algos_online.push_back(new FLMeyerson(seed, fi));
    }
    
    for (const auto& a : algos_offline)
        algos[a->name] = a,
        algos_as_oracle.push_back(a->name);
    for (const auto& a : algos_online)
        algos[a->name] = a;
    assert(algos_offline.size()+algos_online.size() == algos.size());
    // algos.clear(); algos_as_oracle.clear(); // TODO: this
}

Solution induced_solution(FLInstance* instance, const vector<int>& opened) {
    vector<int> connected;
    for (int c = 0; c < instance->c; c++) {
        connected.push_back(instance->ClosestFacility(c, opened).first);
    }
    auto s = Solution(connected, instance, "???");
    return s;
}

void print_solution_log(const string& name, const Solution& sol) {
    cout<<name<<"\t"<<sol.cost<<"\t"<<sol.facility_cost<<"\t"<<sol.service_cost<<"\t"<<sol.opened.size()<<"\t"<<sol.elapsed_ms<<endl;
}
void print_solution_log(const Solution& sol) {
    print_solution_log(sol.name, sol);
}

struct OracleCollection {
    string name;
    vector<Solution> sols;
    Solution merged;
};
vector<OracleCollection> Oracles(FLInstance* instance, const Solution& sol_offline, int seed) {
    vector<OracleCollection> oracles;

    oracles.push_back({"offline", {sol_offline}});

    srand(seed);
    auto sol_random = FLRandom(sol_offline.opened.size()).Run(instance);
    oracles.push_back({"random", {sol_random}});

    // vector<Solution> tmp;
    // for (auto& f : sol_offline.GetOpened())
    //     tmp.push_back(induced_solution(instance, {f}));
    // oracles.push_back({"offline-singleton", std::move(tmp)});
    
    vector<int> mixed, offline = sol_offline.GetOpened(), random = sol_random.GetOpened();
    for (auto& f : sol_offline.GetOpened()) mixed.push_back(f);
    for (auto& f : sol_random.GetOpened()) mixed.push_back(f);

    for (auto& kv : vector<pair<string,vector<int>>>{
        {"offline", offline},
        {"random", random},
        {"mixed", mixed},
    })
        for (int cnt : {1, 2, 3, 4, 5}) {
        // for (int cnt : {1, 2}) {
            random_shuffle(kv.second.begin(), kv.second.end());
            if (kv.second.size() < cnt) continue;
            vector<Solution> tmp;
            for (int i = 0; i < cnt; i++) {
                set<int> opened;
                for (int j = i; j < kv.second.size(); j += cnt)
                    opened.insert(kv.second[j]);
                tmp.push_back(induced_solution(instance, {opened.begin(), opened.end()}));
            }
            oracles.push_back({kv.first+"-k"s+to_string(cnt), std::move(tmp)});
        }

    return oracles;
}

int main(int argc, char* argv[]){
    assert(argc >= 2);
    int seed = time(NULL);
    cost_t f_cost = -1;
    string sorting = "none";
    if (argc > 2)
        f_cost = atof(argv[2]);
    if (argc > 3)
        seed = atoi(argv[3]);
    if (argc > 4)
        sorting = argv[4];
    
    srand(seed);
    InitAlgos(seed);

    #ifdef USE_PTS_GEN
    FLFPts<PointGen> *instance_pts = new FLFPts<PointGen>(argv[1]);
    // Sorting just the clients_pts
    if (sorting == "none") {
        // Done!
    } else {
        cout<<"Unknown sorting crit: "<<sorting<<endl;
        return 1;
    }
    #else
    FLFPts<Point> *instance_pts = new FLFPts<Point>(argv[1]);
    // Sorting just the clients_pts
    if (sorting == "x")
        sort(instance_pts->clients_pts.begin(), instance_pts->clients_pts.end(), [](Point a, Point b) { return a.x < b.x || a.x == b.x && a.y < b.y; });
    else if (sorting == "orig") {
        Point center;
        for (const auto& p : instance_pts->clients_pts)
            center.x += p.x, center.y += p.y;
        center.x /= instance_pts->c;
        center.y /= instance_pts->c;
        sort(instance_pts->clients_pts.begin(), instance_pts->clients_pts.end(), [center](Point a, Point b) { return Point::DistLP2(center,a) < Point::DistLP2(center,b); });
    } else if (sorting == "rand") {
        std::default_random_engine rng;
        rng.seed(seed);
        shuffle(instance_pts->clients_pts.begin(), instance_pts->clients_pts.end(), rng);
    } else if (sorting == "none") {
        // Done!
    } else {
        cout<<"Unknown sorting crit: "<<sorting<<endl;
        return 1;
    }
    #endif
    
    // map<string, Solution> sols;

    instance_pts->PopulateFacs();
    if (f_cost > 0)
        instance_pts->SetFacilityCost(f_cost);
    
    // FLFSimple instance(argv[1]);
    FLInstance *instance = instance_pts;
    vector<OracleCollection> oracles;

    #ifdef SKIP_OFFLINE_SOL
    cerr<<"Skipping offline sol loading"<<endl;
    #else
    // First offline solution, loaded from file
    // As long as the facilities are not scrambled, this should work.
    Solution sol_offline(sol_name(argv[1], f_cost));
    sol_offline.Connect(instance);
    sol_offline.Print();
    print_solution_log(sol_offline);

    oracles = Oracles(instance, sol_offline, seed);

    // sols[sol_offline.name] = sol_offline;
    // algos_as_oracle.push_back(sol_offline.name);
    
    // // Add random with same facility count as the offline
    // auto algo_rnd_offline_count = new FLRandom(sol_offline.opened.size());
    // algo_rnd_offline_count->name = "random(offline_size)";
    // algos[algo_rnd_offline_count->name] = algo_rnd_offline_count,
    // algos_as_oracle.push_back(algo_rnd_offline_count->name);
    #endif

    for (auto& oracle_coll : oracles) {
        set<int> merged;
        for (int i = 0; i < oracle_coll.sols.size(); i++) {
            print_solution_log(oracle_coll.name+"#"s+to_string(i), oracle_coll.sols[i]);
            for (auto& f : oracle_coll.sols[i].GetOpened())
                merged.insert(f);
        }
        oracle_coll.merged = induced_solution(instance, {merged.begin(), merged.end()});
        print_solution_log(oracle_coll.name+"#*"s, oracle_coll.merged);
    }

    Solution sol;
    for (const auto& x : algos) {
        srand(seed);
        sol = x.second->Run(instance);
        sol.Print();
        print_solution_log(sol);
    }
    
    Solution sol_tmp;
    #define RUN_AND_LOG(algo) \
        sol_tmp = algo.Run(instance); \
        sol_tmp.Print(); \
        print_solution_log(sol_tmp);

    
    for (auto& oracle_coll : oracles) {
        if (oracle_coll.merged.fake) continue;
        cerr<<oracle_coll.name<<endl;

        FLAdviceTrust trusty_merged;
        vector<FLAlgo*> trusty_disjoint;
        vector<FLAdviceHST> hsts;
        for (int hst_count : HST_COUNT)
            hsts.emplace_back(seed, hst_count);
            
        for (auto& oracle : oracle_coll.sols) {
            if (oracle.fake) continue;
            AdviceSolution oracle_sol{oracle.GetOpened()};

            trusty_merged.AddAdvice(oracle_sol);
            for (auto& h : hsts)
                h.AddAdvice(oracle_sol);
            
            auto trusty_ptr = new FLAdviceTrust();
            trusty_ptr->AddAdvice(oracle_sol);
            trusty_disjoint.push_back(trusty_ptr);
        }

        trusty_merged.name += "(merge("s+oracle_coll.name+"))";
        RUN_AND_LOG(trusty_merged)

        for (auto& h : hsts) {
            h.Reseed(seed);
            h.name += "("s+oracle_coll.name+")";
            RUN_AND_LOG(h);
        }

        for (float gamma : MAHDIAN_GAMMAS) {
            FLMahdian2 mahdian(new FLMeyerson(seed), &trusty_merged, gamma);
            RUN_AND_LOG(mahdian)

            mahdian = FLMahdian2(&trusty_merged, new FLMeyerson(seed), gamma);
            RUN_AND_LOG(mahdian)


            FLMahdian mahdian_all(trusty_disjoint, oracle_coll.name);
            mahdian = FLMahdian2(new FLMeyerson(seed), &mahdian_all, gamma);
            RUN_AND_LOG(mahdian)

            mahdian_all = FLMahdian(trusty_disjoint, oracle_coll.name);
            mahdian = FLMahdian2(&mahdian_all, new FLMeyerson(seed), gamma);
            RUN_AND_LOG(mahdian)


            for (auto& h : hsts) {
                h.Reseed(seed);
                mahdian = FLMahdian2(new FLMeyerson(seed), &h, gamma);
                RUN_AND_LOG(mahdian)

                h.Reseed(seed);
                mahdian = FLMahdian2(&h, new FLMeyerson(seed), gamma);
                RUN_AND_LOG(mahdian)
            }
        }

        FLMahdian mahdian(trusty_disjoint, oracle_coll.name);
        RUN_AND_LOG(mahdian)

        mahdian = FLMahdian({trusty_disjoint.rbegin(), trusty_disjoint.rend()}, oracle_coll.name + "_reversed");
        RUN_AND_LOG(mahdian)

        auto meyerson = new FLMeyerson(seed);
        trusty_disjoint.push_back(meyerson);
        mahdian = FLMahdian(trusty_disjoint, meyerson->name + "," + oracle_coll.name);
        RUN_AND_LOG(mahdian)
        trusty_disjoint.pop_back();

        meyerson = new FLMeyerson(seed);
        trusty_disjoint.push_back(meyerson);
        mahdian = FLMahdian({trusty_disjoint.rbegin(), trusty_disjoint.rend()}, meyerson->name + "," + oracle_coll.name + "_reversed");
        RUN_AND_LOG(mahdian)
    }

    return 0;
}