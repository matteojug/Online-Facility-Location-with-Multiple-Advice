#include <bits/stdc++.h>
#include "defs.cpp"
#include "mpi.cpp"

using namespace std;

int main(int argc, char* argv[]){
    assert(argc >= 2);
    int seed = time(NULL);
    cost_t f_cost = -1;
    if (argc > 2)
        f_cost = atof(argv[2]);
    if (argc > 3)
        seed = atoi(argv[3]);
    srand(seed);
    // FLAlgo* algo = new FLMpiInterface(UNCAP_FACILITY_LOCATION_TABU,"tabu");
    FLAlgo* algo = new FLMpiInterface(UNCAP_FACILITY_LOCATION_MYZ,"myz");

    #ifdef USE_PTS_GEN
    FLFPts<PointGen> *instance_pts = new FLFPts<PointGen>(argv[1]);
    #else
    FLFPts<Point> *instance_pts = new FLFPts<Point>(argv[1]);
    #endif

    if (instance_pts->f > INSTANCE_LIMIT) {
        cout<<"SKIPPED_TOO_BIG"<<endl;
        return 0;
    }

    instance_pts->PopulateFacs();

    // FLFSimple instance(argv[1]);
    FLInstance *instance = instance_pts;
    if (f_cost > 0)
        instance->SetFacilityCost(f_cost);
    
    string out_name = sol_name(argv[1], f_cost);

    Solution sol;
    sol = algo->Run(instance);
    sol.Print();
    sol.Dump(out_name);
    cout<<sol.name<<"\t"<<sol.cost<<"\t"<<sol.facility_cost<<"\t"<<sol.service_cost<<"\t"<<sol.opened.size()<<"\t"<<sol.elapsed_ms<<endl;

    // Solution sol2(out_name);
    // sol2.Connect(instance);
    // cout<<sol2.name<<"\t"<<sol2.cost<<"\t"<<sol2.facility_cost<<"\t"<<sol2.service_cost<<"\t"<<sol2.opened.size()<<"\t"<<sol2.elapsed_ms<<endl;
    // assert(sol2.connected == sol.connected);
    return 0;
}
