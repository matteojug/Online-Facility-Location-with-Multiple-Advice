#ifndef __MPI__
#define __MPI__

#include <bits/stdc++.h>
#include "mpi-sw/facloc_JMS/facloc_JMS.h"
#include "mpi-sw/facloc_LOCAL/facloc_LOCAL.h"
#include "defs.cpp"

using namespace std;

class FLMpiInterface : public FLAlgo {
    public:
    typedef float mpi_ftype;
    typedef function<bool(mpi_ftype*,mpi_ftype*,int,int,int*,mpi_ftype&)> mpi_fn;
    
    mpi_fn fn;
    FLMpiInterface(const mpi_fn& fn, const string& name_):fn(fn) {
        name = name_;
    };

    protected:
    virtual Solution Solve(const FLInstance* instance) override {
        auto open_cost = new mpi_ftype[instance->f];
        cerr<<"Instantiated open_cost"<<endl;
        auto cost_matrix = new mpi_ftype[instance->f*instance->c];
        cerr<<"Instantiated cost_matrix"<<endl;
        auto connected = new int[instance->c];
        for (int i = 0; i < instance->f; i++) {
            if (i % 50 == 0) {
                cerr<<"\r\tProgress: "<<i<<"/"<<instance->f;
                fflush(stderr);
            }
            open_cost[i] = instance->facilities[i].cost;
            for (int j = 0; j < instance->c; j++)
                cost_matrix[i*instance->c+j] = instance->FCDistance(i,j);
        }
        cerr<<"Populated open_cost"<<endl;
        mpi_ftype cost;
        auto res = fn(open_cost, cost_matrix, instance->f, instance->c, connected, cost);
        cerr<<"Solved by mpi"<<endl;
        delete[] open_cost, cost_matrix;
        vector<int> connected_vec(connected, connected+instance->c);
        delete[] connected;
        auto sol = Solution(connected_vec, instance, name);
        auto acc_err = ((mpi_ftype)abs(sol.cost - cost)) / max(sol.cost,cost);
        if (!(acc_err < 1e-4)) {
            cout<<sol.cost<<" vs "<<cost<<endl;
            cout<<acc_err<<endl;
            assert(acc_err < 1e-4);
        }
        return sol;
    }
};

#endif