#ifndef __MEYERSON__
#define __MEYERSON__

#include <bits/stdc++.h>
#include "defs.cpp"

using namespace std;

class FLMeyerson : public FLAlgo {
    public:
    float alpha;
    FLMeyerson(int seed, float alpha = 1) : alpha(alpha){
        name = "meyerson_rnd_uniform";
        if (alpha != 1)
            name += "("s + to_string(alpha) + ")";
        rng.seed(seed);
    }
    protected:
    virtual Solution Solve(const FLInstance* instance) override {
        assert(instance->IsUniformFL());
        cost_t f_cost = instance->facilities.begin()->cost;
        vector<int> connected;
        vector<int> opened;
        for (int c = 0; c < instance->c; c++) {
            cost_t dist = instance->ClosestFacility(c, opened).second;
            if (GetRandom() < dist/f_cost*alpha) {
                auto f = instance->ClosestFacility(c).first;
                assert(instance->FCDistance(f,c) < 1e-6); // Should open on the same spot
                opened.push_back(f);
            }
            connected.push_back(instance->ClosestFacility(c, opened).first);
        }
        return Solution(connected, instance, name);
    }
};

#endif