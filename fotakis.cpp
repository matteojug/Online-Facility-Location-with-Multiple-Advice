#ifndef __FOTAKIS__
#define __FOTAKIS__

#include <bits/stdc++.h>
#include "defs.cpp"

using namespace std;

class FLFotakis : public FLAlgo {
    public:
    cost_t x;
    FLFotakis(cost_t x): x(x){
        name = "fotakis_det_uniform("s + to_string(x) + ")";
    }
    protected:
    virtual Solution Solve(const FLInstance* instance) override {
        assert(instance->IsUniformFL());
        cost_t f_cost = instance->facilities.begin()->cost;
        vector<int> connected, F;
        set<int> L;
        for (int w = 0; w < instance->c; w++) {
            L.insert(w);
            cost_t d = instance->ClosestFacility(w, F).second;
            cost_t r_w = d / x;
            vector<int> B_w = instance->Ball(w, r_w, L);
            cost_t B_w_pot = instance->Potential(B_w, F);
            if (B_w_pot >= f_cost) {
                int w_cap = w;
                if (d < f_cost) {
                    // Compute the potential of each point in the ball
                    map<int, cost_t> pots;
                    for (const auto& v : B_w)
                        pots[v] = instance->ClosestFacility(v, F).second;
                    // Find the ni value for each u in the ball
                    int ni_min = 0, ni_min_idx = -1;
                    for (const auto& u : B_w) {
                        // Compute the ball centered on u
                        vector<pair<int, cost_t>> u_ball;
                        for (const auto& v : B_w)
                            u_ball.push_back({v, instance->CCDistance(u,v)});
                        // Sort the ball by asc distance from u
                        sort(u_ball.begin(), u_ball.end(), [](pair<int,cost_t> a, pair<int,cost_t> b){return a.second < b.second;});
                        // Compute the potential of the u ball with radius r_w/2^0, and the size of such ball (u_ball_idx)
                        cost_t pot = 0, ni_pow = 1;
                        int ni = 0, u_ball_idx = 0;
                        for (; u_ball_idx < u_ball.size() && u_ball[u_ball_idx].second <= r_w / ni_pow; u_ball_idx++)
                            pot += pots[u_ball[u_ball_idx].first];
                        // Shrink the ball to find the max ni for which the u ball still has enough potential
                        while (pot > B_w_pot / 2 && u_ball_idx) {
                            ni++; ni_pow *= 2;
                            for (; u_ball_idx && u_ball[u_ball_idx-1].second > r_w / ni_pow; u_ball_idx--)
                                pot -= pots[u_ball[u_ball_idx-1].first];
                            // If there are no more points we can remove from the ball (eg multiple demands on the same point)
                            if (u_ball_idx && u_ball[u_ball_idx-1].second == 0)
                                break;
                        }
                        ni--;
                        // ni is the last val st the pot of the ball w/ radius r_w / 2**ni is more than half of B_w pot 
                        if (ni_min_idx == -1)
                            ni_min = ni+1;
                        if (ni < ni_min)
                            ni_min = ni,  ni_min_idx = u;
                    }
                    w_cap = ni_min_idx;
                }
                F.push_back(instance->ClosestFacility(w_cap).first);
                for (const auto& v : B_w)
                    L.erase(v);
            }
            connected.push_back(instance->ClosestFacility(w, F).first);
        }
        return Solution(connected, instance, name);
    }
};

#endif