import argparse, numpy as np, sys
from common import *
from dataclasses import dataclass
from collections import defaultdict
from tqdm import tqdm

fout_dir = "final_days/"
cc = "us"

def keep(row):
    if row.loc_id == "00000000000000000000000000000000": return False
    if not (-90 <= row.latitude <= 90): return False
    if not (-180 <= row.longitude <= 180): return False
    if row.city_country_code.lower() != cc: return False
    return True

def dump_pts(fout, locs, date):
    idx = 0
    indices = []
    users = set()
    while idx < len(locs):
        # if locs[idx].user not in users:
        users.add(locs[idx].user)
        indices.append(idx)
        idx += 1
    indices.sort(key=lambda x: locs[x].time)

    xy = [*map(lambda x: latlon_map(x.latitude, x.longitude), locs)]
    pts = [*map(lambda x: xy[x], indices)]
    
    import types
    args = types.SimpleNamespace()
    args.show_plot = False
    args.save_plot = True
    args.o = fout

    plot_2d(pts, args, swap=True)
    fl = FLFPts.from_pts(pts, pts, MetricNames.GEO.value, 1)
    fl.dump(fout, date)

for fin, start_date, fin_code in [("raw/loc-gowalla_totalCheckins.txt.enriched","2009-10-01","gw"),
                                  ("raw/loc-brightkite_totalCheckins.txt.enriched","2008-05-01","bk")]:

    locs = []
    with open(fin) as fin:
        for l in fin:
            if not l or l.startswith("#"): continue
            if l.split("\t")[1] == '': continue
            row = Chekin(*l.split("\t"))
            if not keep(row): continue
            locs.append(row)
    print("Filtered locs:", len(locs))
    locs.sort(key=lambda x: x.time)

    z = defaultdict(list)
    for c in locs:
        if c.time < start_date: continue
        z[c.time[:10]].append(c)
    k = sorted(z)
    for i,ik in tqdm(enumerate(k), total=len(k)):
        dump_pts("%s/%s-%s-D%s" % (fout_dir, fin_code, cc, i), z[ik], ik)

    z = defaultdict(list)
    day_idx = set()
    for c in locs:
        if c.time < start_date: continue
        z[(c.time[:10],int(c.time[11:13])//6)].append(c) # this?
        day_idx.add(c.time[:10])
    day_idx = {k:i for i,k in enumerate(sorted(day_idx))}
    for d,t in tqdm(z):
        dump_pts("%s/%s-%s-D%s_%s" % (fout_dir, fin_code, cc, day_idx[d], t), z[(d,t)], d)