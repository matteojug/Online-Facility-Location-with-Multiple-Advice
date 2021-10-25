import argparse, numpy as np, sys
from common import *
from dataclasses import dataclass
from collections import defaultdict
from tqdm import tqdm
from dateutil import parser

@dataclass
class Pickup:
    time: str
    latitude: float
    longitude: float
    base: str
    def __post_init__(self):
        self.time = str(parser.parse(self.time.replace("\"","")))
        self.base = self.base.replace("\"","")
        self.latitude = float(self.latitude)
        self.longitude = float(self.longitude)

fout_dir = "final_uber/"

locs = []
print("Loading csvs")
for f in tqdm("uber-raw-data-apr14.csv uber-raw-data-may14.csv uber-raw-data-jun14.csv uber-raw-data-jul14.csv uber-raw-data-aug14.csv uber-raw-data-sep14.csv".split()):
    with open(f"raw/uber-pickups-in-new-york-city/{f}") as fin:
        hdr = True
        for l in fin:
            if hdr:
                hdr = False
                continue
            row = Pickup(*l.split(","))
            locs.append(row)
    # break
print("Filtered locs:", len(locs))
locs.sort(key=lambda x: x.time)
    
def dump_pts(fout, locs, date):
    pts = [*map(lambda x: latlon_map(x.latitude, x.longitude), locs)]
    
    import types
    args = types.SimpleNamespace()
    args.show_plot = False
    args.save_plot = True
    args.o = fout

    plot_2d(pts, args, swap=True)
    fl = FLFPts.from_pts(pts, pts, MetricNames.GEO.value, 1)
    fl.dump(fout, date)

z = defaultdict(list)
for c in tqdm(locs):
#     if c.time < start_date: continue
    z[c.time[:10]].append(c)
k = sorted(z)
print("Dumping pts")
for i,ik in tqdm(enumerate(k), total=len(k)):
    dump_pts(f"{fout_dir}/uber-D{i}", z[ik], ik)

z = defaultdict(list)
day_idx = set()
for c in locs:
    # if c.time < start_date: continue
    z[(c.time[:10],int(c.time[11:13])//6)].append(c) # this?
    day_idx.add(c.time[:10])
day_idx = {k:i for i,k in enumerate(sorted(day_idx))}
for d,t in tqdm(z):
    dump_pts(f"{fout_dir}/uber-D{day_idx[d]}_{t}", z[(d,t)], d)