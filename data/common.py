from dataclasses import dataclass
from collections.abc import Iterable
from functools import partial
from typing import Any
from enum import Enum
import matplotlib.pyplot as plt
import numpy as np

def latlon_map(lat, lon):
    return lat, lon

def plot_2d(pts, args, swap=False):
    if not (args.show_plot or args.save_plot): return
    plt.scatter(*[*zip(*pts)][:2][::-1 if swap else 1])
    if args.save_plot:
        plt.savefig(args.o + ".png")
    if args.show_plot:
        plt.show()
    plt.close()

def fmt_comment(comment):
    return ''.join(map(lambda x: "# %s\n" % x, comment.splitlines()))

def fmt_pt(pt):
    return ' '.join(map(str, pt))

@dataclass
class FLFUncap:
    n: int # number of facilities
    m: int # number of demands
    f_cost: Any # (float[n]) cost of opening
    f_distances: Any # (float[n][m]) service cost
    def __post_init__(self):
        assert(self.n == len(self.f_cost))
        assert(self.n == len(self.f_distances))
        for d in self.f_distances:
            assert(self.m == len(d))
    
    def dump(self, fname, comment=None):
        with open(fname,"w") as fout:
            if comment:
                fout.write(fmt_comment(comment))
            fout.write("%d %d 0\n" % (self.n, self.m))
            for i in range(self.n):
                fout.write("%d %f %s\n" % (i+1, self.f_cost[i], ' '.join(map(str, self.f_distances[i]))))
    
    @staticmethod
    def from_pts(f, d, metric_fn, f_cost):
        if isinstance(f_cost, Iterable):
            fc = [x for x in f_cost]
        else:
            fc = [f_cost]*len(f)
        fd = [[metric_fn(ff,dd) for dd in d] for ff in f]
        return FLFUncap(len(f), len(d), fc, fd)

@dataclass
class FLFPts:
    f: Any # facilities points coord
    d: Any # demands points coord
    metric_name: Any # metric
    f_cost: Any
    dim: int

    def __post_init__(self):
        assert(len(self.f) == len(self.f_cost))
        for x in self.f: assert(len(x) == self.dim)
        for x in self.d: assert(len(x) == self.dim)
    
    def dump(self, fname, comment=None):
        with open(fname,"w") as fout:
            if comment:
                fout.write(fmt_comment(comment))
            fout.write("%d %d 1\n" % (len(self.f), len(self.d)))
            fout.write("%s %d\n" % (self.metric_name, self.dim))
            for c,x in zip(self.f_cost, self.f): fout.write("%d %s\n" % (c, fmt_pt(x)))
            for x in self.d: fout.write("%s\n" % fmt_pt(x))
    
    @staticmethod
    def from_pts(f, d, metric_name, f_cost):
        if isinstance(f_cost, Iterable):
            fc = [x for x in f_cost]
        else:
            fc = [f_cost]*len(f)
        dim = len(f[0])
        return FLFPts(f, d, metric_name, fc, dim)

import pyproj

class MetricNames(Enum):
    GEO = "GEO"
    LP_1 = "LP_1"
    LP_2 = "LP_2"

class Metric:
    G = pyproj.Geod(ellps='WGS84')
    @staticmethod
    def lp_(a, b, p=2):
        s = 0
        for aa,bb in zip(a,b):
            s += abs(aa-bb)**p
        return s**(1.0/p)
    @staticmethod
    def lp(p=2):
        return partial(Metric.lp_, p=p)
    
    def geo(a, b):
        return Metric.G.inv(a[1], a[0], b[1], b[0])[2]


@dataclass
class City:
    geonameid: int
    name: str
    asciiname: str
    alternatenames: str
    latitude: float
    longitude: float
    feature_class: str
    feature_code: str
    country_code: str
    cc2: str
    admin1_code: str
    admin2_code: str
    admin3_code: str
    admin4_code: str
    population: int
    elevation: int
    dem: str
    timezone: str
    modification_date: str
    def __post_init__(self):
        self.geonameid = int(self.geonameid)
        self.latitude = float(self.latitude)
        self.longitude = float(self.longitude)
        self.population = int(self.population)
        self.country_code = self.country_code.lower()

@dataclass
class Chekin:
    user: str
    time: str
    latitude: float
    longitude: float
    loc_id: str
    city_id: int = -1
    city_country_code: str = "??"
    def __post_init__(self):
        self.latitude = float(self.latitude)
        self.longitude = float(self.longitude)
        self.city_country_code = self.city_country_code.strip()