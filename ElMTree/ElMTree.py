"""
An implementation of the list of clusters indexing method for metric spaces
https://doi.org/10.1016/j.patrec.2004.11.014

"""
import os
import random
import bisect 
import pickle as pk

from dataclasses import dataclass
from collections import defaultdict

import numpy as np

from tqdm.contrib.concurrent import process_map
from scipy.spatial.distance import euclidean

from ElMD import ElMD

def main():
    from matminer.datasets import load_dataset
    
    df = load_dataset("matbench_expt_gap").head(2000)

    compositions = [(x, i) for i, x in enumerate(df['composition'])]

    elmtree_lookup = defaultdict(dict)
    db_lookup = defaultdict(dict)

    elmd_compositions = []

    dataset = 'matbench_expt_gap'
    
    for composition, code in compositions:
        elmd_comp = ElMD(composition)
        elmd_compositions.append(elmd_comp)

        if dataset in elmtree_lookup[elmd_comp.pretty_formula]:
            elmtree_lookup[elmd_comp.pretty_formula][dataset].append(code)
        else:
            elmtree_lookup[elmd_comp.pretty_formula][dataset] = [code]

    db_lookup[dataset]['experimental'] = True
    db_lookup[dataset]['structures'] = False

    if not os.path.exists('elmtree_lookup.pk'):
        pk.dump(elmtree_lookup, open('elmtree_lookup.pk', 'wb'))

    if not os.path.exists('db_lookup.pk'):
        pk.dump(db_lookup, open('db_lookup.pk', 'wb'))

    elmtree = ElMTree(elmd_compositions, ElMD("", metric="fast").elmd)

    pk.dump(elmtree, open("indexedElMTree.pk", "wb"))


class ElMTree():
    def __init__(self, points, assigned_metric=euclidean, centroid_ratio=32, on_disk=False):
        self.assigned_metric = assigned_metric
        self.centroid_ratio = centroid_ratio
        self.n = len(points)
        self.m = int(self.n / self.centroid_ratio)
        self.on_disk = on_disk 

        # Used for testing
        self.indexing_metric_counter = 0
        self.querying_metric_counter = 0

        # Select m points to be our centres at random. 
        # Each list item stores the centre point object, a list of children, 
        # and the covering radius
        self.centres = [[point, [], 0] for point in random.sample(points, k=self.m)]

        # Written in separate scripts for the ElMTree
        self.elmtree_lookup = pk.load(open("elmtree_lookup.pk", "rb")) # Returns the dbs a composition string can be found in 
        self.db_lookup = pk.load(open("db_lookup.pk", "rb")) # Gives the metadata about whether the db is experimental or contains structural information

        assignments = process_map(self.get_centroid, points, chunksize=100, max_workers=16)

        for point, ind, distance in assignments:
            try:
                new_entry = self.make_entry(point, distance)
                bisect.insort(self.centres[ind][1], new_entry)
            
                if distance > self.centres[ind][2]:
                    self.centres[ind][2] = distance

            except Exception as e:
                print(e)
                print(point)
                continue

        print()

    def make_entry(self, point, distance):
        db_entries = self.elmtree_lookup[point.pretty_formula]
        experimental = False
        structure = False

        for k, v in db_entries.items():
            if k == "compound_formula":
                continue
            if self.db_lookup[k]["experimental"]:
                experimental = True
            if self.db_lookup[k]["structures"]:
                structure = True

        return Entry(point, distance, experimental, structure)

    def get_centroid(self, point):
        distances = [self.metric(point, centre[0]) for centre in self.centres]
        ind = np.argsort(distances)[0]

        return (point, ind, distances[ind])

    def metric(self, obj1, obj2, advanced_search=None):
        """Overload the metric function to allow greater flexibility"""
        # For large operations assume that the routing object is also a
        # unique filename identifier
        if self.on_disk:
            obj1 = pk.load(open(self.db_folder + str(obj1), "rb"))
            obj2 = pk.load(open(self.db_folder + str(obj2), "rb"))

        else:
            return self.assigned_metric(obj1, obj2)

    def get_matches(self, query, advanced_search=None):
        point, i, query_radius = query 
        centre = self.centres[i]
        distance = self.metric(point, centre[0])
        
        ret = []

        # If the cluster could contain a match
        if distance <= centre[2] + query_radius:
            for leaf in centre[1]:
                if advanced_search is not None:
                    if advanced_search["structures"] and not leaf.structure:
                        continue

                    if advanced_search["experimental"] and not leaf.experimental:
                        continue

                    if len(advanced_search["must_contain"]) > 0 and not set(leaf.entry.normed_composition.keys()).issuperset(advanced_search["must_contain"]):
                        continue

                    if len(advanced_search["must_exclude"]) > 0 and not set(leaf.entry.normed_composition.keys()).isdisjoint(advanced_search["must_exclude"]):
                        continue
                       
                # If the known distances could form a valid triangle
                if abs(distance - leaf.distance) <= query_radius:
                    leaf_distance = self.metric(leaf.entry, point)

                    if leaf_distance <= query_radius:
                        ret.append(self.make_entry(leaf.entry, leaf_distance))

                # Otherwise if query distance to centroid is less than the
                # leaf distance, all subsequent leaves must be even further 
                # away
                elif distance < leaf.distance:
                    break

        return ret

    def range_query(self, query, query_radius=2, advanced_search=None):
        if not isinstance(query):
            query = ElMD(query)

        queries = [(query, i, query_radius) for i in range(len(self.centres))]
        
        cluster_matches = [self.get_matches(c) for c in queries]
            
        res = []

        for match_list in cluster_matches:
            for match in match_list:
                bisect.insort(res, match)

        return [(x.entry, x.distance) for x in res]

    def knn(self, query, k=5, advanced_search=None):
        if not isinstance(query, ElMD):
            query = ElMD(query)

        distances = [self.metric(query, centre[0]) for centre in self.centres]
        inds = np.argsort(distances)

        NN = []
        upper_bound = np.inf

        # This could be parallelisable... Trickier with dynamic NN queue though
        for ind in inds:
            centre = self.centres[ind]

            if ind == 0 or abs(distances[ind] - centre[2]) <= upper_bound:
                for leaf in centre[1]:
                    if advanced_search is not None:
                        if advanced_search["structures"] and not leaf.structure:
                            continue

                        if advanced_search["experimental"] and not leaf.experimental:
                            continue

                        if len(advanced_search["must_contain"]) > 0 and not set(leaf.entry.normed_composition.keys()).issuperset(advanced_search["must_contain"]):
                            continue

                        if len(advanced_search["must_exclude"]) > 0 and not set(leaf.entry.normed_composition.keys()).isdisjoint(advanced_search["must_exclude"]):
                            continue
                 
                    if abs(distances[ind] - leaf.distance) <= upper_bound:
                        leaf_distance = self.metric(leaf.entry, query)

                        if leaf_distance <= upper_bound :
                            bisect.insort(NN, self.make_entry(leaf.entry, leaf_distance))

                            if len(NN) > k:
                                upper_bound = NN[k].distance

                    elif distances[ind] < leaf.distance:
                        break

        return [(x.entry, x.distance) for x in NN[:k]]

@dataclass
class Entry:
    entry : "The indexed object"
    distance: "Numeric distance to the centroid or the query"
    experimental: "Whether the compound is in an experimental DB"
    structure: "Whether the compound has an associated structure"

    def __lt__(self, other):
        return self.distance < other.distance

if __name__=="__main__":
    main()
