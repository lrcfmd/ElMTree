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

from tqdm import tqdm

from tqdm.contrib.concurrent import process_map
from scipy.spatial.distance import euclidean

from ElMD import ElMD

from pymatgen.core import Structure, Composition
from ase import Atoms

import pandas as pd

def main():
    from matminer.datasets import load_dataset
    
    df = load_dataset("matbench_expt_gap").head(2000)

    compositions = [x for i, x in enumerate(df['composition'])]
    # compositions = [x[0] for x in pk.load(open('data/optimade_comps.pk', 'rb'))]

    # for composition, code in compositions:
    #     elmd_comp = ElMD(composition)
    #     elmd_compositions.append(elmd_comp)

    #     if dataset in elmtree_lookup[elmd_comp.pretty_formula]:
    #         elmtree_lookup[elmd_comp.pretty_formula][dataset].append(code)
    #     else:
    #         elmtree_lookup[elmd_comp.pretty_formula][dataset] = [code]

    # db_lookup[dataset]['experimental'] = True
    # db_lookup[dataset]['structures'] = False

    # if not os.path.exists('elmtree_lookup.pk'):
    #     pk.dump(elmtree_lookup, open('elmtree_lookup.pk', 'wb'))

    # if not os.path.exists('db_lookup.pk'):
    #     pk.dump(db_lookup, open('db_lookup.pk', 'wb'))

    df_names  = ['df_outputs_filtout_v1', 'df_outputs_filtout_v6', 'df_outputs_v5', 'df_outputs_v6']
    dfs = [pd.read_pickle('data/' + x + '.pkl') for x in df_names]

    re2fractive_comps, dataset_identifier = zip(*[(str(Structure.from_dict(x).composition), i) for i, df in enumerate(dfs) for x in df['structure']])
    optimade_comps = [x[0] for x in pk.load(open('data/optimade_comps.pk', 'rb'))][:1000]

    compositions = optimade_comps + list(re2fractive_comps)

    elmtree = ElMTree(compositions, verbose=True)

    x = elmtree.knn(compositions[0])
    print(x)

    pk.dump(elmtree, open("indexedElMTree.pk", "wb"))


class ElMTree():
    def __init__(self, 
                 input_compositions, # input_compositions to be indexed, assumed ElMD objects TODO generalize?
                 assigned_metric=ElMD("", metric="fast").elmd, 
                 centroid_ratio=32, 
                 on_disk=False,
                 lookup_tables=False,
                 elmtree_lookup_path=None,
                 db_lookup_path=None,
                 verbose=False,
                 max_workers=None,
                 pre_process=False):
        
        self.assigned_metric = assigned_metric
        self.centroid_ratio = centroid_ratio
        self.n = len(input_compositions)
        self.m = int(self.n / self.centroid_ratio)
        self.on_disk = on_disk 

        self.pre_process = pre_process
        self.verbose = verbose

        # Used for testing
        self.indexing_metric_counter = 0
        self.querying_metric_counter = 0

        if max_workers is None:
            max_workers = os.cpu_count()

        self.max_workers = max_workers

        input_compositions = self.process_input_formula(input_compositions)

        # Select m input_compositions to be our centres at random. 
        # Each list item stores the centre input_composition object, a list of children, 
        # and the covering radius
        if self.verbose: print("Selecting Voronoi centroids")

        self.centres = [[ElMD(input_composition), [], 0] for input_composition in random.sample(input_compositions, k=self.m)]
        # self.centres = [[input_composition, [], 0] for input_composition in random.sample(input_compositions, k=self.m)]
        
        self.lookup_tables = lookup_tables

        if self.verbose: print("Updating composition index lookup")

        if lookup_tables:
            # Written in separate scripts for the ElMTree to check which databases, and which type of db each one is
            self.elmtree_lookup = pk.load(open(elmtree_lookup_path, "rb")) # Returns the dbs a composition string can be found in 
            self.db_lookup = pk.load(open(db_lookup_path, "rb")) # Gives the metadata about whether the db is experimental or contains structural information
        
        elif isinstance(input_compositions[0], ElMD):
            # If there are no separate db lookups, simply index each composition by its position
            self.elmtree_lookup = {input_composition.pretty_formula: i for i, input_composition in enumerate(input_compositions)}

        else:
            if self.verbose: 
                self.elmtree_lookup = {ElMD(input_composition).pretty_formula: i for i, input_composition in tqdm(enumerate(input_compositions), total=len(input_compositions))}
            else:
                self.elmtree_lookup = {ElMD(input_composition).pretty_formula: i for i, input_composition in enumerate(input_compositions)}

        if self.verbose: print("Assigning compositions to centroids")

        assignments = process_map(self.get_centroid, input_compositions, chunksize=1000, max_workers=12)
        
        # Estimated in ~1128 hours single core
        # assignments = []
        # for x in tqdm(points):
        #     assignments.append(self.get_centroid(x)) 

        for input_composition, ind, distance in assignments:
            try:
                new_entry = self.make_entry(input_composition, distance)
                bisect.insort(self.centres[ind][1], new_entry)
            
                if distance > self.centres[ind][2]:
                    self.centres[ind][2] = distance

            except Exception as e:
                print(e)
                print(input_composition)
                continue

        if self.verbose: print("ElMTree Constructed")

    def process_input_formula(self, input_compositions):
        """Preprocess each input into an ElMD object, to save constructing these dynamically.
           This is significantly faster, but can take 100s of GBs RAM for larger datasets"""
        if self.pre_process is None and len(input_compositions) < 20000:
            self.pre_process = True

        if self.pre_process:
            if self.verbose: print("Preprocessing formula into ElMD objects")
            if isinstance(input_compositions[0], ElMD):
                return input_compositions
            elif isinstance(input_compositions[0], str):
                input_compositions = [ElMD(input_composition) for input_composition in input_compositions]
            elif isinstance(input_compositions[0], Structure):
                input_compositions = [ElMD(input_composition.composition) for input_composition in input_compositions]
            elif isinstance(input_compositions[0], Structure.composition):
                input_compositions = [ElMD(input_composition) for input_composition in input_compositions]
            elif isinstance(input_compositions[0], Atoms):
                input_compositions = [ElMD(input_composition.get_chemical_symbols()) for input_composition in input_compositions]

        else:
            if isinstance(input_compositions[0], ElMD):
                return input_compositions
            elif isinstance(input_compositions[0], str):
                return input_compositions
            elif isinstance(input_compositions[0], Structure):
                input_compositions = [input_composition.composition for input_composition in input_compositions]
            elif isinstance(input_compositions[0], Structure.composition):
                input_compositions = [input_composition for input_composition in input_compositions]
            elif isinstance(input_compositions[0], Atoms):
                input_compositions = [input_composition.get_chemical_symbols() for input_composition in input_compositions]

        return input_compositions

    def make_entry(self, input_composition, distance):
        if isinstance(input_composition, ElMD):
            db_entries = self.elmtree_lookup[input_composition.pretty_formula]
        else:
            db_entries = self.elmtree_lookup[ElMD(input_composition).pretty_formula]
        
        experimental = False
        structure = False

        if self.lookup_tables:
            for k, v in db_entries.items():
                if k == "compound_formula":
                    continue
                if self.db_lookup[k]["experimental"]:
                    experimental = True
                if self.db_lookup[k]["structures"]:
                    structure = True

        return Entry(input_composition, distance, experimental, structure)

    def get_centroid(self, input_composition):
        # if self.verbose: print("Computing distances to each centroid")
        input_composition = ElMD(input_composition)

        distances = [self.metric(input_composition, centre[0]) for centre in self.centres]
        ind = np.argsort(distances)[0]

        return (input_composition, ind, distances[ind])

    def metric(self, obj1, obj2, advanced_search=None):
        """Overload assigned metric function to allow more flexibility"""
        # For large operations assume that the routing object is also a
        # unique filename identifier
        if self.on_disk:
            obj1 = pk.load(open(self.db_folder + str(obj1), "rb"))
            obj2 = pk.load(open(self.db_folder + str(obj2), "rb"))

        if isinstance(obj1, str):
            obj1 = ElMD(obj1)

        elif isinstance(obj1, Structure):
            obj1 = ElMD(str(obj1.composition))

        elif isinstance(obj1, Composition):
            obj1 = ElMD(str(obj1))

        elif isinstance(obj1, Atoms):
            obj1 = ElMD(obj1.get_chemical_symbols())

        if isinstance(obj2, str):
            obj2 = ElMD(obj2)

        elif isinstance(obj2, Structure):
            obj2 = ElMD(str(obj2.composition))

        elif isinstance(obj2, Composition):
            obj2 = ElMD(str(obj2))

        elif isinstance(obj2, Atoms):
            obj2 = ElMD(obj2.get_chemical_symbols())   

        return self.assigned_metric(obj1, obj2)

    def get_matches(self, query, advanced_search=None):
        input_composition, i, query_radius = query 
        centre = self.centres[i]
        distance = self.metric(input_composition, centre[0])
        
        ret = []

        # If the cluster could contain a match
        if distance <= centre[2] + query_radius:
            for leaf in centre[1]:
                if advanced_search is not None and self.lookup_tables:
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
                    leaf_distance = self.metric(leaf.entry, input_composition)

                    if leaf_distance <= query_radius:
                        ret.append(self.make_entry(leaf.entry, leaf_distance))

                # Otherwise if query distance to centroid is less than the
                # leaf distance, all subsequent leaves must be even further 
                # away
                elif distance < leaf.distance:
                    break

        return ret

    def range_query(self, query, query_radius=2, advanced_search=None):
        if not isinstance(query, ElMD):
            query = ElMD(query)

        queries = [(query, i, query_radius) for i in range(len(self.centres))]
        
        cluster_matches = [self.get_matches(c) for c in queries]
            
        res = []

        for match_list in cluster_matches:
            for match in match_list:
                bisect.insort(res, match)

        return [(x.entry, x.distance) for x in res]

    def knn(self, query, k=200, advanced_search=None):
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
                    if advanced_search is not None and self.lookup_tables:
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
        
        if isinstance(NN[0], ElMD):
            return [(x.entry, x.distance, self.elmtree_lookup[x.entry.pretty_formula]) for x in NN[:k]]
        else:
            return [(x.entry, x.distance, self.elmtree_lookup[ElMD(x.entry).pretty_formula]) for x in NN[:k]]
        

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
