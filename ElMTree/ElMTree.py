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
    
    df = load_dataset("matbench_expt_gap")

    elmtree = ElMTree(df['composition'], verbose=True, pre_process=True)

    x = elmtree.knn(df['composition'][0])[:5]
    print(x)

    pk.dump(elmtree, open("indexedElMTree.pk", "wb"))


class ElMTree():
    def __init__(self, 
                 input_compositions, # input_compositions to be indexed, assumed ElMD objects TODO generalize?
                 assigned_metric=ElMD("", metric="fast").elmd, 
                 centroid_ratio=32, 
                 k=50,
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
        self.k = k
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
            # Written in separate scripts for the ElMTree to check which databases, and which type of db each one is for advanced searches
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

            return [self.convert_to_ElMD(x) for x in tqdm(input_compositions)]

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
        input_composition = ElMD(input_composition)

        distances = [self.metric(input_composition, centre[0]) for centre in self.centres]
        ind = np.argsort(distances)[0]

        return (input_composition, ind, distances[ind])
    
    def convert_to_ElMD(self, object):
        if isinstance(object, ElMD):
            return object
        elif isinstance(object, str):
            object = ElMD(object)

        elif isinstance(object, Structure):
            object = ElMD(str(object.composition))

        elif isinstance(object, Composition):
            object = ElMD(str(object))

        elif isinstance(object, Atoms):
            object = ElMD(object.get_chemical_symbols())

        else:
            raise Exception("Incorrect input type, must be one of ElMD, str, pymatgen.core.structure.Structure, pymatgen.core.structure.Composition, or ase.Atoms")
        
        return object


    def metric(self, obj1, obj2, advanced_search=None):
        """Overload assigned metric function to allow more flexibility"""
        # For large operations assume that the routing object is a
        # unique filename identifier
        if self.on_disk:
            obj1 = pk.load(open(self.db_folder + str(obj1), "rb"))
            obj2 = pk.load(open(self.db_folder + str(obj2), "rb"))

        obj1 = self.convert_to_ElMD(obj1)
        obj2 = self.convert_to_ElMD(obj2) 

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

    def knn(self, query, k=None, advanced_search=None):
        if k is None:
            k = self.k
        if isinstance(query, list):
            if self.verbose: print("Mapping knn query to processes for list of formalae input")

            return process_map(self.knn, query, chunksize=1000)
        
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
