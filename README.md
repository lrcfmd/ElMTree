# ElMTree
A web application in Flask to publicly share the results of [ElMTree](https://lmds.liverpool.ac.uk/ElMTree) and [ElM2D](https://lmds.liverpool.ac.uk/ElM2D) searches. 

This may be hosted privately, but note that we may not share data contained in the pickled indexing and LC database binaries, which are necessary for the application to work. Please generate these files yourself using the associated script in `ElMTree/ElMTree.py`.

## Installation

Recommended installation is via pip:

```
pip install ElMTreeIndex
```

Which may then be used as so:

```python

from ElMTree import ElMTree

elmtree = ElMTree(YOUR_LIST_OF_COMPOSITIONS_AS_STRINGS)
results = elmtree.knn("NaCl")
```
