# ElMTree
An indexing class that is used in the backend of [ElMTree](https://lmds.liverpool.ac.uk/ElMTree) and [ElM2D](https://lmds.liverpool.ac.uk/ElM2D) searches. A flask application is also provided in `ElMTree/app` for hosting a simplified interface for custom searches.

This may be hosted privately, but note that we may not share data contained in the pickled indexing and ElMTree database binaries, which are necessary for the application to work. Please generate these files yourself using the script provided in `main()` in the file `ElMTree/ElMTree.py`.

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
