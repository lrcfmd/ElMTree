"""Build ElMTree search-index files from a CSV of compositions.

The production index files (indexed_LC.pk + mtree_lookup.pk) are not public.
This builds files in the same format from your own data, as the original
README's ``main()`` did:

    uv run build-index data/sample_compositions.csv \\
        --database matbench_expt_gap --out src/model/artifacts

Input: a CSV with a ``composition`` column and, optionally, an ``id`` column
(record identifiers, shown in search results as ``<database>: <ids>``).

Writes to --out:
  indexed_LC.pk    pickled ElMTree index (the required asset)
  mtree_lookup.pk  {pretty_formula: {database: [ids]}} (the optional asset)
"""

import argparse
import csv
import pickle
import random
from collections import defaultdict
from pathlib import Path

from ElMD import ElMD

from ElMTree import ElMTree


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("csv", type=Path, help="CSV with a 'composition' column")
    parser.add_argument("--database", default="database",
                        help="name to file the ids under in the lookup")
    parser.add_argument("--out", type=Path, default=Path("src/model/artifacts"))
    parser.add_argument("--seed", type=int, default=0,
                        help="seed for the random choice of cluster centres")
    parser.add_argument("--max-workers", type=int, default=None)
    args = parser.parse_args(argv)

    lookup: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    compositions = []
    with args.csv.open(newline="") as f:
        for i, row in enumerate(csv.DictReader(f)):
            formula = (row.get("composition") or "").strip()
            if not formula:
                continue
            pretty = ElMD(formula, metric="fast").pretty_formula
            if not pretty:
                print(f"skipping unparseable composition on row {i}: {formula!r}")
                continue
            if pretty not in lookup:
                compositions.append(formula)  # index each distinct composition once
            lookup[pretty][args.database].append(row.get("id") or str(i))

    print(f"{len(compositions)} distinct compositions")
    random.seed(args.seed)
    tree = ElMTree(compositions, pre_process=True, verbose=True,
                   max_workers=args.max_workers)

    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "indexed_LC.pk").open("wb") as f:
        pickle.dump(tree, f)
    with (args.out / "mtree_lookup.pk").open("wb") as f:
        pickle.dump({k: dict(v) for k, v in lookup.items()}, f)
    print(f"wrote {args.out / 'indexed_LC.pk'} and {args.out / 'mtree_lookup.pk'}")


if __name__ == "__main__":
    main()
