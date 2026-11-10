"""Build an ElMTree search index from a CSV of compositions.

For the production index (all the databases the original web search used),
see ``convert-index``. This builds an index from your own data:

    uv run build-index data/sample_compositions.csv \\
        --database matbench_expt_gap --experimental --out src/model/artifacts

Input: a CSV with a ``composition`` column and, optionally, an ``id`` column
(record identifiers, shown in search results as ``<database>: <ids>``).
Compositions that normalise to the same formula are indexed once, with all
their ids.

Writes ``elmtree_index.db`` (the tool's runtime asset) to --out.
"""

import argparse
import csv
from pathlib import Path

from ElMD import ElMD

from model.index import INDEX_FILE, sparse_vector, write_index


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("csv", type=Path, help="CSV with a 'composition' column")
    parser.add_argument("--database", default="database",
                        help="name to file the ids under in search results")
    parser.add_argument("--description", default="", help="description of the database")
    parser.add_argument("--display", default="", help="name shown in results (default: --database)")
    parser.add_argument("--experimental", action="store_true",
                        help="the database is experimental (for experimental_only searches)")
    parser.add_argument("--structures", action="store_true",
                        help="the database has crystal structures (for structures_only searches)")
    parser.add_argument("--out", type=Path, default=Path("src/model/artifacts"))
    args = parser.parse_args(argv)

    entries: dict[str, dict] = {}
    with args.csv.open(newline="") as f:
        for i, row in enumerate(csv.DictReader(f)):
            formula = (row.get("composition") or "").strip()
            if not formula:
                continue
            elmd = ElMD(formula, metric="fast")
            if not elmd.pretty_formula:
                print(f"skipping unparseable composition on row {i}: {formula!r}")
                continue
            if elmd.pretty_formula not in entries:
                positions, ratios = sparse_vector(elmd)
                entries[elmd.pretty_formula] = {
                    "formula": formula,
                    "normalised_formula": elmd.pretty_formula,
                    "positions": positions,
                    "ratios": ratios,
                    "records": {args.database: []},
                }
            entries[elmd.pretty_formula]["records"][args.database].append(row.get("id") or str(i))

    databases = {args.database: {"display": args.display or args.database,
                                 "description": args.description,
                                 "experimental": args.experimental,
                                 "structures": args.structures, "url": ""}}
    out = args.out / INDEX_FILE
    write_index(out, list(entries.values()), databases, source=f"build-index {args.csv.name}")
    print(f"wrote {out}: {len(entries)} distinct compositions")


if __name__ == "__main__":
    main()
