"""Convert the original ElMTree web search's lookup tables to a search index.

The original server (lmds.liverpool.ac.uk/ElMTree) searched one index of every
composition in its databases, with two lookup tables written offline:

  mtree_lookup.pk  {normalised formula: {database: [ids], ...,
                    "compound_formula": [formulas as given by the databases]}}
  db_lookup.pk     {database: {"description", "experimental", "structures", "url"}}

Its pickled index (indexed_LC.pk) was built from each lookup entry's first
compound formula. Those two tables therefore hold everything needed to rebuild
the same index, without the 5.4 GB pickle:

    uv run convert-index mtree_lookup.pk db_lookup.pk --out src/model/artifacts

Writes ``elmtree_index.db`` (the tool's runtime asset) to --out.
"""

import argparse
import os
import pickle
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from ElMD import ElMD

from model.index import INDEX_FILE, sparse_vector, write_index

MAX_FORMULA_TRIES = 20  # compound formulas tried per entry before falling back


class _DataUnpickler(pickle.Unpickler):
    """Loads plain data only: the lookup tables are dicts/lists of str/int."""

    def find_class(self, module, name):
        raise pickle.UnpicklingError(f"unexpected object {module}.{name} in lookup table")


def _load(path: Path):
    with path.open("rb") as f:
        return _DataUnpickler(f).load()


def _clean(formula) -> str:
    return str(formula).strip().strip("'\"").strip()


def _vector(task):
    """Vector for one lookup entry, from the formula the original index used.

    The lookup keys are rounded to 3 dp, so the compound formula gives the
    exact composition; fall back to the key if none parses to it."""
    key, compound_formulas = task
    for formula in compound_formulas:
        formula = _clean(formula)
        try:
            elmd = ElMD(formula, metric="fast")
        except Exception:
            continue
        if elmd.pretty_formula == key:
            return key, formula, *sparse_vector(elmd), False
    elmd = ElMD(key, metric="fast")
    return key, key, *sparse_vector(elmd), True


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("mtree_lookup", type=Path, help="the original mtree_lookup.pk")
    parser.add_argument("db_lookup", type=Path, help="the original db_lookup.pk")
    parser.add_argument("--out", type=Path, default=Path("src/model/artifacts"))
    parser.add_argument("--max-workers", type=int, default=os.cpu_count())
    args = parser.parse_args(argv)

    t0 = time.time()
    lookup = _load(args.mtree_lookup)
    databases = _load(args.db_lookup)
    print(f"loaded {len(lookup)} lookup entries and {len(databases)} databases "
          f"in {time.time() - t0:.0f} s")

    tasks = [(key, value.get("compound_formula", [])[:MAX_FORMULA_TRIES])
             for key, value in lookup.items() if key.strip()]
    skipped = len(lookup) - len(tasks)

    entries, fallbacks, empty = [], 0, 0
    with ProcessPoolExecutor(args.max_workers) as pool:
        for key, formula, positions, ratios, fallback in pool.map(_vector, tasks, chunksize=2000):
            if len(positions) == 0:
                empty += 1
                continue
            fallbacks += fallback
            records = {db: ids for db, ids in lookup[key].items() if db != "compound_formula"}
            entries.append({"formula": formula, "normalised_formula": key,
                            "positions": positions, "ratios": ratios, "records": records})
    print(f"vectorised {len(entries)} compositions in {time.time() - t0:.0f} s "
          f"({fallbacks} from the rounded key, {skipped + empty} skipped as empty)")

    unknown = sorted({db for e in entries for db in e["records"]} - set(databases))
    if unknown:
        print(f"warning: no metadata for databases {unknown}; they count as "
              "neither experimental nor structural")

    out = args.out / INDEX_FILE
    write_index(out, entries, databases, source=f"convert-index {args.mtree_lookup.name}")
    print(f"wrote {out} ({out.stat().st_size / 1e6:.0f} MB) in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
