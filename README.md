# ElMTree

Nearest-neighbour search over the chemical compositions in many materials
databases, packaged for the LMDS (Local Model Deployment Service) platform.

Compositions are compared with the **Element Movers Distance**
([ElMD](https://github.com/lrcfmd/ElMD), "fast" metric on the modified
Pettifor scale). This is the search behind
[ElMTree](https://lmds.liverpool.ac.uk/ElMTree): the production index holds
1,367,525 distinct compositions from 58 databases (ICSD, CSD, COD, Pearson's,
MPDS, Materials Project, AFLOW, Alexandria, OMDB, JARVIS, the matminer/matbench
datasets and more), each linked to its records in every database it appears
in.

This is a conversion of the original library + Flask app to the LMDS
model-template contract (`schema/model.json` + a single `process` handler).
The `ElMTree` library (list-of-clusters metric index,
[Chávez & Navarro 2005](https://doi.org/10.1016/j.patrec.2004.11.014)) lives
in `src/ElMTree/` with its original import path. The tool itself searches a
compact index built from the same data (see [Search index](#search-index-runtime-asset)).

## Quick Start

```bash
# Production index from the original server's lookup tables (~1 min, 241 MB):
uv run convert-index /path/to/mtree_lookup.pk /path/to/db_lookup.pk

# Or a small demo index from the bundled matbench_expt_gap compositions:
uv run build-index data/sample_compositions.csv --database matbench_expt_gap --experimental

# Nearest neighbours of one composition:
scripts/run_instance.sh LiFePO4 5

# Nearest neighbours for every composition in a CSV (15 trial queries):
scripts/run_dataset.sh data/trial_queries.csv 3     # -> data/output/neighbours.csv
```

Both commands write `src/model/artifacts/elmtree_index.db`. The demo index
covers the 4,601 distinct compositions in the public `matbench_expt_gap`
dataset, and its results list `matbench_expt_gap: <row id>` as their only
database record.

## Modes

| Mode | Response | Inputs | Output |
|------|----------|--------|--------|
| `instance` | immediate | `composition`, `k` (default 100) + filters | summary values (`query`, `normalised_query`, `matches`, `nearest_match`, `nearest_distance`), `elmtree_results.html`, and `neighbours.csv` (`rank`, `formula`, `normalised_formula`, `distance`, `databases`), closest first |
| `dataset` | deferred | `file` (CSV with a `composition` column), `k` (default 10) + filters | `elmtree_results.html` (one section per query, first 500) and `neighbours.csv`: one row per query × neighbour (`query_row`, `query`, `rank`, `formula`, `normalised_formula`, `distance`, `databases`, `error`) |

The defaults of k=100 (single search) and k=10 (batch) match the original web
search and API.

How results are presented:

- `elmtree_results.html` is the readable view. Its header gives the query, how
  it was compared (ElMD's normalised formula) and what was searched. Each
  neighbour is a row with:
  - its formula (`Li1Fe1P1O4` shown as LiFePO₄) and normalised form;
  - its distance, with a bar for scanning and "identical" at 0;
  - one chip per database it is reported in, with its record count.

  Chips are coloured by kind of database: orange for experimental structures,
  blue for experimental properties, green for computational. Hovering a chip
  describes the database, and clicking it opens the database's website in a
  new tab. Record IDs are folded behind "Show N record IDs". Dataset mode adds
  an overview of every query's nearest match, linking to each query's
  results. The page is light-themed and self-contained.
- `neighbours.csv` has the same list for download.
- `formula` is the formula as it was indexed (for the production index, as
  the source database gave it, e.g. `Li1Fe1P1O4`). `normalised_formula` is
  ElMD's fractional form (`Li0.143 Fe0.143 P0.143 O0.571`), which is what the
  distance compares.
- Distances are rounded to 3 dp, as on the original web search.
- `databases` lists every record of that composition, as
  `icsd: 260567 260568; mp: mp-19017; pearsons: 1607459`. CSV cells never
  contain commas, because the frontend's table preview splits on every comma.

`data/trial_queries.csv` is a ready-made dataset-mode input: 15 well-known
materials, with a `name` column that is ignored.

Optional filters, the "advanced search" in both modes:

- `must_contain` / `must_exclude`: lists of element symbols.
- `experimental_only` / `structures_only`: only compositions found in at least
  one experimental database / one database with crystal structures. Which
  databases count is the original server's `db_lookup.pk` metadata (for
  `build-index`, the `--experimental` / `--structures` flags).

Formulas ElMD cannot fully parse are rejected rather than searched. Examples
are unknown symbols (`Xx2O` would otherwise be silently searched as pure `O`)
and lower-case input. In `dataset` mode such rows get an `error` and the rest
of the batch still runs.

## Search index (runtime asset)

The index is **not committed to git** and **not baked into the image**. It is
declared as an asset in `schema/model.json` and provided by the LMDS platform
at run time:

| Asset path | Required | Upload |
|------------|----------|--------|
| `artifacts/elmtree_index.db` | yes | the output of `convert-index` (production) or `build-index` (your own data) |

1. An admin uploads the file for its declared asset path under the version's
   **Manage → Versions → Assets** page (only schema-declared paths are accepted).
2. On each run, the runner downloads it from S3 to `./artifacts/` under the
   container WORKDIR before the handler executes.

**Making the production index.** The original server
(`/home/mis/MIS/app/routes.py`) loaded `indexed_LC.pk` (5.4 GB) and
`mtree_lookup.pk` (198 MB). That pickled index cannot be used as it is: its
1.37 M ElMD objects each carry their own periodic-table dictionaries, so
unpickling it takes more than 10 GB of RAM (it was killed at a 10 GB cap), and
the runner would download all 5.4 GB on every run. It was built from each
lookup entry's first compound formula, so `convert-index` rebuilds the same
index from the two lookup tables alone:

- `mtree_lookup.pk` (`/home/mis/MIS/mtree_lookup.pk` on the old server):
  normalised formula → `{database: [ids], "compound_formula": [...]}`.
- `db_lookup.pk` (`/home/mis/db_lookup.pk`): each database's description and
  whether it is experimental / has structures.

The vectors come from the exact compound formulas (the lookup keys are
rounded to 3 dp, which would shift distances by up to ~0.01). For 449 of the
1,367,525 entries none of the formulas parse to the key, so the key is used.
One empty key (`""`, 76k CSD ids with no formula) is skipped.

**What the file holds.** The file is SQLite and stores each composition's
sparse ElMD ratio vector, its database records, the experimental/structures
flags and the database metadata (see `src/model/index.py`). At load (0.8 s,
about 0.9 GB of RAM for the production index) the tool builds float32
cumulative vectors. The "fast" ElMD is the L1 distance between those, and
entries are stored sorted by a 1-D lower bound on that distance. A search only
scans the window around the query that the bound can't rule out, so a typical
query takes about 10 ms. Results are exact: they match a brute-force scan of
all 1.37 M compositions, and the tests check this.

The index records the ElMD version and vector length it was built with, and
the tool refuses an index whose vectors don't match the installed ElMD.

## Changes from the original code

1. **The search runs on a compact index instead of the pickled `ElMTree`**
   (see [Search index](#search-index-runtime-asset)): one 241 MB file in
   place of 5.4 GB + 198 MB, loaded in under a second instead of needing more
   than 10 GB of RAM. Results are exact and carry all of the original
   lookup's 3,966,299 database records.
2. **Exact k-nearest-neighbour search in the `ElMTree` library.** `ElMTree.knn` skipped a cluster
   unless `abs(d(q, centre) − radius) <= bound`. The correct test is that every
   point in a cluster is at least `d(q, centre) − radius` from the query, so
   the condition should be `d(q, centre) − radius <= bound`. The `abs()` also
   skipped clusters the query sits deep inside. On the bundled data the
   original returned wrong neighbours for about 21% of queries. With the fix,
   results match brute force exactly (the fix is in the query code, not the
   index). The tests check that the library and the tool's index agree.
3. **Advanced-search filters work.** In the original web app they never ran:
   the route used undefined `split()` and dictionary keys the tree didn't read
   (`structural`/`must_include`/`must_exlude` vs
   `structures`/`must_contain`/`must_exclude`). In the library, the element
   filters also work for indexes built without lookup tables, and `knn`
   returns `[]` instead of raising `IndexError` when filters exclude every
   entry.
4. Batch queries run one at a time in-process. The original used `process_map`,
   which copies the whole index into every worker process.
5. The Flask/WTForms/REST layer (`ElMTree/app`, `wsgi.py`) and `setup.py` are
   dropped. The web app's ElM2D (UMAP) and ElMD (single-distance) pages are
   separate LMDS tools and aren't part of this one.

## Project Structure

```
├── schema/
│   ├── model.json              # Contract (instance + dataset modes, index asset)
│   └── model.schema.json       # Shared meta-schema (same as the other converted tools)
├── src/
│   ├── ElMTree/                # The ElMTree library (original import path, kNN fix)
│   ├── cmds/                   # Entrypoints — copied from template, unmodified
│   ├── infrastructure/         # Runner/handler/S3/logging — copied, unmodified
│   └── model/
│       ├── __init__.py         # handler = ElMTreeSearchModel()
│       ├── model.py            # Input validation, search, result formatting
│       ├── index.py            # The index file format and exact kNN search
│       ├── convert_index.py    # `convert-index`: production index from the server's lookups
│       ├── build_index.py      # `build-index`: index from a CSV
│       ├── types.py            # Auto-generated from model.json
│       └── artifacts/          # elmtree_index.db — provided at runtime, NOT in git
├── data/                       # sample_compositions.csv (matbench_expt_gap), sample inputs
├── scripts/                    # run_instance.sh, run_dataset.sh
├── tests/test_model.py
├── Dockerfile
└── pyproject.toml
```

## Local development

```bash
uv sync --extra dev
uv run validate-model schema/model.json
uv run model-types schema/model.json -o src/model/types.py
uv run pytest tests/ -v          # builds its own small index; no assets needed
```

## Dependencies

Runs on the template default **Python 3.12**. `ElMD` is pinned to **0.5.15**
because index and query vectors must come from the same ElMD. ElMD imports
`pkg_resources`, hence `setuptools<81`. `numba` speeds up ElMD. `pymatgen` and
`ase` are needed by the library, which also accepts
`Structure`/`Composition`/`Atoms` objects. The index is SQLite (Python's
standard library).
Docker installs the exact `uv.lock` set (`uv sync --frozen`).

## Docker

```bash
docker build -t elmtree:1.1.0 .

# On the platform the runner downloads the declared asset to ./artifacts before
# running. To run the image by hand, mount the folder holding elmtree_index.db
# where the runner would put it (/app/artifacts):
docker run --rm \
  -v /path/to/artifacts:/app/artifacts:ro \
  elmtree:1.1.0 --local \
  --parameters '{"mode":"instance","inputs":{"composition":{"value":"LiFePO4"},"k":{"value":5}},"parameters":{}}'
```

## Environment Variables

Copy `.env.example` to `.env` and adjust for remote (S3) runs:

| Variable | Default | Description |
|----------|---------|-------------|
| `LOG_LEVEL` | `DEBUG` | Logging verbosity |
| `S3_ENDPOINT` | `http://minio:9000` | S3-compatible storage endpoint |
| `S3_ACCESS_KEY` | `minioadmin` | S3 access key |
| `S3_SECRET_KEY` | `minioadmin` | S3 secret key |
| `S3_BUCKET` | `lmds` | S3 bucket name |

## Licence

GPL-3.0 (see `LICENSE`), as the original ElMTree.
