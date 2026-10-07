# ElMTree

Nearest-neighbour search over a database of chemical compositions, packaged
for the LMDS (Local Model Deployment Service) platform.

Compositions are compared with the **Element Movers Distance**
([ElMD](https://github.com/lrcfmd/ElMD), "fast" metric on the modified
Pettifor scale). The search uses an **ElMTree**, a list-of-clusters metric
index ([Chávez & Navarro 2005](https://doi.org/10.1016/j.patrec.2004.11.014)).
This is the indexing class behind the
[ElMTree](https://lmds.liverpool.ac.uk/ElMTree) search.

This is a conversion of the original library + Flask app to the LMDS
model-template contract (`schema/model.json` + a single `process` handler).
The `ElMTree` library itself lives in `src/ElMTree/` with its original import
path.

## Quick Start

```bash
# Build a demo index from the bundled matbench_expt_gap compositions (~5 s):
uv run build-index data/sample_compositions.csv --database matbench_expt_gap

# Nearest neighbours of one composition:
scripts/run_instance.sh LiFePO4 5

# Nearest neighbours for every composition in a CSV:
scripts/run_dataset.sh data/input.csv 3     # -> data/output/neighbours.csv
```

## Modes

| Mode | Response | Inputs | Output |
|------|----------|--------|--------|
| `instance` | immediate | `composition`, `k` (default 100) + filters | `query` (normalised formula), `neighbours` (`formula`, `distance`, `databases`), closest first |
| `dataset` | deferred | `file` (CSV with a `composition` column), `k` (default 10) + filters | `neighbours.csv`: one row per query × neighbour (`query_row`, `query_composition`, `query_formula`, `rank`, `formula`, `distance`, `databases`, `error`) |

The defaults of k=100 (single search) and k=10 (batch) match the original web
search and API.

Optional filters, the "advanced search" in both modes:

- `must_contain` / `must_exclude`: lists of element symbols. These work with
  any index.
- `experimental_only` / `structures_only`: these only take effect for an index
  built with database lookup tables (`lookup_tables=True`), which record
  per-entry experimental/structure flags. For other indexes they are ignored.

Formulas ElMD cannot fully parse are rejected rather than searched. Examples
are unknown symbols (`Xx2O` would otherwise be silently searched as pure `O`)
and lower-case input. In `dataset` mode such rows get an `error` and the rest
of the batch still runs.

## Search index (runtime assets)

The index files are **not committed to git** and **not baked into the image**.
They are declared as `assets` in `schema/model.json` and provided by the LMDS
platform at run time:

| Asset path | Required | Upload |
|------------|----------|--------|
| `artifacts/indexed_LC.pk` | yes | the original web app's `indexed_LC.pk` (unchanged), or one made by `build-index` |
| `artifacts/mtree_lookup.pk` | no | the original `mtree_lookup.pk` (composition → database ids), or one made by `build-index` |

1. An admin uploads each file for its declared asset path under the version's
   **Manage → Versions → Assets** page (only schema-declared paths are accepted).
2. On each run, the runner downloads them from S3 to `./artifacts/` under the
   container WORKDIR before the handler executes.

See [`src/model/artifacts/README.md`](src/model/artifacts/README.md). The
whole index is unpickled into memory on every run. The demo index is 17 MB for
4.6k compositions, so a large production index takes correspondingly longer to
load and needs matching container memory.

## Changes from the original code

1. **Exact k-nearest-neighbour search.** `ElMTree.knn` skipped a cluster
   unless `abs(d(q, centre) − radius) <= bound`. The correct test is that every
   point in a cluster is at least `d(q, centre) − radius` from the query, so
   the condition should be `d(q, centre) − radius <= bound`. The `abs()` also
   skipped clusters the query sits deep inside. On the bundled data the
   original returned wrong neighbours for about 21% of queries. With the fix,
   results match brute force exactly, including for an index pickled by the
   original code (the fix is in the query code, not the index). The tests check
   this against brute force on every run.
2. **Advanced-search filters work.** In the original web app they never ran:
   the route used undefined `split()` and dictionary keys the tree didn't read
   (`structural`/`must_include`/`must_exlude` vs
   `structures`/`must_contain`/`must_exclude`). The element filters now also
   work for indexes built without lookup tables.
3. `knn` returns `[]` instead of raising `IndexError` when filters exclude
   every entry.
4. Batch queries run one at a time in-process. The original used `process_map`,
   which copies the whole index into every worker process.
5. The Flask/WTForms/REST layer (`ElMTree/app`, `wsgi.py`) and `setup.py` are
   dropped. The web app's ElM2D (UMAP) and ElMD (single-distance) pages are
   separate LMDS tools and aren't part of this one.

## Project Structure

```
├── schema/
│   ├── model.json              # Contract (instance + dataset modes, index assets)
│   └── model.schema.json       # Shared meta-schema (same as the other converted tools)
├── src/
│   ├── ElMTree/                # The ElMTree library (original import path, kNN fix)
│   ├── cmds/                   # Entrypoints — copied from template, unmodified
│   ├── infrastructure/         # Runner/handler/S3/logging — copied, unmodified
│   └── model/
│       ├── __init__.py         # handler = ElMTreeSearchModel()
│       ├── model.py            # Input validation, search, result formatting
│       ├── build_index.py      # `build-index`: make index files from a CSV
│       ├── types.py            # Auto-generated from model.json
│       └── artifacts/          # Index files — provided at runtime, NOT in git
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
because the index pickles embed ElMD objects. ElMD imports `pkg_resources`,
hence `setuptools<81`. `numba` speeds up ElMD. `pymatgen` and `ase` are needed
by the library, which also accepts `Structure`/`Composition`/`Atoms` objects.
Docker installs the exact `uv.lock` set (`uv sync --frozen`).

## Docker

```bash
docker build -t elmtree:1.0.0 .

# On the platform the runner downloads the declared assets to ./artifacts before
# running. To run the image by hand, mount the index files where the runner
# would put them (/app/artifacts):
docker run --rm \
  -v /path/to/artifacts:/app/artifacts:ro \
  elmtree:1.0.0 --local \
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
