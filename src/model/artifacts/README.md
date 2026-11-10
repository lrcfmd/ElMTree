# Search index (runtime asset)

This file is **not committed to git**: the production index contains
identifiers from licensed databases (ICSD, CSD, Pearson's, MPDS). It is
declared as an `asset` in [`schema/model.json`](../../../schema/model.json)
and provided to the tool at run time by the LMDS platform:

1. An admin uploads the file for its declared asset path under the version's
   **Manage → Versions → Assets** page (the upload UI only accepts paths that
   are declared in the schema).
2. When a run starts, the runner downloads it from S3 to its declared path
   under the container WORKDIR (`./artifacts/…`), before the handler runs.

The model reads it from `./artifacts/` (or the package's `artifacts/` dir for
a local source checkout; `MODEL_ARTIFACTS_DIR` overrides both). It fails fast
if the index is missing or was built for a different ElMD vector length.

| Asset path (declared in model.json) | Required | What it is |
|--------------------------------------|----------|------------|
| `artifacts/elmtree_index.db` | yes | SQLite index: every composition's ElMD vector, its records in each database, and the database metadata (see `src/model/index.py`) |

**Production index** (1,367,525 compositions, 58 databases, 241 MB), from the
original web search's lookup tables:

```bash
uv run convert-index mtree_lookup.pk db_lookup.pk
```

**Your own data:** from a CSV with a `composition` column (and an optional
`id` column):

```bash
uv run build-index data/sample_compositions.csv --database matbench_expt_gap --experimental
```

The original server's `indexed_LC.pk` is not used. It needs more than 10 GB
of RAM to unpickle, and `convert-index` rebuilds the same index from the two
lookup tables (see the main README).

**Local development:** both commands write to this directory
(`src/model/artifacts/`) by default, and running from source reads it here.
