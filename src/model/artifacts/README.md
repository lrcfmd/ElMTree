# Search index (runtime assets)

These files are **not committed to git**. The production index contains data
that may not be shared. They are declared as `assets` in
[`schema/model.json`](../../../schema/model.json) and provided to the tool at
run time by the LMDS platform:

1. An admin uploads a file for each declared asset path under the version's
   **Manage → Versions → Assets** page (the upload UI only accepts paths that
   are declared in the schema).
2. When a run starts, the runner downloads each asset from S3 to its declared
   path under the container WORKDIR (`./artifacts/…`), before the handler runs.

The model reads them from `./artifacts/` (or the package's `artifacts/` dir for
a local source checkout; `MODEL_ARTIFACTS_DIR` overrides both). It fails fast if
the required index is missing.

| Asset path (declared in model.json) | Required | What it is |
|--------------------------------------|----------|------------|
| `artifacts/indexed_LC.pk` | yes | Pickled `ElMTree` index of the composition database (the original web app's `indexed_LC.pk`, uploaded unchanged) |
| `artifacts/mtree_lookup.pk` | no | Pickled `{pretty_formula: {database: [ids], ...}}` lookup (the original `mtree_lookup.pk`). Without it, results carry no database ids |

The index pickle references the classes `ElMTree.ElMTree.ElMTree`/`Entry` and
`ElMD.ElMD.ElMD`; this package keeps those import paths, so an index pickled by
the original ElMTree code loads unchanged (verified with an index built by the
original code). `ElMD` is pinned (see `pyproject.toml`).

**Your own data:** build both files from a CSV with a `composition` column (and
an optional `id` column):

```bash
uv run build-index data/sample_compositions.csv --database matbench_expt_gap
```

**Local development:** place the files in this directory
(`src/model/artifacts/`); running from source resolves them here by default.
