import csv
import os
import re
import time
from pathlib import Path

import pandas as pd
from ElMD import ElMD

from infrastructure.logging import newLogger
from model import report
from model.index import INDEX_FILE, CompositionIndex

DEFAULT_K = {"instance": 100, "dataset": 10}  # as the original web search / API
MAX_K = 1000
DISTANCE_DP = 3  # as the original web search displayed distances
RESULTS_HTML = "elmtree_results.html"
MAX_HTML_QUERIES = 500  # dataset-mode queries rendered in the HTML page

# The LMDS frontend previews CSV outputs as a table but splits rows on every
# comma (no quoting), so CSV cells must not contain commas.
CSV_FIELDS = {
    "instance": ["rank", "formula", "normalised_formula", "distance", "databases"],
    "dataset": ["query_row", "query", "rank", "formula", "normalised_formula",
                "distance", "databases", "error"],
}


def _csv_safe(value) -> str:
    return str(value).replace(",", ";") if isinstance(value, str) else value


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows({k: _csv_safe(v) for k, v in row.items()} for row in rows)


# The index file is NOT committed to git; it is declared as an `asset` in
# schema/model.json and provided at runtime. On the platform the runner
# downloads them to their declared paths under the container WORKDIR
# (./artifacts/…), so that is the primary location. For a local source checkout
# the files live in the package's own artifacts/ dir. MODEL_ARTIFACTS_DIR
# overrides both. See artifacts/README.md.
def _artifacts_dir() -> Path:
    env = os.environ.get("MODEL_ARTIFACTS_DIR")
    if env:
        return Path(env)
    workdir_assets = Path("artifacts")  # where the runner downloads assets
    if workdir_assets.is_dir():
        return workdir_assets
    return Path(__file__).resolve().parent / "artifacts"


class ElMTreeSearchModel:
    """Nearest-neighbour search over an indexed database of compositions.

    Ported from the LMDS ElMTree web app. Compositions are compared with the
    Element Movers Distance (ElMD, "fast" 1-D metric on the modified
    Pettifor scale) and searched exactly over a pre-built index of the
    compositions in many databases (see model/index.py). Two Modes share one
    search path:
      ``instance`` — one composition, its k nearest neighbours inline.
      ``dataset``  — a CSV of compositions, a neighbours CSV out.
    """

    def __init__(self):
        self.logger = newLogger("model")
        self._loaded = False

    def _load(self) -> None:
        if self._loaded:
            return
        artifacts = _artifacts_dir()
        index_path = artifacts / INDEX_FILE
        self.logger.info("loading artifacts", dir=str(artifacts))

        if not index_path.is_file():
            raise FileNotFoundError(
                f"Search index missing from {artifacts}: {INDEX_FILE}. "
                "This file is provided to the application at runtime and is "
                "not committed to git — see src/model/artifacts/README.md. Set "
                "MODEL_ARTIFACTS_DIR to the directory that holds it, or mount "
                "it into the default location."
            )
        self.index = CompositionIndex(index_path)
        self.logger.info(
            "loaded index",
            compositions=len(self.index),
            databases=len(self.index.databases),
            source=self.index.meta.get("source"),
        )
        self._loaded = True

    # -- input handling ---------------------------------------------------

    @staticmethod
    def parse(composition) -> ElMD:
        """Parse a formula, rejecting what ElMD would silently drop."""
        if not isinstance(composition, str) or not composition.strip():
            raise ValueError("Composition must be a non-empty string")
        elmd = ElMD(composition, metric="fast")
        symbols = {m[0] for m in re.findall(ElMD.ATOM_REGEX, composition)}
        unknown = sorted(s for s in symbols if s not in elmd.petti_lookup)
        if unknown:
            raise ValueError(
                f"Unknown element(s) {', '.join(unknown)} in {composition!r}"
            )
        if not elmd.normed_composition:
            raise ValueError(f"Could not parse composition {composition!r}")
        return elmd

    @staticmethod
    def _k(values: dict, mode: str) -> int:
        k = values.get("k", DEFAULT_K[mode])
        if isinstance(k, bool) or not isinstance(k, int) or not 1 <= k <= MAX_K:
            raise ValueError(f"k must be an integer between 1 and {MAX_K}")
        return k

    def _filters(self, values: dict) -> dict | None:
        def elements(name):
            raw = values.get(name) or []
            if isinstance(raw, str):  # tolerate "Fe, O"
                raw = raw.split(",")
            out = {str(e).strip() for e in raw if str(e).strip()}
            petti = ElMD("", metric="fast").petti_lookup
            unknown = sorted(e for e in out if e not in petti)
            if unknown:
                raise ValueError(f"{name}: unknown element(s) {', '.join(unknown)}")
            return out

        search = {
            "must_contain": elements("must_contain"),
            "must_exclude": elements("must_exclude"),
            "experimental": bool(values.get("experimental_only", False)),
            "structures": bool(values.get("structures_only", False)),
        }
        if not any(search.values()):
            return None
        if search["must_contain"] & search["must_exclude"]:
            raise ValueError("An element cannot be in both must_contain and must_exclude")
        return search

    # -- search -----------------------------------------------------------

    @staticmethod
    def _ids(db: str, ids) -> list:
        ids = list(ids) if isinstance(ids, (list, tuple, set)) else [ids]
        if db == "mpds":  # one record per MPDS entry, not per phase
            ids = sorted({str(x).split("-")[0] for x in ids})
        return ids

    def _ordered_records(self, records: dict) -> dict:
        """A neighbour's records, ordered as the original site listed them:
        experimental structures, then experimental properties, then
        computational databases."""
        dbs = [db for db in records if db != "compound_formula"]
        dbs.sort(key=lambda db: report.KIND_ORDER[report.kind(self.index.databases.get(db))])
        return {db: self._ids(db, records[db]) for db in dbs}

    @classmethod
    def _databases(cls, records: dict) -> list[str]:
        """Database records for a neighbour, as 'database: id id' strings."""
        return [f"{db}: {' '.join(map(str, cls._ids(db, ids)))}"
                for db, ids in records.items() if db != "compound_formula"]

    def search(self, composition: str, k: int, advanced_search: dict | None = None) -> dict:
        self._load()
        started = time.perf_counter()
        query = self.parse(composition)
        mask = self.index.mask(advanced_search, query.petti_lookup)
        results = self.index.knn(query.ratio_vector, k, mask)
        entries = self.index.entries([row for row, _ in results])
        neighbours = []
        for row, distance in results:
            entry = entries[row]
            records = self._ordered_records(entry["records"])
            # Show the formula as it was indexed (e.g. "LiFePO4"); the
            # normalised ElMD form ("Li0.143 Fe0.143 ...") is what is compared.
            neighbours.append({
                "formula": entry["formula"] or entry["normalised_formula"],
                "normalised_formula": entry["normalised_formula"],
                "distance": distance,
                "databases": self._databases(records),
                "records": list(records.items()),
            })
        return {"input": composition.strip(), "query": query.pretty_formula,
                "neighbours": neighbours, "seconds": time.perf_counter() - started}

    def _stats(self) -> dict:
        return {"compositions": len(self.index), "records": self.index.n_records,
                "databases": len(self.index.databases)}

    @staticmethod
    def _presented(n: dict) -> dict:
        """A neighbour as a CSV row: rounded distance, databases on one line."""
        return {**n, "distance": round(n["distance"], DISTANCE_DP),
                "databases": "; ".join(n["databases"])}

    # -- modes ------------------------------------------------------------

    def process(
        self,
        mode: str,
        values: dict,
        files: dict[str, tuple[Path, str]],
        output_dir: Path,
        parameters: dict,
        logger,
    ) -> dict:
        if mode == "instance":
            return self._process_instance(values, output_dir, logger)
        elif mode == "dataset":
            return self._process_dataset(values, files, output_dir, logger)
        else:
            raise ValueError(f"Unknown mode: {mode!r}")

    def _process_instance(self, values: dict, output_dir: Path, logger) -> dict:
        composition = values.get("composition")
        if composition is None:
            raise ValueError("No composition provided")
        result = self.search(composition, self._k(values, "instance"), self._filters(values))
        neighbours = result["neighbours"]
        logger.info("searched", query=result["query"], neighbours=len(neighbours))

        # The results page shows the neighbours as the original web search did;
        # the CSV has the same list for download, and the inline values are a
        # short summary.
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / RESULTS_HTML).write_text(report.instance_page(
            result, self._filters(values), self._stats(), self.index.databases))
        _write_csv(
            output_dir / "neighbours.csv", CSV_FIELDS["instance"],
            [{"rank": i, **self._presented(n)} for i, n in enumerate(neighbours, start=1)],
        )
        nearest = neighbours[0] if neighbours else None
        return {
            "query": composition.strip(),
            "normalised_query": result["query"],
            "matches": len(neighbours),
            "nearest_match": nearest["formula"] if nearest else "",
            "nearest_distance": round(nearest["distance"], DISTANCE_DP) if nearest else None,
        }

    def _process_dataset(
        self, values: dict, files: dict[str, tuple[Path, str]], output_dir: Path, logger
    ) -> dict:
        if "file" not in files:
            raise ValueError("No input file provided")
        csv_path, _ = files["file"]
        k = self._k(values, "dataset")
        advanced_search = self._filters(values)

        try:
            data = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
        except Exception as e:
            raise ValueError(f"Could not read input CSV: {e}") from e
        if "composition" not in data.columns:
            raise ValueError(
                "Input CSV must have a 'composition' column "
                f"(found: {', '.join(map(str, data.columns))})"
            )
        if data.empty:
            raise ValueError("Input CSV has no rows")

        self._load()
        rows = []
        pages = []
        failed = 0
        # Queries run one at a time in this process (a few ms each); the
        # original list path used process_map, which copied the whole index
        # into every worker.
        for i, composition in enumerate(data["composition"]):
            try:
                result = self.search(composition, k, advanced_search)
            except Exception as e:  # a bad composition shouldn't abort the batch
                logger.warning("search failed", row=i, composition=composition, error=str(e))
                failed += 1
                rows.append({"query_row": i, "query": composition, "error": str(e)})
                pages.append({"row": i, "query": composition, "error": str(e)})
                continue
            pages.append({"row": i, "query": composition, "result": result})
            if not result["neighbours"]:
                rows.append({"query_row": i, "query": composition, "error": "no matches"})
            for rank, n in enumerate(result["neighbours"], start=1):
                rows.append({
                    "query_row": i,
                    "query": composition,
                    "rank": rank,
                    **self._presented(n),
                    "error": "",
                })
        logger.info("searched dataset", queries=len(data), failed=failed, rows=len(rows))

        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "neighbours.csv"
        _write_csv(output_path, CSV_FIELDS["dataset"], rows)
        (output_dir / RESULTS_HTML).write_text(report.dataset_page(
            pages, advanced_search, self._stats(), self.index.databases, MAX_HTML_QUERIES))
        logger.info("saved assets", csv=str(output_path), html=RESULTS_HTML)
        return {}
