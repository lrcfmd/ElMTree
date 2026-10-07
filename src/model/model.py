import csv
import os
import pickle
import re
from pathlib import Path

import pandas as pd
from ElMD import ElMD

from infrastructure.logging import newLogger

INDEX_FILE = "indexed_LC.pk"
LOOKUP_FILE = "mtree_lookup.pk"

DEFAULT_K = {"instance": 100, "dataset": 10}  # as the original web search / API
MAX_K = 1000


# The index files are NOT committed to git; they are declared as `assets` in
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
    Pettifor scale) and searched with a pre-built ElMTree index. Two Modes
    share one search path:
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
        lookup_path = artifacts / LOOKUP_FILE
        self.logger.info("loading artifacts", dir=str(artifacts))

        if not index_path.is_file():
            raise FileNotFoundError(
                f"Search index missing from {artifacts}: {INDEX_FILE}. "
                "This file is provided to the application at runtime and is "
                "not committed to git — see src/model/artifacts/README.md. Set "
                "MODEL_ARTIFACTS_DIR to the directory that holds it, or mount "
                "it into the default location."
            )
        # The index files are trusted platform assets (uploaded by an admin),
        # loaded exactly as the original app did.
        with index_path.open("rb") as f:
            self.tree = pickle.load(f)
        self.lookup = None
        if lookup_path.is_file():
            with lookup_path.open("rb") as f:
                self.lookup = pickle.load(f)
        self.logger.info(
            "loaded index",
            centres=len(self.tree.centres),
            lookup=self.lookup is not None,
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

    def _databases(self, formula: str, tree_value) -> list[str]:
        """Database records for a neighbour, as 'database: id, id' strings."""
        entry = self.lookup.get(formula) if self.lookup is not None else tree_value
        if not isinstance(entry, dict):
            return []
        out = []
        for db, ids in entry.items():
            if db == "compound_formula":
                continue
            if isinstance(ids, (list, tuple, set)):
                ids = list(ids)
                if db == "mpds":  # one record per MPDS entry, not per phase
                    ids = sorted({str(x).split("-")[0] for x in ids})
                out.append(f"{db}: {', '.join(map(str, ids))}")
            else:
                out.append(f"{db}: {ids}")
        return out

    def search(self, composition: str, k: int, advanced_search: dict | None = None) -> dict:
        self._load()
        query = self.parse(composition)
        results = self.tree.knn(query, k, advanced_search=advanced_search)
        neighbours = []
        for entry, distance, tree_value in results:
            formula = ElMD(entry).pretty_formula
            neighbours.append({
                "formula": formula,
                "distance": float(distance),
                "databases": self._databases(formula, tree_value),
            })
        return {"query": query.pretty_formula, "neighbours": neighbours}

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
            return self._process_instance(values, logger)
        elif mode == "dataset":
            return self._process_dataset(values, files, output_dir, logger)
        else:
            raise ValueError(f"Unknown mode: {mode!r}")

    def _process_instance(self, values: dict, logger) -> dict:
        composition = values.get("composition")
        if composition is None:
            raise ValueError("No composition provided")
        result = self.search(composition, self._k(values, "instance"), self._filters(values))
        logger.info("searched", query=result["query"], neighbours=len(result["neighbours"]))
        return result

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
        failed = 0
        # Queries run one at a time in this process: the original list path
        # used process_map, which copies the whole index into every worker.
        for i, composition in enumerate(data["composition"]):
            try:
                result = self.search(composition, k, advanced_search)
            except Exception as e:  # a bad composition shouldn't abort the batch
                logger.warning("search failed", row=i, composition=composition, error=str(e))
                failed += 1
                rows.append({"query_row": i, "query_composition": composition, "error": str(e)})
                continue
            if not result["neighbours"]:
                rows.append({
                    "query_row": i, "query_composition": composition,
                    "query_formula": result["query"], "error": "no matches",
                })
            for rank, n in enumerate(result["neighbours"], start=1):
                rows.append({
                    "query_row": i,
                    "query_composition": composition,
                    "query_formula": result["query"],
                    "rank": rank,
                    "formula": n["formula"],
                    "distance": n["distance"],
                    "databases": "; ".join(n["databases"]),
                    "error": "",
                })
        logger.info("searched dataset", queries=len(data), failed=failed, rows=len(rows))

        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "neighbours.csv"
        fieldnames = [
            "query_row", "query_composition", "query_formula", "rank",
            "formula", "distance", "databases", "error",
        ]
        with output_path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        logger.info("saved asset", path=str(output_path))
        return {}
