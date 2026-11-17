import csv
import os
import re
import time
from pathlib import Path

from ElMD import ElMD

from infrastructure.handler import InputError
from infrastructure.logging import newLogger
from model import report
from model.index import INDEX_FILE, CompositionIndex

DEFAULT_K = {"instance": 100, "dataset": 10}  # as the original web search / API
MAX_K = 1000
MAX_LISTED_ERRORS = 5  # bad CSV rows named in an error message
FORMULA_HINT = ("write each element symbol with a capital first letter, followed by "
                "its amount, e.g. NaCl, Fe2O3 or LiFePO4")


def formula_error(reason: str, hint: str) -> InputError:
    """An InputError for one bad formula; ``reason`` alone goes in a row's error cell."""
    error = InputError(f"{reason}. {hint}.")
    error.reason = reason
    return error
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


def read_compositions_csv(path: Path):
    """The input CSV as strings, with a 'composition' column and at least one row."""
    import pandas as pd  # only needed here; importing it costs ~1 s per run

    try:
        # index_col=False: a row with extra cells must not shift its values
        # into a different column.
        data = pd.read_csv(path, dtype=str, keep_default_na=False, index_col=False)
    except UnicodeDecodeError:
        raise InputError("The input file is not a text CSV file: save it as CSV (UTF-8) and upload it again.")
    except pd.errors.EmptyDataError:
        raise InputError("The input CSV is empty: it needs a header row with a 'composition' "
                         "column, then one chemical formula per row.")
    except (pd.errors.ParserError, ValueError) as e:
        raise InputError(f"Could not read the input CSV ({str(e).strip()}). Check that it is a "
                         "comma-separated file with a header row.")
    if "composition" not in data.columns:
        raise InputError(
            "Input CSV must have a 'composition' column: the first row must name the columns, "
            "with the chemical formulas under 'composition' "
            f"(found: {', '.join(repr(str(c)) for c in data.columns)})."
        )
    if data.empty:
        raise InputError("Input CSV has no rows: add one chemical formula per row under the "
                         "'composition' header.")
    return data


def describe_bad_rows(errors: dict[int, str]) -> str:
    """The first few bad rows as 'row N: reason', rows numbered as in a
    spreadsheet (row 1 is the header)."""
    listed = [f"row {i + 2}: {errors[i]}"
              for i in sorted(errors)[:MAX_LISTED_ERRORS]]
    more = len(errors) - len(listed)
    return "; ".join(listed) + (f"; and {more} more" if more > 0 else "")


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
        """Parse a formula, rejecting what ElMD would silently drop.

        Raises InputError with a ``reason`` attribute (the short reason, for a
        dataset row's error cell) and a message that also says how to fix it."""
        if not isinstance(composition, str) or not composition.strip():
            raise formula_error("The composition is empty", f"Enter a chemical formula: {FORMULA_HINT}")
        try:
            elmd = ElMD(composition, metric="fast")
        except Exception as e:
            problem = ("its brackets do not match" if "bracket" in str(e).lower()
                       else "an amount is not a number")
            raise formula_error(f"Could not parse {composition!r} as a chemical formula: {problem}",
                                FORMULA_HINT[0].upper() + FORMULA_HINT[1:]) from e
        symbols = {m[0] for m in re.findall(ElMD.ATOM_REGEX, composition)}
        unknown = sorted(s for s in symbols if s not in elmd.petti_lookup)
        if unknown:
            raise formula_error(
                f"Unknown element(s) {', '.join(unknown)} in {composition!r}",
                "Check the spelling and capitalisation of each element symbol (e.g. Co, not CO or co)",
            )
        if not elmd.normed_composition:
            raise formula_error(f"Could not parse {composition!r} as a chemical formula",
                                FORMULA_HINT[0].upper() + FORMULA_HINT[1:])
        return elmd

    @staticmethod
    def _k(values: dict, mode: str) -> int:
        k = values.get("k", DEFAULT_K[mode])
        if isinstance(k, bool) or not isinstance(k, int) or not 1 <= k <= MAX_K:
            raise InputError(f"k must be a whole number between 1 and {MAX_K} (got {k!r}).")
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
                raise InputError(
                    f"{name}: unknown element(s) {', '.join(unknown)}. Enter element symbols "
                    "with a capital first letter, e.g. Li, O, Fe."
                )
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
            both = ", ".join(sorted(search["must_contain"] & search["must_exclude"]))
            raise InputError(
                f"{both} cannot be in both must_contain and must_exclude; remove it from one of them."
            )
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
        query = self.parse(composition)  # before loading the index: bad input fails fast
        self._load()
        started = time.perf_counter()
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
        return {"query": query.pretty_formula, "neighbours": neighbours,
                "seconds": time.perf_counter() - started}

    def _stats(self) -> dict:
        return {"compositions": len(self.index), "records": self.index.n_records}

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
            raise InputError("No composition provided: enter a chemical formula, e.g. NaCl.")
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
            raise InputError("No input file provided: upload a CSV with a 'composition' column.")
        csv_path, _ = files["file"]
        k = self._k(values, "dataset")
        advanced_search = self._filters(values)
        data = read_compositions_csv(csv_path)

        # Check every formula before loading the index, so a file with no
        # usable composition fails at once with the reasons.
        errors = {}
        for i, composition in enumerate(data["composition"]):
            try:
                self.parse(composition)
            except InputError as e:
                errors[i] = e.reason
        if len(errors) == len(data):
            raise InputError(
                "None of the compositions in the CSV could be read: "
                + describe_bad_rows(errors)
                + f". In the 'composition' column, {FORMULA_HINT}."
            )

        self._load()
        rows = []
        pages = []
        failed = 0
        # Queries run one at a time in this process (a few ms each); the
        # original list path used process_map, which copied the whole index
        # into every worker.
        for i, composition in enumerate(data["composition"]):
            try:
                if i in errors:  # a bad composition shouldn't abort the batch
                    raise InputError(errors[i])
                result = self.search(composition, k, advanced_search)
            except Exception as e:
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
