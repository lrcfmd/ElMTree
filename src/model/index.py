"""Compact composition index for exact nearest-neighbour search by ElMD.

The tool searches with ElMD's "fast" metric (modified Pettifor scale), which is
the 1-D Earth Mover's Distance between two ratio vectors:

    d(a, b) = sum_j |cumsum(a - b)_j| = || C(a) - C(b) ||_1

so every composition is stored once as its sparse ratio vector, and the search
is over the cumulative vectors C. That is far smaller and faster to load than
the original pickled ElMTree, whose 1.37 M ElMD objects each carried their own
copies of the periodic-table lookups (5.4 GB on disk, >10 GB to unpickle).

The search is exact. Entries are stored sorted by S = sum_j C_j, and
|S(a) - S(q)| <= d(a, q), so scanning outwards from the query's S can stop as
soon as the next unscanned entry's bound exceeds the current k-th distance.
Candidates are screened with float32 cumulative vectors, then re-ranked with
float64 distances computed exactly as ElMD computes them.

Everything lives in one SQLite file (the tool's single runtime asset):
  meta       key/value: format_version, elmd_version, vector_length, ...
  arrays     the sparse vectors (indptr, positions, ratios), the sort key
             (sums) and per-entry experimental/structures flags, as .npy blobs
  entries    id (= row in the arrays), formula, normalised_formula and the
             database records as JSON {database: [ids]}
  databases  name, description, experimental, structures, url
"""

import io
import json
import sqlite3
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import numpy as np
from ElMD import ElMD

FORMAT_VERSION = "1"
INDEX_FILE = "elmtree_index.db"

# The float32 screening distances are within ~1e-4 of the exact ones, so
# keeping every candidate within EPS of the screened k-th distance (and
# scanning until the bound clears it) guarantees the exact top k is kept.
EPS = 1e-3
CHUNK = 65536  # rows per vectorised block (~27 MB of float32 temporaries)
FIRST_BLOCK = 4096


def vector_length() -> int:
    return len(ElMD("", metric="fast").ratio_vector)


def sparse_vector(elmd: ElMD) -> tuple[np.ndarray, np.ndarray]:
    """(positions, ratios) of an ElMD composition's non-zero ratios."""
    positions = np.flatnonzero(elmd.ratio_vector)
    return positions, elmd.ratio_vector[positions]


def _sums(indptr, positions, ratios, width) -> np.ndarray:
    # sum_j cumsum(r)_j = sum_i r_i * (width - p_i)
    weighted = ratios * (width - positions.astype(np.float64))
    return np.add.reduceat(weighted, indptr[:-1]) if len(weighted) else np.zeros(len(indptr) - 1)


def _blob(array: np.ndarray) -> bytes:
    buf = io.BytesIO()
    np.save(buf, array, allow_pickle=False)
    return buf.getvalue()


def write_index(path: Path, entries: list[dict], databases: dict[str, dict], source: str) -> None:
    """Write an index file.

    entries: [{"formula", "normalised_formula", "positions", "ratios",
               "records": {database: [ids]}}]
    databases: {name: {"description", "experimental", "structures", "url"}}
    """
    if not entries:
        raise ValueError("No compositions to index")
    if any(len(e["positions"]) == 0 for e in entries):
        raise ValueError("Every indexed composition needs at least one element")
    width = vector_length()

    def csr(entries):
        lengths = np.array([len(e["positions"]) for e in entries], dtype=np.int64)
        indptr = np.concatenate([[0], np.cumsum(lengths)]).astype(np.int64)
        positions = np.concatenate([e["positions"] for e in entries]).astype(np.uint8)
        ratios = np.concatenate([e["ratios"] for e in entries]).astype(np.float64)
        return indptr, positions, ratios

    # Store the entries sorted by their search key
    order = np.argsort(_sums(*csr(entries), width), kind="stable")
    entries = [entries[i] for i in order]
    indptr, positions, ratios = csr(entries)

    def flag(e, key):
        return any(databases.get(db, {}).get(key, False) for db in e["records"])

    arrays = {
        "indptr": indptr,
        "positions": positions,
        "ratios": ratios,
        "sums": _sums(indptr, positions, ratios, width),
        "experimental": np.array([flag(e, "experimental") for e in entries], dtype=bool),
        "structures": np.array([flag(e, "structures") for e in entries], dtype=bool),
    }
    meta = {
        "format_version": FORMAT_VERSION,
        "elmd_version": version("ElMD"),
        "vector_length": str(width),
        "n_entries": str(len(entries)),
        "source": source,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }

    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.unlink(missing_ok=True)
    con = sqlite3.connect(tmp)
    with con:
        con.execute("CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)")
        con.execute("CREATE TABLE arrays (name TEXT PRIMARY KEY, data BLOB)")
        con.execute(
            "CREATE TABLE entries (id INTEGER PRIMARY KEY, formula TEXT, "
            "normalised_formula TEXT, records TEXT)"
        )
        con.execute(
            "CREATE TABLE databases (name TEXT PRIMARY KEY, description TEXT, "
            "experimental INTEGER, structures INTEGER, url TEXT)"
        )
        con.executemany("INSERT INTO meta VALUES (?, ?)", meta.items())
        con.executemany("INSERT INTO arrays VALUES (?, ?)",
                        ((k, _blob(v)) for k, v in arrays.items()))
        con.executemany(
            "INSERT INTO entries VALUES (?, ?, ?, ?)",
            ((i, e["formula"], e["normalised_formula"], json.dumps(e["records"]))
             for i, e in enumerate(entries)),
        )
        used = {db for e in entries for db in e["records"]}
        con.executemany(
            "INSERT INTO databases VALUES (?, ?, ?, ?, ?)",
            ((name, d.get("description", ""), int(bool(d.get("experimental"))),
              int(bool(d.get("structures"))), d.get("url", ""))
             for name, d in sorted(databases.items()) if name in used),
        )
    con.execute("VACUUM")
    con.close()
    tmp.replace(path)


class CompositionIndex:
    def __init__(self, path: Path):
        self.con = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)
        try:
            meta = dict(self.con.execute("SELECT key, value FROM meta"))
        except sqlite3.DatabaseError as e:
            raise ValueError(f"{path.name} is not an ElMTree index file ({e})") from e
        if meta.get("format_version") != FORMAT_VERSION:
            raise ValueError(
                f"{path.name} has index format {meta.get('format_version')!r}, "
                f"expected {FORMAT_VERSION!r}; rebuild it with build-index/convert-index"
            )
        self.width = int(meta["vector_length"])
        if self.width != vector_length():
            raise ValueError(
                f"{path.name} was built with ElMD {meta.get('elmd_version')} "
                f"({self.width}-element vectors); the installed ElMD uses "
                f"{vector_length()}. Rebuild the index with this version."
            )
        self.meta = meta
        arrays = {name: np.load(io.BytesIO(data), allow_pickle=False)
                  for name, data in self.con.execute("SELECT name, data FROM arrays")}
        self.indptr = arrays["indptr"]
        self.positions = arrays["positions"].astype(np.intp)
        self.ratios = arrays["ratios"]
        self.sums = arrays["sums"]
        self.experimental = arrays["experimental"]
        self.structures = arrays["structures"]
        self.n = len(self.indptr) - 1
        self.databases = {
            name: {"description": desc, "experimental": bool(exp),
                   "structures": bool(struct), "url": url}
            for name, desc, exp, struct, url in self.con.execute("SELECT * FROM databases")
        }

        # float32 cumulative vectors for screening (~0.4 KB per entry)
        self.cdf = np.empty((self.n, self.width), dtype=np.float32)
        for a in range(0, self.n, CHUNK):
            b = min(self.n, a + CHUNK)
            self.cdf[a:b] = np.cumsum(self._dense(np.arange(a, b)), axis=1)
        self._mask_cache: tuple | None = None

    def __len__(self) -> int:
        return self.n

    def _dense(self, rows: np.ndarray) -> np.ndarray:
        """float64 ratio vectors of the given rows."""
        out = np.zeros((len(rows), self.width))
        starts = self.indptr[rows]
        lengths = self.indptr[rows + 1] - starts
        offsets = np.cumsum(lengths) - lengths
        nnz = np.arange(lengths.sum()) - np.repeat(offsets, lengths) + np.repeat(starts, lengths)
        out[np.repeat(np.arange(len(rows)), lengths), self.positions[nnz]] = self.ratios[nnz]
        return out

    def mask(self, filters: dict | None, petti_lookup: dict) -> np.ndarray | None:
        """Boolean mask of entries passing the advanced-search filters."""
        if not filters:
            return None
        key = (frozenset(filters.get("must_contain", ())), frozenset(filters.get("must_exclude", ())),
               bool(filters.get("experimental")), bool(filters.get("structures")))
        if self._mask_cache and self._mask_cache[0] == key:
            return self._mask_cache[1]

        rows = np.repeat(np.arange(self.n), np.diff(self.indptr))
        keep = np.ones(self.n, dtype=bool)
        for el in key[0]:
            has = np.zeros(self.n, dtype=bool)
            has[rows[self.positions == petti_lookup[el]]] = True
            keep &= has
        for el in key[1]:
            keep[rows[self.positions == petti_lookup[el]]] = False
        if key[2]:
            keep &= self.experimental
        if key[3]:
            keep &= self.structures
        self._mask_cache = (key, keep)
        return keep

    def knn(self, vector: np.ndarray, k: int, mask: np.ndarray | None = None) -> list[tuple[int, float]]:
        """The k nearest entries to a ratio vector as (row, distance), closest first."""
        q = np.cumsum(vector)
        q32 = q.astype(np.float32)
        sq = float(np.sum(vector * (self.width - np.arange(self.width))))

        rows = None if mask is None else np.flatnonzero(mask)
        sums = self.sums if rows is None else self.sums[rows]
        m = len(sums)
        if m == 0:
            return []

        found_rows, found_d = [], []

        def screen(a, b):
            for s in range(a, b, CHUNK):
                e = min(b, s + CHUNK)
                sel = slice(s, e) if rows is None else rows[s:e]
                found_rows.append(np.arange(s, e) if rows is None else rows[s:e])
                found_d.append(np.abs(self.cdf[sel] - q32).sum(axis=1))

        lo = hi = int(np.searchsorted(sums, sq))
        block = FIRST_BLOCK
        cand_rows, cand_d = np.empty(0, dtype=np.intp), np.empty(0)
        while True:
            new_lo, new_hi = max(0, lo - block), min(m, hi + block)
            screen(new_lo, lo)
            screen(hi, new_hi)
            lo, hi = new_lo, new_hi
            cand_rows = np.concatenate([cand_rows, *found_rows])
            cand_d = np.concatenate([cand_d, *found_d])
            found_rows.clear()
            found_d.clear()

            kth = np.inf
            if len(cand_d) >= k:
                kth = float(np.partition(cand_d, k - 1)[k - 1])
                keep = cand_d <= kth + EPS
                cand_rows, cand_d = cand_rows[keep], cand_d[keep]
            bound = min(sq - sums[lo - 1] if lo > 0 else np.inf,
                        sums[hi] - sq if hi < m else np.inf)
            if (lo == 0 and hi == m) or bound > kth + EPS:
                break
            block *= 2

        # Exact float64 distances, as ElMD computes them, for the candidates
        exact = np.abs(np.cumsum(self._dense(cand_rows) - vector, axis=1)).sum(axis=1)
        order = np.lexsort((cand_rows, exact))[:k]
        return [(int(cand_rows[i]), float(exact[i])) for i in order]

    def entries(self, rows: list[int]) -> dict[int, dict]:
        """formula, normalised_formula and records for the given rows."""
        out = {}
        for start in range(0, len(rows), 900):  # SQLite variable limit
            batch = rows[start:start + 900]
            marks = ",".join("?" * len(batch))
            for i, formula, normalised, records in self.con.execute(
                f"SELECT id, formula, normalised_formula, records FROM entries WHERE id IN ({marks})",
                batch,
            ):
                out[i] = {"formula": formula, "normalised_formula": normalised,
                          "records": json.loads(records)}
        return out
