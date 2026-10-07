import csv
import pickle
import random
from pathlib import Path

import pytest
from ElMD import ElMD

from infrastructure.logging import newLogger
from model.build_index import main as build_index
from model.model import ElMTreeSearchModel

LOG = newLogger("test")
ROOT = Path(__file__).resolve().parent.parent
SAMPLE = ROOT / "data" / "sample_compositions.csv"  # matbench_expt_gap
INPUT_CSV = ROOT / "data" / "input.csv"
FAST = ElMD("", metric="fast").elmd


@pytest.fixture(scope="session")
def index_dir(tmp_path_factory) -> Path:
    """A real index + lookup built from 1500 sample compositions."""
    tmp = tmp_path_factory.mktemp("index")
    with SAMPLE.open() as f:
        rows = list(csv.DictReader(f))[:1500]
    subset = tmp / "subset.csv"
    with subset.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "composition"])
        w.writeheader()
        w.writerows(rows)
    out = tmp / "artifacts"
    build_index([str(subset), "--database", "matbench", "--out", str(out), "--max-workers", "2"])
    return out


@pytest.fixture
def model(index_dir, monkeypatch) -> ElMTreeSearchModel:
    monkeypatch.setenv("MODEL_ARTIFACTS_DIR", str(index_dir))
    return ElMTreeSearchModel()


@pytest.fixture(scope="session")
def entries(index_dir) -> list[ElMD]:
    with (index_dir / "indexed_LC.pk").open("rb") as f:
        tree = pickle.load(f)
    return [leaf.entry for centre in tree.centres for leaf in centre[1]]


def _brute(query: str, entries, k, keep=lambda e: True) -> list[float]:
    q = ElMD(query, metric="fast")
    return sorted(round(float(FAST(q, e)), 9) for e in entries if keep(e))[:k]


def test_knn_is_exact(model, entries) -> None:
    # The original pruning test (abs(d - r) <= bound) returned wrong
    # neighbours for ~20% of queries; the fixed index must match brute force.
    random.seed(0)
    queries = [e.pretty_formula for e in random.sample(entries, 60)]
    queries += ["NaCl", "LiFePO4", "Fe2O3", "Bi2Te3", "Ca(OH)2"]
    for q in queries:
        got = [round(n["distance"], 9) for n in model.search(q, 10)["neighbours"]]
        assert got == _brute(q, entries, 10), q


def test_element_filters_are_exact(model, entries) -> None:
    filters = {"must_contain": {"O"}, "must_exclude": {"Fe", "Co"},
               "experimental": False, "structures": False}

    def keep(e):
        els = set(e.normed_composition)
        return "O" in els and not els & {"Fe", "Co"}

    for q in ["Fe2O3", "LiCoO2", "NaCl"]:
        result = model.search(q, 8, filters)["neighbours"]
        assert [round(n["distance"], 9) for n in result] == _brute(q, entries, 8, keep)
        for n in result:
            els = set(ElMD(n["formula"]).normed_composition)
            assert "O" in els and not els & {"Fe", "Co"}


def _neighbours_csv(out: Path) -> list[dict]:
    text = (out / "neighbours.csv").read_text()
    # The LMDS frontend splits CSV rows on every comma, so cells must have none.
    widths = {len(line.split(",")) for line in text.strip().splitlines()}
    assert len(widths) == 1, "a CSV cell contains a comma"
    return list(csv.DictReader(text.splitlines()))


def test_instance(model, tmp_path) -> None:
    result = model.process("instance", {"composition": "Fe2O3", "k": 5}, {}, tmp_path, {}, LOG)
    assert result == {
        "query": "Fe2O3", "normalised_query": "Fe0.4 O0.6", "matches": 5,
        "nearest_match": "Fe2O3", "nearest_distance": 0.0,
    }
    rows = _neighbours_csv(tmp_path)
    assert [r["rank"] for r in rows] == ["1", "2", "3", "4", "5"]
    assert rows[0]["formula"] == "Fe2O3" and rows[0]["normalised_formula"] == "Fe0.4 O0.6"
    assert rows[0]["databases"].startswith("matbench: ")
    distances = [float(r["distance"]) for r in rows]
    assert distances == sorted(distances)
    assert all(len(r["distance"].split(".")[-1]) <= 3 for r in rows)  # 3 dp


def test_instance_default_k(model, tmp_path) -> None:
    result = model.process("instance", {"composition": "NaCl"}, {}, tmp_path, {}, LOG)
    assert result["matches"] == 100
    assert len(_neighbours_csv(tmp_path)) == 100


def test_instance_filters_via_process(model, tmp_path) -> None:
    model.process(
        "instance",
        {"composition": "Fe2O3", "k": 3, "must_contain": ["Li"], "must_exclude": "Co, Ni"},
        {}, tmp_path, {}, LOG,
    )
    rows = _neighbours_csv(tmp_path)
    assert len(rows) == 3
    for n in rows:
        els = set(ElMD(n["formula"]).normed_composition)
        assert "Li" in els and not els & {"Co", "Ni"}


def test_filters_excluding_everything(model, tmp_path) -> None:
    result = model.process(
        "instance", {"composition": "NaCl", "must_contain": ["Og"]}, {}, tmp_path, {}, LOG
    )
    assert result["matches"] == 0 and result["nearest_match"] == ""
    assert _neighbours_csv(tmp_path) == []


@pytest.mark.parametrize("bad, match", [
    ("", "non-empty"),
    ("   ", "non-empty"),
    ("nacl", "Could not parse"),
    ("Xx2O", "Unknown element"),     # ElMD would silently search for pure O
    ("NotAnElement", "Unknown element"),
    (42, "non-empty"),
])
def test_invalid_composition(model, bad, match) -> None:
    with pytest.raises(ValueError, match=match):
        model.process("instance", {"composition": bad}, {}, Path("/tmp/elm_out"), {}, LOG)


@pytest.mark.parametrize("values, match", [
    ({"k": 0}, "k must be"),
    ({"k": 1001}, "k must be"),
    ({"k": "5"}, "k must be"),
    ({"must_contain": ["Zz"]}, "unknown element"),
    ({"must_contain": ["O"], "must_exclude": ["O"]}, "both"),
])
def test_invalid_options(model, values, match) -> None:
    with pytest.raises(ValueError, match=match):
        model.process("instance", {"composition": "NaCl", **values}, {}, Path("/tmp/elm_out"), {}, LOG)


def test_instance_missing_composition(model) -> None:
    with pytest.raises(ValueError, match="No composition"):
        model.process("instance", {}, {}, Path("/tmp/elm_out"), {}, LOG)


def test_without_lookup(index_dir, tmp_path, monkeypatch) -> None:
    only_index = tmp_path / "artifacts"
    only_index.mkdir()
    (only_index / "indexed_LC.pk").write_bytes((index_dir / "indexed_LC.pk").read_bytes())
    monkeypatch.setenv("MODEL_ARTIFACTS_DIR", str(only_index))
    result = ElMTreeSearchModel().search("NaCl", 3)
    assert len(result["neighbours"]) == 3
    assert all(n["databases"] == [] for n in result["neighbours"])


def test_missing_index(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("MODEL_ARTIFACTS_DIR", str(tmp_path))
    with pytest.raises(FileNotFoundError, match="indexed_LC.pk"):
        ElMTreeSearchModel().search("NaCl", 3)


def test_databases_formatting(model) -> None:
    model._load()
    model.lookup = {"X": {"icsd": [1, 2], "mpds": ["S1-a", "S1-b", "S2-a"], "compound_formula": "X"}}
    assert model._databases("X", None) == ["icsd: 1 2", "mpds: S1 S2"]
    assert model._databases("missing", None) == []


def test_dataset(model, tmp_path) -> None:
    out = tmp_path / "output"
    result = model.process(
        "dataset", {"k": 3}, {"file": (INPUT_CSV, "text/csv")}, out, {}, LOG
    )
    assert result == {}
    rows = _neighbours_csv(out)
    with INPUT_CSV.open() as f:
        n_queries = sum(1 for _ in csv.DictReader(f))
    assert len(rows) == 3 * n_queries
    assert all(r["error"] == "" for r in rows)
    assert [r["rank"] for r in rows[:3]] == ["1", "2", "3"]
    single = model.search(rows[0]["query"], 3)["neighbours"]
    assert [r["formula"] for r in rows[:3]] == [n["formula"] for n in single]


def test_dataset_bad_rows_do_not_abort(model, tmp_path) -> None:
    src = tmp_path / "in.csv"
    src.write_text("id,composition\n1,NaCl\n2,Xx2O\n3,\n4,GaAs\n")
    out = tmp_path / "output"
    model.process("dataset", {"k": 2}, {"file": (src, "text/csv")}, out, {}, LOG)
    rows = _neighbours_csv(out)
    by_query = {}
    for r in rows:
        by_query.setdefault(r["query_row"], []).append(r)
    assert len(by_query["0"]) == 2 and by_query["0"][0]["error"] == ""
    assert len(by_query["1"]) == 1 and "Unknown element" in by_query["1"][0]["error"]
    assert len(by_query["2"]) == 1 and by_query["2"][0]["error"] != ""
    assert len(by_query["3"]) == 2 and by_query["3"][0]["error"] == ""


def test_dataset_missing_column(model, tmp_path) -> None:
    src = tmp_path / "in.csv"
    src.write_text("formula\nNaCl\n")
    with pytest.raises(ValueError, match="'composition' column"):
        model.process("dataset", {}, {"file": (src, "text/csv")}, tmp_path / "o", {}, LOG)


def test_dataset_no_file(model) -> None:
    with pytest.raises(ValueError, match="No input file"):
        model.process("dataset", {}, {}, Path("/tmp/elm_out"), {}, LOG)


def test_unknown_mode(model) -> None:
    with pytest.raises(ValueError, match="Unknown mode"):
        model.process("bogus", {}, {}, Path("/tmp/elm_out"), {}, LOG)
