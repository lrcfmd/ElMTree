import csv
import pickle
import random
import sqlite3
from pathlib import Path

import pytest
from ElMD import ElMD

from ElMTree import ElMTree
from infrastructure.logging import newLogger
from model.build_index import main as build_index
from model.convert_index import main as convert_index
from model.index import INDEX_FILE
from model import report
from model.model import ElMTreeSearchModel

LOG = newLogger("test")
ROOT = Path(__file__).resolve().parent.parent
SAMPLE = ROOT / "data" / "sample_compositions.csv"  # matbench_expt_gap
INPUT_CSV = ROOT / "data" / "input.csv"
FAST = ElMD("", metric="fast").elmd


@pytest.fixture(scope="session")
def sample_rows() -> list[dict]:
    with SAMPLE.open() as f:
        return list(csv.DictReader(f))[:1500]


@pytest.fixture(scope="session")
def index_dir(tmp_path_factory, sample_rows) -> Path:
    """A real index built from 1500 sample compositions."""
    tmp = tmp_path_factory.mktemp("index")
    subset = tmp / "subset.csv"
    with subset.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "composition"])
        w.writeheader()
        w.writerows(sample_rows)
    out = tmp / "artifacts"
    build_index([str(subset), "--database", "matbench", "--experimental", "--out", str(out)])
    return out


@pytest.fixture
def model(index_dir, monkeypatch) -> ElMTreeSearchModel:
    monkeypatch.setenv("MODEL_ARTIFACTS_DIR", str(index_dir))
    return ElMTreeSearchModel()


@pytest.fixture(scope="session")
def entries(sample_rows) -> list[ElMD]:
    unique = {}
    for row in sample_rows:
        e = ElMD(row["composition"], metric="fast")
        unique.setdefault(e.pretty_formula, e)
    return list(unique.values())


def _brute(query: str, entries, k, keep=lambda e: True) -> list[float]:
    q = ElMD(query, metric="fast")
    return sorted(round(float(FAST(q, e)), 9) for e in entries if keep(e))[:k]


def test_knn_is_exact(model, entries) -> None:
    # The search must match brute force over every indexed composition.
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


def test_instance_results_page(model, tmp_path) -> None:
    model.process("instance", {"composition": "Fe2O3", "k": 12, "must_contain": ["O"]},
                  {}, tmp_path, {}, LOG)
    page = (tmp_path / "elmtree_results.html").read_text()
    assert page.startswith("<!doctype html>")
    assert "for <b>Fe<sub>0.4</sub> O<sub>0.6</sub></b>" in page
    assert "Filters: must contain O" in page
    assert page.count('<td class="rank">') == 12
    assert "Experimental structures reported" in page  # legend
    assert "http" not in page  # self-contained: the frontend sandboxes it


def test_results_page_folds_long_id_lists() -> None:
    n = {"formula": "NaCl", "normalised_formula": "Na0.5 Cl0.5", "distance": 0.0,
         "records": [("icsd", list(range(20)))]}
    page = report.instance_page({"query": "Na0.5 Cl0.5", "neighbours": [n], "seconds": 0.01},
                                None, {"compositions": 1, "records": 20}, {})
    assert "0, 1, 2, 3, 4, 5, 6, 7<details><summary> … +12 more</summary>, 8, 9" in page
    assert report.subscript("Li(Fe0.97Mn0.03)(PO4)") == \
        "Li(Fe<sub>0.97</sub>Mn<sub>0.03</sub>)(PO<sub>4</sub>)"


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


def test_library_elmtree_agrees_with_index(model, entries) -> None:
    # The ElMTree library's list-of-clusters search (with the kNN pruning fix)
    # and the tool's index are both exact, so they must agree.
    random.seed(1)
    tree = ElMTree([e.formula for e in entries[:600]], pre_process=True, max_workers=2)
    for q in ["NaCl", "Fe2O3", "LiCoO2"] + [e.formula for e in random.sample(entries[:600], 15)]:
        got = [round(d, 9) for _, d, _ in tree.knn(q, 5)]
        assert got == _brute(q, entries[:600], 5), q


def test_experimental_and_structures_filters(model) -> None:
    # The sample index's one database was built with --experimental only.
    assert len(model.search("NaCl", 5, {"experimental": True})["neighbours"]) == 5
    assert model.search("NaCl", 5, {"structures": True})["neighbours"] == []


def test_missing_index(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("MODEL_ARTIFACTS_DIR", str(tmp_path))
    with pytest.raises(FileNotFoundError, match=INDEX_FILE):
        ElMTreeSearchModel().search("NaCl", 3)


def test_index_format_is_checked(index_dir, tmp_path, monkeypatch) -> None:
    for key, value, match in [("format_version", "0", "index format"),
                              ("vector_length", "99", "Rebuild the index")]:
        bad = tmp_path / key
        bad.mkdir()
        (bad / INDEX_FILE).write_bytes((index_dir / INDEX_FILE).read_bytes())
        with sqlite3.connect(bad / INDEX_FILE) as con:
            con.execute("UPDATE meta SET value = ? WHERE key = ?", (value, key))
        monkeypatch.setenv("MODEL_ARTIFACTS_DIR", str(bad))
        with pytest.raises(ValueError, match=match):
            ElMTreeSearchModel().search("NaCl", 3)
    not_index = tmp_path / "not_index"
    not_index.mkdir()
    (not_index / INDEX_FILE).write_bytes(b"x" * 4096)
    monkeypatch.setenv("MODEL_ARTIFACTS_DIR", str(not_index))
    with pytest.raises(ValueError, match="not an ElMTree index"):
        ElMTreeSearchModel().search("NaCl", 3)


def test_convert_index(tmp_path, monkeypatch) -> None:
    # A miniature of the original server's mtree_lookup.pk / db_lookup.pk
    lookup = {
        "Na0.5 Cl0.5": {"aflow": ["z"], "mpds_struc": ["S9"], "icsd": ["1", "2"],
                        "compound_formula": ["'NaCl'", "Na1Cl1"]},
        # rounded key: the vector must come from the exact compound formula
        "Cu0.667 Al0.333": {"mp": ["mp-1"], "compound_formula": ["AlCu2\n"]},
        "Fe0.4 O0.6": {"aflow": ["a1"], "compound_formula": ["junk("]},  # falls back to key
        "K0.5 Br0.5": {"mpds_struc": ["S1"], "compound_formula": ["KBr"]},
        "": {"ccdc": ["X"], "compound_formula": [""]},  # skipped
    }
    dbs = {
        "icsd": {"description": "ICSD", "display": "ICSD", "experimental": True,
                 "structures": True, "url": "u"},
        "ccdc": {"description": "CSD", "experimental": True, "structures": True, "url": "u"},
        "mp": {"description": "MP", "experimental": False, "structures": True, "url": "u"},
        "aflow": {"description": "AFLOW", "experimental": False, "structures": False, "url": "u"},
        "mpds_struc": {"description": "MPDS", "display": "MPDS Structural Database",
                       "experimental": True, "structures": False, "url": "u"},
    }
    for name, obj in [("mtree_lookup.pk", lookup), ("db_lookup.pk", dbs)]:
        with (tmp_path / name).open("wb") as f:
            pickle.dump(obj, f)
    out = tmp_path / "artifacts"
    convert_index([str(tmp_path / "mtree_lookup.pk"), str(tmp_path / "db_lookup.pk"),
                   "--out", str(out), "--max-workers", "1"])
    monkeypatch.setenv("MODEL_ARTIFACTS_DIR", str(out))
    model = ElMTreeSearchModel()

    result = model.search("NaCl", 10)["neighbours"]
    assert len(result) == 4
    nacl = result[0]
    assert (nacl["formula"], nacl["normalised_formula"], nacl["distance"]) == ("NaCl", "Na0.5 Cl0.5", 0.0)
    # ordered as the original site: experimental structures, experimental, computational
    assert nacl["databases"] == ["icsd: 1 2", "mpds_struc: S9", "aflow: z"]
    alcu = next(n for n in model.search("AlCu2", 10)["neighbours"] if n["formula"] == "AlCu2")
    assert alcu["distance"] == 0.0  # exact, not the rounded 0.333/0.667
    fallback = next(n for n in result if n["normalised_formula"] == "Fe0.4 O0.6")
    assert fallback["formula"] == "Fe0.4 O0.6"

    def formulas(filters):
        return sorted(n["formula"] for n in model.search("NaCl", 10, filters)["neighbours"])
    assert formulas({"experimental": True}) == ["KBr", "NaCl"]
    assert formulas({"structures": True}) == ["AlCu2", "NaCl"]
    assert formulas({"experimental": True, "structures": True}) == ["NaCl"]

    model.process("instance", {"composition": "NaCl", "k": 4}, {}, tmp_path / "o", {}, LOG)
    page = (tmp_path / "o" / "elmtree_results.html").read_text()
    assert "Searched 4 unique compositions across 7 records" in page
    assert page.index(">ICSD<") < page.index(">MPDS Structural Database<") < page.index(">aflow<")


def test_convert_index_rejects_objects(tmp_path) -> None:
    with (tmp_path / "mtree_lookup.pk").open("wb") as f:
        pickle.dump({"Na0.5 Cl0.5": {"x": Path("/")}}, f)
    with (tmp_path / "db_lookup.pk").open("wb") as f:
        pickle.dump({}, f)
    with pytest.raises(pickle.UnpicklingError, match="unexpected object"):
        convert_index([str(tmp_path / "mtree_lookup.pk"), str(tmp_path / "db_lookup.pk"),
                       "--out", str(tmp_path / "o")])


def test_databases_formatting() -> None:
    records = {"icsd": [1, 2], "mpds": ["S1-a", "S1-b", "S2-a"], "compound_formula": "X"}
    assert ElMTreeSearchModel._databases(records) == ["icsd: 1 2", "mpds: S1 S2"]
    assert ElMTreeSearchModel._databases({}) == []


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
    page = (out / "elmtree_results.html").read_text()
    assert page.count('<section class="query">') == 4
    assert "Unknown element" in page and page.count('<td class="rank">') == 4


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
