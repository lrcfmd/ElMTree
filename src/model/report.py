"""HTML results page for ElMTree searches.

Each neighbour is one row: its composition, its distance (with a bar, for
scanning), and chips for the databases it is reported in. A chip's colour is
the kind of database, hovering it describes the database and clicking it opens
the database's website. Record IDs are folded away per row. The page is
self-contained (no external assets), because the LMDS frontend shows HTML
outputs in a sandboxed iframe.
"""

import json
import re
from html import escape
from urllib.parse import urlparse

# Kinds of database, in the order chips are listed
KINDS = {
    "exp_struct": ("Experimental structures", "Experimental · crystal structures"),
    "exp_prop": ("Experimental properties", "Experimental · measured properties"),
    "computed": ("Computational", "Computational · calculated structures or properties"),
}
KIND_ORDER = {kind: i for i, kind in enumerate(KINDS)}

STYLE = """
:root { color-scheme: light; --text:#1f2937; --muted:#6b7280; --faint:#9ca3af; --line:#e5e7eb;
        --soft:#f9fafb; --accent:#2563eb;
        --exp_struct:#c2410c; --exp_struct-bg:#fff7ed; --exp_struct-line:#fed7aa;
        --exp_prop:#1d4ed8; --exp_prop-bg:#eff6ff; --exp_prop-line:#bfdbfe;
        --computed:#15803d; --computed-bg:#f0fdf4; --computed-line:#bbf7d0; }
* { box-sizing:border-box; }
html, body { background:#ffffff; }
body { margin:0; padding:20px; color:var(--text);
       font:14px/1.5 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
h1 { font-size:18px; font-weight:600; margin:0 0 4px; }
h2 { font-size:16px; font-weight:600; margin:0; }
.meta { color:var(--muted); margin:0; }
.filters { display:flex; flex-wrap:wrap; gap:6px; margin:10px 0 0; }
.filter { font-size:12px; padding:2px 8px; border-radius:999px; background:var(--soft);
          border:1px solid var(--line); color:var(--muted); }
.key { display:flex; flex-wrap:wrap; gap:6px 18px; align-items:center; margin:14px 0 16px;
       padding:10px 12px; background:var(--soft); border:1px solid var(--line); border-radius:8px;
       font-size:12px; color:var(--muted); }
.key .dot { display:inline-block; width:10px; height:10px; border-radius:50%; margin-right:6px;
            vertical-align:-1px; }
table { width:100%; border-collapse:collapse; }
th { text-align:left; font-size:12px; font-weight:600; color:var(--muted); text-transform:uppercase;
     letter-spacing:.04em; padding:8px 10px; border-bottom:1px solid var(--line); }
td { vertical-align:top; padding:12px 10px; border-bottom:1px solid var(--line); }
td.rank { color:var(--faint); width:1%; font-variant-numeric:tabular-nums; }
td.comp { width:26%; }
td.dist { width:15%; }
.formula { font-size:15px; font-weight:600; }
.norm { color:var(--muted); font-size:12px; margin-top:2px; overflow-wrap:anywhere; }
.dval { font-variant-numeric:tabular-nums; font-weight:600; }
.bar { height:4px; background:var(--line); border-radius:2px; margin:6px 0 2px; max-width:120px; }
.bar span { display:block; height:100%; background:var(--accent); border-radius:2px; }
.dnote { color:var(--muted); font-size:12px; }
.chips { display:flex; flex-wrap:wrap; gap:6px; }
.chip { display:inline-flex; align-items:center; gap:6px; padding:3px 9px; border-radius:6px;
        border:1px solid; font-size:13px; text-decoration:none; white-space:nowrap; }
a.chip:hover, a.chip:focus-visible { filter:brightness(.95); text-decoration:underline; outline:none; }
.chip .n { font-size:11px; font-weight:600; padding:0 6px; border-radius:999px; background:#ffffff; }
.chip.exp_struct { color:var(--exp_struct); background:var(--exp_struct-bg); border-color:var(--exp_struct-line); }
.chip.exp_prop { color:var(--exp_prop); background:var(--exp_prop-bg); border-color:var(--exp_prop-line); }
.chip.computed { color:var(--computed); background:var(--computed-bg); border-color:var(--computed-line); }
details.ids { margin-top:8px; font-size:12px; }
details.ids summary { cursor:pointer; color:var(--accent); width:max-content; }
details.ids .idlist { margin:6px 0 0; padding:8px 10px; background:var(--soft); border-radius:6px;
                     max-height:240px; overflow:auto; }
details.ids p { margin:0 0 6px; overflow-wrap:anywhere; }
details.ids p:last-child { margin:0; }
details.ids b { color:var(--muted); font-weight:600; margin-right:6px; }
.empty { color:var(--muted); padding:16px 0; }
#tip { position:fixed; z-index:10; display:none; max-width:340px; padding:10px 12px; background:#ffffff;
       border:1px solid var(--line); border-radius:8px; box-shadow:0 6px 24px rgba(0,0,0,.12);
       font-size:12px; color:var(--text); pointer-events:none; }
#tip b { font-size:13px; }
#tip .kind { margin:2px 0 6px; font-weight:600; }
#tip .desc { color:var(--muted); margin:0 0 6px; }
#tip .foot { color:var(--faint); }
section.query { margin:28px 0 0; padding-top:4px; }
section.query .head { display:flex; flex-wrap:wrap; align-items:baseline; gap:4px 10px; margin-bottom:8px; }
.overview td { padding:8px 10px; }
.overview a { color:var(--accent); text-decoration:none; }
.overview a:hover { text-decoration:underline; }
.error { color:#b42318; }
.note { color:var(--muted); font-size:12px; margin:8px 0 0; }
"""

# A single floating card describes the hovered or focused database chip.
SCRIPT = """
(function () {
  var info = JSON.parse(document.getElementById("db-info").textContent);
  var tip = document.getElementById("tip");
  function esc(s) { var d = document.createElement("div"); d.textContent = s; return d.innerHTML; }
  function show(chip) {
    var db = info[chip.dataset.db]; if (!db) return;
    var n = chip.dataset.n;
    tip.innerHTML = "<b>" + esc(db.name) + "</b>" +
      '<div class="kind" style="color:var(--' + db.kind + ')">' + esc(db.kind_label) + "</div>" +
      (db.description ? '<p class="desc">' + esc(db.description) + "</p>" : "") +
      '<div class="foot">' + n + (n === "1" ? " record" : " records") + " for this composition" +
      (db.host ? " · click to open " + esc(db.host) + " \\u2197" : "") + "</div>";
    tip.style.display = "block";
    var r = chip.getBoundingClientRect(), w = tip.offsetWidth, h = tip.offsetHeight;
    var left = Math.min(Math.max(8, r.left), window.innerWidth - w - 8);
    var top = r.bottom + 6;
    if (top + h > window.innerHeight - 8) top = Math.max(8, r.top - h - 6);
    tip.style.left = left + "px"; tip.style.top = top + "px";
  }
  function hide() { tip.style.display = "none"; }
  document.querySelectorAll(".chip").forEach(function (chip) {
    chip.addEventListener("mouseenter", function () { show(chip); });
    chip.addEventListener("focus", function () { show(chip); });
    chip.addEventListener("mouseleave", hide);
    chip.addEventListener("blur", hide);
  });
  window.addEventListener("scroll", hide, true);
})();
"""


def kind(meta: dict | None) -> str:
    meta = meta or {}
    if meta.get("experimental"):
        return "exp_struct" if meta.get("structures") else "exp_prop"
    return "computed"


def subscript(formula: str) -> str:
    """Escaped formula with its amounts as subscripts."""
    return re.sub(r"(\d+\.?\d*|\.\d+)", r"<sub>\1</sub>", escape(formula, quote=False))


def tidy(formula: str) -> str:
    """Drop unit amounts as chemists write formulas: Li1Fe1P1O4 -> LiFePO4."""
    return re.sub(r"(?<=[A-Za-z)\]])1(?![\d.])", "", formula)


def normalised(formula: str) -> str:
    """ElMD's fractional formula, with amounts it rounds to "0." shown as <0.001."""
    return re.sub(r"(?<=[A-Za-z])0\.(?=\s|$)", "<0.001", formula)


def _db_info(databases: dict[str, dict]) -> str:
    info = {}
    for name, meta in databases.items():
        url = meta.get("url") or ""
        k = kind(meta)
        info[name] = {
            "name": meta.get("display") or name,
            "kind": k,
            "kind_label": KINDS[k][1],
            "description": (meta.get("description") or "").strip(),
            "host": urlparse(url).netloc if url.startswith(("http://", "https://")) else "",
        }
    # Safe inside <script>: no "</" can close the element early
    return json.dumps(info).replace("</", "<\\/")


def _chip(db: str, ids: list, databases: dict[str, dict]) -> str:
    meta = databases.get(db, {})
    url = meta.get("url") or ""
    name = escape(meta.get("display") or db)
    attrs = f'class="chip {kind(meta)}" data-db="{escape(db)}" data-n="{len(ids)}"'
    body = f'{name}<span class="n">{len(ids):,}</span>'
    if url.startswith(("http://", "https://")):
        return (f'<a {attrs} href="{escape(url)}" target="_blank" '
                f'rel="noopener noreferrer">{body}</a>')
    return f'<span {attrs} tabindex="0">{body}</span>'


def _ids(records: list, databases: dict[str, dict]) -> str:
    total = sum(len(ids) for _, ids in records)
    if not total:
        return ""
    rows = "".join(
        f'<p><b>{escape(databases.get(db, {}).get("display") or db)}</b>'
        f'{escape(", ".join(map(str, ids)))}</p>'
        for db, ids in records)
    label = f"Show {total:,} record ID{'s' if total != 1 else ''}"
    return f'<details class="ids"><summary>{label}</summary><div class="idlist">{rows}</div></details>'


def _distance(d: float, longest: float) -> str:
    width = 0 if longest <= 0 else min(100, 100 * d / longest)
    note = '<div class="dnote">identical</div>' if d == 0 else ""
    return (f'<div class="dval">{d:.3f}</div>'
            f'<div class="bar"><span style="width:{width:.1f}%"></span></div>{note}')


def _table(neighbours: list[dict], databases: dict[str, dict]) -> str:
    if not neighbours:
        return '<p class="empty">No compositions match these filters.</p>'
    longest = max(n["distance"] for n in neighbours)
    rows = []
    for rank, n in enumerate(neighbours, start=1):
        formula = tidy(n["formula"])
        norm = ""
        if formula != n["normalised_formula"]:
            norm = f'<div class="norm">{escape(normalised(n["normalised_formula"]))}</div>'
        chips = "".join(_chip(db, ids, databases) for db, ids in n["records"])
        rows.append(
            f'<tr><td class="rank">{rank}</td>'
            f'<td class="comp"><div class="formula">{subscript(formula)}</div>{norm}</td>'
            f'<td class="dist">{_distance(n["distance"], longest)}</td>'
            f'<td><div class="chips">{chips}</div>{_ids(n["records"], databases)}</td></tr>'
        )
    return ('<table><thead><tr><th>#</th><th>Composition</th><th>Distance</th>'
            f'<th>Reported in</th></tr></thead><tbody>{"".join(rows)}</tbody></table>')


def _key() -> str:
    kinds = "".join(f'<span><span class="dot" style="background:var(--{k})"></span>{escape(label)}</span>'
                    for k, (label, _) in KINDS.items())
    return (f'<div class="key">{kinds}<span>Distance: Element Movers Distance, '
            "0 = identical composition · hover a database for details, click to open it</span></div>")


def _filters(filters: dict | None) -> str:
    if not filters:
        return ""
    parts = []
    if filters.get("must_contain"):
        parts.append("contains " + ", ".join(sorted(filters["must_contain"])))
    if filters.get("must_exclude"):
        parts.append("excludes " + ", ".join(sorted(filters["must_exclude"])))
    if filters.get("experimental"):
        parts.append("experimental databases only")
    if filters.get("structures"):
        parts.append("databases with structures only")
    return '<div class="filters">' + "".join(
        f'<span class="filter">{escape(p)}</span>' for p in parts) + "</div>"


def _scope(stats: dict) -> str:
    text = f"{stats['compositions']:,} compositions from {stats['databases']:,} database"
    text += "s" if stats["databases"] != 1 else ""
    if stats.get("records"):
        text += f" ({stats['records']:,} records)"
    return text


def _page(title: str, body: str, databases: dict[str, dict]) -> str:
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width, initial-scale=1">'
            f'<meta name="color-scheme" content="light">'
            f'<title>{escape(title)}</title><style>{STYLE}</style></head><body>{body}'
            f'<div id="tip" role="tooltip"></div>'
            f'<script type="application/json" id="db-info">{_db_info(databases)}</script>'
            f'<script>{SCRIPT}</script></body></html>')


def instance_page(result: dict, filters: dict | None, stats: dict, databases: dict[str, dict]) -> str:
    matches = len(result["neighbours"])
    head = (f'<h1>Compositions most similar to {subscript(result["input"])}</h1>'
            f'<p class="meta">Compared as {escape(result["query"])} · '
            f'{matches:,} closest of {_scope(stats)} · {result["seconds"] * 1000:.0f} ms</p>'
            + _filters(filters))
    return _page("ElMTree results", head + _key() + _table(result["neighbours"], databases), databases)


def dataset_page(queries: list[dict], filters: dict | None, stats: dict,
                 databases: dict[str, dict], max_queries: int) -> str:
    """queries: [{"row", "query", "result" or "error"}]"""
    shown = queries[:max_queries]
    failed = sum("error" in q for q in queries)
    head = (f'<h1>Nearest compositions for {len(queries):,} '
            f'quer{"ies" if len(queries) != 1 else "y"}</h1>'
            f'<p class="meta">Searched {_scope(stats)}'
            + (f' · {failed:,} could not be read' if failed else "") + "</p>" + _filters(filters))
    if len(queries) > len(shown):
        head += (f'<p class="note">Showing the first {len(shown):,} queries; '
                 "neighbours.csv has them all.</p>")

    overview, sections = [], []
    for q in shown:
        anchor = f'q{q["row"] + 1}'
        query = escape(str(q["query"]), quote=False)
        if "error" in q:
            overview.append(f'<tr><td class="rank">{q["row"] + 1}</td><td>{query}</td>'
                            f'<td colspan="2" class="error">{escape(q["error"])}</td></tr>')
            continue
        neighbours = q["result"]["neighbours"]
        nearest = (f'{subscript(tidy(neighbours[0]["formula"]))}</td><td>{neighbours[0]["distance"]:.3f}'
                   if neighbours else '<span class="meta">no matches</span></td><td>')
        overview.append(f'<tr><td class="rank">{q["row"] + 1}</td>'
                        f'<td><a href="#{anchor}">{subscript(str(q["query"]))}</a></td>'
                        f'<td>{nearest}</td></tr>')
        sections.append(
            f'<section class="query" id="{anchor}"><div class="head">'
            f'<h2>{q["row"] + 1}. {subscript(str(q["query"]))}</h2>'
            f'<span class="meta">compared as {escape(q["result"]["query"])}</span></div>'
            f'{_table(neighbours, databases)}</section>')
    table = ('<table class="overview"><thead><tr><th>#</th><th>Query</th><th>Nearest match</th>'
             f'<th>Distance</th></tr></thead><tbody>{"".join(overview)}</tbody></table>')
    return _page("ElMTree results", head + _key() + table + "".join(sections), databases)
