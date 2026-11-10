"""HTML results page, laid out like the original ElMTree web search.

Each neighbour is a row with its composition (subscripted), icons for the kinds
of database it was reported in, its records grouped by database (ordered and
coloured by kind, as on the original site) and its distance. Database names
open the database's website; hovering one describes the database. The page is
self-contained (no external assets), because the LMDS frontend shows HTML
outputs in a sandboxed iframe.
"""

import json
import re
from html import escape
from urllib.parse import urlparse

IDS_SHOWN = 8  # ids per database before the rest fold behind "+N more"

# Database kinds, in the order the original site listed them
KINDS = [
    ("exp_struct", "Experimental structures reported", "#d9773f"),
    ("exp_prop", "Experimentally characterised properties reported", "#2f6fdb"),
    ("computed", "Computationally derived properties or structures reported", "#3f9a4a"),
]
KIND_ORDER = {kind: i for i, (kind, _, _) in enumerate(KINDS)}

ICONS = {
    # crystal unit cell
    "exp_struct": '<path d="M5 8 12 4l7 4v8l-7 4-7-4z" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linejoin="round"/>'
                  '<path d="m5 8 7 4 7-4M12 12v8" fill="none" stroke="currentColor" stroke-width="1.6"/>'
                  '<circle cx="12" cy="4" r="1.6"/><circle cx="5" cy="8" r="1.6"/><circle cx="19" cy="8" r="1.6"/>'
                  '<circle cx="12" cy="12" r="1.6"/><circle cx="5" cy="16" r="1.6"/><circle cx="19" cy="16" r="1.6"/>'
                  '<circle cx="12" cy="20" r="1.6"/>',
    # flask
    "exp_prop": '<path d="M9.5 3h5M10.5 3v6L5 19a1.3 1.3 0 0 0 1.1 2h11.8a1.3 1.3 0 0 0 1.1-2L13.5 9V3" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linejoin="round"/>'
                '<path d="M7.6 15h8.8l2.3 4.4H5.3z"/>',
    # monitor
    "computed": '<rect x="3" y="4" width="18" height="12" rx="1.5" fill="none" stroke="currentColor" stroke-width="1.6"/>'
                '<path d="M9 20h6M12 16v4" stroke="currentColor" stroke-width="1.6"/>'
                '<path d="m8 12 3-3 2 2 3-3" fill="none" stroke="currentColor" stroke-width="1.6"/>',
}

STYLE = """
:root { color-scheme: light; --text:#1d2330; --muted:#667085; --line:#e4e7ec; --stripe:#f5f7fb; --bg:#ffffff; }
* { box-sizing: border-box; }
html { background:var(--bg); }
body { margin:0; padding:16px; background:var(--bg); color:var(--text);
       font:14px/1.45 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
.summary { margin:0 0 4px; }
.filters, .note { margin:0 0 8px; color:var(--muted); font-size:13px; }
.legend { display:flex; flex-wrap:wrap; gap:6px 18px; margin:8px 0 14px; font-size:13px; color:var(--muted); }
.legend span { display:inline-flex; align-items:center; gap:6px; }
table { width:100%; border-collapse:collapse; }
th { text-align:left; font-weight:600; padding:8px 10px; border-bottom:2px solid var(--line); white-space:nowrap; }
td { vertical-align:top; padding:10px; border-bottom:1px solid var(--line); }
tbody tr:nth-child(even) { background:var(--stripe); }
td.rank { color:var(--muted); width:1%; }
td.comp { width:28%; }
td.dist { width:13%; white-space:nowrap; }
.dval { font-variant-numeric:tabular-nums; font-weight:600; }
.bar { height:4px; background:var(--line); border-radius:2px; margin:6px 0 2px; max-width:110px; }
.bar span { display:block; height:100%; background:#2563eb; border-radius:2px; }
.dnote { color:var(--muted); font-size:12px; }
.formula { font-size:15px; }
.reported { color:var(--muted); font-size:12px; margin-top:2px; overflow-wrap:anywhere; }
.icons { display:flex; gap:6px; margin-top:6px; }
svg.icon { width:22px; height:22px; fill:currentColor; }
.db { margin:0 0 3px; overflow-wrap:anywhere; }
.db .dbname { font-weight:600; text-decoration:none; }
.db a.dbname:hover, .db a.dbname:focus-visible { text-decoration:underline; outline:none; }
#tip { position:fixed; z-index:10; display:none; max-width:340px; padding:10px 12px; background:#ffffff;
       border:1px solid var(--line); border-radius:8px; box-shadow:0 6px 24px rgba(0,0,0,.12);
       font-size:12px; line-height:1.45; color:var(--text); pointer-events:none; }
#tip b { font-size:13px; }
#tip .kind { margin:2px 0 6px; font-weight:600; }
#tip .desc { color:var(--muted); margin:0 0 6px; }
#tip .foot { color:#98a2b3; }
details { display:inline; }
summary { display:inline; cursor:pointer; color:var(--muted); list-style:none; }
summary::-webkit-details-marker { display:none; }
details[open] summary { display:none; }
.empty { color:var(--muted); padding:16px 0; }
section.query { margin:0 0 28px; }
section.query h2 { font-size:16px; margin:0 0 6px; }
.error { color:#c4320a; }
"""


# One floating card describes the hovered or focused database name.
SCRIPT = """
(function () {
  var info = JSON.parse(document.getElementById("db-info").textContent);
  var tip = document.getElementById("tip");
  function esc(s) { var d = document.createElement("div"); d.textContent = s; return d.innerHTML; }
  function show(el) {
    var db = info[el.dataset.db]; if (!db) return;
    var n = el.dataset.n;
    tip.innerHTML = "<b>" + esc(db.name) + "</b>" +
      '<div class="kind" style="color:' + db.colour + '">' + esc(db.kind) + "</div>" +
      (db.description ? '<p class="desc">' + esc(db.description) + "</p>" : "") +
      '<div class="foot">' + n + (n === "1" ? " record" : " records") + " for this composition" +
      (db.host ? " \u00b7 click to open " + esc(db.host) + " \u2197" : "") + "</div>";
    tip.style.display = "block";
    var r = el.getBoundingClientRect(), w = tip.offsetWidth, h = tip.offsetHeight;
    tip.style.left = Math.min(Math.max(8, r.left), window.innerWidth - w - 8) + "px";
    var top = r.bottom + 6;
    if (top + h > window.innerHeight - 8) top = Math.max(8, r.top - h - 6);
    tip.style.top = top + "px";
  }
  function hide() { tip.style.display = "none"; }
  document.querySelectorAll(".dbname").forEach(function (el) {
    el.addEventListener("mouseenter", function () { show(el); });
    el.addEventListener("focus", function () { show(el); });
    el.addEventListener("mouseleave", hide);
    el.addEventListener("blur", hide);
  });
  window.addEventListener("scroll", hide, true);
})();
"""


def _url(meta: dict) -> str:
    url = meta.get("url") or ""
    return url if url.startswith(("http://", "https://")) else ""


def _db_info(databases: dict[str, dict]) -> str:
    colours = {k: c for k, _, c in KINDS}
    labels = {k: label for k, label, _ in KINDS}
    info = {
        name: {
            "name": meta.get("display") or name,
            "kind": labels[kind(meta)],
            "colour": colours[kind(meta)],
            "description": (meta.get("description") or "").strip(),
            "host": urlparse(_url(meta)).netloc,
        }
        for name, meta in databases.items()
    }
    # Safe inside <script>: no "</" can close the element early
    return json.dumps(info).replace("</", "<\\/")


def kind(meta: dict | None) -> str:
    meta = meta or {}
    if meta.get("experimental"):
        return "exp_struct" if meta.get("structures") else "exp_prop"
    return "computed"


def subscript(formula: str) -> str:
    """Escaped formula with its amounts as subscripts."""
    return re.sub(r"(\d+\.?\d*|\.\d+)", r"<sub>\1</sub>", escape(formula, quote=False))


def normalised(formula: str) -> str:
    """ElMD's fractional formula, with amounts it rounds to "0." shown as <0.001."""
    return re.sub(r"(?<=[A-Za-z])0\.(?=\s|$)", "<0.001", formula)


def _icon(k: str) -> str:
    title = next(label for key, label, _ in KINDS if key == k)
    colour = next(c for key, _, c in KINDS if key == k)
    return (f'<svg class="icon" viewBox="0 0 24 24" style="color:{colour}" role="img">'
            f'<title>{escape(title)}</title>{ICONS[k]}</svg>')


def _ids(ids: list) -> str:
    ids = [escape(str(i)) for i in ids]
    shown = ", ".join(ids[:IDS_SHOWN])
    if len(ids) <= IDS_SHOWN:
        return shown
    return (f'{shown}<details><summary> … +{len(ids) - IDS_SHOWN} more</summary>, '
            f'{", ".join(ids[IDS_SHOWN:])}</details>')


def _distance(d: float, longest: float) -> str:
    """The distance, with a bar relative to the furthest neighbour shown."""
    width = 0 if longest <= 0 else min(100, 100 * d / longest)
    note = '<div class="dnote">identical</div>' if d == 0 else ""
    return (f'<div class="dval">{d:.3f}</div>'
            f'<div class="bar"><span style="width:{width:.1f}%"></span></div>{note}')


def _rows(neighbours: list[dict], databases: dict[str, dict]) -> str:
    rows = []
    longest = max(n["distance"] for n in neighbours)
    for rank, n in enumerate(neighbours, start=1):
        kinds = sorted({kind(databases.get(db)) for db, _ in n["records"]}, key=KIND_ORDER.get)
        reported = ""
        if n["formula"] != n["normalised_formula"]:
            reported = f'<div class="reported">reported as {subscript(n["formula"])}</div>'
        dbs = []
        for db, ids in n["records"]:
            meta = databases.get(db, {})
            colour = dict((k, c) for k, _, c in KINDS)[kind(meta)]
            name = escape(meta.get("display") or db)
            attrs = f'class="dbname" style="color:{colour}" data-db="{escape(db)}" data-n="{len(ids)}"'
            if _url(meta):
                label = (f'<a {attrs} href="{escape(_url(meta))}" target="_blank" '
                         f'rel="noopener noreferrer">{name}</a>')
            else:
                label = f'<b {attrs} tabindex="0">{name}</b>'
            dbs.append(f'<div class="db">{label}: {_ids(ids)}</div>')
        rows.append(
            f'<tr><td class="rank">{rank}</td>'
            f'<td class="comp"><div class="formula">{subscript(normalised(n["normalised_formula"]))}</div>{reported}'
            f'<div class="icons">{"".join(_icon(k) for k in kinds)}</div></td>'
            f'<td class="dist">{_distance(n["distance"], longest)}</td>'
            f'<td>{"".join(dbs)}</td></tr>'
        )
    return "".join(rows)


def _table(neighbours: list[dict], databases: dict[str, dict]) -> str:
    if not neighbours:
        return '<p class="empty">No compositions match these filters.</p>'
    return ('<table><thead><tr><th>#</th><th>Name</th><th>Distance</th>'
            f'<th>Reported in</th></tr></thead><tbody>{_rows(neighbours, databases)}</tbody></table>')


def _legend() -> str:
    return '<div class="legend">' + "".join(
        f"<span>{_icon(k)}{escape(label)}</span>" for k, label, _ in KINDS) + "</div>"


def _filters(filters: dict | None) -> str:
    if not filters:
        return ""
    parts = []
    if filters.get("must_contain"):
        parts.append("must contain " + ", ".join(sorted(filters["must_contain"])))
    if filters.get("must_exclude"):
        parts.append("must exclude " + ", ".join(sorted(filters["must_exclude"])))
    if filters.get("experimental"):
        parts.append("experimental databases only")
    if filters.get("structures"):
        parts.append("databases with structures only")
    return f'<p class="filters">Filters: {escape("; ".join(parts))}</p>'


def _searched(stats: dict) -> str:
    text = f"Searched {stats['compositions']:,} unique compositions"
    if stats.get("records"):
        text += f" across {stats['records']:,} records"
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
    summary = (f'<p class="summary">{_searched(stats)} for '
               f'<b>{subscript(result["query"])}</b> in {result["seconds"]:.3f} s</p>')
    return _page("ElMTree results",
                 summary + _filters(filters) + _legend() + _table(result["neighbours"], databases),
                 databases)


def dataset_page(queries: list[dict], filters: dict | None, stats: dict,
                 databases: dict[str, dict], max_queries: int) -> str:
    """queries: [{"row", "query", "result" or "error"}]"""
    shown = queries[:max_queries]
    head = f'<p class="summary">{_searched(stats)} for {len(queries):,} queries</p>' + _filters(filters)
    if len(queries) > len(shown):
        head += (f'<p class="note">Showing the first {len(shown):,} queries; '
                 "neighbours.csv has them all.</p>")
    sections = []
    for q in shown:
        title = f'{q["row"] + 1}. {escape(str(q["query"]), quote=False)}'
        if "error" in q:
            sections.append(f'<section class="query"><h2>{title}</h2>'
                            f'<p class="error">{escape(q["error"])}</p></section>')
        else:
            title += f' <span class="filters">→ {subscript(q["result"]["query"])}</span>'
            sections.append(f'<section class="query"><h2>{title}</h2>'
                            f'{_table(q["result"]["neighbours"], databases)}</section>')
    return _page("ElMTree results", head + _legend() + "".join(sections), databases)
