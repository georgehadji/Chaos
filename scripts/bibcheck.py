"""Verify citations against Crossref and arXiv (literature-review/03).

Usage:
  python scripts/bibcheck.py CANDIDATES.md --out DRAFT.bib   fetch every DOI / arXiv ID in the file, write draft BibTeX
  python scripts/bibcheck.py FILE.bib                        check every entry with a DOI against its Crossref record
Set CROSSREF_MAILTO to a contact address for the Crossref polite pool. Standard library only.
Exit code 1 if any record cannot be fetched or any field disagrees.
"""
import concurrent.futures as cf
import html
import json
import os
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET

DOI_RE = re.compile(r"\b(10\.\d{4,9}/[^\s|`,;]+)")
ARXIV_RE = re.compile(r"arXiv:\s*(\d{4}\.\d{4,5}|[a-z-]+(?:\.[A-Z]{2})?/\d{7})", re.I)  # new and old-style IDs
ATOM = "{http://www.w3.org/2005/Atom}"


def ascii_fold(s):
    return unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode()


def norm(s):
    s = html.unescape(re.sub(r"<[^>]+>|[{}\\]", "", s or ""))
    return re.sub(r"[^a-z0-9]", "", ascii_fold(s).lower())


def parse_bib(text):
    """Return [{'key', 'type', field: value}] for each @type{key, ...} entry."""
    out = []
    for m in re.finditer(r"@(\w+)\{([^,]+),(.*?)\n\}", text, re.S):
        entry = {"type": m.group(1).lower(), "key": m.group(2).strip()}
        for f in re.finditer(r"(\w+)\s*=\s*\{(.*?)\},?\s*$", m.group(3), re.M | re.S):
            entry[f.group(1).lower()] = re.sub(r"\s+", " ", f.group(2)).strip()
        out.append(entry)
    return out


def compare(entry, record):
    """Mismatches between a parsed BibTeX entry and a Crossref 'message' record."""
    issues = []
    title = (record.get("title") or [""])[0]
    if norm(title) != norm(entry.get("title", "")):
        issues.append(f"title: record {title!r}")
    fam = [a["family"] for a in record.get("author", []) if a.get("family")]  # skip affiliation-only rows
    mine =[a.strip() for a in entry.get("author", "").split(" and ") if a.strip()]
    if fam and mine and norm(fam[0]) != norm(mine[0].split(",")[0]):
        issues.append(f"first author: record {fam[0]!r}")
    if fam and mine and "others" not in entry.get("author", "") and len(fam) != len(mine):
        issues.append(f"author count: record {len(fam)}, entry {len(mine)}")
    vol = str(record.get("volume", ""))
    if vol and entry.get("volume") and vol != entry["volume"].split("(")[0]:
        issues.append(f"volume: record {vol}")
    pages = (record.get("page") or record.get("article-number") or "").replace("–", "-")
    if pages and entry.get("pages") and pages != entry["pages"].replace("--", "-"):
        issues.append(f"pages: record {pages}")
    return issues


def clean_name(n):
    """Record names as usable text: None or 'None' dropped, ALL-CAPS records title-cased."""
    if not n or n == "None":
        return ""
    return n.title() if n.isupper() else n


def key_base(families, year):
    """Key rule from _config/bibliography/conventions.md: up to three family names, else the first; then year."""
    families = [f for f in map(clean_name, families) if f] or ["Anon"]
    names = families if len(families) <= 3 else families[:1]
    return "".join(re.sub(r"[^A-Za-z]", "", ascii_fold(n)) for n in names) + str(year)


def assign_keys(bases):
    """A base used once stays bare; a clash suffixes every member a, b, ... (Maffione2011a, Maffione2011b)."""
    counts, seen, keys = {}, {}, []
    for b in bases:
        counts[b] = counts.get(b, 0) + 1
    for b in bases:
        if counts[b] == 1:
            keys.append(b)
        else:
            keys.append(b + "abcdefghijklmnopqrstuvwxyz"[seen.get(b, 0)])
            seen[b] = seen.get(b, 0) + 1
    return keys


def http_get(url):
    mail = os.environ.get("CROSSREF_MAILTO", "")
    agent = "chaos-bibcheck/1.0" + (f" (mailto:{mail})" if mail else "")
    req = urllib.request.Request(url, headers={"User-Agent": agent})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=40) as r:
                return r.read()
        except urllib.error.HTTPError as ex:
            if ex.code not in (429, 503) or attempt == 4:
                raise
            time.sleep(float(ex.headers.get("Retry-After") or 2 ** attempt))


def fetch_crossref(doi):
    return json.loads(http_get("https://api.crossref.org/works/" + doi))["message"]


def fetch_arxiv(aid):
    entry = ET.fromstring(http_get("https://export.arxiv.org/api/query?id_list=" + aid)).find(ATOM + "entry")
    names = [a.findtext(ATOM + "name", "") for a in entry.findall(ATOM + "author")]
    return {"title": [" ".join(entry.findtext(ATOM + "title", "").split())],
            "author": [{"family": n.split()[-1], "given": " ".join(n.split()[:-1])} for n in names],
            "issued": {"date-parts": [[int(entry.findtext(ATOM + "published", "0")[:4])]]},
            "eprint": aid}


def rec_year(rec):
    """First non-empty year among the record's date fields (Crossref leaves `issued` empty on some proceedings)."""
    for f in ("issued", "published-print", "published-online", "published"):
        parts = (rec.get(f) or {}).get("date-parts") or [[None]]
        if parts[0] and parts[0][0]:
            return parts[0][0]
    names = [(rec.get("event") or {}).get("name", "")] + (rec.get("container-title") or [])
    for s in names:  # proceedings: the year in the event or volume title, e.g. '1997 1st International Conference'
        m = re.search(r"\b(19|20)\d{2}\b", s or "")
        if m:
            return int(m.group(0))
    parts = (rec.get("created") or {}).get("date-parts") or [[0]]  # metadata date: last resort
    return parts[0][0] or 0


def to_bibtex(key, rec):
    year = rec_year(rec)
    authors = " and ".join(f"{clean_name(a.get('family'))}, {clean_name(a.get('given'))}".strip(", ")
                           for a in rec.get("author", []) if clean_name(a.get("family")))
    fields = {"author": authors, "title": (rec.get("title") or [""])[0], "year": year,
              "journal": (rec.get("container-title") or [""])[0], "volume": rec.get("volume", ""),
              "number": rec.get("issue", ""), "pages": (rec.get("page") or "").replace("–", "--"),
              "doi": rec.get("DOI", ""), "eprint": rec.get("eprint", "")}
    body = ",\n".join(f"  {k} = {{{v}}}" for k, v in fields.items() if v)
    kind = "article" if fields["journal"] else "misc"
    return f"@{kind}{{{key},\n{body}\n}}\n"


def check_bib(path):
    entries = [e for e in parse_bib(open(path, encoding="utf-8").read()) if e.get("doi")]
    bad = 0
    with cf.ThreadPoolExecutor(3) as ex:
        for e, rec in zip(entries, ex.map(lambda e: _safe(fetch_crossref, e["doi"]), entries)):
            issues = [f"fetch failed: {rec}"] if isinstance(rec, str) else compare(e, rec)
            if issues:
                bad += 1
                print(f"{e['key']}: " + "; ".join(issues))
    print(f"checked {len(entries)} entries with a DOI, {bad} with issues")
    return 1 if bad else 0


def fetch_candidates(path, out):
    text = open(path, encoding="utf-8").read()
    ids = [("doi", d.rstrip(".")) for d in dict.fromkeys(DOI_RE.findall(text))]
    ids += [("arxiv", a) for a in dict.fromkeys(ARXIV_RE.findall(text))]
    recs, titles, failed = [], {}, 0
    for kind, ident in ids:
        rec = _safe(fetch_crossref if kind == "doi" else fetch_arxiv, ident)
        if isinstance(rec, str):
            failed += 1
            print(f"{kind}:{ident}: fetch failed: {rec}")
            continue
        year = rec_year(rec)
        t = (norm((rec.get("title") or [""])[0]), year)
        if t in titles:  # same paper under a second identifier (e.g. old and new DOI formats)
            print(f"{kind}:{ident}: duplicate record of {titles[t]}, skipped")
            continue
        titles[t] = f"{kind}:{ident}"
        recs.append(rec)
    bases = [key_base([a.get("family") for a in r.get("author", [])], rec_year(r)) for r in recs]
    chunks = [to_bibtex(k, r) for k, r in zip(assign_keys(bases), recs)]
    open(out, "w", encoding="utf-8", newline="\n").write("\n".join(chunks))
    print(f"wrote {len(chunks)} entries to {out}; {failed} failed")
    return 1 if failed else 0


def _safe(fn, arg):
    try:
        return fn(arg)
    except Exception as ex:  # network or parse failure is reported per record, never silently dropped
        return f"{type(ex).__name__}: {ex}"


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    if argv[1].endswith(".bib"):
        return check_bib(argv[1])
    out = argv[argv.index("--out") + 1] if "--out" in argv else "draft.bib"
    return fetch_candidates(argv[1], out)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
