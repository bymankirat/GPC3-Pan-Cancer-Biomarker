"""Small helpers shared by every step: paths, logging, JSON, HTTP with retries."""
import gzip
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import config  # noqa: E402

RAW = ROOT / "data" / "raw"          # downloads (not committed)
DERIVED = ROOT / "data" / "derived"  # small extracts of public data (committed)
RESULTS = ROOT / "results"
for _d in (RAW, DERIVED, RESULTS):
    _d.mkdir(parents=True, exist_ok=True)


def log(*parts):
    print(time.strftime("%H:%M:%S"), *parts, flush=True)


# ---------------------------------------------------------------- files
def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def save_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True, default=_jsonable)


def load_json(path):
    with open(path) as fh:
        return json.load(fh)


def _jsonable(x):
    try:
        import numpy as np
        if isinstance(x, np.integer):
            return int(x)
        if isinstance(x, np.floating):
            return None if np.isnan(x) else float(x)
        if isinstance(x, np.ndarray):
            return x.tolist()
    except ImportError:
        pass
    raise TypeError(f"not JSON serialisable: {type(x)}")


# ---------------------------------------------------------------- HTTP
_session = None


def session():
    global _session
    if _session is None:
        import requests
        _session = requests.Session()
        _session.headers["User-Agent"] = "GPC3-pan-cancer (github.com/bymankirat)"
    return _session


def http(method, url, *, ok404=False, tries=7, **kw):
    import requests
    kw.setdefault("timeout", 180)
    for attempt in range(tries):
        try:
            r = session().request(method, url, **kw)
        except requests.RequestException as e:
            log(f"  network error ({e.__class__.__name__}), retrying {url}")
            time.sleep(2 ** attempt)
            continue
        if r.status_code == 404 and ok404:
            return None
        if r.status_code == 429 or r.status_code >= 500:
            wait = float(r.headers.get("Retry-After", 2 ** attempt))
            log(f"  HTTP {r.status_code}, waiting {wait:.0f}s: {url}")
            time.sleep(wait)
            continue
        if r.status_code >= 400:
            raise RuntimeError(f"HTTP {r.status_code} for {url}: {r.text[:300]}")
        return r
    raise RuntimeError(f"gave up after {tries} tries: {url}")


def download(url, dest):
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    r = http("GET", url, stream=True)
    with open(tmp, "wb") as fh:
        for chunk in r.iter_content(1 << 20):
            fh.write(chunk)
    tmp.rename(dest)
    return dest
