"""Metabolite display-name cleanup and HMDB cross-table lookup."""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import pandas as pd


MISSING_NAME_TOKENS = {
    "",
    "-",
    "nan",
    "none",
    "null",
    "<na>",
    "unamed",
    "unnamed",
    "unknown",
    "metabolite name unavailable",
    "unspecified metabolite",
}


def _release_dir() -> Path:
    return Path(__file__).resolve().parent.parent.parent / "data" / "databases" / "release"


def _clean_series(series: pd.Series) -> pd.Series:
    cleaned = series.astype("string").fillna("").str.strip()
    return cleaned.mask(cleaned.str.lower().isin(MISSING_NAME_TOKENS), "")


def _is_informative_name(name: str, hmdb_id: str = "") -> bool:
    text = str(name or "").strip()
    hid = str(hmdb_id or "").strip()
    if not text or text.lower() in MISSING_NAME_TOKENS:
        return False
    if hid and text.casefold() == hid.casefold():
        return False
    if text.upper().startswith("HMDB") and text[4:].isdigit():
        return False
    return True


def _first_names_from_frame(df: pd.DataFrame, hmdb_col: str, name_col: str) -> dict[str, str]:
    if hmdb_col not in df.columns or name_col not in df.columns:
        return {}
    ids = _clean_series(df[hmdb_col])
    names = _clean_series(df[name_col])
    valid = ids.str.startswith("HMDB", na=False)
    valid &= pd.Series(
        [bool(_is_informative_name(name, hid)) for hid, name in zip(ids, names)],
        index=df.index,
    )
    if not valid.any():
        return {}
    mapped = pd.DataFrame({"hmdb_id": ids[valid], "name": names[valid]})
    return mapped.drop_duplicates("hmdb_id").set_index("hmdb_id")["name"].to_dict()


@lru_cache(maxsize=128)
def lookup_metabolite_names(hmdb_ids_key: tuple[str, ...]) -> dict[str, str]:
    """Return best available names for HMDB IDs using all local release tables."""
    pending = {str(h).strip() for h in hmdb_ids_key if str(h).strip().startswith("HMDB")}
    if not pending:
        return {}

    found: dict[str, str] = {}
    release_dir = _release_dir()
    files = [
        release_dir / "coremetdb_mdi.csv",
        release_dir / "coremetdb_mdri.csv",
        release_dir / "coremetdb_mei.csv",
        release_dir / "coremetdb_mgwas.csv",
        release_dir / "coremetdb_mmi.csv",
        release_dir / "coremetdb_mpi.csv",
        release_dir / "coremetdb_mgi.csv",
    ]

    for path in files:
        if not path.exists() or not pending:
            continue

        hmdb_col = "HMDB ID" if path.name.endswith("_mpi.csv") else "HMDB_ID"
        name_col = "Metabolite Name" if path.name.endswith("_mpi.csv") else "Metabolite_Name"
        try:
            reader = pd.read_csv(
                path,
                usecols=[hmdb_col, name_col],
                dtype=str,
                keep_default_na=False,
                na_filter=False,
                chunksize=200_000,
            )
            for chunk in reader:
                ids = _clean_series(chunk[hmdb_col])
                sub = chunk[ids.isin(pending)]
                if sub.empty:
                    continue
                for hid, name in zip(_clean_series(sub[hmdb_col]), _clean_series(sub[name_col])):
                    if hid in pending and _is_informative_name(name, hid):
                        found[str(hid)] = str(name)
                pending.difference_update(found)
                if not pending:
                    break
        except Exception:
            continue

    return found


def refine_metabolite_names(
    df: pd.DataFrame,
    name_col: str,
    hmdb_col: str,
    fallback_to_hmdb: bool = True,
) -> pd.DataFrame:
    """Fill bad metabolite names from local mappings; use HMDB ID as last resort."""
    if df.empty or name_col not in df.columns or hmdb_col not in df.columns:
        return df

    out = df.copy()
    out[name_col] = _clean_series(out[name_col])
    out[hmdb_col] = _clean_series(out[hmdb_col])

    ids = _clean_series(out[hmdb_col])
    names = _clean_series(out[name_col])
    needs = ids.str.startswith("HMDB", na=False)
    needs &= pd.Series(
        [not _is_informative_name(name, hid) for hid, name in zip(ids, names)],
        index=out.index,
    )

    if not needs.any():
        blank = names.eq("")
        out.loc[blank, name_col] = "-"
        return out

    local_map = _first_names_from_frame(out, hmdb_col, name_col)
    unresolved = set(ids[needs])

    if local_map:
        fill = ids[needs].map(local_map).fillna("")
        has_fill = fill.ne("")
        target_index = fill[has_fill].index
        out.loc[target_index, name_col] = fill[has_fill]
        unresolved.difference_update(ids.loc[target_index].tolist())

    if unresolved:
        external_map = lookup_metabolite_names(tuple(sorted(unresolved)))
        if external_map:
            current_names = _clean_series(out[name_col])
            current_needs = ids.str.startswith("HMDB", na=False)
            current_needs &= pd.Series(
                [not _is_informative_name(name, hid) for hid, name in zip(ids, current_names)],
                index=out.index,
            )
            fill = ids[current_needs].map(external_map).fillna("")
            has_fill = fill.ne("")
            out.loc[fill[has_fill].index, name_col] = fill[has_fill]

    current_names = _clean_series(out[name_col])
    current_needs = ids.str.startswith("HMDB", na=False)
    current_needs &= pd.Series(
        [not _is_informative_name(name, hid) for hid, name in zip(ids, current_names)],
        index=out.index,
    )
    if fallback_to_hmdb:
        out.loc[current_needs, name_col] = ids[current_needs]

    out.loc[_clean_series(out[name_col]).eq(""), name_col] = "-"
    return out


def refine_metabolite_names_for_known_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Apply metabolite-name cleanup for either MPI or underscore schemas."""
    if "Metabolite Name" in df.columns and "HMDB ID" in df.columns:
        return refine_metabolite_names(df, "Metabolite Name", "HMDB ID")
    if "Metabolite_Name" in df.columns and "HMDB_ID" in df.columns:
        return refine_metabolite_names(df, "Metabolite_Name", "HMDB_ID")
    return df
