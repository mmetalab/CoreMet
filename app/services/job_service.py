"""
Job persistence service — SQLite-based with 7-day retention.
"""

import json
import logging
import sqlite3
import uuid
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional, Dict

import pandas as pd

logger = logging.getLogger(__name__)

DB_PATH = Path(__file__).parent.parent.parent / "data" / "jobs.db"


def _get_connection():
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    try:
        conn = sqlite3.connect(str(DB_PATH), check_same_thread=False)
        conn.execute("PRAGMA journal_mode=WAL")
    except sqlite3.Error:
        logger.exception("Failed to open jobs database at %s", DB_PATH)
        raise
    conn.execute("""
        CREATE TABLE IF NOT EXISTS jobs (
            job_id TEXT PRIMARY KEY,
            created_at TEXT,
            status TEXT,
            input_metabolites TEXT,
            input_proteins TEXT,
            organism TEXT,
            results TEXT,
            expires_at TEXT
        )
    """)
    conn.commit()
    return conn


def create_job(metabolites_json: str, proteins_json: str, organism: str) -> str:
    """Create a new prediction job and return its ID."""
    job_id = str(uuid.uuid4())[:12]
    now = datetime.utcnow()
    expires = now + timedelta(days=7)

    conn = _get_connection()
    conn.execute(
        "INSERT INTO jobs (job_id, created_at, status, input_metabolites, input_proteins, organism, results, expires_at) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        (job_id, now.isoformat(), "running", metabolites_json, proteins_json, organism, "", expires.isoformat()),
    )
    conn.commit()
    conn.close()
    return job_id


def update_job(job_id: str, status: str, results_json: str):
    """Update job status and results."""
    conn = _get_connection()
    conn.execute(
        "UPDATE jobs SET status = ?, results = ? WHERE job_id = ?",
        (status, results_json, job_id),
    )
    conn.commit()
    conn.close()


def get_job(job_id: str) -> Optional[Dict]:
    """Get job by ID."""
    conn = _get_connection()
    cur = conn.execute("SELECT * FROM jobs WHERE job_id = ?", (job_id,))
    row = cur.fetchone()
    conn.close()

    if row is None:
        return None

    cols = ['job_id', 'created_at', 'status', 'input_metabolites', 'input_proteins', 'organism', 'results', 'expires_at']
    return dict(zip(cols, row))


def cleanup_expired():
    """Delete expired jobs."""
    conn = _get_connection()
    now = datetime.utcnow().isoformat()
    conn.execute("DELETE FROM jobs WHERE expires_at < ?", (now,))
    conn.commit()
    conn.close()
