# CoreMet deployment protocol (GitHub + Zenodo + Render)

CoreMet has three parts in three places:

| Part | Where | Notes |
|------|-------|-------|
| **Code** (`CoreMet/`) | GitHub | small (~4 MB); what Render deploys |
| **Citable dataset** (`coremet_dataset_v1.zip`) | Zenodo / Figshare | DOI for the paper |
| **Runtime data** (`coremet_runtime_core.zip`, 29 MB) | Zenodo (fetched at build) | the data the app loads at startup |

The full `data/` tree (~15 GB raw sources) is never committed; `.gitignore` excludes it and
`scripts/fetch_data.py` pulls the 29 MB runtime bundle at deploy time.

---

## Prerequisites

- A GitHub account and the `CoreMet/` folder (this repo).
- A Render account (https://render.com), free to create; the Standard instance is paid.
- The `CoreMetDB/` folder with the two zips (uploaded to Zenodo in Part 2).

---

## Part 1: Push the code to GitHub

```bash
cd CoreMet
git init
git add -A                     # .gitignore keeps data/ out (only data/coremetdb_stats.json is tracked)
git status                     # confirm no data/*.csv, data/models, *.pth are staged
git commit -m "CoreMet web resource"
git branch -M main
git remote add origin https://github.com/<you>/coremet.git
git push -u origin main
```

Sanity check (largest tracked file should be a few hundred KB):
```bash
git ls-files | xargs du -h 2>/dev/null | sort -rh | head
```

---

## Part 2: Upload the data to Zenodo

1. https://zenodo.org → **New upload**.
2. Drag in `CoreMetDB/coremet_dataset_v1.zip` and `CoreMetDB/coremet_runtime_core.zip`.
3. Metadata: Title, Authors, Type = **Dataset**, License = **CC BY 4.0**, short description.
4. **Publish**, then copy:
   - the record **DOI** (for the manuscript), and
   - the **direct download URL** of `coremet_runtime_core.zip`:
     `https://zenodo.org/records/<RECORD_ID>/files/coremet_runtime_core.zip?download=1`

You will paste that runtime URL into Render as `DATA_BUNDLE_URL` (Part 3).

---

## Part 3: Deploy on Render

### Option A, Blueprint (recommended, uses render.yaml)

1. Render dashboard → **New** → **Blueprint**.
2. Connect your GitHub repo. Render reads `render.yaml` and creates the service.
3. When prompted, set the one un-synced variable:
   - **`DATA_BUNDLE_URL`** = the Zenodo runtime URL from Part 2.
4. **Apply** → first build runs `pip install` then `fetch_data.py` (downloads + extracts the 29 MB bundle).

### Option B, Manual web service (field by field)

Render dashboard → **New** → **Web Service** → connect the repo, then set:

| Field | Value |
|-------|-------|
| **Name** | `coremet` |
| **Region** | nearest your users (e.g. Oregon, Frankfurt, Singapore) |
| **Branch** | `main` |
| **Root Directory** | *(blank; repo root is the app)* |
| **Runtime / Language** | `Python 3` |
| **Build Command** | `pip install -r requirements.txt && python scripts/fetch_data.py` |
| **Start Command** | `gunicorn run:server -b 0.0.0.0:$PORT --workers 1 --threads 4 --timeout 300 --max-requests 1000 --max-requests-jitter 50` |
| **Instance Type** | **Standard** (2 GB RAM) for the database app; 4 GB+ if you enable prediction models |
| **Health Check Path** | `/api/v1/health` |
| **Auto-Deploy** | On (deploys on every push to `main`) |

**Environment variables** (Advanced → Add Environment Variable):

| Key | Value | Why |
|-----|-------|-----|
| `DATA_BUNDLE_URL` | Zenodo direct URL to `coremet_runtime_core.zip` | fetched at build |
| `PYTHON_VERSION` | `3.10.13` | pin the interpreter |
| `FLASK_ENV` | `production` | disable debug |
| `CORMET_LAZY_MODELS` | `true` | load optional models only on first predict |
| `PORT` | *(Render sets this automatically; do not override)* | gunicorn binds `$PORT` |

Click **Create Web Service**.

### Build/start parameters explained

- `--workers 1 --threads 4`: one worker keeps the in-memory databases loaded once (~1.5 GB);
  four threads serve concurrent requests. Adding workers multiplies memory and is not needed.
- `--timeout 300`: first request after a cold start loads several CSVs; a long timeout avoids 502s.
- `--max-requests 1000 --max-requests-jitter 50`: periodically recycle the worker to cap memory creep.
- `fetch_data.py` is idempotent: it skips downloading if the data is already present
  (e.g., on a persistent disk), so redeploys are fast.

### Memory and plan sizing

| Configuration | Resident memory | Render plan |
|---------------|----------------:|-------------|
| Database only (default) | ~1.5 GB | Standard (2 GB) |
| Database + prediction models | ~2.7 GB | 4 GB+ (Pro) |

### Optional, persistent disk instead of build-time fetch

To avoid re-downloading on each deploy, add a disk (Render → service → Disks):
mount path `/opt/render/project/src/data`, size `2 GB`. Upload the data once (Render Shell or
a one-off job); `fetch_data.py` then detects it and skips the download.

---

## Part 4: HTTPS and a stable URL

- Render serves every service over HTTPS automatically at `https://coremet.onrender.com`
  (satisfies the NAR HTTPS requirement).
- For a citable, stable address add a **custom domain** (service → Settings → Custom Domains),
  e.g. `coremetdb.org`, and commit to keeping it for at least five years (NAR requirement).
- Replace `[URL]` everywhere once live: manuscript abstract/body, `Cover_Letter.md`,
  `presubmission_query_coremet_db.txt`, and the app `help` page.

---

## Part 5: First-deploy verification

After the deploy goes green, check:
```bash
curl -s https://<your-url>/api/v1/health           # -> {"status": "ok", ...}
```
In a browser: `/home` shows 1,952,688 interactions; `/database` and the seven modules
(`/mpi /mei /mdi /mmi /mdri /mgi /mgwas`) load; `/metabolite?id=HMDB0000039` renders;
`/network` populates on a query; `/downloads` links work.

### Troubleshooting

| Symptom | Cause / fix |
|---------|-------------|
| Build fails at `fetch_data.py` | `DATA_BUNDLE_URL` unset or wrong. Verify the Zenodo direct link ends with `?download=1`. |
| App boots but counts are 0 / pages empty | Runtime bundle not extracted. Check build logs for "fetch_data: done"; confirm `data/databases/release/coremetdb_mgi.csv` exists. |
| 502 on first request | Cold-start data load exceeded timeout. Keep `--timeout 300`; retry once warm. |
| Out-of-memory / worker killed | On Standard with models enabled. Set `CORMET_LAZY_MODELS=true` (default) or move to a 4 GB+ instance. |

---

## Part 6: Updating data later

Re-run the data pipeline locally, rebuild the bundles, publish a new Zenodo version, then
update `DATA_BUNDLE_URL` and redeploy:
```bash
python scripts/compute_db_stats.py
python scripts/build_release_csvs.py
python scripts/build_runtime_bundle.py      # -> ../CoreMetDB/coremet_runtime_core.zip
python scripts/build_data_release_bundle.py # -> ../data_deposit/ (citable)
```
