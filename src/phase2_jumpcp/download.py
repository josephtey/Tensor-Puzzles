"""Download a sample of JUMP-CP TARGET2 control plates across sources.

Selects the first 2 TARGET2 plates from each source that has them (11 sources)
=> 22 plates, ~290 MB total.  Parquets land in data/jump_cp/profiles/.
"""
from __future__ import annotations

import shutil
import urllib.request
from pathlib import Path

import pandas as pd
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[2]
META_PATH = ROOT / "data" / "jump_cp" / "plate.csv"
OUT_DIR = ROOT / "data" / "jump_cp" / "profiles"
OUT_DIR.mkdir(parents=True, exist_ok=True)

S3_URL = "https://cellpainting-gallery.s3.amazonaws.com"
PLATES_PER_SOURCE = 2


def main() -> None:
    plate = pd.read_csv(META_PATH)
    target2 = plate[plate["Metadata_PlateType"] == "TARGET2"].copy()
    selected = target2.groupby("Metadata_Source").head(PLATES_PER_SOURCE).reset_index(drop=True)
    print(f"Selected {len(selected)} plates from {selected['Metadata_Source'].nunique()} sources")

    selected.to_csv(OUT_DIR.parent / "selected_plates.csv", index=False)

    for _, row in tqdm(list(selected.iterrows()), desc="downloading"):
        src = row["Metadata_Source"]
        batch = row["Metadata_Batch"]
        plate_id = row["Metadata_Plate"]
        url = f"{S3_URL}/cpg0016-jump/{src}/workspace/profiles/{batch}/{plate_id}/{plate_id}.parquet"
        out = OUT_DIR / f"{src}__{plate_id}.parquet"
        if out.exists() and out.stat().st_size > 1_000_000:
            continue
        try:
            with urllib.request.urlopen(url, timeout=60) as resp, open(out, "wb") as f:
                shutil.copyfileobj(resp, f)
        except Exception as e:
            print(f"  failed {src}/{plate_id}: {e}")
            if out.exists():
                out.unlink()

    files = sorted(OUT_DIR.glob("*.parquet"))
    total_mb = sum(f.stat().st_size for f in files) / 1e6
    print(f"\nTotal: {len(files)} parquet files ({total_mb:.1f} MB)")


if __name__ == "__main__":
    main()
