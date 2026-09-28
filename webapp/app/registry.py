from __future__ import annotations
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def load_json(name:str):
    with open(ROOT/"benchmark_registry"/name) as fh: return json.load(fh)
