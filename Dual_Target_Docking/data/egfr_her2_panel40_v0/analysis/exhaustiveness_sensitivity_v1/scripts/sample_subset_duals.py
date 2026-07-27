import csv, random
from pathlib import Path
SEED=20260727
mandatory=['EH40_01', 'EH40_02', 'EH40_05', 'EH40_18', 'EH40_23']
root=Path(__file__).resolve().parents[3]
panel=list(csv.DictReader((root/"tables"/"panel_v0_40.csv").open()))
dual_pool=[r["panel_id"] for r in panel if r["class"]=="dual" and r["panel_id"] not in mandatory]
rng=random.Random(SEED)
print(sorted(rng.sample(dual_pool,3)))
