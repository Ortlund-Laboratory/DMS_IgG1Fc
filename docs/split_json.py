import json
import os

# Load your big file
with open("data/heatmap.json") as f:
    obj = json.load(f)

data = obj["x"]["data"]
layout = obj["x"]["layout"]

# Make directory if needed
os.makedirs("data/split", exist_ok=True)

# Save each trace separately
for i, trace in enumerate(data):
    with open(f"data/split/trace_{i}.json", "w") as out:
        json.dump(trace, out)

# Save layout separately
with open("data/split/layout.json", "w") as out:
    json.dump(layout, out)

print(f"✅ Split into {len(data)} trace files")
