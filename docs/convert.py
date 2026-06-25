import re
import json
import os

# ✅ Read HTML
with open("index.html", "r", encoding="utf-8") as f:
    text = f.read()

# ✅ Find ALL JSON script blocks
pattern = r'<script[^>]*type="application/json"[^>]*>(.*?)</script>'
matches = re.findall(pattern, text, re.S)

print(f"✅ Found {len(matches)} JSON blocks")

if not matches:
    raise ValueError("❌ No JSON blocks found")

# ✅ Parse ALL blocks and find the one with real trace data
selected_json = None

for i, m in enumerate(matches):
    try:
        parsed = json.loads(m)

        # ✅ Case 1: list of traces
        if isinstance(parsed, list):
            print(f"✅ Block {i} is a list with {len(parsed)} entries")
            if len(parsed) > 100:  # heuristic: real dataset is large
                selected_json = parsed
                print(f"🎯 Selected block {i}")
                break

        # ✅ Case 2: object with "data"
        elif isinstance(parsed, dict) and "data" in parsed:
            print(f"✅ Block {i} has 'data' with {len(parsed['data'])} traces")
            if len(parsed["data"]) > 100:
                selected_json = parsed["data"]
                print(f"🎯 Selected block {i}")
                break

        # ✅ Case 3: htmlwidgets nested structure
        elif isinstance(parsed, dict) and "x" in parsed and "data" in parsed["x"]:
            print(f"✅ Block {i} has nested 'x.data' with {len(parsed['x']['data'])} traces")
            if len(parsed["x"]["data"]) > 100:
                selected_json = parsed["x"]["data"]
                print(f"🎯 Selected block {i}")
                break

    except Exception:
        continue

# ✅ Fail clearly if nothing found
if selected_json is None:
    raise ValueError("❌ Could not find correct dataset block")

# ✅ Save correct dataset
os.makedirs("data", exist_ok=True)

with open("data/heatmap.json", "w", encoding="utf-8") as out:
    json.dump(selected_json, out)

print("✅ Saved FULL dataset to data/heatmap.json")
