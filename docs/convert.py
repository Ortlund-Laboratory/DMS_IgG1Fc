import re

with open("index.html") as f:
    text = f.read()

pattern = r'<script[^>]*type="application/json"[^>]*>(.*?)</script>'

match = re.search(pattern, text, re.S)

if match:
    # 1. Save JSON to file
    with open("data/heatmap.json", "w") as out:
        out.write(match.group(1))

    # 2. Remove that script block from HTML
    new_html = re.sub(pattern, '', text, count=1, flags=re.S)

    # 3. Write updated HTML back
    with open("index.html", "w") as f:
        f.write(new_html)

    print("✅ Cut JSON from index.html → data/heatmap.json")
else:
    print("❌ No matching script tag found")
