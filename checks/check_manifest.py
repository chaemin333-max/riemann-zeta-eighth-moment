"""Verify the distributed files against SHA256SUMS using the standard library."""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parents[1]
rows = []
for line in (root / 'SHA256SUMS').read_text(encoding='utf-8').splitlines():
    expected, name = line.split(maxsplit=1)
    path = root / name.lstrip('*')
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    rows.append({'file': name, 'matches': actual == expected})
result = {'status': 'PASS' if all(row['matches'] for row in rows) else 'MISMATCH',
          'files': len(rows), 'matching': sum(row['matches'] for row in rows),
          'mismatches': [row['file'] for row in rows if not row['matches']]}
print(json.dumps(result, indent=2))
raise SystemExit(0 if result['status'] == 'PASS' else 1)
