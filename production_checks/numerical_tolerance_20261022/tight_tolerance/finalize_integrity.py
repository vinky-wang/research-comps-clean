from pathlib import Path
import csv, hashlib
root=Path('/tmp/psmelding-tight-tolerance-20261022')
prior=Path('/Users/vinky/research-comps-clean/production_checks/numerical_stability')
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
with (root/'previous_diagnostic_integrity_before.csv').open() as f:
    previous=list(csv.DictReader(f))
verified=[]
for row in previous:
    p=prior/row['path']
    assert p.is_file() and digest(p)==row['sha256'] and p.stat().st_size==int(row['bytes']), p
    verified.append(dict(**row,unchanged=True))
assert {str(p.relative_to(prior)) for p in prior.rglob('*') if p.is_file()} == {x['path'] for x in previous}
with (root/'reports/previous_diagnostic_integrity.csv').open('w') as f:
    w=csv.DictWriter(f,verified[0]);w.writeheader();w.writerows(verified)
(root/'reports/previous_diagnostic_integrity.txt').write_text(f'All {len(verified)} prior diagnostic files unchanged (byte sizes and SHA-256); no added or removed files.\n')
rows=[]
for p in sorted(root.rglob('*')):
    if p.is_file() and p.name!='sha256_manifest.csv':
        rows.append(dict(path=str(p.relative_to(root)),bytes=p.stat().st_size,sha256=digest(p)))
with (root/'sha256_manifest.csv').open('w') as f:
    w=csv.DictWriter(f,rows[0]);w.writeheader();w.writerows(rows)
print(f'Verified all {len(verified)} previous files unchanged and generated {len(rows)}-file output SHA-256 manifest.')
