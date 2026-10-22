from pathlib import Path
import csv, hashlib
root=Path('/tmp/psmelding-fixed-backend-20261022')
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
with (root/'reference_integrity_before.csv').open() as f:
    previous=list(csv.DictReader(f))
verified=[]
for row in previous:
    p=Path(row['reference_root'])/row['path']
    assert p.is_file() and digest(p)==row['sha256'] and p.stat().st_size==int(row['bytes']), p
    verified.append(dict(**row,unchanged=True))
for ref in set(row['reference_root'] for row in previous):
    assert {str(p.relative_to(ref)) for p in Path(ref).rglob('*') if p.is_file()} == {x['path'] for x in previous if x['reference_root']==ref}
with (root/'reports/reference_integrity.csv').open('w') as f:
    w=csv.DictWriter(f,verified[0]);w.writeheader();w.writerows(verified)
(root/'reports/reference_integrity.txt').write_text(f'All {len(verified)} files across both prior diagnostics unchanged (byte sizes and SHA-256); no added or removed files.\n')
rows=[]
for p in sorted(root.rglob('*')):
    if p.is_file() and p.name!='sha256_manifest.csv':
        rows.append(dict(path=str(p.relative_to(root)),bytes=p.stat().st_size,sha256=digest(p)))
with (root/'sha256_manifest.csv').open('w') as f:
    w=csv.DictWriter(f,rows[0]);w.writeheader();w.writerows(rows)
print(f'Verified all {len(verified)} reference files unchanged and generated {len(rows)}-file output SHA-256 manifest.')
