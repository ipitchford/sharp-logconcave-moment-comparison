#!/usr/bin/env python3
"""Verify hashes AND complete payload inventory (including hidden files).

Only .git metadata and the manifest itself are excluded. Run before generating
outputs; put fresh replay artifacts outside the checked checkout.
"""
import argparse
import hashlib
from pathlib import Path, PurePosixPath
import re


def check(root, name='MANIFEST.sha256'):
    root = Path(root).resolve()
    manifest = root/name
    if manifest.is_symlink() or not manifest.is_file():
        raise ValueError('missing/nonregular manifest')
    expected = {}
    for line in manifest.read_text().splitlines():
        match = re.fullmatch(r'([0-9a-f]{64})  (.+)', line)
        if not match:
            raise ValueError('invalid manifest line')
        digest, text = match.groups(); path = PurePosixPath(text)
        if path.is_absolute() or '..' in path.parts or '.' in path.parts or str(path) != text or text == name or '.git' in path.parts or text in expected:
            raise ValueError('unsafe/duplicate/self path: '+text)
        expected[text] = digest
    if not expected:
        raise ValueError('empty manifest')
    actual = {}
    for path in root.rglob('*'):
        rel = path.relative_to(root)
        if '.git' in rel.parts:
            continue
        if path.is_symlink():
            raise ValueError('symlink in payload: '+str(rel))
        if path.is_file() and str(rel) != name:
            actual[rel.as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError('manifest mismatch: missing='+str(sorted(expected.keys()-actual.keys()))+' unlisted='+str(sorted(actual.keys()-expected.keys()))+' changed='+str(sorted(k for k in expected.keys() & actual.keys() if expected[k] != actual[k])))
    return len(actual)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root', default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument('--manifest', default='MANIFEST.sha256')
    args = ap.parse_args()
    if Path(args.manifest).name != args.manifest:
        ap.error('manifest must be a root filename')
    print('PASS:', check(args.root,args.manifest), 'files, complete inventory and SHA-256')
