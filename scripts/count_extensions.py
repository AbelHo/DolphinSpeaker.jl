#!/usr/bin/env python3
"""Count files by extension in a folder (optionally recursive).
Usage examples:
    python3 scripts/count_extensions.py --paths /path/to/dir /other/dir --exts .sud .flac .xml -r
    python3 scripts/count_extensions.py --paths /path/to/dir -r --discover    # auto-discover extensions and count all
"""
import argparse
from pathlib import Path
from collections import Counter
import json


def count_extensions(folder_path, extensions, recursive=True):
    p = Path(folder_path)
    exts_norm = set()
    for e in extensions:
        if not e.startswith('.'):
            e = '.' + e
        exts_norm.add(e.lower())

    counts = Counter()
    if not p.exists():
        raise FileNotFoundError(f"Path does not exist: {folder_path}")

    if recursive:
        iterator = p.rglob('*')
    else:
        iterator = p.iterdir()

    for f in iterator:
        if f.is_file():
            if f.suffix.lower() in exts_norm:
                counts[f.suffix.lower()] += 1

    # Ensure all requested extensions appear in output
    for e in exts_norm:
        counts.setdefault(e, 0)

    return dict(counts)


def discover_extensions(folder_path, recursive=True):
    """Discover all file suffixes in folder and count occurrences.
    Files without a suffix are recorded under the key '<noext>'.
    Returns a dict mapping suffix (including leading dot, lowercased) or '<noext>' to counts.
    """
    p = Path(folder_path)
    if not p.exists():
        raise FileNotFoundError(f"Path does not exist: {folder_path}")

    counts = Counter()
    iterator = p.rglob('*') if recursive else p.iterdir()
    for f in iterator:
        if f.is_file():
            s = f.suffix.lower()
            if s == '':
                s = '<noext>'
            counts[s] += 1

    return dict(counts)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--paths', nargs='+', required=True, help='One or more folder paths to check')
    parser.add_argument('--exts', '--extensions', nargs='*', required=False, help='File extensions to count, e.g. .flac .xml')
    parser.add_argument('--discover', '-a', dest='discover', action='store_true', help='Auto-discover all extensions and count them')
    parser.add_argument('--recursive', '-r', dest='recursive', action='store_true', help='Search recursively')
    parser.add_argument('--no-recursive', dest='recursive', action='store_false', help="Don't search recursively")
    parser.set_defaults(recursive=True)
    args = parser.parse_args()

    results = {}
    if not args.discover and not args.exts:
        parser.error('either --exts or --discover must be provided')

    for path in args.paths:
        try:
            if args.discover:
                results[path] = discover_extensions(path, recursive=args.recursive)
            else:
                # args.exts may be a list or None; we already validated not both missing
                results[path] = count_extensions(path, args.exts, recursive=args.recursive)
        except Exception as e:
            results[path] = {'error': str(e)}

    # Print human readable and JSON
    for path, res in results.items():
        print(f"Path: {path}")
        if 'error' in res:
            print('  ERROR:', res['error'])
        else:
            total = sum(res.values())
            for ext, c in sorted(res.items()):
                print(f"  {ext}: {c}")
            print(f"  total: {total}\n")

    # JSON on last line for programmatic parsing
    print(json.dumps(results))
