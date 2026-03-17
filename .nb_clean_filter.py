#!/usr/bin/env python3
"""Git clean filter: strips volatile cell IDs from Jupyter notebooks.

Reads a notebook from stdin, replaces cell 'id' fields with a stable
sequential value, and writes the result to stdout.  Uses only the
standard library so it works regardless of virtualenv state.
"""
import json
import sys


def clean_notebook(nb):
    for i, cell in enumerate(nb.get("cells", [])):
        cell["id"] = str(i)
    return nb


if __name__ == "__main__":
    try:
        nb = json.load(sys.stdin)
        sys.stdout.write(json.dumps(clean_notebook(nb), indent=1, ensure_ascii=False))
        sys.stdout.write("\n")
    except Exception:
        # If parsing fails, pass the file through unchanged
        sys.stdin.seek(0)
        sys.stdout.write(sys.stdin.read())
