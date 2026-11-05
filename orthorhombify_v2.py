#!/usr/bin/env python3
"""Run without installing: same as the `orthorhombify` command."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
from orthorhombify.cli import orthorhombify_main  # noqa: E402

if __name__ == "__main__":
    sys.exit(orthorhombify_main())
