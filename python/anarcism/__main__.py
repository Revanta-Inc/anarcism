"""Allow ``python -m anarcism`` to invoke the compatibility CLI."""

from .cli import main

raise SystemExit(main())
