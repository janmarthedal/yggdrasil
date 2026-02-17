# Yggdrasil

A Python library for finite element analysis, built on NumPy and SciPy.

## Setup

This project uses [uv](https://docs.astral.sh/uv/) for package management.

```bash
uv sync
```

## Running tests

```bash
uv run pytest
```

## Running examples

Examples have additional dependencies (e.g. matplotlib) managed via the `examples` dependency group:

```bash
uv run --group examples python examples/elements/line3.py
```
