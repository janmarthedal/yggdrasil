# Yggdrasil

A Python library for finite element analysis, built on NumPy and SciPy.

## Setup

This project uses [uv](https://docs.astral.sh/uv/) for package management.

```bash
uv sync
```

## Running tests

```bash
uv run pytest                # all tests
uv run pytest tests/unit/    # unit tests only
uv run pytest tests/system/  # system tests only
```

## Linting

```bash
uv run ruff check
```

## Running examples

Examples have additional dependencies (e.g. matplotlib) managed via the `examples` dependency group:

```bash
uv run --group examples python examples/elements/line3.py
```
