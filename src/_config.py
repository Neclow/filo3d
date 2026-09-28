"""Shared configuration helpers for Filo3D."""

from pathlib import Path

import yaml


def load_params(path="params.yaml"):
    """Load the shared params.yaml into a dict."""
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def load_json(path):
    import json

    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)
    return Path(path)
