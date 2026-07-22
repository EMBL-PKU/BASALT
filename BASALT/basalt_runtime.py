#!/usr/bin/env python3

"""Small runtime checks shared by BASALT command modules."""

from __future__ import annotations

import os
from pathlib import Path


EXPECTED_MODEL_ENSEMBLES = 5


def require_model_directory() -> Path:
    """Return a validated BASALT model directory or raise an actionable error."""
    configured = os.environ.get("BASALT_WEIGHT")
    if not configured:
        raise RuntimeError(
            "BASALT_WEIGHT is not set. Download the five BASALT model ensembles "
            "with BASALT_models_download.py, then export BASALT_WEIGHT to that "
            "absolute directory."
        )

    model_dir = Path(configured).expanduser().resolve()
    if not model_dir.is_dir():
        raise RuntimeError(
            "BASALT_WEIGHT does not name a directory: {}".format(model_dir)
        )

    descriptors = sorted(model_dir.glob("*_ensemble.csv"))
    if len(descriptors) != EXPECTED_MODEL_ENSEMBLES:
        raise RuntimeError(
            "BASALT_WEIGHT must contain {} top-level *_ensemble.csv files; "
            "found {} in {}.".format(
                EXPECTED_MODEL_ENSEMBLES, len(descriptors), model_dir
            )
        )

    missing_directories = [
        descriptor.with_suffix("").name
        for descriptor in descriptors
        if not descriptor.with_suffix("").is_dir()
    ]
    if missing_directories:
        raise RuntimeError(
            "BASALT_WEIGHT is missing checkpoint directories: {}".format(
                ", ".join(missing_directories)
            )
        )

    missing_checkpoints = []
    unsafe_checkpoints = []
    for descriptor in descriptors:
        checkpoint_dir = descriptor.with_suffix("").resolve()
        entries = [
            line.strip()
            for line in descriptor.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        if not entries:
            missing_checkpoints.append("{} (empty descriptor)".format(descriptor.name))
            continue
        for entry in entries:
            checkpoint = (checkpoint_dir / entry).resolve()
            if checkpoint_dir not in checkpoint.parents:
                unsafe_checkpoints.append("{}/{}".format(checkpoint_dir.name, entry))
            elif not checkpoint.is_file():
                missing_checkpoints.append("{}/{}".format(checkpoint_dir.name, entry))
    if unsafe_checkpoints:
        raise RuntimeError(
            "BASALT_WEIGHT descriptors contain unsafe checkpoint paths: {}".format(
                ", ".join(unsafe_checkpoints[:10])
            )
        )
    if missing_checkpoints:
        raise RuntimeError(
            "BASALT_WEIGHT has missing checkpoint files: {}".format(
                ", ".join(missing_checkpoints[:10])
                + (" ..." if len(missing_checkpoints) > 10 else "")
            )
        )

    return model_dir
