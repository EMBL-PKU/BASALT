#!/usr/bin/env python3

"""Download and validate the trained BASALT model weights."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Optional
from zipfile import BadZipFile, ZipFile


HUGGING_FACE_REPO = "PKU-EMBL/BASALT_WEIGHT"
FIGSHARE_URL = "https://figshare.com/ndownloader/files/41093033"
BAIDU_URL = "https://pan.baidu.com/s/1ouKqabxHYr1GmvpquQCzqw?pwd=embl"
EXPECTED_MODELS = 5


def resolve_target_dir(path: Optional[str] = None) -> Path:
    """Resolve the model directory from CLI, environment, or local default."""
    configured = path or os.environ.get("BASALT_WEIGHT")
    return Path(configured).expanduser() if configured else Path.cwd() / "BASALT_WEIGHT"


def validate_models(target_dir: Path) -> list[Path]:
    """Validate the five top-level ensemble descriptors and checkpoint folders."""
    descriptors = sorted(target_dir.glob("*_ensemble.csv"))
    if len(descriptors) != EXPECTED_MODELS:
        raise RuntimeError(
            "Expected {} top-level *_ensemble.csv descriptors in {}, found {}.".format(
                EXPECTED_MODELS, target_dir, len(descriptors)
            )
        )
    missing_directories = [
        descriptor.with_suffix("").name
        for descriptor in descriptors
        if not descriptor.with_suffix("").is_dir()
    ]
    if missing_directories:
        raise RuntimeError(
            "Missing model checkpoint directories: {}".format(
                ", ".join(missing_directories)
            )
        )
    return descriptors


def models_are_ready(target_dir: Path) -> bool:
    """Return whether a complete model set is already present."""
    try:
        validate_models(target_dir)
    except RuntimeError:
        return False
    return True


def download_huggingface(
    target_dir: Path,
    endpoint: Optional[str] = None,
    revision: str = "main",
) -> None:
    """Download the model repository concurrently with resumable HF metadata."""
    try:
        from huggingface_hub import snapshot_download
    except ImportError as exc:
        raise RuntimeError(
            "Hugging Face download requires huggingface_hub. Install it with "
            "`python -m pip install huggingface_hub`."
        ) from exc

    target_dir.mkdir(parents=True, exist_ok=True)
    print("Downloading BASALT models from Hugging Face: {}".format(HUGGING_FACE_REPO))
    snapshot_download(
        repo_id=HUGGING_FACE_REPO,
        revision=revision,
        local_dir=str(target_dir),
        allow_patterns=["*_ensemble.csv", "*_ensemble/**"],
        endpoint=endpoint,
        max_workers=8,
    )


def download_archive(url: str, target_dir: Path) -> Path:
    """Download a ZIP archive atomically and return its local path."""
    try:
        import requests
        from tqdm import tqdm
    except ImportError as exc:
        raise RuntimeError(
            "Archive downloads require requests and tqdm. Install them with "
            "`python -m pip install requests tqdm`."
        ) from exc

    target_dir.mkdir(parents=True, exist_ok=True)
    archive = target_dir / "BASALT.zip"
    if archive.is_file() and archive.stat().st_size > 0:
        print("Using existing archive: {}".format(archive))
        return archive

    partial_archive = archive.with_suffix(".zip.part")
    print("Downloading BASALT model archive: {}".format(url))
    print("Saving archive to: {}".format(archive))
    try:
        with requests.get(url, stream=True, timeout=(30, 300)) as response:
            response.raise_for_status()
            total_size = int(response.headers.get("content-length", 0))
            with partial_archive.open("wb") as output, tqdm(
                total=total_size or None,
                unit="iB",
                unit_scale=True,
                desc="BASALT models",
            ) as progress:
                for block in response.iter_content(chunk_size=1024 * 1024):
                    if block:
                        output.write(block)
                        progress.update(len(block))
        partial_archive.replace(archive)
    except Exception:
        partial_archive.unlink(missing_ok=True)
        raise
    return archive


def unpack_model(archive: Path, target_dir: Path) -> None:
    """Extract a BASALT model archive after checking every member path."""
    if not archive.is_file():
        raise FileNotFoundError("Model archive not found: {}".format(archive))
    target_dir.mkdir(parents=True, exist_ok=True)
    target_root = target_dir.resolve()
    try:
        with ZipFile(archive) as zipped:
            for member in zipped.infolist():
                member_path = (target_root / member.filename).resolve()
                if target_root not in member_path.parents and member_path != target_root:
                    raise ValueError(
                        "Unsafe path in model archive: {}".format(member.filename)
                    )
            zipped.extractall(target_root)
    except BadZipFile as exc:
        raise RuntimeError("Invalid model ZIP archive: {}".format(archive)) from exc


def use_local_archive(archive: Path, target_dir: Path) -> None:
    """Extract an existing archive without copying or deleting the source file."""
    print("Extracting local model archive: {}".format(archive))
    unpack_model(archive.expanduser().resolve(), target_dir)


def obtain_models(args: argparse.Namespace, target_dir: Path) -> str:
    """Obtain weights from the selected source and return its display name."""
    if models_are_ready(target_dir) and not args.force:
        print("A complete BASALT model set already exists; no download needed.")
        return "existing files"

    if args.archive:
        use_local_archive(Path(args.archive), target_dir)
        return "local archive"
    if args.url:
        archive = download_archive(args.url, target_dir)
        unpack_model(archive, target_dir)
        return "custom URL"

    if args.source == "archive":
        raise RuntimeError("--source archive requires --archive.")
    if args.source == "url":
        raise RuntimeError("--source url requires --url.")
    if args.source == "huggingface":
        download_huggingface(target_dir, args.hf_endpoint, args.revision)
        return "Hugging Face"
    if args.source == "figshare":
        archive = download_archive(FIGSHARE_URL, target_dir)
        unpack_model(archive, target_dir)
        return "Figshare"

    # Automatic policy: the official Hugging Face repository is concurrent and
    # resumable, while the legacy Figshare ZIP remains a compatibility fallback.
    try:
        download_huggingface(target_dir, args.hf_endpoint, args.revision)
        validate_models(target_dir)
        return "Hugging Face"
    except Exception as exc:
        print("Hugging Face download failed: {}".format(exc))
        print("Falling back to the Figshare archive ...")
        archive = download_archive(FIGSHARE_URL, target_dir)
        unpack_model(archive, target_dir)
        return "Figshare fallback"


def build_parser() -> argparse.ArgumentParser:
    """Define the command-line interface."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Download trained BASALT weights. Automatic mode tries the official "
            "Hugging Face repository before the Figshare archive."
        ),
    )
    parser.add_argument(
        "--path",
        help="Target directory. Defaults to BASALT_WEIGHT or ./BASALT_WEIGHT.",
    )
    parser.add_argument(
        "--source",
        choices=("auto", "huggingface", "figshare", "archive", "url"),
        default="auto",
        help="Model source.",
    )
    parser.add_argument("--archive", help="Existing BASALT model ZIP archive.")
    parser.add_argument("--url", help="Custom BASALT model ZIP URL.")
    parser.add_argument(
        "--hf-endpoint",
        help="Optional Hugging Face-compatible endpoint; official HF is the default.",
    )
    parser.add_argument(
        "--revision", default="main", help="Hugging Face branch, tag, or commit."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Refresh the selected source even if five descriptors already exist.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.archive and args.source not in ("auto", "archive"):
        parser.error("--archive is only compatible with --source auto/archive")
    if args.url and args.source not in ("auto", "url"):
        parser.error("--url is only compatible with --source auto/url")
    if args.archive and args.url:
        parser.error("--archive and --url are mutually exclusive")

    target_dir = resolve_target_dir(args.path).resolve()
    try:
        source = obtain_models(args, target_dir)
        descriptors = validate_models(target_dir)
    except Exception as exc:
        raise RuntimeError(
            "Model installation failed: {}\n"
            "Manual China-mainland fallback: {} (extraction code: embl). "
            "After downloading the ZIP, rerun with --source archive --archive FILE."
            .format(exc, BAIDU_URL)
        ) from exc

    print("Models ready in: {}".format(target_dir))
    print("Source: {}".format(source))
    print("Ensemble descriptors found: {}".format(len(descriptors)))
    print('export BASALT_WEIGHT="{}"'.format(target_dir))


if __name__ == "__main__":
    main()
