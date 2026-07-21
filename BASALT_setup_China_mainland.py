#!/usr/bin/env python3

"""Install BASALT with mainland-China-friendly, non-persistent mirrors.

The installer prefers micromamba, then mamba, then conda.  Mirror settings
apply only to the child installation processes; the user's global Conda and
pip configuration files are never modified.
"""

from __future__ import print_function

import argparse
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


REPOSITORY_URL = "https://github.com/PKU-EMBL/BASALT.git"
MICROMAMBA_URL = "https://micro.mamba.pm/api/micromamba/{platform}/latest"

MIRRORS = {
    "tuna": {
        "label": "TUNA (Tsinghua University)",
        "channels": [
            "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge",
            "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda",
        ],
        "pip": "https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple",
    },
    "bfsu": {
        "label": "BFSU (Beijing Foreign Studies University)",
        "channels": [
            "https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge",
            "https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda",
        ],
        "pip": "https://mirrors.bfsu.edu.cn/pypi/web/simple",
    },
    "ustc": {
        "label": "USTC (University of Science and Technology of China)",
        "channels": [
            "https://mirrors.ustc.edu.cn/anaconda/cloud/conda-forge",
            "https://mirrors.ustc.edu.cn/anaconda/cloud/bioconda",
        ],
        "pip": "https://mirrors.ustc.edu.cn/pypi/simple",
    },
    "upstream": {
        "label": "upstream conda-forge and Bioconda",
        "channels": [
            "https://conda.anaconda.org/conda-forge",
            "https://conda.anaconda.org/bioconda",
        ],
        "pip": None,
    },
}


class InstallationError(RuntimeError):
    """A user-facing installation failure."""


def quote_command(command: Sequence[str]) -> str:
    """Return a shell-readable command without executing it through a shell."""
    return shlex.join([str(item) for item in command])


def run_command(
    command: Sequence[str],
    cwd: Optional[Path] = None,
    env: Optional[Dict[str, str]] = None,
    dry_run: bool = False,
    capture: bool = False,
) -> subprocess.CompletedProcess:
    """Run a command safely, or print it in dry-run mode."""
    location = " (in {})".format(cwd) if cwd else ""
    print("+ {}{}".format(quote_command(command), location))
    if dry_run:
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")
    try:
        return subprocess.run(
            [str(item) for item in command],
            cwd=str(cwd) if cwd else None,
            env=env,
            check=True,
            text=True,
            capture_output=capture,
        )
    except FileNotFoundError as exc:
        raise InstallationError("Command not found: {}".format(command[0])) from exc
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or exc.stdout or "").strip()
        message = "Command failed (exit {}): {}".format(
            exc.returncode, quote_command(command)
        )
        if detail:
            message += "\n{}".format(detail)
        raise InstallationError(message) from exc


def is_basalt_source(path: Path) -> bool:
    """Return whether *path* looks like the BASALT source repository."""
    return (
        (path / "basalt_environment.yml").is_file()
        and (path / "install.sh").is_file()
        and (path / "BASALT" / "BASALT.py").is_file()
    )


def resolve_source(args: argparse.Namespace) -> Path:
    """Use an existing checkout when possible, otherwise clone BASALT."""
    if args.source_dir:
        source = Path(args.source_dir).expanduser().resolve()
        if not is_basalt_source(source):
            raise InstallationError(
                "--source-dir is not a BASALT checkout: {}".format(source)
            )
        return source

    candidates = [Path(__file__).resolve().parent, Path.cwd().resolve()]
    for candidate in candidates:
        if is_basalt_source(candidate):
            print("Using existing BASALT source: {}".format(candidate))
            return candidate

    destination = Path(args.install_root).expanduser().resolve() / "BASALT"
    if destination.exists():
        raise InstallationError(
            "Clone destination exists but is not a BASALT checkout: {}".format(
                destination
            )
        )
    run_command(
        ["git", "clone", "--depth", "1", args.repo_url, str(destination)],
        dry_run=args.dry_run,
    )
    if args.dry_run:
        raise InstallationError(
            "Dry-run cannot inspect a checkout that has not been cloned yet; "
            "clone it first or pass --source-dir."
        )
    if not is_basalt_source(destination):
        raise InstallationError("The cloned repository is missing BASALT files.")
    return destination


def micromamba_platform() -> str:
    """Map the current Linux architecture to a micromamba artifact name."""
    machine = platform.machine().lower()
    mapping = {
        "x86_64": "linux-64",
        "amd64": "linux-64",
        "aarch64": "linux-aarch64",
        "arm64": "linux-aarch64",
        "ppc64le": "linux-ppc64le",
    }
    if sys.platform != "linux" or machine not in mapping:
        raise InstallationError(
            "Automatic micromamba bootstrap supports Linux x86-64, ARM64, "
            "and ppc64le only. Install micromamba manually and pass "
            "--manager-path."
        )
    return mapping[machine]


def bootstrap_micromamba(destination: Path, url_template: str) -> Path:
    """Download the official self-contained micromamba binary safely."""
    artifact = micromamba_platform()
    url = url_template.format(platform=artifact)
    destination = destination.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    print("Downloading micromamba from: {}".format(url))

    request = Request(url, headers={"User-Agent": "BASALT-installer/1.2"})
    with tempfile.TemporaryDirectory(prefix="basalt-micromamba-") as temp_name:
        archive = Path(temp_name) / "micromamba.tar.bz2"
        try:
            with urlopen(request, timeout=120) as response, archive.open("wb") as output:
                shutil.copyfileobj(response, output)
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            raise InstallationError(
                "Could not download micromamba from {}: {}".format(url, exc)
            ) from exc

        try:
            with tarfile.open(str(archive), "r:bz2") as bundle:
                member = bundle.getmember("bin/micromamba")
                source = bundle.extractfile(member)
                if source is None or not member.isfile():
                    raise InstallationError("Invalid micromamba binary in archive.")
                temporary_binary = destination.with_suffix(".part")
                with source, temporary_binary.open("wb") as output:
                    shutil.copyfileobj(source, output)
                temporary_binary.chmod(0o755)
                temporary_binary.replace(destination)
        except (KeyError, tarfile.TarError, OSError) as exc:
            raise InstallationError(
                "The downloaded archive does not contain a valid bin/micromamba."
            ) from exc

    print("Installed micromamba: {}".format(destination))
    return destination


def infer_manager(executable: str) -> str:
    """Infer a supported manager name from an executable path."""
    name = Path(executable).name.lower()
    if "micromamba" in name:
        return "micromamba"
    if name == "mamba":
        return "mamba"
    if name == "conda":
        return "conda"
    raise InstallationError(
        "Cannot infer manager type from {!r}; also pass --manager.".format(executable)
    )


def resolve_manager(args: argparse.Namespace) -> Tuple[str, str]:
    """Choose micromamba, mamba, or conda in that order."""
    if args.manager_path:
        executable = str(Path(args.manager_path).expanduser().resolve())
        if not Path(executable).is_file():
            raise InstallationError("--manager-path does not exist: {}".format(executable))
        kind = infer_manager(executable) if args.manager == "auto" else args.manager
        return kind, executable

    requested = ["micromamba", "mamba", "conda"] if args.manager == "auto" else [args.manager]
    for name in requested:
        executable = shutil.which(name)
        if executable:
            return name, executable
        if (
            name == "micromamba"
            and args.bootstrap_micromamba
            and args.manager in ("auto", "micromamba")
        ):
            destination = Path(args.micromamba_destination).expanduser().resolve()
            if destination.is_file():
                print("Using bootstrapped micromamba: {}".format(destination))
                return "micromamba", str(destination)
            if args.dry_run:
                print("Dry-run: micromamba would be bootstrapped to {}".format(destination))
                return "micromamba", str(destination)
            executable = bootstrap_micromamba(destination, args.micromamba_url)
            return "micromamba", str(executable)

    if args.dry_run and args.manager != "auto":
        print("Dry-run: using {!r} as a command placeholder.".format(args.manager))
        return args.manager, args.manager

    raise InstallationError(
        "No supported environment manager was found. Install micromamba, mamba, "
        "or conda; or rerun with --bootstrap-micromamba."
    )


def mirror_is_reachable(channels: Sequence[str], timeout: float = 4.0) -> bool:
    """Probe noarch repodata without changing package-manager state."""
    for channel in channels:
        url = channel.rstrip("/") + "/noarch/repodata.json.zst"
        request = Request(
            url,
            headers={"User-Agent": "BASALT-installer/1.2", "Range": "bytes=0-0"},
        )
        try:
            with urlopen(request, timeout=timeout) as response:
                if response.status not in (200, 206):
                    return False
        except (HTTPError, URLError, TimeoutError, OSError):
            return False
    return True


def resolve_mirror(
    args: argparse.Namespace,
) -> Tuple[str, List[str], Optional[str]]:
    """Resolve Conda channels and a matching PyPI endpoint."""
    if args.channel:
        pip_index = None if args.no_pip_mirror else args.pip_index_url
        return "custom", list(args.channel), pip_index

    name = args.mirror
    if name == "auto":
        if args.dry_run:
            name = "tuna"
            print("Dry-run: using TUNA as the representative automatic mirror.")
        else:
            for candidate in ("tuna", "bfsu", "ustc", "upstream"):
                preset = MIRRORS[candidate]
                print("Checking {} ...".format(preset["label"]))
                if mirror_is_reachable(preset["channels"]):
                    name = candidate
                    break
            else:
                raise InstallationError(
                    "No configured Conda mirror is reachable. Check proxy/DNS "
                    "settings or pass two explicit --channel URLs."
                )

    preset = MIRRORS[name]
    pip_index = args.pip_index_url or preset["pip"]
    if args.no_pip_mirror:
        pip_index = None
    print("Selected mirror: {}".format(preset["label"]))
    return name, list(preset["channels"]), pip_index


def render_environment(source: Path, destination: Path, channels: Sequence[str]) -> None:
    """Render a temporary environment file containing explicit channels."""
    lines = source.read_text(encoding="utf-8").splitlines()
    try:
        channel_start = next(
            index for index, line in enumerate(lines) if line.strip() == "channels:"
        )
        dependency_start = next(
            index
            for index, line in enumerate(lines[channel_start + 1 :], channel_start + 1)
            if line.strip() == "dependencies:"
        )
    except StopIteration as exc:
        raise InstallationError(
            "Invalid basalt_environment.yml: channels/dependencies section missing."
        ) from exc

    rendered = lines[:channel_start] + ["channels:"]
    rendered.extend("  - {}".format(channel.rstrip("/")) for channel in channels)
    rendered.append("  - nodefaults")
    rendered.extend(lines[dependency_start:])
    destination.write_text("\n".join(rendered) + "\n", encoding="utf-8")


def write_temporary_condarc(destination: Path, channels: Sequence[str]) -> None:
    """Create an isolated strict-priority Conda configuration."""
    lines = ["channels:"]
    lines.extend("  - {}".format(channel.rstrip("/")) for channel in channels)
    lines.extend(
        [
            "channel_priority: strict",
            "show_channel_urls: true",
            "remote_connect_timeout_secs: 20",
            "remote_read_timeout_secs: 120",
            "remote_max_retries: 5",
        ]
    )
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def installation_environment(
    manager: str,
    condarc: Path,
    pip_index: Optional[str],
    mamba_root_prefix: Path,
) -> Dict[str, str]:
    """Build child-only environment variables for package downloads."""
    environment = os.environ.copy()
    environment["CONDARC"] = str(condarc)
    environment["CONDA_CHANNEL_PRIORITY"] = "strict"
    environment["MAMBA_CHANNEL_PRIORITY"] = "strict"
    environment["PIP_DEFAULT_TIMEOUT"] = "120"
    environment["PIP_RETRIES"] = "10"
    if manager == "micromamba":
        environment["MAMBA_ROOT_PREFIX"] = str(mamba_root_prefix)
    if pip_index:
        environment["PIP_INDEX_URL"] = pip_index
        environment["UV_DEFAULT_INDEX"] = pip_index
    return environment


def environment_exists(
    executable: str,
    env_name: str,
    environment: Dict[str, str],
    dry_run: bool,
) -> bool:
    """Check for a named environment through the selected manager."""
    if dry_run:
        return False
    try:
        result = run_command(
            [executable, "env", "list", "--json"],
            env=environment,
            capture=True,
        )
        paths = json.loads(result.stdout).get("envs", [])
    except (InstallationError, ValueError, TypeError):
        return False
    return any(Path(path).name == env_name for path in paths)


def environment_command(
    manager: str,
    executable: str,
    env_name: str,
    environment_file: Path,
    exists: bool,
    update: bool,
    offline: bool,
) -> List[str]:
    """Build a create/update command for the selected manager."""
    if exists and not update:
        return []

    if exists:
        command = [
            executable,
            "env",
            "update",
            "--name",
            env_name,
            "--file",
            str(environment_file),
        ]
        if manager == "conda":
            command.extend(["--prune", "--solver", "libmamba"])
        else:
            command.append("--prune")
    elif manager == "micromamba":
        command = [
            executable,
            "create",
            "--name",
            env_name,
            "--file",
            str(environment_file),
        ]
    else:
        command = [
            executable,
            "env",
            "create",
            "--name",
            env_name,
            "--file",
            str(environment_file),
        ]
        if manager == "conda":
            command.extend(["--solver", "libmamba"])

    command.append("--yes")
    if offline:
        command.append("--offline")
    return command


def manager_run(executable: str, env_name: str, command: Sequence[str]) -> List[str]:
    """Run a command without requiring shell activation."""
    return [executable, "run", "--name", env_name] + [str(item) for item in command]


def install_models(
    args: argparse.Namespace,
    executable: str,
    env_name: str,
    source: Path,
    environment: Dict[str, str],
) -> None:
    """Invoke the canonical model downloader inside the BASALT environment."""
    if args.model_source == "none":
        print("Skipping model weights (--model-source none).")
        return

    model_dir = Path(args.model_dir).expanduser().resolve()
    downloader = source / "BASALT" / "BASALT_models_download.py"
    command = [
        "python",
        str(downloader),
        "--path",
        str(model_dir),
        "--source",
        args.model_source,
    ]
    if args.model_archive:
        command.extend(["--archive", str(Path(args.model_archive).expanduser().resolve())])
    if args.model_url:
        command.extend(["--url", args.model_url])
    if args.hf_endpoint:
        command.extend(["--hf-endpoint", args.hf_endpoint])
    run_command(
        manager_run(executable, env_name, command),
        env=environment,
        dry_run=args.dry_run,
    )


def build_parser() -> argparse.ArgumentParser:
    """Define the command-line interface."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Install BASALT on Linux with selectable China-mainland mirrors. "
            "No global Conda or pip configuration is modified."
        ),
    )
    parser.add_argument("--source-dir", help="Use an existing BASALT checkout.")
    parser.add_argument(
        "--install-root",
        default=str(Path.cwd()),
        help="Parent directory used when the repository must be cloned.",
    )
    parser.add_argument("--repo-url", default=REPOSITORY_URL, help="Git clone URL.")
    parser.add_argument(
        "--manager",
        choices=("auto", "micromamba", "mamba", "conda"),
        default="auto",
        help="Environment manager; auto prefers micromamba.",
    )
    parser.add_argument(
        "--manager-path", help="Explicit micromamba, mamba, or conda executable."
    )
    parser.add_argument(
        "--bootstrap-micromamba",
        action="store_true",
        help="Download the official standalone micromamba when none is installed.",
    )
    parser.add_argument(
        "--micromamba-destination",
        default=str(Path.home() / ".local" / "bin" / "micromamba"),
        help="Destination used by --bootstrap-micromamba.",
    )
    parser.add_argument(
        "--micromamba-url",
        default=MICROMAMBA_URL,
        help="Bootstrap URL template; {platform} is replaced automatically.",
    )
    parser.add_argument(
        "--mamba-root-prefix",
        default=os.environ.get("MAMBA_ROOT_PREFIX", str(Path.home() / "micromamba")),
        help="Root prefix for micromamba environments and its shared package cache.",
    )
    parser.add_argument(
        "--mirror",
        choices=("auto", "tuna", "bfsu", "ustc", "upstream"),
        default="auto",
        help="Conda/PyPI mirror preset; auto probes several endpoints.",
    )
    parser.add_argument(
        "--channel",
        action="append",
        metavar="URL",
        help="Explicit Conda channel URL; repeat for conda-forge and Bioconda.",
    )
    parser.add_argument("--pip-index-url", help="Explicit temporary PyPI index URL.")
    parser.add_argument(
        "--no-pip-mirror",
        action="store_true",
        help="Do not override the inherited/default PyPI index.",
    )
    parser.add_argument("--env-name", default="basalt", help="Environment name.")
    parser.add_argument(
        "--update",
        action="store_true",
        help="Update and prune an existing environment instead of reusing it.",
    )
    parser.add_argument(
        "--offline", action="store_true", help="Use only existing package caches."
    )
    parser.add_argument(
        "--model-source",
        choices=("auto", "huggingface", "figshare", "archive", "url", "none"),
        default="auto",
        help="Weight source; auto tries Hugging Face before Figshare.",
    )
    parser.add_argument(
        "--model-dir",
        default=str(Path.cwd() / "BASALT_WEIGHT"),
        help="Persistent model-weight directory.",
    )
    parser.add_argument("--model-archive", help="Existing BASALT model ZIP archive.")
    parser.add_argument("--model-url", help="Custom BASALT model ZIP URL.")
    parser.add_argument(
        "--hf-endpoint",
        help="Optional Hugging Face-compatible endpoint; official HF is the default.",
    )
    parser.add_argument(
        "--skip-verify", action="store_true", help="Skip the BASALT --help smoke test."
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="Print commands without changing files."
    )
    return parser


def validate_arguments(args: argparse.Namespace) -> None:
    """Reject conflicting or unsafe option combinations."""
    if not args.dry_run and sys.platform != "linux":
        raise InstallationError(
            "BASALT's Conda installation is supported on Linux. Use --dry-run "
            "to inspect commands on another platform."
        )
    if args.channel and len(args.channel) < 2:
        raise InstallationError(
            "Repeat --channel for both conda-forge and Bioconda, in that order."
        )
    if args.model_source == "archive" and not args.model_archive:
        raise InstallationError("--model-source archive requires --model-archive.")
    if args.model_source == "url" and not args.model_url:
        raise InstallationError("--model-source url requires --model-url.")
    if args.model_archive and args.model_source not in ("auto", "archive"):
        raise InstallationError(
            "--model-archive is only compatible with --model-source archive/auto."
        )
    if args.model_url and args.model_source not in ("auto", "url"):
        raise InstallationError(
            "--model-url is only compatible with --model-source url/auto."
        )
    if args.model_archive and args.model_url:
        raise InstallationError(
            "--model-archive and --model-url are mutually exclusive."
        )


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Install the environment, BASALT scripts, and optional model weights."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        validate_arguments(args)
        source = resolve_source(args)
        manager, executable = resolve_manager(args)
        print("Environment manager: {} ({})".format(manager, executable))
        if not args.dry_run:
            run_command([executable, "--version"])

        _, channels, pip_index = resolve_mirror(args)
        print("Conda channels (highest priority first):")
        for channel in channels:
            print("  - {}".format(channel))
        print("PyPI index: {}".format(pip_index or "inherited/default"))

        with tempfile.TemporaryDirectory(prefix="basalt-install-") as temp_name:
            temporary = Path(temp_name)
            environment_file = temporary / "basalt_environment.yml"
            condarc = temporary / "condarc.yml"
            render_environment(
                source / "basalt_environment.yml", environment_file, channels
            )
            write_temporary_condarc(condarc, channels)
            child_environment = installation_environment(
                manager,
                condarc,
                pip_index,
                Path(args.mamba_root_prefix).expanduser().resolve(),
            )

            exists = environment_exists(
                executable,
                args.env_name,
                child_environment,
                args.dry_run,
            )
            if exists and not args.update:
                print(
                    "Environment {!r} already exists; reusing it. Pass --update "
                    "to synchronize dependencies.".format(args.env_name)
                )
            command = environment_command(
                manager,
                executable,
                args.env_name,
                environment_file,
                exists,
                args.update,
                args.offline,
            )
            if command:
                run_command(
                    command,
                    env=child_environment,
                    dry_run=args.dry_run,
                )

            run_command(
                manager_run(executable, args.env_name, ["bash", "install.sh"]),
                cwd=source,
                env=child_environment,
                dry_run=args.dry_run,
            )
            install_models(
                args,
                executable,
                args.env_name,
                source,
                child_environment,
            )
            if not args.skip_verify:
                run_command(
                    manager_run(executable, args.env_name, ["BASALT", "--help"]),
                    env=child_environment,
                    dry_run=args.dry_run,
                )

        model_dir = Path(args.model_dir).expanduser().resolve()
        print("\nBASALT installation completed.")
        if manager == "micromamba":
            root_prefix = Path(args.mamba_root_prefix).expanduser().resolve()
            print(
                "Run without activation: MAMBA_ROOT_PREFIX={} {}".format(
                    shlex.quote(str(root_prefix)),
                    quote_command(
                        [
                            executable,
                            "run",
                            "--name",
                            args.env_name,
                            "BASALT",
                            "--help",
                        ]
                    ),
                )
            )
            print("For bash activation:")
            print('  export MAMBA_ROOT_PREFIX="{}"'.format(root_prefix))
            print(
                '  eval "$({})"'.format(
                    quote_command([executable, "shell", "hook", "--shell", "bash"])
                )
            )
            print("  micromamba activate {}".format(args.env_name))
        else:
            print("Activate: conda activate {}".format(args.env_name))
        if args.model_source != "none":
            print('Set models: export BASALT_WEIGHT="{}"'.format(model_dir))
        print("Download CheckM2 database: checkm2 database --download")
        return 0
    except InstallationError as exc:
        print("\nERROR: {}".format(exc), file=sys.stderr)
        print(
            "Inspect the selected route with --dry-run, or retry with an explicit "
            "--mirror/--channel after checking the error above.",
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    sys.exit(main())
