#!/usr/bin/env python3

"""Bounded release-smoke tests for BASALT's dependency-light interfaces."""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from zipfile import ZipFile
from pathlib import Path
from unittest.mock import patch


REPOSITORY = Path(__file__).resolve().parents[1]
MODULES = REPOSITORY / "BASALT"
sys.path.insert(0, str(MODULES))

from basalt_runtime import require_model_directory  # noqa: E402
from BASALT_models_download import unpack_model, validate_models  # noqa: E402
from S1e_extra_binners import (  # noqa: E402
    _metabinner_executable,
    _read_vamb_assignments,
    _write_assigned_bins,
    lorbin,
    vamb,
)


class ModelDirectoryTests(unittest.TestCase):
    def test_accepts_five_descriptor_and_checkpoint_pairs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            model_dir = Path(temporary)
            for index in range(5):
                descriptor = model_dir / f"model_{index}_ensemble.csv"
                descriptor.write_text("checkpoint.pth\n", encoding="utf-8")
                checkpoint_dir = descriptor.with_suffix("")
                checkpoint_dir.mkdir()
                (checkpoint_dir / "checkpoint.pth").write_bytes(b"test checkpoint")

            with patch.dict(os.environ, {"BASALT_WEIGHT": str(model_dir)}):
                self.assertEqual(require_model_directory(), model_dir.resolve())
                self.assertEqual(len(validate_models(model_dir)), 5)

    def test_rejects_an_unset_model_directory(self) -> None:
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaisesRegex(RuntimeError, "BASALT_WEIGHT is not set"):
                require_model_directory()

    def test_rejects_a_descriptor_with_a_missing_checkpoint(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            model_dir = Path(temporary)
            for index in range(5):
                descriptor = model_dir / f"model_{index}_ensemble.csv"
                descriptor.write_text("missing.pth\n", encoding="utf-8")
                descriptor.with_suffix("").mkdir()

            with self.assertRaisesRegex(RuntimeError, "Missing model checkpoint files"):
                validate_models(model_dir)

    def test_rejects_archive_path_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = root / "unsafe.zip"
            with ZipFile(archive, "w") as zipped:
                zipped.writestr("../outside.txt", "unsafe")

            with self.assertRaisesRegex(ValueError, "Unsafe path"):
                unpack_model(archive, root / "models")
            self.assertFalse((root / "outside.txt").exists())

    def test_rejects_checkpoint_path_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            model_dir = Path(temporary)
            for index in range(5):
                descriptor = model_dir / f"model_{index}_ensemble.csv"
                descriptor.write_text("../outside.pth\n", encoding="utf-8")
                descriptor.with_suffix("").mkdir()
            (model_dir / "outside.pth").write_bytes(b"not a model")

            with self.assertRaisesRegex(RuntimeError, "Unsafe checkpoint paths"):
                validate_models(model_dir)


class OptionalBinnerTests(unittest.TestCase):
    def test_metabinner_home_requires_an_executable_script(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            home = Path(temporary)
            script = home / "run_metabinner.sh"
            script.write_text("#!/bin/sh\n", encoding="utf-8")
            with patch.dict(
                os.environ,
                {"METABINNER_HOME": str(home), "PATH": ""},
                clear=False,
            ):
                with self.assertRaisesRegex(RuntimeError, "executable"):
                    _metabinner_executable()
                script.chmod(0o755)
                self.assertEqual(_metabinner_executable(), script.resolve())

    def test_reads_vamb_5_unsplit_cluster_table(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_dir = Path(temporary)
            (output_dir / "vae_clusters_unsplit.tsv").write_text(
                "clustername\tcontigname\n1\tcontig_a\n2\tcontig_b\n",
                encoding="utf-8",
            )
            self.assertEqual(
                _read_vamb_assignments(output_dir),
                {"contig_a": "1", "contig_b": "2"},
            )

    def test_materializes_cluster_assignments_as_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            assembly = root / "assembly.fa"
            assembly.write_text(
                ">contig_a\nAAAA\n>contig_b\nCCCCCC\n", encoding="utf-8"
            )
            output_dir = root / "bins"
            paths = _write_assigned_bins(
                assembly,
                {"contig_a": "cluster 1", "contig_b": "cluster 2"},
                output_dir,
                "candidate",
                minimum_size=4,
            )
            self.assertEqual([path.name for path in paths], [
                "candidate.cluster_1.fa",
                "candidate.cluster_2.fa",
            ])
            self.assertTrue(all(path.stat().st_size > 0 for path in paths))

    def test_vamb_adapter_uses_v5_command_and_sorted_bam(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "1_assembly.fa").write_text(
                ">contig_a\n" + "A" * 500_000 + "\n", encoding="utf-8"
            )
            (root / "1_DNA-1_sorted.bam").write_bytes(b"bam")
            commands = []

            def fake_run(command, cwd=None):
                commands.append([str(item) for item in command])
                output_dir = Path(command[command.index("--outdir") + 1])
                output_dir.mkdir()
                (output_dir / "vae_clusters_unsplit.tsv").write_text(
                    "clustername\tcontigname\n1\tcontig_a\n", encoding="utf-8"
                )

            with patch("S1e_extra_binners._run", side_effect=fake_run), patch(
                "S1e_extra_binners._run_quality_check"
            ):
                vamb(
                    "1_assembly.fa",
                    {"1": ["reads_R1.fastq", "reads_R2.fastq"]},
                    8,
                    str(root),
                    "checkm2",
                )

            self.assertEqual(commands[0][:3], ["vamb", "bin", "default"])
            self.assertIn("--bamfiles", commands[0])
            self.assertIn("--minfasta", commands[0])
            self.assertIn("-p", commands[0])
            self.assertIn(str((root / "1_DNA-1_sorted.bam").resolve()), commands[0])
            self.assertTrue((root / "1_assembly.fa_100_vamb_genomes").is_dir())

    def test_lorbin_adapter_uses_sorted_bams_and_multi_mode(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "1_assembly.fa").write_text(">contig_a\nAAAA\n", encoding="utf-8")
            for index in (1, 2):
                (root / f"1_DNA-{index}_sorted.bam").write_bytes(b"bam")
            commands = []

            def fake_run(command, cwd=None):
                commands.append([str(item) for item in command])
                output_dir = Path(command[command.index("-o") + 1])
                output_dir.mkdir()
                (output_dir / "lorbin.1.fa").write_text(
                    ">contig_a\nAAAA\n", encoding="utf-8"
                )

            with patch("S1e_extra_binners._run", side_effect=fake_run), patch(
                "S1e_extra_binners._run_quality_check"
            ):
                lorbin(
                    "1_assembly.fa",
                    {
                        "1": ["sample_1_R1.fastq", "sample_1_R2.fastq"],
                        "2": ["sample_2_R1.fastq", "sample_2_R2.fastq"],
                    },
                    8,
                    str(root),
                    "checkm2",
                )

            self.assertEqual(commands[0][:2], ["LorBin", "bin"])
            self.assertIn("--multi", commands[0])
            self.assertIn("--num_process", commands[0])
            self.assertIn(str((root / "1_DNA-1_sorted.bam").resolve()), commands[0])
            self.assertIn(str((root / "1_DNA-2_sorted.bam").resolve()), commands[0])
            self.assertTrue((root / "1_assembly.fa_100_lorbin_genomes").is_dir())


class CommandLineTests(unittest.TestCase):
    def run_cli(self, *arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [sys.executable, str(MODULES / "BASALT.py"), *arguments],
            cwd=REPOSITORY,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def test_help_is_available_without_pipeline_imports(self) -> None:
        result = self.run_cli("--help")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--assemblies", result.stdout)

    def test_rejects_nonpositive_threads(self) -> None:
        result = self.run_cli("--threads", "0")
        self.assertEqual(result.returncode, 2)
        self.assertIn("must be a positive integer", result.stderr)

    def test_rejects_malformed_paired_reads(self) -> None:
        result = self.run_cli("-a", "assembly.fa", "-s", "reads_R1.fastq")
        self.assertEqual(result.returncode, 2)
        self.assertIn("requires R1,R2", result.stderr)

    def test_rejects_a_normal_run_without_an_assembly(self) -> None:
        result = self.run_cli()
        self.assertEqual(result.returncode, 2)
        self.assertIn("requires at least one assembly", result.stderr)


if __name__ == "__main__":
    unittest.main()
