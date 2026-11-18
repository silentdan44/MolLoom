"""Optional Cargo build hook for the native polymer-growth core."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.dist import Distribution


ROOT = Path(__file__).parent.resolve()
RUST_BUILD_ENV = "MOLLOOM_BUILD_RUST"
LEGACY_RUST_BUILD_ENV = "POLYMER_GENERATOR_BUILD_RUST"


def rust_build_enabled() -> bool:
    enabled = {"1", "true", "yes"}
    return any(
        os.environ.get(name, "").strip().lower() in enabled
        for name in (RUST_BUILD_ENV, LEGACY_RUST_BUILD_ENV)
    )


def library_filename() -> str:
    if os.name == "nt":
        return "polymer_growth_core.dll"
    if sys.platform == "darwin":
        return "libpolymer_growth_core.dylib"
    return "libpolymer_growth_core.so"


def cargo_executable() -> Path | None:
    found = shutil.which("cargo")
    if found:
        return Path(found)

    candidates = [Path.home() / ".cargo" / "bin" / "cargo"]
    toolchains = Path.home() / ".rustup" / "toolchains"
    if toolchains.is_dir():
        candidates.extend(sorted(toolchains.glob("*/bin/cargo")))
    return next((path for path in candidates if path.is_file()), None)


class BuildPyWithRust(build_py):
    def run(self):
        super().run()
        if not rust_build_enabled():
            return

        cargo = cargo_executable()
        if cargo is None:
            raise RuntimeError(
                f"{RUST_BUILD_ENV}=1 requires a Rust toolchain; install Rust or unset the flag"
            )
        manifest = ROOT / "polymer_lib" / "rust_growth" / "Cargo.toml"
        environment = os.environ.copy()
        rustc = cargo.parent / "rustc"
        if rustc.is_file():
            environment["RUSTC"] = str(rustc)
        elif shutil.which("rustc") is None:
            raise RuntimeError("Cargo was found, but rustc is missing from the Rust toolchain")
        subprocess.run(
            [cargo, "build", "--release", "--manifest-path", str(manifest)],
            cwd=ROOT,
            env=environment,
            check=True,
        )
        source = manifest.parent / "target" / "release" / library_filename()
        destination = Path(self.build_lib) / "polymer_lib" / library_filename()
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)

    def get_outputs(self, include_bytecode=1):
        outputs = super().get_outputs(include_bytecode)
        if rust_build_enabled():
            outputs.append(
                str(Path(self.build_lib) / "polymer_lib" / library_filename())
            )
        return outputs


class OptionalRustDistribution(Distribution):
    def has_ext_modules(self):
        return rust_build_enabled()


setup(
    cmdclass={"build_py": BuildPyWithRust},
    distclass=OptionalRustDistribution,
)
