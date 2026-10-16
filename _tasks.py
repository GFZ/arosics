"""Cross-platform task helpers for pixi."""

import glob
import os
import shutil
import subprocess
import sys
import webbrowser
from pathlib import Path


def clean_build():
    for d in ["build", "dist", ".eggs"]:
        shutil.rmtree(d, ignore_errors=True)
    for p in glob.glob("**/*.egg-info", recursive=True):
        shutil.rmtree(p, ignore_errors=True)
    for p in glob.glob("**/*.egg", recursive=True):
        os.remove(p)


def clean_pyc():
    for pat in ["**/*.pyc", "**/*.pyo", "**/*~"]:
        for p in glob.glob(pat, recursive=True):
            os.remove(p)
    for p in glob.glob("**/__pycache__", recursive=True):
        shutil.rmtree(p, ignore_errors=True)


def clean_test():
    for f in [".coverage", "report.html", "report.xml", "coverage.xml"]:
        if os.path.exists(f):
            os.remove(f)
    for d in ["htmlcov", ".pytest_cache"]:
        shutil.rmtree(d, ignore_errors=True)
    for p in glob.glob(".coverage.*"):
        os.remove(p)


def docs(pkg):
    for f in [os.path.join("docs", pkg + ".rst"), os.path.join("docs", "modules.rst")]:
        if os.path.exists(f):
            os.remove(f)
    subprocess.run(
        ["sphinx-apidoc", pkg, "-o", "docs/", "--private",
         "--doc-project", "Python API reference"],
        check=True,
    )
    shutil.rmtree(os.path.join("docs", "_build"), ignore_errors=True)
    subprocess.run(
        ["sphinx-build", "-b", "html", "docs", os.path.join("docs", "_build", "html")],
        check=True,
    )
    webbrowser.open(os.path.join("docs", "_build", "html", "index.html"))


def dist():
    subprocess.run([sys.executable, "-m", "build"], check=True)
    for f in os.listdir("dist"):
        print(f)


def pytest(pkg):
    """Run pytest.

    - puts the coverage results in the folder 'htmlcov'
    - generates cobertura 'coverage.xml' (needed to show coverage in GitLab MR changes)
    - generates 'report.html' based on pytest-reporter-html1
    - generates JUnit 'report.xml' to show the test report as a new tab in a GitLab MR

    NOTE: Additional options pytest and coverage (plugin pytest-cov) are defined in .pytest.ini and .coveragerc.
    """
    import pytest
    return pytest.main(
        [
            "tests",
            "--verbosity=3",
            "--color=yes",
            "--tb=short",
            f"--cov={pkg}",
            "--cov-report=html:htmlcov",
            "--cov-report=term-missing",
            "--cov-report=xml:coverage.xml",
            "--template=html1/index.html",
            "--report=report.html",
            "--junitxml=report.xml",
        ]
    )


def lint(pkg):
    def _run_check(command):
        result = subprocess.run(command)
        if result.returncode != 0:
            sys.exit(result.returncode)

    _run_check(["flake8", "--max-line-length=120", pkg])
    _run_check(["pycodestyle", pkg, "--max-line-length=120"])

    # Report docstring violations, but don't fail the CI job.
    subprocess.run(
        ["pydocstyle", pkg],
        check=False,
    )


def urlcheck():
    try:
        result = subprocess.run(
            [
                "lychee",
                "**/*.md",
                "**/*.rst",
                "**/*.py",
                "**/*.json",
                "--no-progress",
                "--timeout", "2",
                "--verbose",
                "--exclude-path", ".pixi",
                "--exclude-path", ".git",
            ],
            timeout=120,
        )
    except subprocess.TimeoutExpired:
        print("ERROR: lychee timed out after 120 seconds.")
        raise SystemExit(1)

    if result.returncode != 0:
        raise SystemExit(result.returncode)


def create_env_info(env_name: str):
    """Capture information about the current Pixi environment."""
    out = Path(f"env_info__{env_name}")
    out.mkdir(exist_ok=True)

    (out / "pixi-list.txt").write_text(
        subprocess.check_output(["pixi", "list", "-e", env_name], text=True)
    )
    (out / "pixi-info.txt").write_text(
        subprocess.check_output(["pixi", "info"], text=True)
    )
    (out / "pixi-list.json").write_text(
        subprocess.check_output(
            ["pixi", "list", "-e", env_name, "--json"], text=True
        )
    )


if __name__ == "__main__":
    task = sys.argv[1] if len(sys.argv) > 1 else None
    extra = sys.argv[2] if len(sys.argv) > 2 else None
    tasks = {
        "clean-build": clean_build,
        "clean-pyc": clean_pyc,
        "clean-test": clean_test,
        "docs": lambda: docs(extra) if extra else sys.exit("docs requires package name"),
        "dist": dist,
        "pytest": lambda: pytest(extra) if extra else sys.exit("pytest requires package name"),
        "lint": lambda: lint(extra) if extra else sys.exit("lint requires package name"),
        "urlcheck": urlcheck,
        "create-env-info": lambda: (create_env_info(extra) if extra
                                    else sys.exit("env-info requires environment name"))
    }

    if task not in tasks:
        print("Unknown task: {}. Available: {}".format(task, ", ".join(tasks)))
        sys.exit(1)

    result = tasks[task]()

    if isinstance(result, int):
        sys.exit(result)
