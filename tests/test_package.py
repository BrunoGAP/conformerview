"""Smoke tests for the initial project structure."""

from importlib import import_module

import conformer_analyzer


def test_package_version_is_available() -> None:
    assert conformer_analyzer.__version__ == "0.1.0"


def test_project_modules_import_cleanly() -> None:
    module_names = (
        "alignment",
        "interface",
        "models",
        "parsing",
        "rmsd",
        "validation",
        "visualization",
        "web_interface",
    )

    for module_name in module_names:
        import_module(f"conformer_analyzer.{module_name}")
