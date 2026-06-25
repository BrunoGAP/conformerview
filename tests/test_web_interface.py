"""Tests for the dependency-free local graphical interface."""

from pathlib import Path
from threading import Thread
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import pytest

from conformer_analyzer.elements import ATOMIC_NUMBER_TO_SYMBOL
from conformer_analyzer.parsing import parse_gaussian_log
from conformer_analyzer.web_interface import (
    ConformerHTTPServer,
    UploadedFile,
    WebApplication,
    create_server,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
COOH_LOGS = (
    PROJECT_ROOT / "conformers" / "COOH" / "conf1.log",
    PROJECT_ROOT / "conformers" / "COOH" / "conf2.log",
)


@pytest.fixture
def gui_server(tmp_path: Path):
    server = create_server(port=0, runs_directory=tmp_path)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield server
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=2)


def _uploads(*paths: Path) -> tuple[UploadedFile, ...]:
    return tuple(UploadedFile(path.name, path.read_bytes()) for path in paths)


def _xyz_upload_from_gaussian(path: Path, filename: str) -> UploadedFile:
    conformer = parse_gaussian_log(path).conformers[0]
    rows = "\n".join(
        f"{ATOMIC_NUMBER_TO_SYMBOL[atom.atomic_number]} "
        f"{atom.x:.8f} {atom.y:.8f} {atom.z:.8f}"
        for atom in conformer.geometry.atoms
    )
    content = f"{len(conformer.geometry.atoms)}\nconverted\n{rows}\n"
    return UploadedFile(filename, content.encode("utf-8"))


def _multipart_request(
    url: str,
    files: tuple[UploadedFile, ...],
    fields: dict[str, str] | None = None,
) -> Request:
    boundary = "----conformer-analyzer-test"
    chunks: list[bytes] = []
    for name, value in (fields or {}).items():
        chunks.extend(
            (
                f"--{boundary}\r\n".encode(),
                f'Content-Disposition: form-data; name="{name}"\r\n\r\n'.encode(),
                value.encode(),
                b"\r\n",
            )
        )
    for uploaded in files:
        chunks.extend(
            (
                f"--{boundary}\r\n".encode(),
                (
                    'Content-Disposition: form-data; name="log_files"; '
                    f'filename="{uploaded.filename}"\r\n'
                ).encode(),
                b"Content-Type: application/octet-stream\r\n\r\n",
                uploaded.content,
                b"\r\n",
            )
        )
    chunks.append(f"--{boundary}--\r\n".encode())
    return Request(
        url,
        data=b"".join(chunks),
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
        method="POST",
    )


def test_gui_defaults_are_independent_of_launch_directory() -> None:
    application = WebApplication()

    assert application.runs_directory == PROJECT_ROOT / "outputs" / "gui_runs"
    assert ConformerHTTPServer.allow_reuse_address is False


def test_gui_refuses_a_second_server_on_the_same_port(tmp_path: Path) -> None:
    first_server = create_server(port=0, runs_directory=tmp_path / "first")
    port = first_server.server_address[1]
    try:
        with pytest.raises(OSError):
            create_server(port=port, runs_directory=tmp_path / "second")
    finally:
        first_server.server_close()


def test_home_page_contains_gaussian_upload_controls(gui_server) -> None:
    with urlopen(gui_server.base_url + "/") as response:
        page = response.read().decode()

    assert response.status == 200
    assert "Choose conformer files" in page
    assert 'accept=".log,.xyz"' in page
    assert "Run RMSD analysis" in page


def test_gui_analysis_displays_results_and_opens_viewer(
    gui_server,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    opened_urls: list[str] = []
    monkeypatch.setattr(
        "conformer_analyzer.web_interface.webbrowser.open_new_tab",
        lambda url: opened_urls.append(url) or True,
    )
    request = _multipart_request(
        gui_server.base_url + "/analyze",
        _uploads(*COOH_LOGS),
        {"viewer_fit_atoms": "heavy", "show_hydrogens": "on"},
    )

    with urlopen(request) as response:
        page = response.read().decode()

    assert response.status == 200
    assert "Fitted RMSD results" in page
    assert "2.1200" in page
    assert "1.5507" in page
    assert "Scientific warnings" in page
    assert "TD optimization job 1" in page
    assert "Connectivity was inferred" in page
    assert len(opened_urls) == 1
    assert opened_urls[0].endswith("/viewer")

    with urlopen(opened_urls[0]) as viewer_response:
        viewer_html = viewer_response.read().decode()
    assert viewer_response.status == 200
    assert "Aligned conformer overlay" in viewer_html
    assert viewer_html.count('class="color-picker" type="color"') == 2
    assert "function setConformerColor" in viewer_html
    assert 'id="hydrogens" type="checkbox" checked' in viewer_html

    report_url = opened_urls[0].removesuffix("viewer") + "report"
    with urlopen(report_url) as report_response:
        report = report_response.read().decode()
    assert report_response.status == 200
    assert "2.1199556775" in report


def test_gui_requires_at_least_two_files(gui_server) -> None:
    request = _multipart_request(
        gui_server.base_url + "/analyze",
        _uploads(COOH_LOGS[0]),
    )

    with pytest.raises(HTTPError) as captured:
        urlopen(request)

    assert captured.value.code == 400
    assert "Select at least two Gaussian .log or XYZ .xyz files" in (
        captured.value.read().decode()
    )


def test_gui_accepts_mixed_gaussian_and_xyz_uploads(
    gui_server,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    opened_urls: list[str] = []
    monkeypatch.setattr(
        "conformer_analyzer.web_interface.webbrowser.open_new_tab",
        lambda url: opened_urls.append(url) or True,
    )
    request = _multipart_request(
        gui_server.base_url + "/analyze",
        (
            UploadedFile("conf2.log", COOH_LOGS[1].read_bytes()),
            _xyz_upload_from_gaussian(COOH_LOGS[1], "conf2_copy.xyz"),
        ),
    )

    with urlopen(request) as response:
        page = response.read().decode()

    assert response.status == 200
    assert "conf2_copy" in page
    assert "<td>0.0000</td>" in page
    assert len(opened_urls) == 1


def test_gui_rejects_unsupported_uploads(gui_server) -> None:
    request = _multipart_request(
        gui_server.base_url + "/analyze",
        (
            UploadedFile("geometry.sdf", b"not supported yet"),
            UploadedFile("conf2.log", COOH_LOGS[1].read_bytes()),
        ),
    )

    with pytest.raises(HTTPError) as captured:
        urlopen(request)

    assert captured.value.code == 400
    assert "Only Gaussian .log and XYZ .xyz files are accepted" in (
        captured.value.read().decode()
    )
