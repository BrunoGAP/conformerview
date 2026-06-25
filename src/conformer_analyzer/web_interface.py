"""Dependency-free local browser interface for conformer comparison."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from email import policy
from email.parser import BytesParser
from html import escape
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import re
import socket
from threading import Timer
from typing import Mapping, Sequence
from urllib.parse import urlsplit
from uuid import UUID, uuid4
import webbrowser

from conformer_analyzer.interface import (
    SUPPORTED_INPUT_SUFFIXES,
    AnalysisRun,
    analyze_log_files,
)


DEFAULT_HOST = "127.0.0.1"
DEFAULT_PORT = 5000
MAX_UPLOAD_BYTES = 100 * 1024 * 1024
PACKAGE_DIRECTORY = Path(__file__).resolve().parent
PROJECT_DIRECTORY = PACKAGE_DIRECTORY.parent.parent


@dataclass(frozen=True, slots=True)
class UploadedFile:
    """One file received from the local multipart form."""

    filename: str
    content: bytes


@dataclass(frozen=True, slots=True)
class AnalysisPage:
    """Successful analysis data needed to render the result page."""

    run_id: str
    result: AnalysisRun
    viewer_url: str
    viewer_opened: bool


class WebApplication:
    """State and rendering for the local HTTP interface."""

    def __init__(self, runs_directory: str | Path | None = None) -> None:
        self.runs_directory = (
            Path(runs_directory)
            if runs_directory is not None
            else PROJECT_DIRECTORY / "outputs" / "gui_runs"
        )
        self.completed_runs: dict[str, AnalysisRun] = {}

    def render_page(
        self,
        *,
        error: str | None = None,
        analysis: AnalysisPage | None = None,
    ) -> str:
        template = (PACKAGE_DIRECTORY / "templates" / "index.html").read_text(
            encoding="utf-8"
        )
        error_block = ""
        if error:
            error_block = (
                '<div class="alert error" role="alert">'
                "<strong>Analysis could not run.</strong> "
                f"{escape(error)}</div>"
            )
        results_block = _render_results(analysis) if analysis else ""
        return template.replace("{{ERROR_BLOCK}}", error_block).replace(
            "{{RESULTS_BLOCK}}", results_block
        )

    def analyze(
        self,
        uploads: Sequence[UploadedFile],
        form: Mapping[str, str],
        base_url: str,
    ) -> AnalysisPage:
        upload_tuple = tuple(upload for upload in uploads if upload.filename)
        validation_error = _validate_uploads(upload_tuple)
        if validation_error:
            raise ValueError(validation_error)

        viewer_fit_atoms = form.get("viewer_fit_atoms", "all")
        show_hydrogens = form.get("show_hydrogens") == "on"
        run_id = uuid4().hex
        run_directory = self.runs_directory / run_id
        log_paths = _save_uploads(upload_tuple, run_directory / "uploads")
        result = analyze_log_files(
            log_paths,
            run_directory / "results",
            viewer_fit_atoms=viewer_fit_atoms,
            show_hydrogens=show_hydrogens,
            open_viewer=False,
        )
        self.completed_runs[run_id] = result
        viewer_url = f"{base_url.rstrip('/')}/runs/{run_id}/viewer"
        viewer_opened = webbrowser.open_new_tab(viewer_url)
        return AnalysisPage(run_id, result, viewer_url, viewer_opened)

    def get_run(self, run_id: str) -> AnalysisRun | None:
        try:
            UUID(hex=run_id)
        except ValueError:
            return None
        return self.completed_runs.get(run_id)


class ConformerHTTPServer(ThreadingHTTPServer):
    """HTTP server carrying application state and its public local URL."""

    # On Windows, SO_REUSEADDR can allow two local servers to listen on the
    # same port. Refusing that prevents a new GUI from competing with a stale
    # process that still has older application code in memory.
    allow_reuse_address = False

    def server_bind(self) -> None:
        """Bind exclusively so Windows cannot run two GUIs on one port."""

        if hasattr(socket, "SO_EXCLUSIVEADDRUSE"):
            self.socket.setsockopt(
                socket.SOL_SOCKET,
                socket.SO_EXCLUSIVEADDRUSE,
                1,
            )
        super().server_bind()

    def __init__(
        self,
        server_address: tuple[str, int],
        application: WebApplication,
    ) -> None:
        super().__init__(server_address, ConformerRequestHandler)
        self.application = application

    @property
    def base_url(self) -> str:
        host, port = self.server_address[:2]
        display_host = "127.0.0.1" if host in {"", "0.0.0.0"} else host
        return f"http://{display_host}:{port}"


class ConformerRequestHandler(BaseHTTPRequestHandler):
    """Serve the upload form, generated results, and static assets."""

    server: ConformerHTTPServer

    def do_GET(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler API
        path = urlsplit(self.path).path
        if path == "/":
            self._send_html(self.server.application.render_page())
            return
        if path == "/static/app.css":
            self._send_file(PACKAGE_DIRECTORY / "static" / "app.css", "text/css")
            return
        if path == "/static/app.js":
            self._send_file(
                PACKAGE_DIRECTORY / "static" / "app.js",
                "text/javascript",
            )
            return

        match = re.fullmatch(r"/runs/([0-9a-f]{32})/(viewer|report)", path)
        if match is None:
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        run_id, artifact = match.groups()
        result = self.server.application.get_run(run_id)
        if result is None:
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        if artifact == "viewer":
            self._send_file(result.viewer_file, "text/html")
        else:
            self._send_file(
                result.matrix_file,
                "text/plain; charset=utf-8",
                download_name="rmsd_matrix.txt",
            )

    def do_POST(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler API
        if urlsplit(self.path).path != "/analyze":
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        try:
            uploads, form = self._read_multipart_form()
            analysis = self.server.application.analyze(
                uploads,
                form,
                self.server.base_url,
            )
        except (OSError, ValueError) as error:
            self._send_html(
                self.server.application.render_page(error=str(error)),
                status=HTTPStatus.BAD_REQUEST,
            )
            return
        self._send_html(self.server.application.render_page(analysis=analysis))

    def _read_multipart_form(
        self,
    ) -> tuple[tuple[UploadedFile, ...], dict[str, str]]:
        try:
            content_length = int(self.headers.get("Content-Length", "0"))
        except ValueError as error:
            raise ValueError("The upload request has an invalid size.") from error
        if content_length <= 0:
            raise ValueError("No upload data was received.")
        if content_length > MAX_UPLOAD_BYTES:
            raise ValueError("The combined upload is larger than the 100 MB limit.")
        content_type = self.headers.get("Content-Type", "")
        if not content_type.lower().startswith("multipart/form-data"):
            raise ValueError("The upload form must use multipart data.")

        body = self.rfile.read(content_length)
        message = BytesParser(policy=policy.default).parsebytes(
            b"Content-Type: "
            + content_type.encode("ascii", errors="replace")
            + b"\r\nMIME-Version: 1.0\r\n\r\n"
            + body
        )
        if not message.is_multipart():
            raise ValueError("The uploaded form data could not be read.")

        uploads: list[UploadedFile] = []
        fields: dict[str, str] = {}
        for part in message.iter_parts():
            field_name = part.get_param("name", header="content-disposition")
            if not field_name:
                continue
            filename = part.get_filename()
            payload = part.get_payload(decode=True) or b""
            if filename is not None:
                uploads.append(UploadedFile(filename, payload))
            else:
                charset = part.get_content_charset() or "utf-8"
                fields[field_name] = payload.decode(charset, errors="replace")
        return tuple(uploads), fields

    def _send_html(
        self,
        content: str,
        *,
        status: HTTPStatus = HTTPStatus.OK,
    ) -> None:
        payload = content.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def _send_file(
        self,
        path: Path,
        content_type: str,
        *,
        download_name: str | None = None,
    ) -> None:
        try:
            payload = path.resolve().read_bytes()
        except OSError:
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(payload)))
        if download_name:
            self.send_header(
                "Content-Disposition",
                f'attachment; filename="{download_name}"',
            )
        self.end_headers()
        self.wfile.write(payload)

    def log_message(self, format: str, *args: object) -> None:
        """Keep routine local requests quiet; errors still reach the page."""


def create_server(
    *,
    host: str = DEFAULT_HOST,
    port: int = DEFAULT_PORT,
    runs_directory: str | Path | None = None,
) -> ConformerHTTPServer:
    """Create a configured local server without starting its event loop."""

    return ConformerHTTPServer((host, port), WebApplication(runs_directory))


def _validate_uploads(uploads: tuple[UploadedFile, ...]) -> str | None:
    if len(uploads) < 2:
        return "Select at least two Gaussian .log or XYZ .xyz files."
    invalid = [
        upload.filename
        for upload in uploads
        if Path(upload.filename).suffix.lower() not in SUPPORTED_INPUT_SUFFIXES
    ]
    if invalid:
        return "Only Gaussian .log and XYZ .xyz files are accepted: " + ", ".join(
            invalid
        )
    return None


def _safe_filename(filename: str, fallback: str) -> str:
    basename = Path(filename.replace("\\", "/")).name
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", basename).strip("._")
    return sanitized or fallback


def _save_uploads(
    uploads: tuple[UploadedFile, ...],
    upload_directory: Path,
) -> tuple[Path, ...]:
    paths: list[Path] = []
    for index, upload in enumerate(uploads, start=1):
        filename = _safe_filename(upload.filename, f"conformer_{index}.dat")
        destination = upload_directory / f"{index:03d}" / filename
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(upload.content)
        paths.append(destination)
    return tuple(paths)


def _render_results(analysis: AnalysisPage) -> str:
    result = analysis.result
    labels = tuple(f"C{index}" for index in range(1, len(result.conformers) + 1))
    conformer_cards = "".join(
        '<article class="conformer-card">'
        f'<span class="conformer-label">{label}</span><div>'
        f"<strong>{escape(conformer.name)}</strong><small>"
        f"{len(conformer.geometry.atoms)} atoms · charge "
        f"{conformer.charge if conformer.charge is not None else 'unknown'} · "
        "multiplicity "
        f"{conformer.multiplicity if conformer.multiplicity is not None else 'unknown'}"
        "</small></div></article>"
        for label, conformer in zip(labels, result.conformers)
    )
    warning_panel = ""
    if result.warnings:
        warning_items = "".join(
            f"<li>{escape(warning)}</li>" for warning in result.warnings
        )
        warning_panel = (
            '<div class="warning-panel"><h3>Scientific warnings</h3>'
            f"<ul>{warning_items}</ul></div>"
        )
    auto_open_notice = ""
    if not analysis.viewer_opened:
        auto_open_notice = (
            '<div class="alert info">The visualization tab could not be opened '
            "automatically. Use “Open 3D overlay” above.</div>"
        )
    report_url = f"/runs/{analysis.run_id}/report"
    viewer_url = escape(analysis.viewer_url, quote=True)
    return f"""
      <section class="results" aria-labelledby="results-title">
        <div class="results-heading"><div><p class="eyebrow dark">Analysis complete</p>
          <h2 id="results-title">Fitted RMSD results</h2></div>
          <div class="actions"><a class="secondary-button" href="{report_url}">Download report</a>
            <a class="primary-button compact" href="{viewer_url}" target="_blank" rel="noopener">Open 3D overlay ↗</a></div></div>
        {auto_open_notice}
        <div class="conformer-grid">{conformer_cards}</div>
        {warning_panel}
        <div class="matrix-grid">
          {_render_matrix_card('All atoms', 'Alignment and RMSD include hydrogen atoms.', labels, result.all_atom_matrix)}
          {_render_matrix_card('Heavy atoms', 'Hydrogen atoms are excluded from alignment and RMSD.', labels, result.heavy_atom_matrix)}
        </div>
        <p class="units">RMSD values are in ångström (Å). Every comparison uses optimal rigid-body alignment.</p>
      </section>
    """


def _render_matrix_card(
    title: str,
    description: str,
    labels: tuple[str, ...],
    matrix: tuple[tuple[float, ...], ...],
) -> str:
    headings = "".join(f"<th>{label}</th>" for label in labels)
    rows = "".join(
        f"<tr><th>{label}</th>"
        + "".join(f"<td>{value:.4f}</td>" for value in row)
        + "</tr>"
        for label, row in zip(labels, matrix)
    )
    return (
        '<article class="matrix-card">'
        f"<div><h3>{title}</h3><p>{description}</p></div>"
        '<div class="table-wrap"><table><thead><tr><th></th>'
        f"{headings}</tr></thead><tbody>{rows}</tbody></table></div></article>"
    )


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="conformer-analyzer-gui",
        description="Open the local graphical interface for conformer analysis.",
    )
    parser.add_argument("--host", default=DEFAULT_HOST)
    parser.add_argument("--port", type=int, default=DEFAULT_PORT)
    parser.add_argument("--no-open", action="store_true")
    return parser


def main() -> None:
    """Start the local graphical interface."""

    parser = _argument_parser()
    arguments = parser.parse_args()
    try:
        server = create_server(host=arguments.host, port=arguments.port)
    except OSError:
        parser.error(
            f"Cannot start on {arguments.host}:{arguments.port}. Close the "
            "existing Conformer Analyzer GUI before launching another one."
        )
    if not arguments.no_open:
        Timer(0.4, webbrowser.open_new_tab, args=(server.base_url + "/",)).start()
    print(f"Conformer Analyzer GUI: {server.base_url}/")
    print("Press Ctrl+C to stop the local server.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
