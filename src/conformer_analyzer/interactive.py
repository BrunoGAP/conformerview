"""Generate a lightweight interactive py3Dmol conformer overlay."""

import argparse
from dataclasses import dataclass
from html import escape
from pathlib import Path
import re
from typing import Sequence

import py3Dmol

from conformer_analyzer.elements import ATOMIC_NUMBER_TO_SYMBOL
from conformer_analyzer.models import Conformer, MoleculeGeometry
from conformer_analyzer.parsing import parse_conformer_file
from conformer_analyzer.visualization import (
    DEFAULT_COLORS,
    Bond,
    align_conformers_for_overlay,
    infer_bonds,
)


@dataclass(frozen=True, slots=True)
class InteractiveExportResult:
    """Metadata and warnings from an interactive viewer export."""

    output_file: Path
    bonds: tuple[Bond, ...]
    warnings: tuple[str, ...]


def export_interactive_overlay_html(
    conformers: Sequence[Conformer],
    output_file: str | Path,
    *,
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
    exclude_hydrogens_from_fit: bool = False,
    show_hydrogens: bool = False,
    width: int = 1200,
    height: int = 760,
) -> InteractiveExportResult:
    """Export aligned conformers as a self-contained viewer page shell.

    The molecular data and controls are embedded in the HTML. 3Dmol.js itself
    is loaded by py3Dmol from its version-matched CDN URL when the page opens.
    """

    conformer_tuple = tuple(conformers)
    if not conformer_tuple:
        raise ValueError("Interactive viewing requires at least one conformer.")
    if width <= 0 or height <= 0:
        raise ValueError("Viewer width and height must be positive.")

    display_labels = (
        tuple(labels)
        if labels is not None
        else tuple(conformer.name for conformer in conformer_tuple)
    )
    if len(display_labels) != len(conformer_tuple):
        raise ValueError("The number of labels must match the conformer count.")

    if colors is None:
        repeats = (len(conformer_tuple) + len(DEFAULT_COLORS) - 1) // len(
            DEFAULT_COLORS
        )
        display_colors = (DEFAULT_COLORS * repeats)[: len(conformer_tuple)]
    else:
        display_colors = tuple(colors)
    if len(display_colors) != len(conformer_tuple):
        raise ValueError("The number of colors must match the conformer count.")

    aligned_geometries = align_conformers_for_overlay(
        conformer_tuple,
        exclude_hydrogens=exclude_hydrogens_from_fit,
    )
    bonds = infer_bonds(aligned_geometries[0])
    if not bonds:
        raise ValueError("No bonds could be inferred for interactive viewing.")

    viewer = py3Dmol.view(
        width=width,
        height=height,
        options={"backgroundColor": "white"},
    )
    for index, (geometry, color) in enumerate(
        zip(aligned_geometries, display_colors)
    ):
        viewer.addModel(_mol_block(geometry, bonds, display_labels[index]), "mol")
        viewer.setStyle(
            {"model": index},
            {"stick": {"color": color, "radius": 0.18}},
        )
        if not show_hydrogens:
            viewer.setStyle({"model": index, "elem": "H"}, {})

    viewer.setProjection("orthographic")
    viewer.setViewStyle({"style": "outline", "color": "black", "width": 0.08})
    viewer.zoomTo()

    viewer_html = viewer._make_html()
    viewer_html = re.sub(
        r"width:\s*\d+px;\s*height:\s*\d+px;",
        "width: 100%; height: 100%;",
        viewer_html,
        count=1,
    )
    page = _viewer_page(
        viewer_html,
        viewer.uniqueid,
        display_labels,
        display_colors,
        show_hydrogens,
    )

    path = Path(output_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(page, encoding="utf-8")
    warnings = (
        "Connectivity was inferred from covalent radii and should be verified "
        "when explicit bond data becomes available.",
        "The viewer loads 3Dmol.js from the internet when opened.",
    )
    return InteractiveExportResult(path, bonds, warnings)


def _mol_block(
    geometry: MoleculeGeometry,
    bonds: tuple[Bond, ...],
    title: str,
) -> str:
    atom_lines: list[str] = []
    for atom in geometry.atoms:
        try:
            symbol = ATOMIC_NUMBER_TO_SYMBOL[atom.atomic_number]
        except KeyError as error:
            raise ValueError(
                "Cannot create interactive model: no element symbol is "
                f"configured for atomic number {atom.atomic_number}."
            ) from error
        atom_lines.append(
            f"{atom.x:10.4f}{atom.y:10.4f}{atom.z:10.4f} "
            f"{symbol:<3} 0  0  0  0  0  0  0  0  0  0  0  0"
        )
    bond_lines = [
        f"{first + 1:3d}{second + 1:3d}{1:3d}  0  0  0  0"
        for first, second in bonds
    ]
    counts = f"{len(geometry.atoms):3d}{len(bonds):3d}  0  0  0  0            999 V2000"
    return "\n".join(
        (
            title,
            "Conformer Analyzer / py3Dmol",
            "",
            counts,
            *atom_lines,
            *bond_lines,
            "M  END",
            "",
        )
    )


def _viewer_page(
    viewer_html: str,
    unique_id: str,
    labels: tuple[str, ...],
    colors: tuple[str, ...],
    show_hydrogens: bool,
) -> str:
    viewer_name = f"viewer_{unique_id}"
    model_controls = "\n".join(
        f'''<div class="model-control">
          <label class="model-toggle">
            <input type="checkbox" checked onchange="toggleModel({index}, this.checked)">
            <span id="swatch-{index}" class="swatch" style="background:{escape(color)}"></span>
            {escape(label)}
          </label>
          <label class="color-control">
            Color
            <input class="color-picker" type="color" value="{escape(color)}"
              aria-label="Choose color for {escape(label)}"
              title="Choose color for {escape(label)}"
              oninput="setConformerColor({index}, this.value)">
          </label>
        </div>'''
        for index, (label, color) in enumerate(zip(labels, colors))
    )
    color_array = ", ".join(f'"{escape(color)}"' for color in colors)
    hydrogen_checked = " checked" if show_hydrogens else ""
    return f'''<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Aligned conformer overlay</title>
  <style>
    * {{ box-sizing: border-box; }}
    html, body {{ width: 100%; height: 100%; margin: 0; background: #f5f6f8; }}
    body {{ display: grid; grid-template-rows: auto 1fr auto; font: 15px system-ui, sans-serif; color: #17191c; }}
    header {{ display: flex; flex-wrap: wrap; align-items: center; gap: 14px; padding: 12px 18px; background: white; border-bottom: 1px solid #dfe2e7; }}
    h1 {{ margin: 0 18px 0 0; font-size: 18px; }}
    button, .model-control, .hydrogen-toggle {{ border: 1px solid #c9ced6; border-radius: 7px; background: white; padding: 7px 10px; cursor: pointer; }}
    button:hover, .model-toggle:hover, .hydrogen-toggle:hover {{ background: #f0f3f7; }}
    .model-control {{ display: inline-flex; align-items: center; gap: 8px; padding: 3px 4px 3px 10px; }}
    .model-toggle, .hydrogen-toggle {{ display: inline-flex; align-items: center; gap: 7px; }}
    .model-toggle {{ cursor: pointer; }}
    .color-control {{ display: inline-flex; align-items: center; gap: 4px; color: #555c66; font-size: 12px; cursor: pointer; }}
    .swatch {{ width: 24px; height: 5px; border-radius: 5px; display: inline-block; }}
    .color-picker {{ width: 34px; height: 30px; padding: 2px; border: 0; border-radius: 5px; background: transparent; cursor: pointer; }}
    main {{ min-height: 0; padding: 12px; }}
    .viewer-shell {{ width: 100%; height: 100%; min-height: 480px; overflow: hidden; border-radius: 10px; background: white; box-shadow: 0 1px 5px #0002; }}
    footer {{ padding: 8px 18px 12px; color: #555c66; background: white; border-top: 1px solid #dfe2e7; }}
  </style>
</head>
<body>
  <header>
    <h1>Aligned conformer overlay</h1>
    {model_controls}
    <label class="hydrogen-toggle"><input id="hydrogens" type="checkbox"{hydrogen_checked} onchange="toggleHydrogens(this.checked)">Hydrogens</label>
    <button type="button" onclick="resetView()">Reset view</button>
    <button type="button" onclick="saveSnapshot()">Save PNG</button>
  </header>
  <main><div class="viewer-shell">{viewer_html}</div></main>
  <footer>Drag to rotate · Scroll to zoom · Right-drag to translate · Shift-drag to zoom</footer>
  <script>
    const conformerColors = [{color_array}];
    function getViewer() {{ return window.{viewer_name}; }}
    function toggleModel(index, visible) {{
      const viewer = getViewer();
      if (!viewer) return;
      const model = viewer.getModel(index);
      visible ? model.show() : model.hide();
      viewer.render();
    }}
    function setConformerColor(index, color) {{
      conformerColors[index] = color;
      document.getElementById(`swatch-${{index}}`).style.background = color;
      applyConformerStyle(index);
    }}
    function applyConformerStyle(index) {{
      const viewer = getViewer();
      if (!viewer) return;
      const color = conformerColors[index];
      viewer.setStyle(
        {{model: index}},
        {{stick: {{color: color, radius: 0.18}}}}
      );
      viewer.setStyle(
        {{model: index, elem: "H"}},
        document.getElementById("hydrogens").checked
          ? {{stick: {{color: color, radius: 0.11}}}}
          : {{}}
      );
      viewer.render();
    }}
    function toggleHydrogens(visible) {{
      const viewer = getViewer();
      if (!viewer) return;
      conformerColors.forEach((color, index) => {{
        viewer.setStyle(
          {{model: index, elem: "H"}},
          visible ? {{stick: {{color: color, radius: 0.11}}}} : {{}}
        );
      }});
      viewer.render();
    }}
    function resetView() {{
      const viewer = getViewer();
      if (!viewer) return;
      viewer.zoomTo();
      viewer.render();
    }}
    function saveSnapshot() {{
      const viewer = getViewer();
      if (!viewer) return;
      Promise.resolve(viewer.pngURI()).then((uri) => {{
        const link = document.createElement("a");
        link.href = uri;
        link.download = "aligned_conformer_view.png";
        link.click();
      }});
    }}
  </script>
</body>
</html>
'''


def main() -> None:
    """Generate the interactive viewer from supported conformer input paths."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="*", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/cooh_interactive_overlay.html"),
    )
    parser.add_argument("--show-hydrogens", action="store_true")
    parser.add_argument("--fit-heavy-atoms", action="store_true")
    arguments = parser.parse_args()

    input_paths = arguments.logs or sorted(Path("conformers/COOH").glob("*.log"))
    if not input_paths:
        parser.error("No conformer input files were provided or found.")
    conformers = [
        conformer
        for path in input_paths
        for conformer in parse_conformer_file(path).conformers
    ]
    result = export_interactive_overlay_html(
        conformers,
        arguments.output,
        labels=tuple(f"C{index}" for index in range(1, len(conformers) + 1)),
        exclude_hydrogens_from_fit=arguments.fit_heavy_atoms,
        show_hydrogens=arguments.show_hydrogens,
    )
    print(f"Interactive viewer written to {result.output_file.resolve()}")
    for warning in result.warnings:
        print(f"Warning: {warning}")


if __name__ == "__main__":
    main()
