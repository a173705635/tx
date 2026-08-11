import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import os
from pathlib import Path
import subprocess
import sys

from antenna_diversity.analysis.plotting import (
    plot_acquisition,
    plot_results,
    save_figures,
)
from antenna_diversity.config import default_config
from antenna_diversity.pipeline import run_end_to_end


def test_acquisition_plot_contains_maps_slices_and_peak_markers(tmp_path):
    """Catch loss of the visual acquisition evidence requested by the user."""
    cfg = default_config()
    cfg.plot.enable = False
    results = run_end_to_end(cfg)

    figure = plot_acquisition(results)
    paths = save_figures({"acquisition_search": figure}, tmp_path)

    assert len(figure.axes) >= 6
    assert sum(bool(axis.images) for axis in figure.axes) >= 3
    assert sum(len(axis.lines) for axis in figure.axes) >= 9
    assert paths["acquisition_search"].is_file()
    assert paths["acquisition_search"].stat().st_size > 0
    plt.close(figure)


def test_plot_results_produces_four_named_figures(tmp_path):
    """Catch omission of a stage-level diagnostic figure from the full report."""
    cfg = default_config()
    cfg.plot.enable = False
    results = run_end_to_end(cfg)

    figures = plot_results(results)
    paths = save_figures(figures, tmp_path)

    assert set(figures) == {
        "rf_dfe_spectra",
        "array_response",
        "acquisition_search",
        "tracking_comparison",
    }
    assert all(path.is_file() and path.stat().st_size > 0 for path in paths.values())
    for figure in figures.values():
        plt.close(figure)


def test_cli_avoids_show_warning_with_noninteractive_backend():
    """Catch calls to plt.show() when the selected backend cannot display a GUI."""
    project_dir = Path(__file__).resolve().parents[1]
    environment = os.environ.copy()
    environment["MPLBACKEND"] = "Agg"

    completed = subprocess.run(
        [sys.executable, "-W", "error", "run_end_to_end.py"],
        cwd=project_dir,
        env=environment,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
