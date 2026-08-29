"""Render a seamlessly looping turntable animation of the RL10 thrust chamber.

Loads ``validation/RL10/regen_results.pkl``, rebuilds the cooling-channel
geometry with :func:`pyskyfire.viz.make_engine_3d`, and orbits the camera one
full revolution about the vertical axis with the engine lying on its side.

Frames are rendered at ``360 / n_frames`` degree increments and the closing
frame is omitted, so the last frame is exactly one increment away from the
first and the loop is seamless.

The PNG frames are assembled into an animated WebP (or a GIF, if the output
path ends in .gif) with ffmpeg.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

import numpy as np
import pyvista as pv

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import pyskyfire as psf  # noqa: E402

REGEN_RESULTS = REPO_ROOT / "validation" / "RL10" / "regen_results.pkl"

# Matches the circuit colouring used by the published RL10 validation report.
PYSKYFIRE_RED = "#d62728"
CIRCUIT_COLORS = ("#2b2d31", PYSKYFIRE_RED)


def build_plotter(stride: int) -> pv.Plotter:
    if not REGEN_RESULTS.is_file():
        raise FileNotFoundError(
            f"Missing {REGEN_RESULTS}. Run validation/RL10/regen_sim.py first."
        )

    results = psf.common.Results.load(str(REGEN_RESULTS))
    viewer = psf.viz.make_engine_3d(
        results.thrust_chamber,
        show=False,
        stride=stride,
        circuit_colors=CIRCUIT_COLORS,
    )
    return viewer.plotter


def frame_paths(plotter, frame_dir, n_frames, size, elevation, zoom):
    """Render one full revolution, one frame per step.

    The engine lies on its side in the horizontal test-stand pose: its axis (x)
    runs across the frame and the camera orbits the world vertical (z). Because
    the orbit axis is perpendicular to the engine's symmetry axis, the full
    revolution is visible -- orbiting about the symmetry axis itself would
    produce a pixel-identical image at every angle.
    """

    plotter.window_size = list(size)
    plotter.set_background("white")
    # 360 tubes at this scale alias badly into moire fringes without it.
    plotter.enable_anti_aliasing("ssaa")

    plotter.camera_position = [
        (0.0, -1.0, 0.0),  # position, normalised below by reset_camera
        (0.0, 0.0, 0.0),   # focal point, likewise
        (0.0, 0.0, 1.0),   # view up = world vertical = the orbit axis
    ]
    plotter.reset_camera()
    # Raise the camera so it looks down on the engine. View up deliberately
    # stays on z (no OrthogonalizeViewUp) so the orbit remains a true turntable.
    plotter.camera.elevation = elevation
    plotter.camera.zoom(zoom)

    step = 360.0 / n_frames
    paths = []
    for i in range(n_frames):
        path = frame_dir / f"frame_{i:04d}.png"
        plotter.screenshot(str(path))
        paths.append(path)
        # pyvista's azimuth setter is absolute, not incremental: it undoes the
        # previous rotation before applying the new angle. Assigning a constant
        # step here would pin the camera at that one angle for every frame.
        plotter.camera.azimuth = (i + 1) * step
        plotter.render()
    return paths


def assemble_animation(
    frame_dir: Path, output: Path, fps: int, max_colors: int, quality: int
) -> None:
    """Mux the frames into an animated WebP or GIF, chosen by file extension.

    WebP is preferred and roughly halves the size of an equivalent GIF while
    keeping full 24-bit colour, because GIF has to quantise every frame to a
    256-entry palette. Note that ffmpeg can encode animated WebP but cannot
    decode it, so inspect the result with Pillow rather than ffprobe.
    """

    if shutil.which("ffmpeg") is None:
        raise RuntimeError("ffmpeg is required to assemble the animation")

    pattern = str(frame_dir / "frame_%04d.png")
    common = ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps),
              "-i", pattern]

    if output.suffix == ".webp":
        subprocess.run(
            common + ["-c:v", "libwebp_anim", "-lossless", "0",
                      "-q:v", str(quality), "-compression_level", "4",
                      "-loop", "0", "-an", str(output)],
            check=True,
        )
        return

    # GIF fallback. Error-diffusion dithering is deliberately disabled: it adds
    # per-pixel noise that differs between frames, which defeats GIF's
    # inter-frame compression and roughly quadruples the file size here.
    palette = frame_dir / "palette.png"
    subprocess.run(
        common + ["-vf", f"palettegen=max_colors={max_colors}:stats_mode=diff",
                  str(palette)],
        check=True,
    )
    subprocess.run(
        common + ["-i", str(palette), "-lavfi", "paletteuse=dither=none",
                  "-loop", "0", str(output)],
        check=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=REPO_ROOT / "images" / "RL10_turntable.webp",
                        help="Output path; .webp or .gif selects the encoder.")
    parser.add_argument("--frames", type=int, default=192,
                        help="Frames per revolution.")
    parser.add_argument("--fps", type=int, default=24,
                        help="Playback rate; frames/fps sets the loop period.")
    parser.add_argument("--width", type=int, default=880)
    parser.add_argument("--height", type=int, default=690)
    parser.add_argument("--stride", type=int, default=1,
                        help="Mesh subsampling passed to make_engine_3d.")
    parser.add_argument("--elevation", type=float, default=18.0,
                        help="Degrees the camera is raised above the engine.")
    parser.add_argument("--zoom", type=float, default=1.18)
    parser.add_argument("--max-colors", type=int, default=64,
                        help="Palette size for GIF output only.")
    parser.add_argument("--quality", type=int, default=60,
                        help="libwebp quality for WebP output only.")
    parser.add_argument("--keep-frames", type=Path, default=None,
                        help="Directory to retain the rendered PNG frames in.")
    args = parser.parse_args()

    plotter = build_plotter(args.stride)

    with tempfile.TemporaryDirectory() as td:
        frame_dir = Path(args.keep_frames) if args.keep_frames else Path(td)
        frame_dir.mkdir(parents=True, exist_ok=True)

        print(f"Rendering {args.frames} frames at {args.width}x{args.height}")
        frame_paths(plotter, frame_dir, args.frames,
                    (args.width, args.height), args.elevation, args.zoom)

        args.output.parent.mkdir(parents=True, exist_ok=True)
        assemble_animation(frame_dir, args.output, args.fps,
                           args.max_colors, args.quality)

    plotter.close()
    size_mb = args.output.stat().st_size / 1e6
    period = args.frames / args.fps
    print(f"Wrote {args.output} ({size_mb:.1f} MB, {period:.1f} s per revolution)")


if __name__ == "__main__":
    main()
