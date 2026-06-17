# Watercolor: physically-simulated digital watercolor painter

A single-file C/SDL2 program that simulates real watercolor behavior, with flowing
water, pigment that suspends and settles, edge-darkening blooms as washes
dry, and color mixing that genuinely follows pigment absorption (so opposing
colors go muddy, just like on paper).

## Build

Requires SDL2 and SDL2_ttf packages.

**Debian/Ubuntu:**
```
sudo apt-get install libsdl2-dev libsdl2-ttf-dev
gcc watercolor.c -O2 -o watercolor $(pkg-config --cflags --libs sdl2 SDL2_ttf) -lm
./watercolor
```

**macOS (Homebrew):**
```
brew install sdl2 sdl2_ttf
gcc watercolor.c -O2 -o watercolor $(pkg-config --cflags --libs sdl2 SDL2_ttf) -lm
./watercolor
```

**Windows:** easiest via MSYS2 — install the `mingw-w64-x86_64-SDL2` and
`mingw-w64-x86_64-SDL2_ttf` packages, then build with the same gcc command
inside the MSYS2 shell.

The program looks for a DejaVu/Liberation TrueType font on common system
paths to draw the UI panel's text labels. If none is found it still runs, just without text labels.

## What's actually simulated

Each layer is a real (if simplified) fluid + pigment grid, not a bitmap with
a brush stamped onto it:

- **Water flows.** Wet areas spread across the paper based on local water
  height, the same way a wash pools and creeps on real paper.
- **Pigment suspends, then fixes.** Pigment riding in the water moves with
  it; as water evaporates, some of that pigment permanently stains the paper while the rest stays mobile.
- **Color mixing uses pigment absorption, not RGB blending.** Layer colors
  glaze on top of one another the way transparent paint does, so e.g. a near
  complementary pair turns appropriately muddy/dark rather than just
  averaging.
- **Paper texture matters.** A procedurally generated fibrous paper texture
  biases how water flows and where pigment settles (granulation), and shows
  through every wash.
- **Brushes run dry.** Pulling a long stroke depletes the brush's water
  load, so strokes naturally fade out as you drag, matching how a loaded
  watercolor brush behaves.

## Controls

- **Left-click + drag** on the canvas (left ~⅔ of the window) to paint.
- **Right panel:** click a pigment swatch to pick a color, or drag the R/G/B
  sliders for a custom mix. Adjust brush **size** and **wetness** with their
  sliders.
- **Layers:** "+ Add" adds a layer, click a layer row to make it active,
  click the dot to toggle visibility, "Delete" removes the active layer
  (at least one always remains).
- **New Paper** clears everything and regenerates a fresh paper texture.
  **Clear Layer** wipes just the active layer. **Save Painting** writes a
  timestamped `.bmp` of the composited painting to the current directory.

**Keyboard shortcuts:** `[` / `]` brush size down/up (mouse wheel also
works), `1`–`6` select layer by number, `v` toggle active layer visibility,
`c` clear active layer, `n` new paper, `s` save, `Esc` quit.

## Notes

This deliberately leaves out most traditional digital painting tooling
(undo history, selections, blend modes, opacity sliders, eraser, etc.). The
goal is an emulation of sitting down with real watercolors, where the medium
itself does a lot of the work and "happy accidents" from the physics are the
point.
