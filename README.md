# Pendulum Cascade website

Static page that runs the N-link pendulum live in the browser. No build step and no dependencies: plain HTML, CSS and ES modules.

- `sim.js`: integrator, a port of `sim_code/src/f90`, using the same IMEX-SSP3(4,3,3) scheme and the same 1:2 length and mass ratio. The O(N²) sums are rewritten as O(N) prefix sums. After 2000 steps it agrees with the Fortran to about 1e-9 relative error.
- `app.js`: drawing, controls, flick gesture, energy spectrum.
- `index.html`, `style.css`: page and layout.
- `og-dark.png`, `og-light.png`: 1200 × 630 link previews; `index.html` points to the dark one.
- `favicon.svg`, `apple-touch-icon.png`: tab and home-screen icons.

## Link preview and icons

Both are generated from a real frame of the simulation:

```bash
node tools/snapshot.mjs
python3.12 tools/og.py
```

`snapshot.mjs` runs `sim.js` with the default parameters and a fixed seed, picks an extended, hooked frame, and writes `tools/snapshot.json`. `og.py` draws that frame as `app.js` does, with the colours read from `style.css` and the text in IBM Plex Sans (`tools/fonts/`, SIL Open Font License). To restyle the images without picking a new frame, run only `og.py`; it needs matplotlib and not Node.

## Run locally

ES modules do not load from `file://`, so serve the folder:

```bash
python3 -m http.server 8000 --directory website
```

Then open http://localhost:8000.

## Deploy on GitHub Pages

1. Create a public repo (for example `pendulum-web`) and copy the contents of this folder to its root.
2. On GitHub, under Settings → Pages, choose "Deploy from a branch", branch `main`, folder `/`.
3. The site goes live at `https://<user>.github.io/pendulum-web/`. Every push updates it.
