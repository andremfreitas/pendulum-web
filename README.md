# Pendulum Cascade website

Static page that runs the N-link pendulum live in the browser. No build step and no dependencies: plain HTML, CSS and ES modules.

- `sim.js`: integrator, a port of `sim_code/src/f90`, using the same IMEX-SSP3(4,3,3) scheme and the same 1:2 length and mass ratio. The O(N²) sums are rewritten as O(N) prefix sums. After 2000 steps it agrees with the Fortran to about 1e-9 relative error.
- `app.js`: drawing, controls, flick gesture, energy spectrum.
- `index.html`, `style.css`: page and layout.

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
