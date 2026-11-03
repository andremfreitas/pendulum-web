"""Render the link preview image and the icons.

Draws the frame chosen by snapshot.mjs the way app.js draws the page,
with the colours read from style.css, and writes

    og-dark.png, og-light.png   1200 x 630 link previews
    favicon.svg                 browser tab icon (follows the tab theme)
    apple-touch-icon.png        180 x 180 home-screen icon

into the website folder.

    python3.12 tools/og.py
"""

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

TOOLS = Path(__file__).resolve().parent
SITE = TOOLS.parent
FONTS = TOOLS / "fonts"

W, H = 1200, 630
DPI = 100
PT = 72 / DPI   # points per pixel
DPR = 2.0       # line widths as the page draws them on a 2x screen
CHAIN = 1.3     # the chain a little heavier than on the page, for thumbnails
TRAIL = 0.45    # fraction of the page's trail to show

URL = "andremfreitas.github.io/pendulum-web"
TITLE = "Pendulum Cascade"
SUBTITLE = "Energy cascading down a forced\nN-link pendulum, live in your browser"

REGULAR = font_manager.FontProperties(fname=FONTS / "IBMPlexSans-Regular.ttf")
SEMIBOLD = font_manager.FontProperties(fname=FONTS / "IBMPlexSans-SemiBold.ttf")

# tick labels such as 10^-3 go through mathtext, set in Plex as well
for f in FONTS.glob("*.ttf"):
    font_manager.fontManager.addfont(f)
plt.rcParams.update({
    "mathtext.fontset": "custom",
    "mathtext.rm": "IBM Plex Sans",
    "mathtext.it": "IBM Plex Sans",
    "mathtext.bf": "IBM Plex Sans:semibold",
})


# ------------------------------------------------------------ |
# colours from style.css                                       |

def palettes():
    css = (SITE / "style.css").read_text()

    def tokens(block):
        return dict(re.findall(r"--([a-z-]+):\s*(#[0-9a-fA-F]{6})", block))

    light = tokens(re.search(r":root\s*\{(.*?)\}", css, re.S).group(1))
    dark = tokens(re.search(r':root\[data-theme="dark"\]\s*\{(.*?)\}', css, re.S).group(1))
    return {"light": light, "dark": dark}


# ------------------------------------------------------------ |
# link preview                                                 |

def chain_joints(snap, cx, cy, R):
    th = np.array(snap["theta"])
    d = np.array(snap["disp"])
    x = cx + R * np.concatenate([[0], np.cumsum(d * np.sin(th))])
    y = cy + R * np.concatenate([[0], np.cumsum(d * np.cos(th))])
    return x, y


def draw_chain(ax, snap, c, cx, cy, R):
    n = len(snap["theta"])

    # reach circle
    ax.add_patch(Circle((cx, cy), R, fill=False, ec=c["line"], lw=DPR * PT,
                        ls=(0, (2 * DPR, 4 * DPR))))

    # tip trail, fading with age
    tr = np.array(snap["trail"])
    pts = np.stack([cx + R * tr[:, 1], cy + R * tr[:, 2]], axis=1)
    segs = np.stack([pts[:-1], pts[1:]], axis=1)
    keep = tr[:, 0] <= TRAIL
    tr = tr[keep]
    pts = pts[keep]
    segs = np.stack([pts[:-1], pts[1:]], axis=1)
    age = tr[:-1, 0] / TRAIL
    rgba = np.array([to_rgba(c["accent"])] * len(segs))
    rgba[:, 3] = 0.55 * np.clip(1 - age, 0, 1) ** 2
    ax.add_collection(LineCollection(segs, colors=rgba, linewidths=1.5 * DPR * PT,
                                     capstyle="round"))

    # links and masses
    x, y = chain_joints(snap, cx, cy, R)
    for i in range(n):
        f = i / (n - 1)
        ax.plot(x[i:i + 2], y[i:i + 2], color=c["ink"], lw=(3.4 - 2.4 * f) * DPR * CHAIN * PT,
                solid_capstyle="round")
    for i in range(n):
        f = i / (n - 1)
        ax.add_patch(Circle((x[i + 1], y[i + 1]), (3.6 - 2.2 * f) * DPR * CHAIN, color=c["ink"],
                            lw=0))

    # pivot
    ax.add_patch(Circle((cx, cy), 4.5 * DPR * CHAIN, fc=c["bg"], ec=c["ink"],
                        lw=1.5 * DPR * CHAIN * PT, zorder=3))


def draw_spectrum(ax, snap, c, x0, y0, w, h):
    avg = np.log10(np.array(snap["avg"]))
    now = np.log10(np.array(snap["spec"]))
    n = len(avg)
    lo, hi = np.floor(avg.min() - 0.3), np.ceil(avg.max() + 0.3)

    pad_l = 44
    X = lambda k: x0 + pad_l + k / (n - 1) * (w - pad_l)  # noqa: E731
    Y = lambda v: y0 + (1 - (v - lo) / (hi - lo)) * h  # noqa: E731

    every = max(1, int(np.ceil((hi - lo + 1) / 4)))
    for d in range(int(lo), int(hi) + 1):
        if (d - int(lo)) % every:
            continue
        ax.plot([x0 + pad_l, x0 + w], [Y(d)] * 2, color=c["line"], lw=1.2 * PT)
        label = "1" if d == 0 else rf"$10^{{{d}}}$"
        ax.text(x0 + pad_l - 10, Y(d), label, color=c["muted"], fontproperties=REGULAR,
                fontsize=15 * PT, ha="right", va="center")

    k = np.arange(n)
    ok = (now > lo - 0.5) & (now < hi + 0.5)
    ax.scatter(X(k[ok]), Y(np.clip(now[ok], lo, hi)), s=(4 * PT * 2) ** 2, color=c["muted"],
               alpha=0.45, lw=0)
    ax.plot(X(k), Y(avg), color=c["accent"], lw=2.4 * PT)
    ax.scatter(X(k), Y(avg), s=(5.2 * PT * 2) ** 2, color=c["accent"], lw=0, zorder=3)

    ax.text(x0, y0 - 26, "Energy per link", color=c["ink"], fontproperties=SEMIBOLD,
            fontsize=17 * PT, ha="left", va="baseline")
    ax.text(x0 + pad_l, y0 + h + 26, "large links", color=c["muted"], fontproperties=REGULAR,
            fontsize=14 * PT, ha="left", va="center")
    ax.text(x0 + w, y0 + h + 26, "small links", color=c["muted"], fontproperties=REGULAR,
            fontsize=14 * PT, ha="right", va="center")


def preview(snap, c, path):
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    fig.patch.set_facecolor(c["bg"])
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.set_aspect("equal")
    ax.axis("off")

    # centre what is drawn (chain and trail), not the pivot
    R = 245
    x, y = chain_joints(snap, 0, 0, R)
    tr = np.array(snap["trail"])
    tr = tr[tr[:, 0] <= TRAIL]
    xs = np.concatenate([x, R * tr[:, 1]])
    ys = np.concatenate([y, R * tr[:, 2]])
    cx = 345 - (xs.min() + xs.max()) / 2
    cy = H / 2 - (ys.min() + ys.max()) / 2
    draw_chain(ax, snap, c, cx=cx, cy=cy, R=R)

    tx = 700
    ax.text(tx, 132, TITLE, color=c["ink"], fontproperties=SEMIBOLD, fontsize=50 * PT,
            ha="left", va="baseline")
    ax.text(tx, 176, SUBTITLE, color=c["muted"], fontproperties=REGULAR, fontsize=22 * PT,
            ha="left", va="top", linespacing=1.45)

    draw_spectrum(ax, snap, c, x0=tx, y0=330, w=430, h=170)

    ax.text(tx, H - 52, URL, color=c["accent"], fontproperties=REGULAR, fontsize=17 * PT,
            ha="left", va="baseline")

    fig.savefig(path, dpi=DPI, facecolor=c["bg"])
    plt.close(fig)


# ------------------------------------------------------------ |
# icons: three tapered links bending away from a pivot         |

ICON_LINKS = [(11.0, 0.45, 2.6), (8.0, 1.35, 2.0), (6.0, 2.3, 1.5)]  # length, angle, width
ICON_DOTS = (2.0, 1.6, 1.3)
ICON_PIVOT_R, ICON_PIVOT_W = 2.2, 1.2


def icon_points():
    # chain from the origin, then centred in the 32 x 32 box
    x, y = 0.0, 0.0
    pts = [(x, y)]
    for L, a, _ in ICON_LINKS:
        x, y = x + L * np.sin(a), y + L * np.cos(a)
        pts.append((x, y))
    pts = np.array(pts)
    lo = pts.min(axis=0) - ICON_PIVOT_R
    hi = pts.max(axis=0) + ICON_PIVOT_R
    return [tuple(q) for q in pts + 16 - (lo + hi) / 2]


def favicon(pal, path):
    pts = icon_points()
    # the pivot ring is open, so start the first link at its edge
    a0 = ICON_LINKS[0][1]
    start = (pts[0][0] + ICON_PIVOT_R * np.sin(a0), pts[0][1] + ICON_PIVOT_R * np.cos(a0))
    pts = [start] + pts[1:]
    lines = "\n".join(
        f'  <line x1="{pts[i][0]:.2f}" y1="{pts[i][1]:.2f}" x2="{pts[i + 1][0]:.2f}" '
        f'y2="{pts[i + 1][1]:.2f}" stroke-width="{w}"/>'
        for i, (_, _, w) in enumerate(ICON_LINKS))
    dots = "\n".join(f'  <circle cx="{x:.2f}" cy="{y:.2f}" r="{r}"/>'
                     for (x, y), r in zip(pts[1:], ICON_DOTS))
    px, py = icon_points()[0]
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 32 32">
<style>
  @media (prefers-color-scheme: dark) {{ g {{ stroke: {pal['dark']['accent']}; fill: {pal['dark']['accent']}; }} }}
</style>
<g stroke="{pal['light']['accent']}" fill="{pal['light']['accent']}" stroke-linecap="round">
{lines}
{dots}
  <circle cx="{px:.2f}" cy="{py:.2f}" r="{ICON_PIVOT_R}" fill="none" stroke-width="{ICON_PIVOT_W}"/>
</g>
</svg>
"""
    path.write_text(svg)


def touch_icon(c, path):
    S = 180
    fig = plt.figure(figsize=(S / DPI, S / DPI), dpi=DPI)
    fig.patch.set_facecolor(c["bg"])
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 32)
    ax.set_ylim(32, 0)
    ax.axis("off")
    k = S / 32 * PT   # icon units to points
    pts = icon_points()
    for i, (_, _, w) in enumerate(ICON_LINKS):
        ax.plot([pts[i][0], pts[i + 1][0]], [pts[i][1], pts[i + 1][1]], color=c["accent"],
                lw=w * k, solid_capstyle="round")
    for (x, y), r in zip(pts[1:], ICON_DOTS):
        ax.add_patch(Circle((x, y), r, color=c["accent"], lw=0))
    ax.add_patch(Circle(pts[0], ICON_PIVOT_R, fc=c["bg"], ec=c["accent"], lw=ICON_PIVOT_W * k,
                        zorder=3))
    fig.savefig(path, dpi=DPI, facecolor=c["bg"])
    plt.close(fig)


def main():
    snap = json.loads((TOOLS / "snapshot.json").read_text())
    pal = palettes()
    for theme in ("dark", "light"):
        preview(snap, pal[theme], SITE / f"og-{theme}.png")
    favicon(pal, SITE / "favicon.svg")
    touch_icon(pal["dark"], SITE / "apple-touch-icon.png")


if __name__ == "__main__":
    main()
