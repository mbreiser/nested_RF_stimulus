#!/usr/bin/env python3
"""
Generate bar flash stimulus position cartoons.

Creates a row of 11 schematics, each showing a 30×30 arena box with
a 4-unit-wide vertical bar at one of 11 positions (stepping by 2 units).
The active bar is filled solid; the other 10 are shown as faint outlines.

Two color variants:
  1. Bright yellow-green bar on dim green background  (LED ON)
  2. Black bar on dim green background                (LED OFF / dark bar)

Outputs:
  stim_cartoons_bright.pdf / .svg / .png
  stim_cartoons_dark.pdf   / .svg / .png

Usage:
  /usr/bin/python3 scripts/generate_stim_cartoons.py
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path

# --- Geometry (internal units) ---
BOX_W = 30          # arena crop region width
BOX_H = 30          # arena crop region height
BAR_W = 4           # bar width
N_POS = 11          # number of positions
STEP  = 2           # step between bar centers

# Bar centers: position 6 (middle) at box center (x=15)
center_pos6 = BOX_W / 2                        # = 15
bar_centers = [center_pos6 + (i - 5) * STEP for i in range(N_POS)]
# → [5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25]

# --- Layout ---
GAP = 8             # gap between boxes (internal units)
MARGIN_X = 4        # horizontal margin
MARGIN_TOP = 4      # top margin
MARGIN_BOT = 7      # bottom margin (room for labels)
total_w = N_POS * BOX_W + (N_POS - 1) * GAP + 2 * MARGIN_X
total_h = BOX_H + MARGIN_TOP + MARGIN_BOT

# Scale to figure inches (aim for ~14 inches wide for clarity)
SCALE = 14.0 / total_w   # inches per internal unit

# --- Color palettes ---
PALETTES = {
    'bright': {
        'name': 'bright',
        'bg':             '#0a3a0a',     # dim green background (outer)
        'box_fill':       '#0d4a0d',     # slightly lighter green inside box
        'box_outline':    '#cccccc',     # white/light box outline
        'bar_fill':       '#b8e600',     # bright yellow-green (LED ON)
        'bar_outline':    '#cccccc',     # white outline on active bar
        'ghost_outline':  '#cccccc',     # white outlines for inactive
        'ghost_style':    '-',           # solid lines
        'ghost_lw':       0.4,
        'label_color':    '#bbbbbb',     # position labels
    },
    'dark': {
        'name': 'dark',
        'bg':             '#0a3a0a',     # dim green background (outer)
        'box_fill':       '#0d4a0d',     # slightly lighter green inside box
        'box_outline':    '#cccccc',     # white/light box outline
        'bar_fill':       '#050505',     # near-black bar
        'bar_outline':    '#cccccc',     # white outline on active bar
        'ghost_outline':  '#cccccc',     # white outlines for inactive
        'ghost_style':    '-',           # solid lines
        'ghost_lw':       0.4,
        'label_color':    '#bbbbbb',     # position labels
    },
}

output_dir = Path('/Users/reiserm/Documents/GitHub/nested_RF_stimulus/figures/stim_cartoons')
output_dir.mkdir(parents=True, exist_ok=True)


def draw_row(palette, ax):
    """Draw a row of 11 bar-flash cartoons on the given axes."""
    p = PALETTES[palette]

    # Background
    ax.set_facecolor(p['bg'])

    for pos_idx in range(N_POS):
        # Box origin (lower-left)
        box_x = MARGIN_X + pos_idx * (BOX_W + GAP)
        box_y = MARGIN_BOT

        # --- Box fill (slightly lighter green "screen") ---
        fill_rect = patches.Rectangle(
            (box_x, box_y), BOX_W, BOX_H,
            linewidth=0, facecolor=p['box_fill'], zorder=1
        )
        ax.add_patch(fill_rect)

        # --- Ghost outlines for all 11 bar positions ---
        for gi in range(N_POS):
            cx = bar_centers[gi]
            # Bar left/right edges, clipped to box
            left  = max(0, cx - BAR_W / 2)
            right = min(BOX_W, cx + BAR_W / 2)
            bx = box_x + left
            bw = right - left
            if bw <= 0:
                continue
            ghost_rect = patches.Rectangle(
                (bx, box_y), bw, BOX_H,
                linewidth=p['ghost_lw'], edgecolor=p['ghost_outline'],
                facecolor='none', linestyle=p['ghost_style'],
                zorder=2
            )
            ax.add_patch(ghost_rect)

        # --- Active bar (filled, clipped to box) ---
        cx = bar_centers[pos_idx]
        left  = max(0, cx - BAR_W / 2)
        right = min(BOX_W, cx + BAR_W / 2)
        active_x = box_x + left
        active_w = right - left
        active_rect = patches.Rectangle(
            (active_x, box_y), active_w, BOX_H,
            linewidth=0.8, edgecolor=p['bar_outline'],
            facecolor=p['bar_fill'], zorder=4
        )
        ax.add_patch(active_rect)

        # --- Box outline (on top) ---
        box_rect = patches.Rectangle(
            (box_x, box_y), BOX_W, BOX_H,
            linewidth=1.4, edgecolor=p['box_outline'],
            facecolor='none', zorder=5
        )
        ax.add_patch(box_rect)

        # --- Position label below box ---
        ax.text(
            box_x + BOX_W / 2, box_y - 2.0,
            str(pos_idx + 1),
            ha='center', va='top', fontsize=7, fontweight='normal',
            color=p['label_color'], fontfamily='sans-serif'
        )

    # Axis setup
    ax.set_xlim(0, total_w)
    ax.set_ylim(0, total_h)
    ax.set_aspect('equal')
    ax.axis('off')


def draw_single(palette, ax, pos_idx):
    """Draw a single bar-flash cartoon (one position) on the given axes.
    Minimal version: box outline + filled active bar only (no ghost lines).
    """
    p = PALETTES[palette]

    box_x, box_y = 0, 0

    # --- Box fill ---
    fill_rect = patches.Rectangle(
        (box_x, box_y), BOX_W, BOX_H,
        linewidth=0, facecolor=p['box_fill'], zorder=1
    )
    ax.add_patch(fill_rect)

    # --- Active bar only ---
    cx = bar_centers[pos_idx]
    left  = max(0, cx - BAR_W / 2)
    right = min(BOX_W, cx + BAR_W / 2)
    active_rect = patches.Rectangle(
        (box_x + left, box_y), right - left, BOX_H,
        linewidth=0, facecolor=p['bar_fill'], zorder=4
    )
    ax.add_patch(active_rect)

    # --- Box outline ---
    box_rect = patches.Rectangle(
        (box_x, box_y), BOX_W, BOX_H,
        linewidth=1.4, edgecolor=p['box_outline'],
        facecolor='none', zorder=5
    )
    ax.add_patch(box_rect)

    ax.set_xlim(-0.5, BOX_W + 0.5)
    ax.set_ylim(-0.5, BOX_H + 0.5)
    ax.set_aspect('equal')
    ax.axis('off')


def generate_variant(palette_name):
    """Generate combined row (PDF/SVG/PNG) and individual PNGs."""
    p = PALETTES[palette_name]

    # --- Combined row ---
    fig_w = total_w * SCALE
    fig_h = total_h * SCALE
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h))
    fig.patch.set_facecolor(p['bg'])

    draw_row(palette_name, ax)

    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)

    stem = f'stim_cartoons_{palette_name}'
    for ext in ['pdf', 'svg', 'png']:
        out_path = output_dir / f'{stem}.{ext}'
        dpi = 300 if ext == 'png' else None
        fig.savefig(
            out_path, format=ext, facecolor=fig.get_facecolor(),
            edgecolor='none', bbox_inches='tight', pad_inches=0.05,
            dpi=dpi, transparent=False
        )
        print(f'  Saved: {out_path}')
    plt.close(fig)

    # --- Individual cartoons (transparent background) ---
    indiv_dir = output_dir / palette_name
    indiv_dir.mkdir(parents=True, exist_ok=True)

    SINGLE_INCHES = 2.0   # each square is 2×2 inches
    for pos_idx in range(N_POS):
        fig_s, ax_s = plt.subplots(1, 1, figsize=(SINGLE_INCHES, SINGLE_INCHES))
        fig_s.patch.set_alpha(0.0)   # transparent figure background

        draw_single(palette_name, ax_s, pos_idx)

        for ext in ['png', 'pdf']:
            fname = f'pos{pos_idx + 1:02d}_{palette_name}.{ext}'
            out_path = indiv_dir / fname
            fig_s.savefig(
                out_path, format=ext,
                facecolor='none', edgecolor='none',
                bbox_inches='tight', pad_inches=0.02,
                dpi=300, transparent=True
            )
        plt.close(fig_s)

    print(f'  Saved: {indiv_dir}/pos01–pos11_{palette_name}.png + .pdf (transparent)')


if __name__ == '__main__':
    print('Generating bar flash stimulus cartoons...\n')
    print(f'  Box: {BOX_W}x{BOX_H} units')
    print(f'  Bar: {BAR_W} units wide, full height')
    print(f'  Bar centers: {[f"{c:.0f}" for c in bar_centers]}')
    print(f'  Bar left edges: {[f"{c - BAR_W/2:.0f}" for c in bar_centers]}')
    print(f'  Output: {output_dir}\n')

    for variant in ['bright', 'dark']:
        print(f'--- {variant} variant ---')
        generate_variant(variant)
        print()

    print('Done.')
