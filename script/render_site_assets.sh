#!/usr/bin/env bash
# Regenerate every PNG the workflowr site shows *from* the panel PDFs in
# figures/final-selected/, so the site can never drift from the published
# panel set.
#
# Run from the repo root, after any script/FigureN.R re-run and before
# workflowr::wflow_build():
#
#     bash script/render_site_assets.sh
#
# Why convert instead of re-plotting: the site PNGs used to be produced by
# their own ggsave() calls (see the retired render_new_panels.R). Any panel
# whose layout is not fully determined by the data -- ggrepel labels, jitter,
# graph layouts -- then rendered *differently* in the PNG than in the PDF, and
# the two silently diverged as scripts were re-run at different times.
# Converting the PDF guarantees the site shows exactly the published panel.
#
# Resolution: 72 dpi, i.e. 1 px per PDF point, matching the existing assets.
# Bump DENSITY to 144 for 2x (sharper on HiDPI screens, ~4x the file size).
set -euo pipefail

DENSITY=${DENSITY:-72}
SRC="figures/final-selected"
DST="analysis/assets"

# Panels deliberately NOT on the site:
#   Figure 1/1B, Figure 7/7a  -- hand-drawn schematics, no code, not shown
#   Figure 8/8A, 8C-8G        -- out of scope (Figma + Matplotlib, no source)
# Not handled here:
#   Figure 8/8B -- script/Figure8b_plot.py saves the published PDF and
#   analysis/assets/Figure8/8B.png from the SAME matplotlib figure in one call,
#   so unlike every other panel there is no second plotting path that could
#   drift from the first, and converting the PDF here would only lower the
#   preview's resolution.
skip() {
  case "$1" in
    "Figure 1/1B" | "Figure 7/7a") return 0 ;;
    "Figure 8/"*) return 0 ;;
    *) return 1 ;;
  esac
}

n=0
for pdf in "$SRC"/Figure*/*.pdf; do
  dir=$(basename "$(dirname "$pdf")")   # e.g. "Figure S1"
  base=$(basename "$pdf" .pdf)          # e.g. "S1A"
  skip "$dir/$base" && continue
  outdir="$DST/${dir// /}"              # "Figure S1" -> "FigureS1"
  mkdir -p "$outdir"
  out="$outdir/$base.png"
  # -strip / exclude-chunk: ImageMagick otherwise writes creation-time text and
  # tIME chunks, so re-running would rewrite all 60 PNGs with new bytes and show
  # up as a 60-file diff even when nothing changed.
  magick -density "$DENSITY" "${pdf}[0]" -background white -alpha remove -alpha off \
    -strip -define png:exclude-chunk=date,time "$out"
  printf '  %-46s -> %s\n' "$pdf" "$out"
  n=$((n + 1))
done
echo "converted $n panel PDFs at ${DENSITY} dpi"

# workflowr builds each page with rmarkdown::render() in its own session, which
# does NOT do rmarkdown site-resource copying, so docs/assets is never refreshed
# by wflow_build -- not even with republish = TRUE. Mirror it here, or the built
# site keeps serving whatever PNGs were there before.
if [ -d docs ]; then
  rsync -a --delete --exclude '.DS_Store' "$DST"/ docs/assets/
  echo "mirrored $DST -> docs/assets"
fi

# Fail loudly if the site references a PNG we did not produce.
echo "checking every include_graphics() reference resolves..."
missing=0
grep -rho 'assets/[A-Za-z0-9]*/[^"]*\.png' analysis/*.Rmd | sort -u | while read -r ref; do
  [ -f "analysis/$ref" ] || { echo "  MISSING: analysis/$ref"; missing=1; }
done
[ "$missing" -eq 0 ] && echo "all references resolve"
