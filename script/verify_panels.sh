#!/usr/bin/env bash
# Pixel-compare every regenerated panel against its published counterpart in
# figures/Previous/bits/, using the old->new panel mapping (Figure 1 was split
# into 1+2 and Figures 2/3/4 shifted up to 3/4/5, so same-letter filenames are
# NOT counterparts -- see script/README.md).
#
#     bash script/verify_panels.sh            # all panels
#     bash script/verify_panels.sh 5d 6d 3M   # just these
#
# Prints one line per panel: RMSE (0 = pixel-identical) and a verdict. Panels
# that differ on purpose are listed in script/README.md's "Verification status"
# -- check any *new* non-zero entry against that list before shipping.
#
# Do NOT verify by comparing file sizes: a panel can plot the wrong variable at
# almost exactly the same size (this is how a Level-1/Level-2 mix-up in 5d, and
# a permuted 6d-6f gallery, both survived an earlier review). md5 never matches
# either, because each PDF embeds a fresh /CreationDate.
set -uo pipefail

NEW_ROOT="figures/final-selected"
OLD_ROOT="figures/Previous/bits"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
DENSITY=${DENSITY:-110}

# new_panel <TAB> published_panel   (paths relative to the two roots above)
#
# Pairs are matched by CONTENT, not by letter. That matters for Figure 3's
# J and L: we letter the gdT/CD8aa/DN loading panels GP3 = J, GP29 = K,
# GP22 = L, where the published 2J/2K/2L were GP22/GP29/GP3. So 3J's
# counterpart is 2L and 3L's is 2J -- comparing them letter-to-letter would
# compare different GPs and look like a huge regression. It also matters for
# Figures 6 / S3 / S6 after the 2026-07-28 re-lettering: Figure 6 was reordered
# (its KLRG1 and CD69 panels moved to the front, the protein-program heatmap and
# the two CD69 mean-activity heatmaps moved into Figure S6), and Figure S3's
# panels each dropped one letter when its s3a/s3b were merged into a single s3a.
# And for Figure 4 after 2026-07-29: its TF-GP network was dropped, so the two
# heatmaps that followed moved up a letter and the published 3e is now unpaired.
MAP=$(
  cat <<'EOF'
Figure 1/1A	Figure 1/1A
Figure 1/1B	Figure 1/1B
Figure 1/1D	Figure 1/1C
Figure 2/2A	Figure 1/1D
Figure 2/2B	Figure 1/1E
Figure 2/2C	Figure 1/1F
Figure 2/2D	Figure 1/1G
Figure 2/2E	Figure 1/1H
Figure 2/2F	Figure 1/1I
Figure 3/3A	Figure 2/2A
Figure 3/3B	Figure 2/2B
Figure 3/3C	Figure 2/2C
Figure 3/3D	Figure 2/2D
Figure 3/3E	Figure 2/2E
Figure 3/3F	Figure 2/2F
Figure 3/3G	Figure 2/2G
Figure 3/3H	Figure 2/2H
Figure 3/3I	Figure 2/2I
Figure 3/3J	Figure 2/2L
Figure 3/3K	Figure 2/2K
Figure 3/3L	Figure 2/2J
Figure 3/3M	Figure 2/2M
Figure 4/4a	Figure 3/3c
Figure 4/4b	Figure 3/3d
Figure 4/4c	Figure 3/3f
Figure 4/4d	Figure 3/3g
Figure 5/5a	Figure 4/4a
Figure 5/5b	Figure 4/4b
Figure 5/5c	Figure 4/4c
Figure 5/5d	Figure 4/4d
Figure 5/5e	Figure 4/4e
Figure 6/6a	Figure 6/6a
Figure 6/6b	Figure 6/6g
Figure 6/6c	Figure 6/6h
Figure 6/6d	Figure 6/6i
Figure 6/6e	Figure 6/6c
Figure 6/6f	Figure 6/6d
Figure 6/6g	Figure 6/6e
Figure 6/6h	Figure 6/6f
Figure S6/s6a	Figure 6/6b
Figure S6/s6b	Figure 6/6j
Figure S6/s6c	Figure 6/6k
Figure S1/S1A	Figure S1/S1A
Figure S1/S1B	Figure S1/S1B
Figure S2/S2A	Figure S2/S2A
Figure S2/S2B	Figure S2/S2B
Figure S2/S2C	Figure S2/S2C
Figure S3/s3b	Figure S3/s3c
Figure S3/s3c	Figure S3/s3d
Figure S3/s3d	Figure S3/s3e
Figure S3/s3e	Figure S3/s3f
Figure S3/s3f	Figure S3/s3g
Figure S3/s3g	Figure S3/s3h
EOF
)
# Published panels with no counterpart HERE any more: Figure 3/3e, the TF-GP
# network dropped from Figure 4 on 2026-07-29, and the published Figure S1/S1C
# (mean vs. variance of per-IGT mean loading) and S1D (top-variance heatmap),
# both replaced on 2026-07-30 by a different analysis -- see the header of
# script/FigureS1.R.
# Panels with no published counterpart: Figure 1/1C, Figure S1/S1C, S1D, S1E and
# S1F (S1C/S1D are the 2026-07-30 replacements above and are new analyses, not
# reproductions; S1E is new; S1F is the former S1E, which never had one),
# Figure S4/*, Figure 6/6i (GP77) and 6j (GP8), Figure S6/s6d-s6g (their GPs
# were only ever drawn inside the retired s6-1/s6-2 gallery pages, not as
# standalone panels), and Figure 7/7B (ours since 2026-07-28 -- it is the former
# Extended Data Figure 5, assembled into one panel, and it replaced a different
# collaborator panel, so the published 7B is NOT its counterpart).
# The rest of Figure 7 and all of S7 are out of scope (straight copies).

render() { # $1 = pdf, $2 = out png -- normalise onto one square canvas so a
           # changed aspect ratio shows up as a difference rather than an error
  magick -density "$DENSITY" "${1}[0]" -background white -alpha remove \
    -resize 700x700 -background white -gravity center -extent 700x700 \
    -colorspace Gray "$2" 2>/dev/null
}

printf '%-16s %-16s %-12s %s\n' PANEL PUBLISHED RMSE VERDICT
printf '%s\n' "-------------------------------------------------------------"
n=0
nz=0
while IFS=$'\t' read -r new old; do
  [ -n "${new:-}" ] || continue
  if [ $# -gt 0 ]; then
    want=0
    for a in "$@"; do [ "$(basename "$new")" = "$a" ] && want=1; done
    [ "$want" = 1 ] || continue
  fi
  np="$NEW_ROOT/$new.pdf"
  op="$OLD_ROOT/$old.pdf"
  if [ ! -f "$np" ] || [ ! -f "$op" ]; then
    printf '%-16s %-16s %-12s %s\n' "$new" "$old" "-" "MISSING FILE"
    continue
  fi
  render "$np" "$TMP/a.png"
  render "$op" "$TMP/b.png"
  r=$(magick compare -metric RMSE "$TMP/a.png" "$TMP/b.png" null: 2>&1 >/dev/null |
    sed -n 's/.*(\(.*\)).*/\1/p')
  v=$(awk -v x="${r:-1}" 'BEGIN {
        if (x + 0 == 0)      print "identical"
        else if (x + 0 < 0.02) print "~ label/AA jitter"
        else if (x + 0 < 0.06) print "differs -- check README"
        else                 print "*** DIFFERS ***"
      }')
  printf '%-16s %-16s %-12s %s\n' "$(basename "$new")" "$(basename "$old")" "$r" "$v"
  n=$((n + 1))
  [ "${r:-1}" != "0" ] && nz=$((nz + 1))
done <<<"$MAP"
printf '%s\n' "-------------------------------------------------------------"
echo "$n panels compared, $nz not pixel-identical"
echo "Expected non-identical panels are documented in script/README.md."
