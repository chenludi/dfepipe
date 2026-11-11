#!/bin/bash
# 2026-10-27  Two discrete-DFE figures, drawn with the manuscript's own plotter.
#
# Calls the published script unchanged:
#   $OLD_SCRIPT/plot_DFE_mutload2_fromcsv.py -format dfe -listcsv ... -colors ...
# It reads each csv with index_col=0, so one row per csv, the index being the
# legend label and the columns the six |s| bins. One colour per row, in order.
#
# figS2   Fig. S2 (old numbering figS1_dfe_enhancer_promoterall) with two bars
#         inserted after EnhA_ELS, from the old-ENCODE ELS run *_261025:
#         EnhA, ELS, EnhA_ELS, [ELS_x_Prom, ELS_not_Prom], dELS, pELS,
#         PromF, PLS, PromF_PLS, nonsyn
#         The nine original bars are the published csvs, reused untouched.
#
# oldnew  old ENCODE vs ENCODE v4, for ELS and its two Prom subsets, plus PLS:
#         ELS, ELS_x_Prom, ELS_not_Prom (run *_261025)
#         ELSv4, ELSv4_x_Prom, ELSv4_not_Prom, PLSv4 (run *_261026)
#         PLS (old) is the published csv: that bed was not re-run, so it is the
#         one bar here whose provenance differs. Everything else is a fresh run.
#
# Bars are the best-fit discrete DFE (the delta_50_* columns), not gamma by
# assumption; dfe_csv_from_bestpdf.py prints which pdf each bar came from.
#
# Output. On the cluster: $OUT/figs/261027/. The copies pulled down on
# 2026-10-27 were moved on 2026-11-10 to
#   DFE_non-coding/a_results_weekly/mv_261027_dfepipe_els_v4_phylop/figs/
# They were in dfepipe/exchange/ before that; that directory no longer exists.
# $OUT/figs/261027/residuals/ was never copied down.
#
#   bash records/record_261027_dfe_figs.sh            both
#   FIGS=figS2 bash records/record_261027_dfe_figs.sh
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
source config/profiles/hoffman_replicate.conf
source config/pipeline.conf
source config/env.sh
activate_conda            # the cluster's bare `python` is Python 2

FIGS="${FIGS:-figS2 oldnew}"
NC="$(dirname "$REFS_DIR")"
OLD_SCRIPT="${OLD_SCRIPT:-$REFS_DIR/script}"
PUB="${PUB:-$REFS_DIR/plots_dfe_mutload/manuscript_0830/figS1_dfe_enhancer_promoterall}"
FIGDIR="${FIGDIR:-$OUT/figs/261027}"
B25="$OUT/summary/sum_261025/bestpdf_sumdfe.tsv"
B26="$OUT/summary/sum_261026/bestpdf_sumdfe.tsv"

for f in "$OLD_SCRIPT/plot_DFE_mutload2_fromcsv.py" "$PUB" "$B25" "$B26"; do
  [ -e "$f" ] || { echo "missing: $f" >&2; exit 1; }
done
mkdir -p "$FIGDIR"

run_plot() {        # <outname> <csv list> <colour list>
  local out="$1"; shift
  local csvs="$1"; shift
  local cols="$1"
  echo "  plotting $out"
  # shellcheck disable=SC2086
  python "$OLD_SCRIPT/plot_DFE_mutload2_fromcsv.py" \
    -format dfe -listcsv $csvs -colors $cols -o "$FIGDIR/$out" >/dev/null
  echo "  wrote $FIGDIR/$out"
  echo "        $FIGDIR/${out%.pdf}_plot.csv"
}

for fig in $FIGS; do
case "$fig" in

figS2)
  echo "== figS2: Fig. S2 with the two ELS/Prom bars inserted after EnhA_ELS"
  D="$FIGDIR/figS2_csv"; mkdir -p "$D"
  python records/dfe_csv_from_bestpdf.py -i "$B25" -o "$D" \
    --group exHET_noq.GRCh38-ELS_x_Prom_qui2   --label ELS_x_Prom \
    --group exHET_noq.GRCh38-ELS_not_Prom_qui2 --label ELS_not_Prom
  # the nine published bars, in their published order, with the two inserted 3rd
  csvs="$PUB/dfe_EnhA_qui2.chrs_1t22_gamma.csv
$PUB/dfe_exHET_noq.GRCh38-ELS_qui2.chrs_1t22_gamma.csv
$PUB/dfe_EnhA_ELS_qui2.chrs_1t22_gamma.csv
$D/dfe_ELS_x_Prom.csv
$D/dfe_ELS_not_Prom.csv
$PUB/dfe_exHET_noq.GRCh38-cCREs.dELS_qui2.chrs_1t22_gamma.csv
$PUB/dfe_exHET_noq.GRCh38-cCREs.pELS_qui2.chrs_1t22_gamma.csv
$PUB/dfe_PromF_qui2.chrs_1t22_gamma.csv
$PUB/dfe_exHET_noq.GRCh38-PLS_qui2.chrs_1t22_gamma.csv
$PUB/dfe_PromF_PLS_qui2.chrs_1t22_gamma.csv
$PUB/dfe_nonsyn_Kim.csv"
  # published colours, with two blues added for the inserted ELS subsets
  cols="Blue LightSteelBlue DarkBlue CornflowerBlue SteelBlue Turquoise DarkTurquoise Red LightSalmon DarkRed MediumPurple"
  run_plot "figS2_dfe_enh_prom_withELSprom.pdf" "$csvs" "$cols"
  ;;

oldnew)
  echo "== oldnew: old ENCODE vs v4, ELS subsets and PLS"
  D="$FIGDIR/oldnew_csv"; mkdir -p "$D"
  python records/dfe_csv_from_bestpdf.py -i "$B25" "$B26" -o "$D" \
    --group exHET_noq.GRCh38-ELS_qui2            --label ELS_old \
    --group exHET_noq.GRCh38-ELS_x_Prom_qui2     --label ELS_x_Prom_old \
    --group exHET_noq.GRCh38-ELS_not_Prom_qui2   --label ELS_notProm_old \
    --group exHET_noq.GRCh38-ELSv4_qui2          --label ELS_v4 \
    --group exHET_noq.GRCh38-ELSv4_x_Prom_qui2   --label ELS_x_Prom_v4 \
    --group exHET_noq.GRCh38-ELSv4_not_Prom_qui2 --label ELS_notProm_v4 \
    --group exHET_noq.GRCh38-PLSv4_qui2          --label PLS_v4
  csvs="$D/dfe_ELS_old.csv
$D/dfe_ELS_v4.csv
$D/dfe_ELS_x_Prom_old.csv
$D/dfe_ELS_x_Prom_v4.csv
$D/dfe_ELS_notProm_old.csv
$D/dfe_ELS_notProm_v4.csv
$PUB/dfe_exHET_noq.GRCh38-PLS_qui2.chrs_1t22_gamma.csv
$D/dfe_PLS_v4.csv"
  # old/new pairs share a hue, the v4 member lighter
  cols="DarkBlue CornflowerBlue DarkGreen MediumSeaGreen Sienna SandyBrown DarkRed LightSalmon"
  run_plot "dfe_oldENCODE_vs_v4.pdf" "$csvs" "$cols"
  ;;

*) echo "unknown figure key: $fig" >&2; exit 1 ;;
esac
echo
done
echo "figures in $FIGDIR"
