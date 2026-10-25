#!/bin/bash
# ==========================================================
# QC combinado: SeqKit + BUSCO
# Version adaptada para carpeta PLANA (todos los .fna.gz
# directamente en el directorio actual, sin subcarpetas)
#
# IMPORTANTE: este script debe correrse desde una carpeta que
# sea HIJA DIRECTA de ~/project (ej. ~/project/00_genomas o
# ~/project/prueba_qc), porque guarda los resultados en
# ../01_qc/ (un nivel arriba, en la carpeta de QC del proyecto)
#
# Requiere (instalar una sola vez):
#   conda install -c bioconda seqkit
#   conda create -n busco_env -c bioconda -c conda-forge busco=5.7.1
#   conda activate busco_env
#   busco --download eudicots_odb10
#   conda deactivate
#
# Uso:
#   cd ~/project/00_genomas   (o ~/project/prueba_qc)
#   bash qc_seqkit_busco_flat.sh
# ==========================================================

set -uo pipefail

LINEAGE="eudicots_odb10"      # cambia si alguna especie no es eudicotiledonea
THREADS=32                     # ajusta segun 'nproc'
BUSCO_ENV="busco_env"

QC_DIR="../01_qc"
BUSCO_RESULTS_DIR="$QC_DIR/busco"
mkdir -p "$BUSCO_RESULTS_DIR"

SUMMARY="$QC_DIR/qc_summary.csv"
echo "assembly_file,num_seqs,sum_len_bp,N50,GC_percent,busco_complete_pct,busco_single_pct,busco_duplicated_pct,busco_fragmented_pct,busco_missing_pct" > "$SUMMARY"

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

shopt -s nullglob
FILES=(*_genomic.fna.gz)

if [ ${#FILES[@]} -eq 0 ]; then
    echo "ERROR: no encuentro archivos *_genomic.fna.gz en esta carpeta."
    exit 1
fi

if [ ! -d "$QC_DIR" ]; then
    echo "ERROR: no encuentro la carpeta $QC_DIR"
    echo "Verifica que estas corriendo el script desde una carpeta"
    echo "hija directa de ~/project (ej: 00_genomas o prueba_qc)."
    exit 1
fi

echo "Se encontraron ${#FILES[@]} genomas."
echo "Los resultados se guardaran en: $QC_DIR"

get_col() {
    local header="$1"
    local data="$2"
    local colname="$3"
    awk -F'\t' -v col="$colname" '
        NR==1 { for (i=1; i<=NF; i++) if ($i==col) idx=i }
        NR==2 { print $idx }
    ' <<< "$header"$'\n'"$data"
}

for FNA_GZ in "${FILES[@]}"; do
    BASENAME="${FNA_GZ%.fna.gz}"
    echo "=============================================="
    echo ">> Procesando $BASENAME"
    echo "=============================================="

    FNA="${BASENAME}.fna"
    if [ ! -f "$FNA" ]; then
        gunzip -k "$FNA_GZ"
    fi

    echo "  -> SeqKit stats..."
    STATS_FULL=$(seqkit stats -a -T "$FNA")
    HEADER=$(echo "$STATS_FULL" | head -n1)
    DATA=$(echo "$STATS_FULL" | tail -n1)

    NUM_SEQS=$(get_col "$HEADER" "$DATA" "num_seqs")
    SUM_LEN=$(get_col "$HEADER" "$DATA" "sum_len")
    N50=$(get_col "$HEADER" "$DATA" "N50")
    GC=$(get_col "$HEADER" "$DATA" "GC(%)")

    echo "  -> BUSCO (puede tardar)..."
    conda activate "$BUSCO_ENV"

    BUSCO_OUT="busco_${BASENAME}"
    rm -rf "$BUSCO_RESULTS_DIR/$BUSCO_OUT"

    busco -i "$FNA" \
          -o "$BUSCO_OUT" \
          -l "$LINEAGE" \
          -m genome \
          -c "$THREADS" \
          --out_path "$BUSCO_RESULTS_DIR" \
          -f > "$BUSCO_RESULTS_DIR/${BASENAME}_busco_log.txt" 2>&1

    conda deactivate

    SUMMARY_FILE=$(find "$BUSCO_RESULTS_DIR/$BUSCO_OUT" -name "short_summary*.txt" 2>/dev/null | head -n1)

    if [ -n "$SUMMARY_FILE" ]; then
        C_PCT=$(grep -oP 'C:\K[0-9.]+(?=%)' "$SUMMARY_FILE")
        S_PCT=$(grep -oP 'S:\K[0-9.]+(?=%)' "$SUMMARY_FILE")
        D_PCT=$(grep -oP 'D:\K[0-9.]+(?=%)' "$SUMMARY_FILE")
        F_PCT=$(grep -oP 'F:\K[0-9.]+(?=%)' "$SUMMARY_FILE")
        M_PCT=$(grep -oP 'M:\K[0-9.]+(?=%)' "$SUMMARY_FILE")
    else
        echo "  [!] BUSCO no genero resumen para $BASENAME. Revisa $BUSCO_RESULTS_DIR/${BASENAME}_busco_log.txt"
        C_PCT="NA"; S_PCT="NA"; D_PCT="NA"; F_PCT="NA"; M_PCT="NA"
    fi

    echo "$BASENAME,$NUM_SEQS,$SUM_LEN,$N50,$GC,$C_PCT,$S_PCT,$D_PCT,$F_PCT,$M_PCT" >> "$SUMMARY"
    echo "  -> Listo: $BASENAME (BUSCO completo: ${C_PCT}%)"
    echo ""
done

echo "=============================================="
echo "QC terminado. Resumen en: $SUMMARY"
echo "Logs y resultados detallados de BUSCO en: $BUSCO_RESULTS_DIR/"
echo "=============================================="
column -s',' -t "$SUMMARY"
