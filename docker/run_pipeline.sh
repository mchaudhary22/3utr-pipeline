#!/bin/bash
set -e

# ============================================================
# URSA Pipeline v13
# Changes from v12:
#   (1) StringTie -c lowered from 3 to 2
#       Allows assembly of transcripts with lower coverage,
#       recovering more lowly-expressed genes (rare cell types).
#       The Python --min-coverage 2.0 gate remains the operative
#       biological filter downstream.
#   (2) --min-coverage lowered from 3.0 to 2.0
#       Matches the new StringTie -c 2 threshold. Recovers
#       genuine extensions on lowly-expressed genes that were
#       silently dropped at 3.0x. Splice junction filter is
#       still the real biological quality gate.
#   (3) --max-extension raised from 2000 to 3000 bp
#       Plant UTR literature supports genuine 3' UTRs up to
#       2-3kb for specific gene families. Neighbor overlap check
#       and chromosome boundary cap still protect against false
#       extensions regardless of this parameter.
#   All other logic, filters, and steps unchanged from v12.
# ============================================================

# Input/Output/Reference paths
SAMPLE_DIR="/app/input"
OUTPUT_DIR="/app/output"
REFERENCE_DIR="/app/output/reference"

# Initialize variables
SRR_ID=""
SRR_LIST_FILE=""
EXPECTED_CELLS=""       # Empty = let Cell Ranger auto-detect
FORCE_CELLS=""          # If set, overrides Cell Ranger knee detection entirely
GENOME_NAME=""

# -----------------------------------
# Parse command-line arguments
# -----------------------------------
while getopts "s:l:c:f:g:h" opt; do
    case $opt in
        s) SRR_ID="$OPTARG" ;;
        l) SRR_LIST_FILE="$OPTARG" ;;
        c) EXPECTED_CELLS="$OPTARG" ;;
        f) FORCE_CELLS="$OPTARG" ;;
        g) GENOME_NAME="$OPTARG" ;;
        h)
            echo "Usage: $0 [-s SRR_ID | -l SRR_LIST_FILE] [-c EXPECTED_CELLS] [-f FORCE_CELLS] -g GENOME_NAME"
            echo "  -c  Hint Cell Ranger toward expected cell count (soft guidance, e.g. 10000)"
            echo "  -f  Force exact cell count, bypasses knee detection (e.g. 10000)"
            echo "  If neither -c nor -f is given, Cell Ranger auto-detects cell count"
            echo "  Input files: GENOME_NAME.fa and GENOME_NAME.gff3 under /app/input/"
            exit 0
            ;;
        *) echo "Invalid option"; exit 1 ;;
    esac
done

if [[ -z "$SRR_ID" && -z "$SRR_LIST_FILE" ]]; then
    echo "ERROR: Provide either -s (single SRR) or -l (list of SRRs)"
    exit 1
fi

if [[ -z "$GENOME_NAME" ]]; then
    echo "ERROR: Provide species name with -g"
    exit 1
fi

# -----------------------------------
# Locate genome and annotation files
# -----------------------------------
for ext in fa fasta fna; do
    for variant in "$GENOME_NAME" "${GENOME_NAME,,}" "${GENOME_NAME^}" "${GENOME_NAME^^}"; do
        if [[ -f "/app/input/${variant}.${ext}" ]]; then
            GENOME_FILE="/app/input/${variant}.${ext}"
            break 2
        fi
    done
done

for ext in gff3 gff; do
    for variant in "$GENOME_NAME" "${GENOME_NAME,,}" "${GENOME_NAME^}" "${GENOME_NAME^^}"; do
        if [[ -f "/app/input/${variant}.${ext}" ]]; then
            GFF3_FILE="/app/input/${variant}.${ext}"
            break 2
        fi
    done
done

if [[ ! -f "$GENOME_FILE" || ! -f "$GFF3_FILE" ]]; then
    echo "ERROR: Genome or GFF3 not found for '$GENOME_NAME'"
    ls -la /app/input/
    exit 1
fi

echo "Genome file: $GENOME_FILE"
echo "Annotation file: $GFF3_FILE"

# -----------------------------------
# Core processing function
# -----------------------------------
process_single_srr() {
    local current_srr="$1"
    local expected_cells="$2"
    local force_cells="$3"
    echo "========== Processing sample: $current_srr =========="

    local srr_path="$SAMPLE_DIR/$current_srr"
    if [[ ! -d "$srr_path" ]]; then
        echo "ERROR: SRR folder not found: $srr_path"
        exit 1
    fi

    mkdir -p "$OUTPUT_DIR/$current_srr" "$REFERENCE_DIR"
    cd "$OUTPUT_DIR/$current_srr"

    # Decompress FASTQ files
    echo "[INFO] Preparing FASTQ files..."
    for i in 1 2 3; do
        fq_gz="$srr_path/${current_srr}_${i}.fastq.gz"
        fq_raw="$srr_path/${current_srr}_${i}.fastq"
        out_file="${current_srr}_S1_L001_R${i}_001.fastq"

        if [[ -f "$fq_gz" ]]; then
            gunzip -c "$fq_gz" > "$out_file"
        elif [[ -f "$fq_raw" ]]; then
            cp "$fq_raw" "$out_file"
        fi
    done

    [[ -f "${current_srr}_S1_L001_R3_001.fastq" ]] && mv "${current_srr}_S1_L001_R3_001.fastq" "${current_srr}_S1_L001_I1_001.fastq"

    # Convert GFF3 → GTF
    cd "$REFERENCE_DIR"
    if [[ ! -f "${GENOME_NAME}.gtf" ]]; then
        echo "[INFO] Converting GFF3 to GTF..."
        gffread "$GFF3_FILE" -T -E -o "${GENOME_NAME}.gtf"

        # ============================================================
        # V10 FIX: Strip gene: prefix for consistent gene IDs
        # gffread keeps gene: prefix from GFF3 but StringTie does not,
        # causing the same gene to appear as both "gene:SolycXXX" and
        # "SolycXXX" — doubling genes and breaking overlap detection.
        # ============================================================
        echo "[INFO] Stripping gene: prefix from gene IDs..."
        sed -i 's/gene_id "gene:/gene_id "/g; s/transcript_id "transcript:/transcript_id "/g' "${GENOME_NAME}.gtf"
        echo "[INFO] Gene ID normalization completed!"

        # ============================================================
        # V6 NEW: Remove overlapping genes
        # ============================================================
        echo "[INFO] Removing same-strand overlapping genes..."
        python3 /app/remove_overlap.py \
            "${GENOME_NAME}.gtf" \
            "${GENOME_NAME}_filtered.gtf" \
            --report "${GENOME_NAME}_overlap_report.tsv" \
            --policy prefer_protein_coding

        # Replace original with filtered version
        mv "${GENOME_NAME}_filtered.gtf" "${GENOME_NAME}.gtf"
        echo "[INFO] Overlap removal completed!"
    fi

    # ============================================================
    # V12: Generate genome index (.fai) for chromosome boundary
    # checks in extend_utr_boundaries.py. samtools faidx creates
    # GENOME_NAME.fa.fai which maps chrom -> length.
    # This prevents extensions from exceeding contig boundaries
    # (confirmed root cause of Poplar chromosome overflow failure).
    # ============================================================
    if [[ ! -f "${GENOME_FILE}.fai" ]]; then
        echo "[INFO] Indexing genome FASTA for chromosome size lookup..."
        samtools faidx "$GENOME_FILE"
        echo "[INFO] Genome index created: ${GENOME_FILE}.fai"
    fi

    # Create base genome reference
    if [[ ! -d "${GENOME_NAME}_genome" ]]; then
        echo "[INFO] Creating base genome reference..."
        cellranger mkref --genome="${GENOME_NAME}_genome" --fasta="$GENOME_FILE" --genes="${GENOME_NAME}.gtf"
    fi

    # Build Cell Ranger cell-count argument
    local cellcount_arg=""
    if [[ -n "$FORCE_CELLS" ]]; then
        cellcount_arg="--force-cells=$FORCE_CELLS"
        echo "[INFO] Cell count mode: forced to $FORCE_CELLS cells"
    elif [[ -n "$expected_cells" ]]; then
        cellcount_arg="--expect-cells=$expected_cells"
        echo "[INFO] Cell count mode: expecting ~$expected_cells cells"
    else
        echo "[INFO] Cell count mode: auto-detect (no -c or -f given)"
    fi

    # Cell Ranger count (base genome)
    cd "$OUTPUT_DIR"
    local matrix_id="${current_srr}_matrix"
    if [[ ! -d "$matrix_id" ]]; then
        echo "[INFO] Running Cell Ranger count on base genome..."
        cellranger count --id="$matrix_id" \
            --transcriptome="$REFERENCE_DIR/${GENOME_NAME}_genome" \
            --fastqs="$OUTPUT_DIR/$current_srr" \
            --sample="$current_srr" \
            $cellcount_arg --create-bam=true
    fi

    # ============================================================
    # V13: StringTie with further relaxed parameters
    # -c 2   minimum 2x coverage — lowered from v12's -c 3
    #        Matches the new --min-coverage 2.0 threshold below.
    #        Recovers real extensions on lowly-expressed genes
    #        (rare cell types) that were dropped at -c 3.
    #        Splice junction filter remains the real quality gate.
    # -f 0.1 minimum 10% of major isoform (unchanged)
    # -g 50  merge gaps under 50bp (unchanged)
    # ============================================================
    cd "$OUTPUT_DIR"
    local bam_file="${matrix_id}/outs/possorted_genome_bam.bam"
    local stringtie_gtf="${current_srr}_stringtie.gtf"
    local formatted_gtf="${current_srr}_stringtie_v2.gtf"
    local filtered_gtf="${current_srr}_stringtie_filtered.gtf"
    local filter_report="${current_srr}_splice_filter_report.tsv"

    if [[ -f "$bam_file" && ! -f "$stringtie_gtf" ]]; then
        echo "[INFO] Running StringTie (v13: -c 2 -f 0.1 -g 50)..."
        stringtie "$bam_file" \
            -G "$REFERENCE_DIR/${GENOME_NAME}.gtf" \
            -c 2 -f 0.1 -g 50 \
            -o "$stringtie_gtf"
    fi

    # ============================================================
    # V10 FIX: Splice filter BEFORE ID remapping (fixes ? class bug)
    # In v9, update_ref_gene_format.py ran first and renamed STRG IDs
    # to Solyc IDs. The filter then looked up Solyc IDs in the annotated
    # GTF which only has STRG IDs — causing 13,836 transcripts to get
    # ? class and be kept by default as noise.
    # Fix: run gffcompare and filter on original STRG IDs first,
    # then remap IDs on the already-filtered output.
    # ============================================================
    local annotated_gtf="$OUTPUT_DIR/${current_srr}_stringtie_srr_annotation.annotated.gtf"
    local extended_ref_gtf="$OUTPUT_DIR/${current_srr}_reference_extended.gtf"
    local extension_report="$OUTPUT_DIR/${current_srr}_extension_report.tsv"

    # ============================================================
    # V11 STEP 1: Splice filter on original StringTie STRG IDs
    # (same as v10 — filter before ID remapping)
    # ============================================================
    if [[ ! -f "$filtered_gtf" && -f "$stringtie_gtf" ]]; then
        echo "[INFO] Running gffcompare on original StringTie output..."
        cd "$OUTPUT_DIR"
        gffcompare -r "$REFERENCE_DIR/${GENOME_NAME}.gtf" \
            -o "${current_srr}_stringtie_srr_annotation" \
            "$stringtie_gtf"

        echo "[INFO] Running splice junction filter on STRG IDs..."
        if python3 /app/filter_stringtie_splicing.py \
            --stringtie "$stringtie_gtf" \
            --reference "$REFERENCE_DIR/${GENOME_NAME}.gtf" \
            --annotated "$annotated_gtf" \
            --output "$filtered_gtf" \
            --report "$filter_report"; then
            echo "[INFO] Splice filter completed. Report: $filter_report"
        else
            echo "[ERROR] Splice filter failed!"
            exit 1
        fi
    fi

    # ============================================================
    # V13 STEP 2: Extend reference GTF terminal exon coordinates
    # Changes from v12:
    #   --max-extension raised from 2000 to 3000 bp
    #     Plant UTR literature supports genuine 3' UTRs up to 2-3kb.
    #     Neighbor overlap check and chromosome boundary cap still
    #     protect against false extensions.
    #   --min-coverage lowered from 3.0 to 2.0
    #     Matches new StringTie -c 2 threshold. Splice junction
    #     filter is the real biological quality gate.
    #   --fai unchanged: chromosome boundary hard cap still active.
    # ============================================================
    if [[ ! -f "$extended_ref_gtf" && -f "$filtered_gtf" ]]; then
        echo "[INFO] Extending reference UTR boundaries (v13: max 3000bp, min-coverage 2.0)..."
        if python3 /app/extend_utr_boundaries.py \
            --stringtie "$filtered_gtf" \
            --annotated "$annotated_gtf" \
            --reference "$REFERENCE_DIR/${GENOME_NAME}.gtf" \
            --output "$extended_ref_gtf" \
            --report "$extension_report" \
            --fai "${GENOME_FILE}.fai" \
            --max-extension 3000 \
            --min-coverage 2.0; then
            echo "[INFO] UTR extension completed. Report: $extension_report"
        else
            echo "[ERROR] UTR extension failed!"
            exit 1
        fi
    fi

    # Create final Cell Ranger reference using extended REFERENCE GTF
    # (not the StringTie GTF — this is the key v11 change, unchanged in v13)
    cd "$REFERENCE_DIR"
    local final_genome="${current_srr}_stringtie_genome"
    if [[ ! -d "$final_genome" ]]; then
        echo "[INFO] Creating final Cell Ranger reference (extended coordinates)..."
        cellranger mkref --genome="$final_genome" --fasta="$GENOME_FILE" --genes="$extended_ref_gtf"
    fi

    # Final Cell Ranger count
    cd "$OUTPUT_DIR"
    local final_matrix_id="${current_srr}_stringtie_matrix"
    if [[ ! -d "$final_matrix_id" ]]; then
        echo "[INFO] Running final Cell Ranger count..."
        cellranger count --id="$final_matrix_id" \
            --transcriptome="$REFERENCE_DIR/$final_genome" \
            --fastqs="$OUTPUT_DIR/$current_srr" \
            --sample="$current_srr" \
            $cellcount_arg --create-bam=true
    fi

    echo "✅ Completed $current_srr → $final_matrix_id"
}

# -----------------------------------
# Main
# -----------------------------------
if [[ -n "$SRR_LIST_FILE" ]]; then
    while IFS= read -r line; do
        [[ -z "$line" || "$line" =~ ^# ]] && continue
        process_single_srr "$(echo "$line" | tr -d '[:space:]')" "$EXPECTED_CELLS"
    done < "$SRR_LIST_FILE"
else
    process_single_srr "$SRR_ID" "$EXPECTED_CELLS" "$FORCE_CELLS"
fi
