#!/bin/bash
set -e

# Input/Output/Reference paths 
SAMPLE_DIR="/app/input"
OUTPUT_DIR="/app/output"
REFERENCE_DIR="/app/output/reference"  

# Initialize variables
SRR_ID=""
SRR_LIST_FILE=""
EXPECTED_CELLS=1800000
GENOME_NAME=""

# Parse command-line options
while getopts "s:l:c:g:h" opt; do
    case $opt in
        s) SRR_ID="$OPTARG" ;;
        l) SRR_LIST_FILE="$OPTARG" ;;
        c) EXPECTED_CELLS="$OPTARG" ;;
        g) GENOME_NAME="$OPTARG" ;;
        h)
            echo "Usage: $0 [-s SRR_ID | -l SRR_LIST_FILE] [-c EXPECTED_CELLS] -g GENOME_NAME"
            echo "  - Input files: GENOME_NAME.fa and GENOME_NAME.gff3 under /app/input/"
            echo "  - Use -s for a single SRR, or -l for a list (in srr_list.txt)"
            exit 0
            ;;
        *)
            echo "Invalid option"
            exit 1
            ;;
    esac
done

# Validate arguments
if [[ -z "$SRR_ID" && -z "$SRR_LIST_FILE" ]]; then
    echo "ERROR: Must provide either a single SRR (-s) or a list file (-l)"
    exit 1
fi

if [[ -z "$GENOME_NAME" ]]; then
    echo "ERROR: Must provide species name with -g"
    exit 1
fi

# Locate genome and GFF3 file based on species name (case-insensitive)
for ext in fa fasta; do
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

# Confirm files found
if [[ ! -f "$GENOME_FILE" ]]; then
    echo "ERROR: Genome file not found for species '$GENOME_NAME'"
    ls -la /app/input/
    exit 1
fi

if [[ ! -f "$GFF3_FILE" ]]; then
    echo "ERROR: GFF3 file not found for species '$GENOME_NAME'"
    ls -la /app/input/
    exit 1
fi

echo "Genome: $GENOME_FILE"
echo "Annotation: $GFF3_FILE"

# Core function
process_single_srr() {
    local current_srr="$1"
    local expected_cells="$2"

    # Check input files
    local srr_path="$SAMPLE_DIR/$current_srr"
    if [[ ! -d "$srr_path" ]]; then
        echo "ERROR: SRR folder not found: $srr_path"
        exit 1
    fi

    mkdir -p "$OUTPUT_DIR/$current_srr" "$REFERENCE_DIR"

    # Prepare FASTQ files
    cd "$OUTPUT_DIR/$current_srr"
    for i in 1 2 3; do
        fq_gz="$srr_path/${current_srr}_${i}.fastq.gz"
        fq_raw="$srr_path/${current_srr}_${i}.fastq"
        out_file="${current_srr}_S1_L001_R${i}_001.fastq"

        if [[ -f "$fq_gz" ]]; then
            echo "Decompressing: $fq_gz → $out_file"
            gunzip -c "$fq_gz" > "$out_file"
        elif [[ -f "$fq_raw" ]]; then
            echo "Copying: $fq_raw → $out_file"
            cp "$fq_raw" "$out_file"
        else
            echo "Warning: FASTQ file not found for read R$i"
        fi
    done

    # Rename R3 to I1 if exists
    if [[ -f "${current_srr}_S1_L001_R3_001.fastq" ]]; then
        mv "${current_srr}_S1_L001_R3_001.fastq" "${current_srr}_S1_L001_I1_001.fastq"
    fi

    # Convert GFF3 to GTF
    cd "$REFERENCE_DIR"
    if [[ ! -f "${GENOME_NAME}.gtf" ]]; then
        echo "Converting GFF3 to GTF..."

        # Conversion
        gffread "$GFF3_FILE" -T -E -o "${GENOME_NAME}.gtf"

        exon_count=$(grep -c "exon" "${GENOME_NAME}.gtf" || echo "0")
        echo "GTF conversion completed. Exon features found: $exon_count"

        if [[ "$exon_count" -eq 0 ]]; then
            echo "WARNING: No exon features found in GTF file"
            echo "First 10 lines of converted GTF:"
            head -10 "${GENOME_NAME}.gtf"
        fi
    fi

    # Make initial reference genome
    if [[ ! -d "${GENOME_NAME}_genome" ]]; then
        cellranger mkref --genome="${GENOME_NAME}_genome" --fasta="$GENOME_FILE" --genes="${GENOME_NAME}.gtf"
    fi

    # Run Cell Ranger count (initial)
    cd "$OUTPUT_DIR"
    matrix_id="${current_srr}_matrix"
    if [[ ! -d "$matrix_id" ]]; then
        cellranger count --id="$matrix_id" --transcriptome="$REFERENCE_DIR/${GENOME_NAME}_genome" --fastqs="$OUTPUT_DIR/$current_srr" --sample="$current_srr" --expect-cells="$expected_cells" --create-bam=true
    fi

    # Create CDS-only GTF and CDS genome
    cd "$REFERENCE_DIR"
    if [[ ! -f "${GENOME_NAME}.cds2exon.gtf" ]]; then
        grep 'CDS' "${GENOME_NAME}.gtf" | sed 's/CDS/exon/g' > "${GENOME_NAME}.cds2exon.gtf"
    fi

    if [[ ! -d "${GENOME_NAME}_cds2exon_genome" ]]; then
        cellranger mkref --genome="${GENOME_NAME}_cds2exon_genome" --fasta="$GENOME_FILE" --genes="${GENOME_NAME}.cds2exon.gtf"
    fi

    # Run Cell Ranger count on CDS genome
    cd "$OUTPUT_DIR"
    cds_matrix_id="${current_srr}_cds2exon_matrix"
    if [[ ! -d "$cds_matrix_id" ]]; then
        cellranger count --id="$cds_matrix_id" --transcriptome="$REFERENCE_DIR/${GENOME_NAME}_cds2exon_genome" --fastqs="$OUTPUT_DIR/$current_srr" --sample="$current_srr" --expect-cells="$expected_cells" --create-bam=true
    fi

    # Run StringTie on BAM
    cd "$OUTPUT_DIR"
    local stringtie_gtf="${current_srr}_stringtie.gtf"
    local bam_file="${cds_matrix_id}/outs/possorted_genome_bam.bam"
    if [[ ! -f "$stringtie_gtf" && -f "$bam_file" ]]; then
        stringtie "$bam_file" -G "$REFERENCE_DIR/${GENOME_NAME}.cds2exon.gtf" -o "$stringtie_gtf"
    fi

    # Format final StringTie GTF
    local formatted_gtf="${current_srr}_stringtie_v2.gtf"
    if [[ ! -f "$formatted_gtf" && -f "$stringtie_gtf" ]]; then
        python3 /app/update_ref_gene_format.py -i "$stringtie_gtf" -o "$formatted_gtf"
    fi

    # Create final reference genome 
    cd "$REFERENCE_DIR"
    local final_genome_name="${current_srr}_stringtie_genome"
    if [[ ! -d "$final_genome_name" && -f "$OUTPUT_DIR/$formatted_gtf" ]]; then
        cellranger mkref --genome="$final_genome_name" --fasta="$GENOME_FILE" --genes="$OUTPUT_DIR/$formatted_gtf"
    fi

    # Final Cell Ranger count
    cd "$OUTPUT_DIR"
    local final_matrix_id="${current_srr}_stringtie_matrix"
    if [[ ! -d "$final_matrix_id" ]]; then
        cellranger count --id="$final_matrix_id" --transcriptome="$REFERENCE_DIR/$final_genome_name" --fastqs="$OUTPUT_DIR/$current_srr" --sample="$current_srr" --expect-cells="$expected_cells" --create-bam=true
    fi

    echo "Done with $current_srr → $final_matrix_id"
}

# Main execution
if [[ -n "$SRR_LIST_FILE" ]]; then
    while IFS= read -r line; do
        [[ -z "$line" || "$line" =~ ^# ]] && continue
        current_srr=$(echo "$line" | tr -d '[:space:]')
        [[ -n "$current_srr" ]] && process_single_srr "$current_srr" "$EXPECTED_CELLS"
    done < "$SRR_LIST_FILE"
else
    process_single_srr "$SRR_ID" "$EXPECTED_CELLS"
fi
