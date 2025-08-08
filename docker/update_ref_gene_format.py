#!/usr/bin/env python3
import argparse

def update_gtf(input_gtf, output_gtf):
    with open(input_gtf, "r") as infile, open(output_gtf, "w") as outfile:
        for line in infile:
            if line.startswith("#") or not line.strip():
                outfile.write(line)
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                outfile.write(line)
                continue

            if fields[6] == ".":
                fields[6] = "+"

            attr_text = fields[8]
            attrs = []
            for attr in attr_text.strip().split(";"):
                attr = attr.strip()
                if attr:
                    try:
                        key, value = attr.split(" ", 1)
                        value = value.strip().strip('"')
                        attrs.append((key, value))
                    except ValueError:
                        continue

            attr_dict = dict(attrs)

            if "reference_id" in attr_dict and "ref_gene_id" in attr_dict:
                new_gene = attr_dict["ref_gene_id"]
                if new_gene.startswith("gene:"):
                    new_gene = new_gene[len("gene:"):]
                attr_dict["gene_id"] = new_gene

                new_transcript = attr_dict["reference_id"]
                if new_transcript.startswith("transcript:"):
                    new_transcript = new_transcript[len("transcript:"):]
                attr_dict["transcript_id"] = new_transcript

                attr_dict.pop("ref_gene_id", None)
                attr_dict.pop("reference_id", None)

                new_attrs = []
                seen_keys = set()
                for key, _ in attrs:
                    if key in ("ref_gene_id", "reference_id"):
                        continue
                    if key in attr_dict and key not in seen_keys:
                        new_attrs.append(f'{key} "{attr_dict[key]}"')
                        seen_keys.add(key)
                for key in attr_dict:
                    if key not in seen_keys:
                        new_attrs.append(f'{key} "{attr_dict[key]}"')
                fields[8] = "; ".join(new_attrs) + ";"

            outfile.write("\t".join(fields) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update gene_id and transcript_id using reference_id and ref_gene_id when present, and update strand (column 7) if it is '.' in a GTF file.")
    parser.add_argument("-i", "--input", required=True, help="Input GTF file")
    parser.add_argument("-o", "--output", required=True, help="Output (updated) GTF file")
    args = parser.parse_args()
    
    update_gtf(args.input, args.output)