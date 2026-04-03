def parse_vcf(file):
    """
    Robust VCF parser for real-world files.
    """

    variants = []

    for line in file:

        # Convert bytes → string
        if isinstance(line, bytes):
            line = line.decode("utf-8")

        # Skip headers
        if line.startswith("#"):
            continue

        line = line.strip()

        if not line:
            continue

        # Split columns
        parts = line.split("\t") if "\t" in line else line.split()

        if len(parts) < 5:
            continue

        try:
            chrom = parts[0].replace("chr", "")
            pos   = int(parts[1])   # ✅ fixed
            ref   = parts[3]
            alt_field = parts[4]

            # Handle multiple ALT alleles
            alt_list = alt_field.split(",")

            # INFO field
            info = parts[7] if len(parts) >= 8 else ""

            gene = "Unknown"

            # Try extracting gene from common formats
            if "GENE=" in info:
                gene = info.split("GENE=")[1].split(";")[0]

            elif "ANN=" in info:
                # SnpEff format
                gene = info.split("ANN=")[1].split("|")[3]

            elif "CSQ=" in info:
                # VEP format (simplified)
                gene = info.split("CSQ=")[1].split("|")[3]

            # Create one entry per ALT allele
            for alt in alt_list:

                if alt in [".", "*", ""]:
                    continue

                variants.append({
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "gene": gene
                })

        except Exception:
            continue  # skip bad lines safely

    return variants[:10]