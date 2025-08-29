import gzip

def count_nonzeros_tsv(input_file):
    open_func = gzip.open if input_file.endswith(".gz") else open
    count = 0

    with open_func(input_file, "rt") as f:
        _ = f.readline()  # skip header
        for line_num, line in enumerate(f, start=1):
            values = line.rstrip("\n").split("\t")[1:]  # skip sample ID
            for v in values:
                if v != "0" and v != "0.0" and v != "":
                    count += 1

            if line_num % 1000 == 0:
                print(f"Processed {line_num} rows...")

    print(f"✅ Total non-zero values: {count:,}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python3 count_nonzeros.py <input.tsv[.gz]>")
    else:
        count_nonzeros_tsv(sys.argv[1])