import sys
import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
import time
import gzip
import os
import tempfile

def open_file_auto(path, mode='rt', buffering=None):
    if path.endswith('.gz'):
        return gzip.open(path, mode, encoding='utf-8')
    else:
        if buffering is not None:
            return open(path, mode, buffering=buffering, encoding='utf-8')
        else:
            return open(path, mode, encoding='utf-8')

def format_eta(seconds_left):
    if seconds_left < 60:
        return f"{int(seconds_left)}s"
    elif seconds_left < 3600:
        return f"{int(seconds_left // 60)}m {int(seconds_left % 60)}s"
    else:
        hours = int(seconds_left // 3600)
        minutes = int((seconds_left % 3600) // 60)
        return f"{hours}h {minutes}m"

def transpose_sparse_matrix_with_triplet_file(input_file, output_file=None, keep_triplets=False):
    global_start = time.time()
    print(f"Input file: {input_file}")

    if not output_file:
        output_file = input_file.replace('.tsv.gz', '_transposed.tsv.gz') \
                                .replace('.tsv', '_transposed.tsv')

    # --- Step 1: Parse header to get OTU IDs ---
    print("Reading header row...")
    with open_file_auto(input_file, 'rt') as f:
        header = [h.strip('"') for h in f.readline().strip().split('\t')]
        if header[0].strip() == "" or header[0].lower() in ["sampleid", "sample_id", "run"]:
            otu_ids = header[1:]
        else:
            otu_ids = header
        num_otus = len(otu_ids)

    sample_ids = []

    # --- Step 2: Create temporary triplet file ---
    triplet_fd, triplet_path = tempfile.mkstemp(suffix=".tsv", prefix="triplets_", text=True)
    os.close(triplet_fd)  # Will use open() not os.fdopen()

    print(f"Streaming input and writing triplets to: {triplet_path}")
    read_start = time.time()
    row_count = 0
    total_rows_estimate = None

    with open_file_auto(input_file, 'rt') as f_in, open(triplet_path, 'wt') as f_triplets:
        _ = f_in.readline()  # skip header
        for row_num, line in enumerate(f_in):
            parts = [p.strip('"') for p in line.rstrip('\n').split('\t')]
            if len(parts) == num_otus + 1:
                sample_id = parts[0]
                values = parts[1:]
            elif len(parts) == num_otus:
                sample_id = f"Sample_{row_num+1}"
                values = parts
            else:
                raise ValueError(f"Row {row_num+1} has unexpected number of columns: {len(parts)}")

            sample_ids.append(sample_id)

            for col_num, val in enumerate(values):
                if val != "0" and val != "0.0":
                    try:
                        f_triplets.write(f"{row_num}\t{col_num}\t{val}\n")
                    except ValueError:
                        raise ValueError(f"Invalid numeric value at row {row_num+1}, column {col_num+1}: '{val}'")

            row_count += 1
            if row_num % 1000 == 0 and row_num > 0:
                elapsed = time.time() - read_start
                print(f"  Processed ~{row_num} rows... Elapsed: {format_eta(elapsed)}")

    num_samples = len(sample_ids)
    print(f"Matrix shape: {num_samples} rows × {num_otus} columns")
    print(f"Finished triplet writing. Reading back into sparse matrix...")

    # --- Step 3: Read triplets from file and build matrix ---
    data = []
    row_idx = []
    col_idx = []

    with open(triplet_path, 'rt') as f_triplets:
        for line in f_triplets:
            r, c, v = line.strip().split('\t')
            row_idx.append(int(r))
            col_idx.append(int(c))
            data.append(float(v))

    print(f"Creating sparse matrix from {len(data)} non-zero entries...")
    mat = coo_matrix((data, (row_idx, col_idx)), shape=(num_samples, num_otus), dtype=np.float32)

    # --- Step 4: Transpose ---
    print("Transposing...")
    mat_T = mat.transpose().tocsc()

    # --- Step 5: Write transposed matrix ---
    print(f"Writing transposed matrix to: {output_file}")
    write_start = time.time()
    with open_file_auto(output_file, 'wt', buffering=16*1024*1024) as out:
        out.write("\t" + "\t".join(sample_ids) + "\n")

        total_columns = mat_T.shape[0]
        for i in range(total_columns):
            row = mat_T.getrow(i).toarray().flatten()
            row_str = [otu_ids[i]] + [str(int(x)) if x != 0 else "0" for x in row]
            out.write("\t".join(row_str) + "\n")

            if i % 1000 == 0 and i > 0:
                elapsed = time.time() - write_start
                percent = (i / total_columns) * 100
                eta = format_eta(elapsed / (i / total_columns) - elapsed)
                print(f"  Wrote {percent:.1f}% of columns... ETA: {eta}")

    if not keep_triplets:
        os.remove(triplet_path)
    else:
        print(f"✅ Kept intermediate triplet file: {triplet_path}")

    print(f"✅ Done in {round(time.time() - global_start, 2)} seconds.")

# --- Usage ---
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python transpose_sparse_triplets.py <input.tsv[.gz]>")
        print("Optional: add '--keep-triplets' to retain intermediate file")
    else:
        keep = "--keep-triplets" in sys.argv
        input_file = sys.argv[1]
        transpose_sparse_matrix_with_triplet_file(input_file, keep_triplets=keep)

    total_rows_estimate = None
