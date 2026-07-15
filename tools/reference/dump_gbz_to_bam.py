"""GBZ → BAM extractor for pangenome-aware DV testing.

Pangenome-aware DV's --pangenome flag accepts BAM/CRAM at runtime
(or GBZ via the upstream load_gbz_into_shared_memory preprocessor).
v2's native binary supports BAM-only — to test against Docker output
without GBZ runtime support, we pre-extract the synthetic reads from
the GBZ once via this script.

Critical: iterate in SMALL CHUNKS (1 kb), not one big query() call.
The GbzReader's Query() returns the haplotype segments overlapping
the queried sub-range; a single big call only returns one chunk's
worth (~94 reads at the same start position). Chunked iteration with
deduplication recovers the full set of synthetic reads (~8.7k for
chr20:10M-10.1M from hprc-v1.1-mc-grch38.gbz).

Run inside the upstream pangenome Docker image (extracted zipapp on
PYTHONPATH for sam.SamReader → GbzReader dispatch):

  cd /tmp && rm -rf dv_src && mkdir dv_src && cd dv_src && \\
    unzip -o -q /opt/deepvariant/bin/make_examples_pangenome_aware_dv.zip
  pip install --quiet pysam
  PYTHONPATH=/tmp/dv_src/runfiles/com_google_protobuf/python:\\
    /tmp/dv_src/runfiles/com_google_deepvariant:\\
    /tmp/dv_src/runfiles \\
    python3 dump_gbz_to_bam.py
"""
import os, sys
RF = "/tmp/dv_src/runfiles"
sys.path.insert(0, f"{RF}/com_google_protobuf/python")
sys.path.insert(0, f"{RF}/com_google_deepvariant")
sys.path.insert(0, RF)

from third_party.nucleus.io import sam
from third_party.nucleus.protos import range_pb2
import pysam

GBZ      = os.environ.get("GBZ",      "/data/hprc.gbz")
OUT      = os.environ.get("OUT",      "/out/pangenome.bam")
REFS_BAM = os.environ.get("REFS_BAM", "/refs/HG003.chr20.bam")
CHR      = os.environ.get("CHR",      "chr20")
START    = int(os.environ.get("START", "10000000"))
END      = int(os.environ.get("END",   "10100000"))
CHUNK    = int(os.environ.get("CHUNK", "1000"))

reader = sam.SamReader(GBZ, ref_name='GRCh38', use_loaded_shared_memory=False)

src = pysam.AlignmentFile(REFS_BAM, "rb")
hdr_dict = src.header.to_dict()
hdr_dict["RG"] = [{"ID": "pangenome", "SM": "pangenome", "PL": "ILLUMINA"}]
out = pysam.AlignmentFile(OUT, "wb", header=hdr_dict)
src.close()

nuc2pysam = {1:0, 2:1, 3:2, 4:3, 5:4, 6:5, 7:6, 8:7, 9:8}

seen = set()
n = 0
for chunk_start in range(START, END, CHUNK):
    chunk_end = min(chunk_start + CHUNK, END)
    r = range_pb2.Range(reference_name="chr20",
                         start=chunk_start, end=chunk_end)
    try:
        reads = list(reader.query(r))
    except Exception as e:
        print(f"  chunk {chunk_start}-{chunk_end}: error {e}")
        continue
    for read in reads:
        key = (read.fragment_name,
               read.alignment.position.position,
               len(read.aligned_sequence))
        if key in seen:
            continue
        seen.add(key)
        aln = pysam.AlignedSegment(out.header)
        aln.query_name = read.fragment_name or f"r{n}"
        aln.query_sequence = read.aligned_sequence or "N"
        if read.aligned_quality:
            try:
                aln.query_qualities = pysam.qualitystring_to_array(
                    "".join(chr(q + 33) for q in read.aligned_quality))
            except Exception:
                pass
        aln.flag = 16 if read.alignment.position.reverse_strand else 0
        aln.reference_id = out.header.get_tid(read.alignment.position.reference_name)
        aln.reference_start = read.alignment.position.position
        aln.mapping_quality = max(read.alignment.mapping_quality, 0)
        aln.cigartuples = [(nuc2pysam.get(cu.operation, 0), cu.operation_length)
                            for cu in read.alignment.cigar]
        out.write(aln)
        n += 1

out.close()
pysam.sort("-o", OUT + ".sorted", OUT)
os.rename(OUT + ".sorted", OUT)
pysam.index(OUT)
print(f"Done. Wrote {n} reads to {OUT}")
