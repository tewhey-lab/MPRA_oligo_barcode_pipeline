### A guide to output files for `MPRAmatch` and `MPRAcount`

`MPRAmatch`

  - `.merged.match`
    1. Sequence ID
    2. Found Barcode Sequence
    3. Found Oligo Sequence
    4. Found Oligo Length
    5. Flashed Sequence Length


  - `.merged.match.enh.mapped`
    1. Sequence ID
    2. Barcode
    3. Mapping flag from MiniMap2
    4. Best mapped ID
    5. All mapped IDs
    6. count
    7. mapped length
    8. CIGAR
    9. Adjusted score ( (mismatched bp + inserted bp + deleted bp + non-aligned length) / total oligo length )
    10. oligo sequence
    11. PASS/FAIL (PASS if adjusted score <= 0.05)
    12. Raw error
    13. MD/cs tag
    14. start/stop from mapping within reference oligo sequence


  - `.merged.match.enh.mapped.barcode.ct`
    1. Barcode
    2. Oligo(s) (comma separated if multiple oligos)
    3. individual seen (comma separated for each oligo)
    4. total seen
    5. overall flag (0/1/2)
      - 0: pass, no conflict with mapping; single oligo associated with barcode
      - 1: conflict, barcode associated with multiple oligos
      - 2: fail, barcode not associated with any oligo or failed mapping to any
    6. individual flags (comma separated for each oligo)
    7. error (best for barcode-oligo combo)
    8. CIGAR (best for barcode-oligo combo)
    9. MD/cs tag (best for barcode-oligo combo)
    10. start/stop from mapping within reference oligo sequence (best for barcode-oligo combo)


  - `.merged.match.enh.mapped.barcode.ct.parsed` <br>
  These columns are the same as for `.merged.match.enh.mapped.barcode.ct`. The barcodes that were identified as "conflict" (column 5 == 1) are checked for obvious pass/fail ratios.


`MPRAcount`

  - `*.match`
    1. Sequence ID
    2. Found Barcode
    3. (matching oligo)


  - `*.tag`
    1. Barcode
    2. Total Seen
    3. Tag Flag
      -  0: Barcode found in parsed file and no issue mapping
      - -9: Barcode not found in parsed file
      - -1: Barcode found in parsed file, collision accepted
      - -4: Barcode found in parsed file, collision not accepted
      - -5: Barcode found in parsed file, mapping failed
      - -6: Barcode found in parsed file, mapped reverse complement
    4. Oligo ID (from parsed)
    5. Mapping Flag (from parsed; `-` if barcode not found in dictionary)
    6. Barcode
    7. error (from parsed)
    8. CIGAR (from parsed)
    9. MD/cs tag (from parsed)
    10. start/stop (from parsed)


  - `.stats`
    1. replicate ID
    2. good barcodes (mapping flag == 0)
    3. good reads (mapping flag == 0)
    4. total matched reads (mapping flag in [0,1,2])
    5. percent (good reads)/(total matched reads)
    6. all barcodes (mapping flag in [-,0,1,2])
    7. all reads (mapping flag in [-,0,1,2])
    8. percent (good reads)/(all reads)


  - `_condition.txt`
    1. replicate ID (row names when read into R for `MPRAmodel`)
    2. condition
