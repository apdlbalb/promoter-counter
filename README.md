# promoter-counter

## Dependencies
- Requires SeqIO from Biopython

## Usage
```
promoterCounter.py genesOfInterest.txt promoterSeqs.txt FastaGFF/
```

Where the `FastaGFF/` directory contains the Ensembl GFF and Fasta files of a species,

where `promoterSeqs.txt` looks like:
```
AAACCA
AAAG
AATAAA
AATCTGATCG
ACATAAAATAAAAAAAGGCA
acatgTGTAAAGgtatt
acatgTGTAAAGgtgaa
```

and where `genesOfInterest.txt` looks like:
```
Zm00001d002384
Zm00001d038084
Zm00001d011797
Zm00001d011709
Zm00001d021620
Zm00001d003756
Zm00001d008559
```
