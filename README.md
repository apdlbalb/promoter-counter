# promoter-counter
This program reads in a text file of gene IDs, a text file of promoter motif sequences, and a directory of GFF3 and fasta files and creates a report of the number of times each motif was found. 

## Dependencies
- Written in Python 3
- Requires SeqIO from Biopython
- Program assumes all GFF3 files in `FastaGFF/` have a matching fasta file in the directory with chromosome number as the key

## Usage
Assuming all are in the same directory, please run in command line as follows:
```
./promoterCounter.py genesOfInterest.txt promoterSeqs.txt FastaGFF/
```

where `genesOfInterest.txt` looks like:
```
Zm00001d002384
Zm00001d038084
Zm00001d011797
Zm00001d011709
Zm00001d021620
Zm00001d003756
Zm00001d008559
```

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

and where the `FastaGFF/` directory contains the Ensembl GFF3 and Fasta files of a species.



