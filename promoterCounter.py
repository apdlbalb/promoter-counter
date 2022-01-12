# IMPORT MODULE ---------------------------------------------------------
import sys
import os
from Bio import SeqIO
import re

# DEFINE CONSTANT -------------------------------------------------------
PROMOTER_REGION_LEN = 500

# DEFINE CLASS ----------------------------------------------------------
class Feature:

    # This Feature constructor takes a feature line and breaks it up into
    # variables that correspond to the fields of a GFF file
    def __init__(self, GFFLine):
        fields = GFFLine.split("\t")
        self.seqname = fields[0]
        self.source = fields[1]
        self.feature = fields[2]
        self.start = fields[3]
        self.end = fields[4]
        self.score = fields[5]
        self.strand = fields[6]
        self.frame = fields[7]
        self.attribute = fields[8]

# DEFINE FUNCTIONS -------------------------------------------------------

# reverseComplement takes a DNA sequence and returns its reverse complemet
# String <- String
def reverseComplement(seq):
    pos = len(seq) - 1
    newSeq = ""
    while pos >= 0:
        letter = seq[pos]
        if letter == "A":
            newSeq = newSeq + "T"
        elif letter == "T":
            newSeq = newSeq + "A"
        elif letter == "G":
            newSeq = newSeq + "C"
        elif letter == "C":
            newSeq = newSeq + "G"
        pos -= 1
    return newSeq

# GFFToObjects takes the lines of a GFF file and converts them into objects
# (listof Feature) <- (listof String)
def GFFToObjects(GFF):
    GFFObjects = []
    for line in GFF:
        if line[0] != "#":
            GFFObjects.append(Feature(line))
    return GFFObjects

# getGeneObjects takes a list of Features and filters for only genes
# (listof Feature) <- (listof Feature)
def getGeneObjects(GFFObjects):
    geneObjects = []
    for object in GFFObjects:
        if object.feature == "gene":
            geneObjects.append(object)
    return geneObjects

# chrGenesOfInterest filters a list of Feature objects that are genes for only
# those that correspond to the genes listed in the gene file
# (listof Feature) <- (listof Feature) (listof String)
def chrGenesOfInterest(geneObjects, geneIDs):
    objectsInGeneID = []
    for gene in geneIDs:
        for geneObject in geneObjects:
            if geneObject.attribute.find(gene) != -1:
                objectsInGeneID.append(geneObject)
    return objectsInGeneID

# getUpstreamSeqs takes a list of Feature objects that are genes in a
# chromosome and the corresponding DNA sequence and finds the promoter regions
# (listof String) <- String (listof Feature)
def getUpstreamSeqs(chrFastaSeq, chrGeneObjectList):
    upstreamSeqs = []
    for geneObject in chrGeneObjectList:
        promoterSeqLen = 0
        promoterSeq = ""
        if geneObject.strand == "+":
            i = int(geneObject.start) - 2
            while (chrFastaSeq[i] != "N") and (promoterSeqLen < PROMOTER_REGION_LEN) and (i >= 0):
                promoterSeq = chrFastaSeq[i] + promoterSeq
                i -= 1
                promoterSeqLen += 1
        else:
            i = int(geneObject.end)
            while (chrFastaSeq[i] != "N") and (promoterSeqLen < PROMOTER_REGION_LEN) and (i < len(chrFastaSeq)):
                promoterSeq = promoterSeq + chrFastaSeq[i]
                i += 1
                promoterSeqLen += 1
            promoterSeq = reverseComplement(promoterSeq)
        upstreamSeqs.append(promoterSeq)
    return upstreamSeqs

# countMotifs counts the number of times a motif appears in a list of promoter
# region sequences
# (dictof String : Int) <- (listof String) (listof String)
def countMotifs(promoters, upstreamSeqs):
    promoterDict = {}
    for promoterSeq in promoters:
        count = 0
        promoterSeqUp = promoterSeq.upper()
        for seq in upstreamSeqs:
            count += len(re.findall("(?=" + promoterSeqUp + ")", seq))
        promoterDict[promoterSeq] = count
    return promoterDict

# GET ARGUMENTS ---------------------------------------------------------
geneFilename = sys.argv[1]
promoterFilename = sys.argv[2]
GFFFastaFolderFilename = sys.argv[3]

# UNPACK PROMOTERS FILE -------------------------------------------------
promoterFile = open(promoterFilename, "r")
promoters = promoterFile.read().splitlines()
promoterFile.close()

# UNPACK GENES FILE ------------------------------------------------------
geneFile = open(geneFilename, "r")
genes = geneFile.read().splitlines()
geneFile.close()

# UNPACK FOLDER -----------------------------------------------------------
fastaFilenames = []
GFFFilenames = []

# Organize fasta files and GFF files into lists
for name in os.listdir(GFFFastaFolderFilename):
    if name.endswith(".fa"):
        fastaFilenames.append(name)
    elif name.endswith(".gff3"):
        GFFFilenames.append(name)

# UNPACK FOLDER FOR FASTA -------------------------------------------------
seqRecords = {}
for fastaFile in fastaFilenames:
    temp = next(SeqIO.parse(GFFFastaFolderFilename + fastaFile, "fasta"))
    seqRecords[int(temp.id)] = temp

# UNPACK FOLDER FOR GFFs ---------------------------------------------------
GFFs = {}

for GFFFilename in GFFFilenames:
    GFFFile = open(GFFFastaFolderFilename + GFFFilename, "r")
    GFFLines = GFFFile.readlines()
    GFFs[int(GFFLines[7].split("\t")[0])] = GFFLines
    GFFFile.close()

# COMBINE ALL UPSTREAM SEQUENCES ------------------------------------------
upSeqs = []

for key in GFFs:
    seq = seqRecords[key].seq
    GFFLines = GFFs[key]
    allFeatures = GFFToObjects(GFFLines)
    allGenes = getGeneObjects(allFeatures)
    genesInChr = chrGenesOfInterest(allGenes, genes)
    upSeqs.extend(getUpstreamSeqs(seq, genesInChr))

# COUNT --------------------------------------------------------------------
promoterDict = countMotifs(promoters, upSeqs)

# PRINT REPORT ------------------------------------------------------------
reportFile = open("report.txt", "a")
for i in range(0, len(promoters)):
    reportFile.write(promoters[i])
    reportFile.write("\t")
    reportFile.write(str(promoterDict[promoters[i]]))
    reportFile.write("\n")
reportFile.close()
