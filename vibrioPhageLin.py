

# create linear vibrio phage map using biopython
# imported required modules

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors
from reportlab.lib.units import cm
import os


# first read through PHAST results and extract information
# dictionary structure is: {phage : {geneNum : {start : 1, stop : 100, name : hypoth}, totalLength : 10000 } }

phage_regions_dir = '/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/1.phageRegions/'
phage_regions_file = '/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/1.phageRegions/regions.txt'

max_len = 0
phageDict = {}

for phage_regions in open(phage_regions_file):
    phage_regions = phage_regions.strip()
    print "processing file:", phage_regions
    phage_handle = open(phage_regions_dir + phage_regions)
    phageDict[phage_regions] = {}
    lineCount = 0
    total_len = 0
    for line in phage_handle:
        lineCount += 1
        line = line.strip()
        cols = line.split()
        geneNum = lineCount
        region = cols[0]
        if region.startswith("complement"):
            region = region.lstrip("complement(").rstrip(")").split("..")
            if lineCount == 1:
                initPosition = int(region[0])
            startRegion = int(region[0]) - initPosition
            stopRegion = int(region[1]) - initPosition
            total_len = max(total_len, stopRegion)
            phageDict[phage_regions][geneNum] = {"start" : startRegion, "stop" : stopRegion, "rev" : True}

        else:
            region = region.split("..")
            if lineCount == 1:
                initPosition = int(region[0])
            startRegion = int(region[0]) - initPosition
            stopRegion = int(region[1]) - initPosition
            total_len = max(total_len, stopRegion)
            phageDict[phage_regions][geneNum] = {"start" : startRegion, "stop" : stopRegion}
    phageDict[phage_regions]["totalLength"] = total_len


print phageDict

# ok create the genome diagram and tracks for each phage

genomeDiag = GenomeDiagram.Diagram("vibrioPhages")
genomeTrack = genomeDiag.new_track(1, name=phage_regions, greytrack=True, start=0, end=total_len)
genomeSet = genomeTrack.new_set()


feature = SeqFeature(FeatureLocation(startRegion, stopRegion), strand=-1)
feature = SeqFeature(FeatureLocation(startRegion, stopRegion), strand=+1)

genomeSet.add_feature(feature, name=name, label=False, label_position="start", sigil="ARROW")
