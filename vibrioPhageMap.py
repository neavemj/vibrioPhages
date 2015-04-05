
# create linear vibrio phage map using biopython

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm
import os


genomeDiag = GenomeDiagram.Diagram("vibrioPhages")
phage_region_dir = '/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/phageRegions/'
max_len = 0


for phage_regions in os.listdir(phage_region_dir):
    print "processing file: ", phage_regions
    phage_handle = open(phage_region_dir + phage_regions)
    lineCount = 0
    total_len = 0
    for line in phage_handle:
        lineCount += 1
        line = line.strip()
        cols = line.split()
        region = cols[0]
        if region.startswith("complement"):
            region = region.lstrip("complement(").rstrip(")").split("..")
            if lineCount == 1:
                initPosition = int(region[0])
            startRegion = int(region[0]) - initPosition
            stopRegion = int(region[1]) - initPosition
            total_len = max(total_len, stopRegion)

        else:
            region = region.split("..")
            if lineCount == 1:
                initPosition = int(region[0])
            startRegion = int(region[0]) - initPosition
            stopRegion = int(region[1]) - initPosition
            total_len = max(total_len, stopRegion)


    genomeTrack = genomeDiag.new_track(1, name=phage_regions, greytrack=True, start=0, end=total_len)
    genomeSet = genomeTrack.new_set()
    lineCount = 0

    phage_handle.seek(0)
    for line in phage_handle:

        lineCount += 1
        line = line.strip()
        cols = line.split()
        region = cols[0]
        name = cols[1]
        if region.startswith("complement"):
            region = region.lstrip("complement(").rstrip(")").split("..")
            if lineCount == 1:
                initPosition = int(region[0])
            startRegion = int(region[0]) - initPosition
            stopRegion = int(region[1]) - initPosition
            feature = SeqFeature(FeatureLocation(startRegion, stopRegion), strand=-1)
            max_len = max(max_len, stopRegion)

        else:
            region = region.split("..")
            if lineCount == 1:
                initPosition = int(region[0])
            startRegion = int(region[0]) - initPosition
            stopRegion = int(region[1]) - initPosition
            feature = SeqFeature(FeatureLocation(startRegion, stopRegion), strand=+1)
            max_len = max(max_len, stopRegion)

        genomeSet.add_feature(feature, name=name, label=True, label_position="start", sigil="ARROW")



genomeDiag.draw(format='linear', pagesize=("A4"), fragments=1, start=0, end=max_len)
genomeDiag.write("BAA450_VCY.pdf", "pdf")
genomeDiag.write("BAA450_VCY.svg", "svg")



