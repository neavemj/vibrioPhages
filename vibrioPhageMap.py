
# create linear vibrio phage map using biopython

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm
import os

genomeDiag = GenomeDiagram.Diagram("vibrioPhages")


phage_region_dir = '/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/phageRegions/'

for phage_regions in os.listdir(phage_region_dir):

    phage_handle = open(phage_region_dir + phage_regions)
    genomeTrack = genomeDiag.new_track(1, name=phage_regions, greytrack=False)
    genomeSet = genomeTrack.new_set()
    lineCount = 0

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

        else:
            region = region.split("..")
            if lineCount == 1:
                initPosition = int(region[0])
            startRegion = int(region[0]) - initPosition
            stopRegion = int(region[1]) - initPosition
            feature = SeqFeature(FeatureLocation(startRegion, stopRegion), strand=+1)

        genomeSet.add_feature(feature, name=name, label=True, sigil="ARROW")



genomeDiag.draw(format='linear', pagesize=("A4"), fragments=1, start=0, end=10000)
genomeDiag.write("BAA450_VCY.pdf", "pdf")
genomeDiag.write("BAA450_VCY.svg", "svg")



