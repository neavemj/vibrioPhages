
# create linear vibrio phage map using biopython

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors
from reportlab.lib.units import cm
import os


genomeDiag = GenomeDiagram.Diagram("vibrioPhages")
phage_regions_dir = '/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/1.phageRegions/'
phage_regions_file = '/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/1.phageRegions/regions.txt'

max_len = 0


for phage_regions in open(phage_regions_file):
    phage_regions = phage_regions.strip()
    print "processing file: ", phage_regions
    phage_handle = open(phage_regions_dir + phage_regions)
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

        genomeSet.add_feature(feature, name=name, label=False, label_position="start", sigil="ARROW")
    print 'initialPosition: ', initPosition
    print 'max_len: ', max_len

## add cross links between tracks to indicate blast similarity


BAAvsP1_4_handle = open("/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/3.allVsall/BAA450vsP1_4.edit.txt")

track_BAA450 = genomeDiag.tracks[9]
track_P1_4 = genomeDiag.tracks[8]

for connection in BAAvsP1_4_handle:
    connection = connection.strip()
    values = connection.split("\t")
    track_BAA450_start = float(values[0])
    track_BAA450_end = float(values[1])
    track_P1_4_start = float(values[2])
    track_P1_4_end = float(values[3])
    score = float(values[4])

    color = colors.linearlyInterpolatedColor(colors.white, colors.firebrick, 0, 100, score)
    link_xy = CrossLink((track_BAA450, track_BAA450_start, track_BAA450_end), (track_P1_4, track_P1_4_start, track_P1_4_end), color=color)

    feature = SeqFeature(FeatureLocation(int(track_BAA450_start), int(track_BAA450_end), strand=0))
    genomeSet = track_BAA450.new_set()
    genomeSet.add_feature(feature, name=name, label=False, label_position="start", sigil="BOX", color=color)


    genomeDiag.cross_track_links.append(link_xy)


genomeDiag.draw(format='linear', pagesize=("A4"), fragments=1, start=0, end=max_len)
genomeDiag.write("BAA450_VCY.pdf", "pdf")
genomeDiag.write("BAA450_VCY.svg", "svg")



