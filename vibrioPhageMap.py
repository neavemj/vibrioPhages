
# create linear vibrio phage map using biopython

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm


genomeDiag = GenomeDiagram.Diagram("vibrioPhages")
genomeTrack = genomeDiag.new_track(1, greytrack=False)
genomeSet = genomeTrack.new_set()


BAA450_VCY = open("/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/phageRegions/BAA450_VCY_region.fa")

for line in BAA450_VCY:
    line = line.strip()
    cols = line.split()
    region = cols[0]
    name = cols[1]
    region = region.split("..")
    startRegion = int(region[0])
    stopRegion = int(region[1])

    feature = SeqFeature(FeatureLocation(startRegion, stopRegion), strand=+1)
    genomeSet.add_feature(feature, name=name, label=True, sigil="ARROW")


genomeDiag.draw(format='linear', pagesize=(15*cm, 4*cm), fragments=1, start=4431042, end=4438916)
genomeDiag.write("BAA450_VCY.pdf", "pdf")




