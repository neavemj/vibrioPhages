

# create linear vibrio phage map using biopython
# Matthew J. Neave 2.4.2015
# imported required modules

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors
from reportlab.lib.units import cm


# first read through PHAST results and extract information
# dictionary structure is: {phage : {geneNum : {start : 1, stop : 100, name : hypoth}, totalLength : 10000 } }

phage_regions_dir = '/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/1.phageRegions/'
phage_regions_file = '/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/1.phageRegions/regions.txt'

max_len = 0
phageDict = {}
phageList = []

for phage_regions in open(phage_regions_file):
    phage_regions = phage_regions.strip()
    phageList.append(phage_regions)
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

for genome in phageList:
    genomeTrack = genomeDiag.new_track(1, name=genome, greytrack=True, start=0, end=phageDict[genome]["totalLength"])


# dictionary to translate phage name into track number

phageTrack = {"BAA450_VCY_region.fa" : 9, "P1_4_VCY.fa" : 8, "RE98_web_1_kappa.fa" : 7, "OCN008_K139_region.fa" : 6, "P1_1_KS14.fa" : 5, "P1_3_Yersin.fa" : 4, "RE98_web_2_ep3.fa" : 3, "P1_2_CTX.fa" : 2, "OCN014_VP882_region.fa" : 1}


# read in all-vs-all blast results to create links between similar genes
# file should be first gene name, second gene name then percent similarity
# want to add these links first so that the arrows go over the top

links_file_handle = open("/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/3.allVsall/links.txt")
links_dir = '/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/3.allVsall/'


for link_tracks in open(links_file_handle):
    link_tracks = link_tracks.strip()
    links_handle = open(links_dir + link_tracks)


    track_1 = genomeDiag.tracks[phageTrack[link_tracks]]

    for connection in links_handle:
        connection = connection.strip()
        values = connection.split("\t")
        track_1_start = float(values[0])
        track_1_end = float(values[1])
        track_2_start = float(values[2])
        track_2_end = float(values[3])
        score = float(values[4])



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
    genomeSet.add_feature(feature, label=False, label_position="start", sigil="BOX", color=color)


    genomeDiag.cross_track_links.append(link_xy)



# feature = SeqFeature(FeatureLocation(startRegion, stopRegion), strand=-1)
# feature = SeqFeature(FeatureLocation(startRegion, stopRegion), strand=+1)
#
# genomeSet = genomeTrack.new_set()
# genomeSet.add_feature(feature, name=name, label=False, label_position="start", sigil="ARROW")
#
# genomeDiag.draw(format='linear', pagesize=("A4"), fragments=1, start=0, end=max_len)
# genomeDiag.write("phageLinear.pdf", "pdf")
# genomeDiag.write("phageLinear.svg", "svg")
