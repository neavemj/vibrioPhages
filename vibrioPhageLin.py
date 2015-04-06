

# create linear vibrio phage map using biopython
# Matthew J. Neave 2.4.2015
# imported required modules

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors
from reportlab.lib.units import cm
import brewer2mpl

# first read through PHAST results and extract information
# dictionary structure is: {phage : {geneNum : {start : 1, stop : 100, name : hypoth}, totalLength : 10000 } }

phage_regions_dir = '/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/1.phageRegions/'
phage_regions_file = '/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/1.phageRegions/regions.txt'

max_len = 0
phageDict = {}
phageMaxLengths = {}
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
        geneType = cols[2]
        region = cols[0]
        if region.startswith("complement"):
            region = region.lstrip("complement(").rstrip(")").split("..")
            if lineCount == 1:
                initPosition = int(region[0])
            startRegion = int(region[0]) - initPosition
            stopRegion = int(region[1]) - initPosition
            total_len = max(total_len, stopRegion)
            phageDict[phage_regions][geneNum] = {"start" : startRegion, "stop" : stopRegion, "rev" : True, "geneType" : geneType}

        else:
            region = region.split("..")
            if lineCount == 1:
                initPosition = int(region[0])
            startRegion = int(region[0]) - initPosition
            stopRegion = int(region[1]) - initPosition
            total_len = max(total_len, stopRegion)
            phageDict[phage_regions][geneNum] = {"start" : startRegion, "stop" : stopRegion, "geneType" : geneType}
    phageMaxLengths[phage_regions] = total_len


## take the reverse complement of a couple of phages so it looks better

phageToRevComplement = ["OCN008_K139_region.fa", "RE98_web_2_ep3.fa"]

for phageRev in phageToRevComplement:
    for genes in phageDict[phageRev]:
        phageDict[phageRev][genes]["start"] = phageMaxLengths[phageRev] - phageDict[phageRev][genes]["start"]
        phageDict[phageRev][genes]["stop"] = phageMaxLengths[phageRev] - phageDict[phageRev][genes]["stop"]
        if "rev" in phageDict[phageRev][genes]:
            del phageDict[phageRev][genes]["rev"]
        else:
            phageDict[phageRev][genes]["rev"] = "True"


# ok create the genome diagram and tracks for each phage

genomeDiag = GenomeDiagram.Diagram("vibrioPhages")

for genome in phageList:
    genomeTrack = genomeDiag.new_track(1, name=genome, greytrack=False, start=0, end=phageMaxLengths[genome])

# dictionary to translate phage name into track number

phageTrack = {"BAA450" : 9, "P1_4" : 8, "RE98_1" : 7, "OCN008" : 6, "P1_1" : 5, "P1_3" : 4, "RE98_2" : 3, "P1_2" : 2, "OCN014" : 1}
phageFullNameConvert = {"BAA450" : "BAA450_VCY_region.fa", "P1_4" : "P1_4_VCY.fa", "RE98_1" : "RE98_web_1_kappa.fa", "OCN008": "OCN008_K139_region.fa", "P1_1" : "P1_1_KS14.fa", "P1_3" : "P1_3_Yersin.fa", "RE98_2" : "RE98_web_2_ep3.fa", "P1_2" : "P1_2_CTX.fa", "OCN014" : "OCN014_VP882_region.fa"}
fullNameToNum = {"BAA450_VCY_region.fa" : 9, "P1_4_VCY.fa" : 8, "RE98_web_1_kappa.fa" : 7, "OCN008_K139_region.fa" : 6, "P1_1_KS14.fa" : 5, "P1_3_Yersin.fa" : 4, "RE98_web_2_ep3.fa" : 3, "P1_2_CTX.fa" : 2, "OCN014_VP882_region.fa" : 1}


# read in all-vs-all blast results to create links between similar genes
# file should be first gene name, second gene name then percent similarity
# i'll add these first so the arrows come out on top

links_file_handle = open("/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/3.allVsall/links.txt")
links_dir = '/Users/neavemj/otherProjects/vibrioCoral/4.17.2.15/PHAST_phage_finder/3.allVsall/'


for link_tracks in links_file_handle:
    link_tracks = link_tracks.strip()
    links_handle = open(links_dir + link_tracks)

    for connection in links_handle:
        connection = connection.strip()
        values = connection.split("\t")

        phage1 = values[0]
        phage1Full = phageFullNameConvert[phage1]
        phage1_gene = int(values[1])
        phage2 = values[2]
        phage2Full = phageFullNameConvert[phage2]
        phage2_gene = int(values[3])
        score = float(values[4])

        track1 = genomeDiag.tracks[phageTrack[phage1]]
        track2 = genomeDiag.tracks[phageTrack[phage2]]

        color = colors.linearlyInterpolatedColor(colors.white, colors.firebrick, 0, 100, score)

        if phage1 == "P1_1":
            link_xy = CrossLink((track1, phageDict[phage1Full][phage1_gene]["start"], phageDict[phage1Full][phage1_gene]["stop"]), (track2, phageDict[phage2Full][phage2_gene]["start"], phageDict[phage2Full][phage2_gene]["stop"]), color=color, flip=False)
        else:
            link_xy = CrossLink((track1, phageDict[phage1Full][phage1_gene]["start"], phageDict[phage1Full][phage1_gene]["stop"]), (track2, phageDict[phage2Full][phage2_gene]["start"], phageDict[phage2Full][phage2_gene]["stop"]), color=color, flip=True)

        # add link features to first track

        BoxFeatureTrack1 = SeqFeature(FeatureLocation(phageDict[phage1Full][phage1_gene]["start"], phageDict[phage1Full][phage1_gene]["stop"], strand=0))

        genomeSet = track1.new_set()
        genomeSet.add_feature(BoxFeatureTrack1, label=False, label_position="start", sigil="BOX", color=color)
        genomeDiag.cross_track_links.append(link_xy)

        # add link features to second track

        BoxFeatureTrack2 = SeqFeature(FeatureLocation(phageDict[phage2Full][phage2_gene]["start"], phageDict[phage2Full][phage2_gene]["stop"], strand=0))

        genomeSet = track2.new_set()
        genomeSet.add_feature(BoxFeatureTrack2, label=False, label_position="start", sigil="BOX", color=color)
        genomeDiag.cross_track_links.append(link_xy)




# colors for the arrows - currently similar to the PHAST website

dark2 = brewer2mpl.get_map('Dark2', 'qualitative', 8).mpl_colors


arrowColorDict = {'Oth':"white", 'fib':dark2[1], 'PLP':dark2[2], 'Int':"red", 'RNA':"pink", 'Por':dark2[3], 'Hyp':dark2[4],
                  'sha':dark2[5], 'Att':"black", 'Coa':dark2[6], 'Pla':dark2[7], 'Ter':dark2[0]}

# typeSet = set()
#
# for phage in phageDict:
#     for geneRegion in phageDict[phage]:
#         typeSet.add(phageDict[phage][geneRegion]["geneType"])
#
# print typeSet


# add the gene arrows
# if else statement to separate out genes in the reverse direction

for phage in phageDict:
    genomeSet = genomeDiag.tracks[fullNameToNum[phage]].new_set()
    for geneRegion in phageDict[phage]:
        if "rev" in phageDict[phage][geneRegion]:
            feature = SeqFeature(FeatureLocation(phageDict[phage][geneRegion]["start"], phageDict[phage][geneRegion]["stop"]), strand=-1)
        else:
            feature = SeqFeature(FeatureLocation(phageDict[phage][geneRegion]["start"], phageDict[phage][geneRegion]["stop"]), strand=+1)

        genomeSet.add_feature(feature, label=True, name=str(geneRegion), label_position="start", sigil="ARROW", color=arrowColorDict[phageDict[phage][geneRegion]["geneType"]])




# ok now draw diagram

genomeDiag.draw(format='linear', pagesize=(60*cm, 40*cm), fragments=1, start=0, end=max_len)
genomeDiag.write("phageLinear.pdf", "pdf")
genomeDiag.write("phageLinear.svg", "svg")
