#!/usr/bin/env python
import sys
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from pylab import *
import random

def main():
    if len(sys.argv) < 2 or "-h" in sys.argv or "--help" in sys.argv:
        print >> sys.stderr, "Calculates some statistics on single exon transcripts that we can use for figuring out if they are real or not."
        print >> sys.stderr, "single_exon_stats.py\t<transcripts.gtf>"
        sys.exit()


    filename = sys.argv[1]

    """
    scaffold1       StringTie       transcript      4845    5565    1000    .       .       gene_id "Hl.SW.STRG.1"; transcript_id "Hl.SW.STRG.1.1"; cov "3.993066"; FPKM "0.041226";
    scaffold1       StringTie       exon    4845    5565    1000    .       .       gene_id "Hl.SW.STRG.1"; transcript_id "Hl.SW.STRG.1.1"; exon_number "1"; cov "3.993066";

    maybe get coverage fpkm and closest transcript?
    """
    current_buffer = [] #contains exons in current transcript ? 
    scaffold_buffer = [] #contains list of exons
    ref,source,type,start,stop,length,strand = range(7)
    cur_scaffold = ""
    stats = {}
    for line in open(filename, 'r'):
        if line[0] == "#":
            continue
        line_arr = line.strip().split("\t")
        if line_arr[ref] != cur_scaffold:
            if cur_scaffold != "":
                stats[cur_scaffold] = compute_scaffold_stats(scaffold_buffer)
            cur_scaffold = line_arr[ref]
            scaffold_buffer = []    
        if line_arr[type] == "transcript":
            if current_buffer != []:
                scaffold_buffer.append(current_buffer)
            current_buffer = []
        current_buffer.append(line)
    num_exons = []
    for scaffold in stats.values():
        for transcript in scaffold:
            num_exons.append(transcript["num_exons"])
    max_exon = max(num_exons)

    pure_data = []
    for scaffold in stats.values():
        pure_data += scaffold
    
    run_stats(stats, max_exon)
    classifier(pure_data)

def compute_scaffold_stats(scaffold):
    ref,source,type,start,stop,length,strand,mystery,info = range(9)
    if scaffold == []:
        return
    else:
        scaffold_stats = []
        for features in scaffold:
            num_exons = len(features) - 1
            feature_arr = features[0].split("\t")
            cov, fpkm = get_cov_fpkm(feature_arr[info])
            trans_id = get_trans_id(feature_arr[info])
            start_p, stop_p = feature_arr[start], feature_arr[stop]
            length = int(stop_p) - int(start_p)
            transcript_stats = {"coverage": float(cov), "fpkm": float(fpkm), "start": int(start_p), "stop":int(stop_p), "length": length, "num_exons":num_exons, \
                                "distance_to_next":0, "id": trans_id}
            scaffold_stats.append(transcript_stats)
        sorted_stats = sorted(scaffold_stats, key=lambda k: k['start']) 
    
        for i in range(0, len(sorted_stats) - 1):
            cur = sorted_stats[i]
            next = sorted_stats[i+1]
            cur["distance_to_next"] = 0
            if next["start"] - cur["stop"] > 0:
                cur["distance_to_next"] = next["start"] - cur["stop"] 
        return sorted_stats

def get_trans_id(info):
    #gene_id "Hl.SW.STRG.1"; transcript_id "Hl.SW.STRG.1.1"; exon_number "1"; cov "3.993066";
    transcript_id = info.split("; ")[1].split("transcript_id ")[1].strip('"')
    return transcript_id
 
def get_cov_fpkm(info):
    arr = info.split(";")
    cov = arr[-3].strip().split(" ")[-1].strip('"')
    fpkm = arr[-2].strip().split(" ")[-1].strip('"')
    return cov, fpkm

def run_stats(data_dic, max_exon):
    distances = [[] for i in range(max_exon)] #array of num exons
    fpkms = [[] for i in range(max_exon)] #array of num exons
    covs = [[] for i in range(max_exon)] #array of num exons
    lengths = [[] for i in range(max_exon)] #array of num exons
    
    for scaffold in data_dic.values():
        for data_point in scaffold:
            if "distance_to_next" in data_point:
                distances[data_point["num_exons"] - 1].append(data_point["distance_to_next"])
            if "fpkm" in data_point:
                fpkms[data_point["num_exons"] - 1].append(data_point["fpkm"])
            if "coverage" in data_point:
                covs[data_point["num_exons"] - 1].append(data_point["coverage"])
            if "length" in data_point:
                lengths[data_point["num_exons"] - 1].append(data_point["length"])
        
    fig, axes = subplots()
    axes.set_title("Distance to next transcript (bp)")
    # basic plot
    axes.boxplot(distances[:10])
    ylim([50, 5000])
    ylabel("Distance (bp)")
    xlabel("Number of exons")
    savefig("distance.pdf")
    
    # notched plot
    fig, axes = subplots()
    axes.set_title("Transcript FPKM")
    # basic plot
    axes.boxplot(fpkms[:10])
    ylabel("FPKM")
    xlabel("Number of exons")
    ylim([0, 5])
    
    savefig("fpkm.pdf")

    fig, axes = subplots()
    axes.set_title("Transcript coverage")
    # basic plot
    axes.boxplot(covs[:10])
    ylim([0, 200])
    ylabel("Transcript Coverage")
    xlabel("Number of exons")
    savefig("coverage.pdf")

    fig, axes = subplots()
    axes.set_title("Transcript Lengths")
    # basic plot
    axes.boxplot(lengths[:10])
    ylim([0, 20000])
    ylabel("Distance (bp)")
    xlabel("Number of exons")
    savefig("lengths.pdf")
    
        
def classifier(data):
    from sklearn.covariance import EllipticEnvelope
    from sklearn.svm import OneClassSVM
    from sklearn.datasets import load_boston
    from sklearn import preprocessing
    # Get data

    # Define "classifiers" to be used
    legend1 = {}
    legend2 = {}
    evaluation = [[val["coverage"], val["num_exons"], val["distance_to_next"]] for val in data] 
    X = [[val["coverage"], val["num_exons"], val["distance_to_next"]] for val in data]  
    X = preprocessing.scale(X)
    evaluation = preprocessing.scale(evaluation)
    # Learn a frontier for outlier detection with several classifiers
    sample = random.sample(X, 20000)
    clf = OneClassSVM(nu=.1, kernel='rbf')
    test = random.sample(evaluation, 2000)
    print >> sys.stderr, "fitting data"    
    clf.fit(sample)
    print >> sys.stderr, "predicting data"
    Y = clf.predict(test)
    print >> sys.stderr, "plotting data"
    fig, axes = subplots()
    
    for i in range(len(test)):
        if Y[i] == 1:
            color = 'blue'
        else:
            color = 'red'
        axes.scatter(test[i][2], test[i][1], c=color)
    #ylim([50,2000]) #num exons
    ylabel("distance")
    #xlim([3,10])
    xlabel("coverage")
    savefig("DistanceVCoverage.pdf")

    fig, axes = subplots()
    """
    for i in range(len(test)):
        if Y[i] == 1:
            color = 'blue'
        else:
            color = 'red'
        axes.scatter(test[i][1], test[i][0], c=color)
    #xlim([0,10]) #num exons
    xlabel("number of exons")
    #ylim([3,15])
    ylabel("coverage")
    savefig("ExonsvsCoverage.pdf")
    """
    full_test = clf.predict(evaluation)
    novel, regular = [],[]
    for i in range(len(full_test)):
        result = full_test[i]
        if result == -1:
            print data[i]["id"]
            novel.append(data[i]["num_exons"])
        else:
            regular.append(data[i]["num_exons"])
    multi_exon_novel = [val for val in novel if val > 1]
    multi_exon_regular = [val for val in regular if val > 1]
    print >> sys.stderr, "novel, regular"
    print >> sys.stderr, len(novel), len(regular)
    print >> sys.stderr, mean(multi_exon_novel), mean(multi_exon_regular), len(multi_exon_novel), len(multi_exon_regular)
    
     

if __name__ == "__main__":
    main()
