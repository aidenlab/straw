import straw_ext

bins = straw_ext.straw("NONE /Users/nchernia/Downloads/intra_nofrag.hic 1:1000000:7400000 1:1000000:7400000 BP 25000");
# the values returned are in bins.x / bins.y / bins.counts
for i in range(len(bins.x)):
    print bins.x[i], bins.y[i], bins.counts[i]
