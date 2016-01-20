#OneClass SVM

The script here takes a gtf and trains a one-class SVM on length, coverage, and number of exons. The script is really rough since it was a onetime use thing. Right now it's configued to print out the outliers rather than the inliers since the dataset used had > 1 million transcripts that were mostly single exonic low coverage transcripts. By changing the printing clause from -1 to 1 (or zero, not sure what it outputs) it should get inliers. It also produces some plots, the exon plot is commented out, some work should be done to make this a different function.
