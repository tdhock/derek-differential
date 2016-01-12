figure-RNAseq/index.html: figure-RNAseq.R RNAseq.RData
	R --no-save < $<
RNAseq.RData: RNAseq.R
	R --no-save < $<