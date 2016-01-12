works_with_R("3.2.3",
             "tdhock/namedCapture@05175927a45c301a18e8c6ebae67ea39a842d264",
             data.table="1.9.6")

TableS1 <- fread("TableS1.csv")

pattern1 <- paste0(
  "(?<before>[^/])",
  "/L ",
  "(?<treatment>[^ ]+)",
  " R",
  "(?<replicate>[12])",
  " +",
  "(?<hours>[^h]+)")
header1 <- str_match_named(names(TableS1)[-1], pattern1, list(
  hours=as.numeric,
  replicate=as.integer))

TableS2 <- fread("TableS2.csv")

pattern2 <- paste0(
  "R",
  "(?<replicate>[12])",
  " ",
  "(?<treatment>[^ ]+)",
  " ",
  "(?<hours>[^h]+)")
header2.short <- str_match_named(names(TableS2)[-1], pattern2, list(
  hours=as.numeric,
  replicate=as.integer))
header2 <- data.frame(before=NA, header2.short)

wide.data.list <- list(
  list(data=TableS1, header=header1),
  list(data=TableS2, header=header2))

RNAseq.list <- list()
for(file.i in seq_along(wide.data.list)){
  wide.data <- wide.data.list[[file.i]]
  for(col.i in 1:nrow(wide.data$header)){
    sub.dt <- wide.data$data[, c(1, col.i+1), with=FALSE]
    setnames(sub.dt, c("genes", "expression"))
    RNAseq.list[[paste(file.i, col.i)]] <- data.table(
      file.i, wide.data$header[col.i, ], sub.dt)
  }
}

RNAseq <- do.call(rbind, RNAseq.list)

save(RNAseq, file="RNAseq.RData")
