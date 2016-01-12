works_with_R("3.2.3",
             "tdhock/animint@3b1f84ec926ffbd765f0aa004596e43203750fd4",
             "tdhock/ggplot2@a8b06ddb680acdcdbd927773b1011c562134e4d2",
             data.table="1.9.6")

load("RNAseq.RData")

RNAseq[, table(genes, file.i)]

RNAseq[, code := ifelse(is.na(before),
                "NA (TableS2)",
                paste0(before, "/L (TableS1)"))]

hist(RNAseq$expression)
min(RNAseq$expression, na.rm=TRUE)
RNAseq[expression==0, ]

hist(log10(RNAseq$expression), breaks=100)
q.vec <- quantile(RNAseq$expression, seq(0, 1, by=0.0002), na.rm=TRUE)
hist(RNAseq$expression, breaks=q.vec)
hist(log(RNAseq$expression), breaks=log(q.vec)[-1])

RNAseq[, log.expr := log10(ifelse(expression < q.vec[2], 0, expression))]
hist(RNAseq$log.expr, breaks=100)

all.genes.vec <- unique(RNAseq$genes)
some.genes.vec <- head(all.genes.vec)

some.RNAseq <- RNAseq[genes %in% some.genes.vec,]

breaks.vec <- unique(RNAseq$hours, na.rm=TRUE)

ggplot()+
  geom_histogram(aes(log.expr), data=RNAseq)

fmean <- function(x, fun=mean){
  f <- x[is.finite(x)]
  fun(f)
}

gene.ranges <- RNAseq[, {
  data.table(min.expr=fmean(log.expr, min),
             max.expr=fmean(log.expr, max))
}, by=genes]

ggplot()+
  geom_segment(aes(genes, min.expr,
                   xend=genes, yend=max.expr),
               data=gene.ranges)

ggplot()+
  geom_point(aes(min.expr, max.expr),
               data=gene.ranges)

ggplot()+
  geom_line(aes(hours, log.expr, color=code,
                linetype=treatment,
                group=interaction(replicate, treatment, code)),
            data=some.RNAseq)+
  geom_point(aes(hours, log.expr, color=code),
             data=some.RNAseq)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  scale_x_log10(breaks=breaks.vec[-1])+
  facet_grid(genes ~ file.i)

scatter.dt <- RNAseq[, {
  LPS <- .SD[treatment=="LPS",]
  MOCK <- .SD[treatment=="MOCK",]
  stopifnot(nrow(LPS) == nrow(MOCK))
  first.LPS <- LPS[hours==1,]
  first.MOCK <- MOCK[hours==1,]
  last.LPS <- LPS[hours==12,]
  last.MOCK <- MOCK[hours==12,]
  data.table(h1.diff=fmean(first.LPS$log.expr)-fmean(first.MOCK$log.expr),
             h12.diff=fmean(last.LPS$log.expr)-fmean(last.MOCK$log.expr))
}, by=.(genes, code)]

equality.df <- data.frame(intercept=0, slope=1)

ggplot()+
  geom_abline(aes(intercept=intercept, slope=slope),
              data=equality.df)+
  coord_equal()+
  geom_point(aes(h1.diff, h12.diff, color=code, clickSelects=genes),
             data=scatter.dt)

scatter.tall.dt <- RNAseq[hours %in% c(1, 12), {
  LPS <- .SD[treatment=="LPS",]
  MOCK <- .SD[treatment=="MOCK",]
  stopifnot(nrow(LPS) == nrow(MOCK))
  data.table(LPS=fmean(LPS$log.expr),
             MOCK=fmean(MOCK$log.expr))
}, by=.(genes, code, hours)]

ggplot()+
  geom_abline(aes(intercept=intercept, slope=slope),
              data=equality.df)+
  coord_equal()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(. ~ hours)+
  geom_point(aes(LPS, MOCK, color=code, clickSelects=genes),
             data=scatter.tall.dt)

not.na <- RNAseq[!is.na(log.expr), ]
text.dt <- not.na[, {
  h <- max(hours)
  sub.dt <- .SD[hours==h,]
  data.table(hours=h, log.expr=mean(sub.dt$log.expr))
}, by=.(code, treatment, genes)]

viz <- list(
  title="RNA-seq data viz",
  ranges=ggplot()+
    ggtitle("select genes")+
    theme_bw()+
    xlab("min(log expression)")+
    ylab("max(log expression)")+
    geom_point(aes(min.expr, max.expr, clickSelects=genes),
               alpha=0.6,
               data=gene.ranges)+
    geom_text(aes(min.expr, max.expr,
                  showSelected=genes,
                  clickSelects=genes,
                  label=genes),
              data=gene.ranges),
  timeseries=ggplot()+
    ggtitle("time series for selected genes")+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "cm"))+
    theme_animint(width=500)+
    geom_line(aes(hours, log.expr, color=code,
                  linetype=treatment,
                  clickSelects=genes,
                  showSelected=genes,
                  group=interaction(replicate, treatment, code, genes)),
              data=not.na)+
    geom_point(aes(hours, log.expr, color=code,
                   clickSelects=genes,
                   showSelected=genes),
               data=not.na)+
    geom_text(aes(hours, log.expr, color=code, label=genes,
                  clickSelects=genes,
                  showSelected=genes),
               hjust=0,
               data=text.dt)+
    scale_x_log10(breaks=breaks.vec[-1]),
  first=list(genes=c("Ccl5", "Tmsb4x")),
  selector.types=list(genes="multiple"))

animint2dir(viz, "figure-RNAseq")
