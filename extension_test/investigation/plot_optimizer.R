# Summarize retained controlled fits without selecting or replacing any fit
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else 'investigation')
paths <- list.files(file.path(root,'optimizer/checkpoints'),pattern='\\.rds$',full.names=TRUE)
records <- lapply(paths,readRDS)
records <- Filter(function(x)isTRUE(x$success),records)
ids <- c('strength_case_0003_Melding','clustering_case_0003_Melding','strength_case_0003_PSmelding')
settings <- c('original','default_repeat_1','default_repeat_2','tight_saved',
              'tight_area_0p001','tight_area_0p1','tight_saved_repeat')
labels <- c('Original','Default repeat 1','Default repeat 2','Tight: saved mode',
            'Tight: area = 0.001','Tight: area = 0.1','Tight: saved repeat')
rows <- list()
for (phase in c('pilot','finer_mesh')) for(id in ids) {
  r <- Filter(function(x)x$job$phase==phase && x$job$job_id==id,records)
  if(!length(r)) next
  original <- r[[1]]$original_parameters
  for(setting in c('original',vapply(r,function(x)x$job$setting,''))) {
    item <- if(setting=='original') r[[1]] else r[[which(vapply(r,function(x)x$job$setting==setting,logical(1)))]]
    p <- if(setting=='original')original else item$parameters
    p <- p[p$parameter=='tau2_B',]
    rows[[length(rows)+1L]] <- data.frame(phase=phase,job_id=id,setting=setting,
      mean=p$posterior_mean,median=p$median,lower=p$lower,upper=p$upper,
      mode_status=if(setting=='original')item$original_mode$mode.status else item$mode$mode.status)
  }
}
table <- do.call(rbind,rows)
write.csv(table,file.path(root,'reports/areal_variance_comparison.csv'),row.names=FALSE)
png(file.path(root,'reports/optimizer_areal_variance.png'),width=2100,height=1500,res=180)
par(mfrow=c(2,3),mar=c(4.4,7.9,3,1),oma=c(0,0,3,0))
for(phase in c('pilot','finer_mesh')) for(id in ids) {
  p <- table[table$phase==phase & table$job_id==id,]
  p <- p[order(match(p$setting,settings)),]
  y <- rev(seq_len(nrow(p)))
  colour <- ifelse(p$setting=='original','#2166ac',ifelse(p$mode_status!=0,'#b2182b','#222222'))
  plot(p$median,y,log='x',xlim=c(1e-5,2),ylim=c(.5,nrow(p)+.5),yaxt='n',
    xlab='Areal variance (log scale)',ylab='',pch=16,col=colour,
    main=paste(if(phase=='pilot')'Production mesh'else'Finer mesh',
    c('Baseline Melding','Delta=2 Melding','Baseline PSmelding')[match(id,ids)],sep='\n'))
  axis(2,at=y,labels=labels[match(p$setting,settings)],las=1,cex.axis=.75)
  segments(p$lower,y,p$upper,y,col=colour,lwd=1.5)
  points(p$median,y,pch=16,col=colour)
  points(p$mean,y,pch=4,col=colour)
  abline(v=.1,lty=2,col='grey50')
}
mtext('Dot: posterior median   Cross: mean   Line: 95% interval   Dashed: generating variance 0.1',
      outer=TRUE,line=1.6,cex=.85)
mtext('Blue: original   Red: nonzero mode status', outer=TRUE,line=.3,cex=.85)
dev.off()
