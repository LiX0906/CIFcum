
# PLOT ESTIMATED WEIGHT FUNCTION

par(mfrow=c(3,2),mar=c(5,5,3,2))
for (q in 1:4){
  plot(elapsed[1:lag],weight_mat[q,][1:lag],type="n",ylim=c(-0.025,0.038),
       ylab = 'Estimated weight function', xlab = 'Days before last follow-up', cex.lab=1.5, cex.axis=1.5,frame.plot=FALSE)
  for (i in 1:50){
    lines(output_rc_unconstrained[i,(120*(q-1)+1):(120*q)],lwd=1.5,col=grey(0.8))
  }
  lines(elapsed[1:lag],weight_mat[q,][1:lag],lwd=3,pch=19)
  legend("top",c("True","Estimated"),lty=c(1,1),cex=1.5,
         col=c(1,grey(0.8)),horiz=T,bty="n",x.intersp=0.1,inset=0.03)
}

for (q in 1:2){
  plot(elapsed[1:lag],weight_mat[(q+1),][1:lag],type="n",ylim=c(-0.02,0.05),
       ylab = 'Estimated weight function', xlab = 'Days before last follow-up', cex.lab=1.5, cex.axis=1.5,frame.plot=FALSE)
  for (i in 1:50){
    lines(output_rc_constrained[i,(120*(q-1)+1):(120*q)],lwd=1.5,col=grey(0.8))
  }
  lines(elapsed[1:lag],weight_mat[(q+1),][1:lag],lwd=3,pch=19)
  legend("top",c("True","Estimated"),lty=c(1,1),cex=1.5,
         col=c(1,grey(0.8)),horiz=T,bty="n",x.intersp=0.1,inset=0.03)
}







