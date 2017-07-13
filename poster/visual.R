devtools::install_github("franciscorichter/dmea")
library(dmea)
s = sim_phyl()
par(mfrow=(c(1,3))) # 1 row but 2 plot panels
color = 'darkgreen'
plot(s$newick,direction = 'upwards',show.tip.label=FALSE,edge.width=7,edge.color = color,type="fan", no.margin = TRUE)
# this arrows are not working yet, just spend 5 minuts and solve it..
plot(1:10,type="n", yaxt="n",xaxt="n",xlab="", ylab="", bty="n")
arrows(x0=8, y0=7, x1 = 3, y1 = 7, lwd = 5,col = 'darkblue')
arrows(x0=3, y0=4, x1 = 8, y1 = 4, lwd = 5,col = 'darkblue')
dropex <- drop.fossil(s$newick)
plot(dropex,direction = 'upwards',show.tip.label=FALSE,edge.width=7,edge.color = color,type="fan", no.margin = TRUE)

