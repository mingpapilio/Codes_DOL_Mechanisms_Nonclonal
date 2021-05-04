library(lattice)
library(readr)
###
label_unequal<- c("#FFFFFF", "#E2E2E2", "#E2E2E2", "#E2E2E2", "#C6C6C6", "#C6C6C6", "#C6C6C6", "#C6C6C6", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#909090", "#909090", "#909090", "#909090", "#767676", "#767676", "#767676", "#767676", "#5E5E5E", "#5E5E5E", "#5E5E5E", "#5E5E5E", "#464646", "#464646", "#464646", "#464646", "#303030", "#303030", "#303030", "#303030", "#1B1B1B", "#1B1B1B", "#1B1B1B", "#1B1B1B", "#000000", "#000000", "#000000", "#000000")
#
label_unequal_c<- c("#26548A", "#4364A5", "#4364A5", "#4364A5", "#677195", "#677195", "#677195", "#677195", "#8B7E7F", "#8B7E7F", "#8B7E7F", "#8B7E7F", "#AF8B6A", "#AF8B6A", "#AF8B6A", "#AF8B6A", "#D39855", "#D39855", "#D39855", "#D39855", "#E8A553", "#E8A553", "#E8A553", "#E8A553", "#EEB465", "#EEB465", "#EEB465", "#EEB465", "#F4C277", "#F4C277", "#F4C277", "#F4C277", "#FAD189", "#FAD189", "#FAD189", "#FAD189", "#FFE09C", "#FFE09C", "#FFE09C", "#FFE09C")
#
raw<- read_csv(" level_summary.csv")
y_grid<- seq(4, 40, by= 4)
x_grid<- seq(0.5, 1.0, length.out= 10)
grid<- expand.grid(X= x_grid, Y= y_grid)
#
dev.new()
levelplot(raw$avg_porp~raw$essentiality*raw$rp_size, grid, col.regions=label_unequal, at= seq(0, 1, by=0.025), 
          xlab="Trait essentiality, e", ylab="Group size, n")
#
dev.new()
levelplot(raw$avg_s~raw$essentiality*raw$grp_size, grid, col.regions=label_unequal_c, at= seq(0, 1, by=0.025), 
          xlab="Trait essentiality, e", ylab="Group size, n")
#
