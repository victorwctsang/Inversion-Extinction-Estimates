# library(tikzDevice)
library(ggplot2)
#For some reason, Rstudio needs to know the time zone...
options(tz="AEST")

#Create a .tex file that will contain your plot as vectors
#You need to set the size of your plot here, if you do it in LaTeX, font consistency with the rest of the document will be lost
# tikz(file = "plot_test.tex", width = 5, height = 5)
#Simple plot of the dummy data using LaTeX elements
plot <- ggplot() + 
  geom_function(fun = function(e) 10 - e) +
  xlim(-10, 10) +
  theme_classic()
#This line is only necessary if you want to preview the plot right after compiling
print(plot)
#Necessary to close or the tikxDevice .tex file will not be written
# dev.off()
