library(tikzDevice)
library(ggplot2)
#For some reason, Rstudio needs to know the time zone...
options(tz="CA")
#Dummy data for the plot
y <- exp(seq(1,10,.1))
x <- 1:length(y)
data <- data.frame(x = x, y = y)

#Create a .tex file that will contain your plot as vectors
#You need to set the size of your plot here, if you do it in LaTeX, font consistency with the rest of the document will be lost
tikz(file = "plot_test.tex", width = 5, height = 5)
#Simple plot of the dummy data using LaTeX elements
plot <- ggplot(data, aes(x = x, y = y)) + 
  geom_line() +
  #Space does not appear after Latex
  ggtitle( paste("Fancy \\LaTeX ", "\\hspace{0.01cm} title")) +
  labs( x = "$x$ = Time", y = "$\\Phi$ = Innovation output") +
  theme_bw()
#This line is only necessary if you want to preview the plot right after compiling
print(plot)
#Necessary to close or the tikxDevice .tex file will not be written
dev.off()
