###Below is the code for the soma


#Required packages
library(spatstat)


###CHANGE THE BELOW TO THE CORRECT PATH#######
#read in the mask file - read in every mask file, give it a new name
neuron_mask_1 <- read.table(file="E:/Libraries/Documents/Test/MaskFiles/VAMP2pHluorin_488_wt_4_mask_file.csv",header=F, sep=",")
neuron_mask_2 <- read.table(file="E:/Libraries/Documents/SCIENCE_GUPTON/ANALYSE_THESE/5_15_16/WT/WT/Stream_3_mask_file_soma.csv",header=F, sep=",")


#read in the x,y,t file - different variable for each
neuron_datapoints_1 <- read.table(file = "E:/Libraries/Documents/Test/DataFiles/Stream_4_fusion_stats.csv", header = T, sep=",")
neuron_datapoints_2 <- read.table(file = "E:/Libraries/Documents/SCIENCE_GUPTON/ANALYSE_THESE/5_15_16/WT/WT/stream_3_pointsxyt.csv", header = T, sep=",")




#Convert the mask file to logic - do this for each mask.
neuron_temp_1 <- (neuron_mask_1 > 0.5)
neuron_mask_log_1<- unname(neuron_temp_1, force = TRUE)

neuron_temp_2 <- (neuron_mask_2 > 0.5)
neuron_mask_log_2<- unname(neuron_temp_2, force = TRUE)

#extract the x, y from the original data file - feel free to give it whatever name you want. Must do for each.
#Cell 1

x_pos_1 = na.omit(as.numeric(neuron_datapoints_1$x_pos))
y_pos_1 = na.omit(as.numeric(neuron_datapoints_1$y_pos))

#Cell 2

x_pos_2 = na.omit(as.numeric(neuron_datapoints_2$x))
y_pos_2 = na.omit(as.numeric(neuron_datapoints_2$y))

#Convert our x,y positions into a point process format
pp_neuron_1<- ppp(x_pos_1,y_pos_1, mask = neuron_mask_log_1)
pp_neuron_2<- ppp(x_pos_2,y_pos_2, mask = neuron_mask_log_2)

#pp_neuron_2<- ppx(t_pos_1)
#Lest(pp_neuron_2)

#Plot a heat map of the density of points. 0.4 is what I use, feel free to change it BUT keep it CONSISTENT if comparing 
#between two density plots
plot(density(pp_neuron_1,20))
plot(density(pp_neuron_2,10))

#Run the Ripley's K function for the envelope functions to get an idea of CSR
#NOTE: Pool function now calculates the error
#E1 <- envelope(ppX, Lest, nsim = 39, savefuns = TRUE)
#E2 <- envelope(ppY, Lest, nsim = 39, savefuns = TRUE)

#Here, run Ripley's K function for each individual dataset
K1 <- Lest(pp_neuron_1, ratio = TRUE)
K2 <- Lest(pp_neuron_2, ratio = TRUE)



#We can plot the curves to look at them. Not necessary step
#plot(K1, cbind(trans,theo)~r)
#plot(K2, cbind(trans,theo)~r)

#Next, we'll pool all of the samples together into one Ripley's per group. Pool as many groups as you have,
#for example, if you have 10 cells WT and 12 cells TRIM67-/-, pool the WT 10 cells into one variable,
#and the TRIM67-/- in another.
wt_set_pool <- pool(K1,K2)

#Plot the graph. Here, pooltrans is the mean dataline, pooltheo is the theoretically random,
#and hitrans/lotrans are the regions of error around the mean
plot(wt_set_pool, cbind(pooltrans,pooltheo) ~ r, shade=c("hitrans", "lotrans"), xlim=c(0,120))
