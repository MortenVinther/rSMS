
library(RTMB)
library(tidyverse)
library(ggpubr)

library(FishStomachs)

rsms.root.prog<-file.path(root.prog,'r_prog','rsms')
rsms.root.prog.area<-file.path(root.prog,'r_prog','rsms',modelArea)
# 
# if (modelArea==c('NorthSea','Baltic')[1]) {
# 
#   stomcon_model<-"stomcon_list_Boots_0500_haul_as_observed.dat"
# } else if (modelArea==c('NorthSea','Baltic')[2]) {
#   stom.input<-file.path(root,"rsms_Baltic","stom_input")
# 
# }
# 


source(file.path(rsms.root.prog,"various_functions.R"))
source(file.path(rsms.root.prog,"rsms.control.R"))
source(file.path(rsms.root.prog,"parameters.R"))
source(file.path(rsms.root.prog,"plotCompareRuns.R"))
source(file.path(rsms.root.prog,"make_rsms_data_function.R"))
source(file.path(rsms.root.prog,"pick_species.R"))
source(file.path(rsms.root.prog,"quarter2annual.R"))
source(file.path(rsms.root.prog,"rsms_function.R"))
#source(file.path(rsms.root.prog,"summary_plot.R"))
source(file.path(rsms.root.prog,"map_param.R"))
source(file.path(rsms.root.prog,"lowerUpper.R"))
source(file.path(rsms.root.prog,"residuals_catch_survey.R"))

source(file.path(rsms.root.prog.area,"batch_control.R"))

