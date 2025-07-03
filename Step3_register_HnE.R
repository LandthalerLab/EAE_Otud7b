#Step 3

options(rgl.useNULL=TRUE)

#Install Voltron from github
devtools::install_github("Artur-man/VoltRon")

library(Seurat)
library(VoltRon)


setwd("/path/to/working/directory/")


#Register H&E images on d0 object

all_d0 <- readRDS("all_d0_mols_celltypeann.rds")

#Register H&E images
#Check assay and sample names
vrImageChannelNames(all_d0)
all_d0@sample.metadata

#All images were registered manually

#Register d0_WT1
imgdata <- importImageData("/Volumes/Storage/MoreApplications/stuff/Xenium/Stainings/EAE_HE/d0_WT1.jpg")
d0_WT1 <- subset(all_d0, samples="d0_WT1")
#Invert/negate DAPI image
xen_reg <- registerSpatialData(object_list = list(d0_WT1, imgdata))
imgdata_reg <- xen_reg$registered_spat[[2]]
vrImages(all_d0[["Assay1"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
vrImages(all_d0[["Assay2"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
rm(d0_WT1)

#Register d0_WT2
imgdata <- importImageData("/Volumes/Storage/MoreApplications/stuff/Xenium/Stainings/EAE_HE/d0_WT2.jpg")
d0_WT2 <- subset(all_d0, samples="d0_WT2")
#Invert/negate DAPI image
xen_reg <- registerSpatialData(object_list = list(d0_WT2, imgdata))
imgdata_reg <- xen_reg$registered_spat[[2]]
vrImages(all_d0[["Assay3"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
vrImages(all_d0[["Assay4"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
rm(d0_WT2)

#Register d0_WT3
imgdata <- importImageData("/Volumes/Storage/MoreApplications/stuff/Xenium/Stainings/EAE_HE/d0_WT3.jpg")
d0_WT3 <- subset(all_d0, samples="d0_WT3")
#Invert/negate DAPI image
xen_reg <- registerSpatialData(object_list = list(d0_WT3, imgdata))
imgdata_reg <- xen_reg$registered_spat[[2]]
vrImages(all_d0[["Assay5"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
vrImages(all_d0[["Assay6"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
rm(d0_WT3)
saveRDS(all_d0, "all_d0_mols_reg.rds")

#Register d0_KO1
imgdata <- importImageData("/Volumes/Storage/MoreApplications/stuff/Xenium/Stainings/EAE_HE/d0_KO1_slice1.jpg")
d0_KO1_slice1 <- subset(all_d0, samples="d0_KO1_slice1")
#Invert/negate DAPI image
xen_reg <- registerSpatialData(object_list = list(d0_KO1_slice1, imgdata))
imgdata_reg <- xen_reg$registered_spat[[2]]
vrImages(all_d0[["Assay7"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
vrImages(all_d0[["Assay8"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
rm(d0_KO1_slice1)

#Register d0_KO2
imgdata <- importImageData("/Volumes/Storage/MoreApplications/stuff/Xenium/Stainings/EAE_HE/d0_KO2_slice2.jpg")
d0_KO2_slice2 <- subset(all_d0, samples="d0_KO2_slice2")
#Invert/negate DAPI image
xen_reg <- registerSpatialData(object_list = list(d0_KO2_slice2, imgdata))
imgdata_reg <- xen_reg$registered_spat[[2]]
vrImages(all_d0[["Assay9"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
vrImages(all_d0[["Assay10"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
rm(d0_KO2_slice2)

#Register d0_KO3
imgdata <- importImageData("/Volumes/Storage/MoreApplications/stuff/Xenium/Stainings/EAE_HE/d0_KO3.jpg")
d0_KO3 <- subset(all_d0, samples="d0_KO3")
#Invert/negate DAPI image
xen_reg <- registerSpatialData(object_list = list(d0_KO3, imgdata))
imgdata_reg <- xen_reg$registered_spat[[2]]
vrImages(all_d0[["Assay11"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
vrImages(all_d0[["Assay12"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
rm(d0_KO3)


saveRDS(all_d0, "all_d0_mols_celltypeann_reg.rds")


#Export images in lower resolution for having a quick look
for (the_assay in unique(vrImageChannelNames(all_d0)$Assay)) {
  magick::image_write(vrImages(all_d0[[the_assay]], name = "main", channel = "H&E", scale.perc = 50), paste0(the_assay, "_HE.jpg"))
  magick::image_write(vrImages(all_d0[[the_assay]], name = "main", channel = "DAPI", scale.perc = 50), paste0(the_assay, "_DAPI.jpg"))
}



all_d15 <- readRDS("all_d15_mols_celltypeann.rds")

#Register H&E images done on the slides following the Xenium run
#Check assay and sample names
vrImageChannelNames(all_d15)
all_d15@sample.metadata

#All images were registered manually, modulation was not done but would be possible with
#d15_WT2 <- modulateImage(d15_WT2, brightness = 300)

#Register d15_WT1
imgdata <- importImageData("/Volumes/Storage/MoreApplications/stuff/Xenium/Stainings/EAE_HE/d15_WT1.jpg")
d15_WT1 <- subset(all_d15, samples="d15_WT1")
#Invert/negate DAPI image
xen_reg <- registerSpatialData(object_list = list(d15_WT1, imgdata))
imgdata_reg <- xen_reg$registered_spat[[2]]
vrImages(all_d15[["Assay1"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
vrImages(all_d15[["Assay2"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
rm(d15_WT1)

#Register d15_WT2 - this had to be done manually
imgdata <- importImageData("/Volumes/Storage/MoreApplications/stuff/Xenium/Stainings/EAE_HE/d15_WT2.jpg")
d15_WT2 <- subset(all_d15, samples="d15_WT2")
#Invert/negate DAPI image
xen_reg <- registerSpatialData(object_list = list(d15_WT2, imgdata))
imgdata_reg <- xen_reg$registered_spat[[2]]
vrImages(all_d15[["Assay3"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
vrImages(all_d15[["Assay4"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
rm(d15_WT2)

#Register d15_WT3
imgdata <- importImageData("/Volumes/Storage/MoreApplications/stuff/Xenium/Stainings/EAE_HE/d15_WT3.jpg")
d15_WT3 <- subset(all_d15, samples="d15_WT3")
#Invert/negate DAPI image
xen_reg <- registerSpatialData(object_list = list(d15_WT3, imgdata))
imgdata_reg <- xen_reg$registered_spat[[2]]
vrImages(all_d15[["Assay5"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
vrImages(all_d15[["Assay6"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
rm(d15_WT3)

#Register d15_KO1
imgdata <- importImageData("/Volumes/Storage/MoreApplications/stuff/Xenium/Stainings/EAE_HE/d15_KO1.jpg")
d15_KO1 <- subset(all_d15, samples="d15_KO1")
#d15_WT3 <- modulateImage(d15_KO1, brightness = 300)
#Invert/negate DAPI image
xen_reg <- registerSpatialData(object_list = list(d15_KO1, imgdata))
imgdata_reg <- xen_reg$registered_spat[[2]]
vrImages(all_d15[["Assay7"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
vrImages(all_d15[["Assay8"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
rm(d15_KO1)

#Register d15_KO2
imgdata <- importImageData("/Volumes/Storage/MoreApplications/stuff/Xenium/Stainings/EAE_HE/d15_KO2.jpg")
d15_KO2 <- subset(all_d15, samples="d15_KO2")
#Invert/negate DAPI image
xen_reg <- registerSpatialData(object_list = list(d15_KO2, imgdata))
imgdata_reg <- xen_reg$registered_spat[[2]]
vrImages(all_d15[["Assay9"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
vrImages(all_d15[["Assay10"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
rm(d15_KO2)

#Register d15_KO3_slice2
imgdata <- importImageData("/Volumes/Storage/MoreApplications/stuff/Xenium/Stainings/EAE_HE/d15_KO3_slice2.jpg")
d15_KO3_slice2 <- subset(all_d15, samples="d15_KO3_slice2")
#Invert/negate DAPI image
xen_reg <- registerSpatialData(object_list = list(d15_KO3_slice2, imgdata))
imgdata_reg <- xen_reg$registered_spat[[2]]
vrImages(all_d15[["Assay11"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
vrImages(all_d15[["Assay12"]], name = "main", channel = "H&E") <-  vrImages(imgdata_reg, assay = "Assay1", name = "main_reg")
rm(d15_KO3_slice2)

saveRDS(all_d15, "all_d15_mols_celltypeann_reg.rds")

#Export images in lower resolution for having a quick look
for (the_assay in c("Assay1", "Assay2", "Assay3", "Assay4", "Assay5", "Assay6", "Assay7", "Assay8", "Assay9", "Assay10", "Assay11", "Assay12")) {
  magick::image_write(vrImages(all_d15[[the_assay]], name = "main", channel = "H&E", scale.perc = 20), paste0(the_assay, "_HE.jpg"))
  magick::image_write(vrImages(all_d15[[the_assay]], name = "main", channel = "DAPI", scale.perc = 20), paste0(the_assay, "_DAPI.jpg"))
}

