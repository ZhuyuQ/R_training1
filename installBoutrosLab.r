path <- "/cloud/project/"
input <- c("BoutrosLab.utilities_1.9.10.tar.gz",
           "BoutrosLab.dist.overload_1.0.2.tar.gz",
           "BoutrosLab.statistics.general_2.1.3.tar.gz",
           "BoutrosLab.statistics.survival_0.4.20.tar.gz",
           "BoutrosLab.prognosticsignature.general_1.3.13.tar.gz",
           "BoutrosLab.plotting.general_5.9.8.tar.gz",
           "BoutrosLab.plotting.survival_3.0.10.tar.gz"
);

for (i in 1:length(input)){
   install.packages(paste(path, input[i], sep="/"), repo=NULL, type="source", dependencies=TRUE);
   library(unlist(strsplit(input[i], "_"))[1],character.only=TRUE);
}
