
# First, let us do the analysis manually, from scratch
# ----------------------------------------------------
#  
library(CellNOptR)
args <- commandArgs(TRUE)
warnings()

# The one step version --------------------------------
#dataToy<-readMIDAS("MD-p53.csv")
dataToy<-readMIDAS(args[1])
CNOlistToy<-makeCNOlist(dataToy,subfield=FALSE)
ToyModel<-readSIF(args[2])
#ToyModel<-readSIF("all_pathways.sif")
#ToyModel<-readSIF("pathway_egfr_signaling.sif")
#ToyModel<-readSIF("pathway_non_small_cell_lung.sif")

res <- CNORwrap(paramsList=NA, name="Toy", namesData=list(CNOlist="ToyData",model="ToyModel"),data=CNOlistToy, model=ToyModel)

pList = defaultParameters(CNOlistToy, ToyModel)
res <- CNORwrap( paramsList=pList,    name="Toy1Step",    namesData=list(CNOlist="ToyData",model="ToyModel"),    data=CNOlistToy,    model=ToyModel)

