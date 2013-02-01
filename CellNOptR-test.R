
# First, let us do the analysis manually, from scratch
# ----------------------------------------------------
#  
library(CellNOptR)
warnings()

# The one step version --------------------------------
dataToy<-readMIDAS("MD-p53.csv")
CNOlistToy<-makeCNOlist(dataToy,subfield=FALSE)
ToyModel<-readSIF("p53_feedback_loops.sif")

res <- CNORwrap(paramsList=NA, name="Toy", namesData=list(CNOlist="ToyData",model="ToyModel"),data=CNOlistToy, model=ToyModel)

pList = defaultParameters(CNOlistToy, ToyModel)
res <- CNORwrap( paramsList=pList,    name="Toy1Step",    namesData=list(CNOlist="ToyData",model="ToyModel"),    data=CNOlistToy,    model=ToyModel)

