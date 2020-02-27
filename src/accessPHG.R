library("rJava")
library("rTASSEL")
.jaddClassPath("/Users/peterbradbury/git/practicalhaplotypegraph/dist/phg.jar")

#use MethodTableReportPlugin to get the methods table from a database
config = "/Users/peterbradbury/temp/phgSmallSeq/data/configSQLite.txt"
plugin <- new (J("net/maizegenetics/pangenome/api/MethodTableReportPlugin"))
plugin <- plugin$configFile(config)
plugin$configFile()

#retrieve the methods table as a TableReport
ds <- plugin$performFunction(.jnull("net/maizegenetics/plugindef/DataSet"))
datum <- ds$getData(0L)
tr <- .jcast(datum$getData(), new.class = "net/maizegenetics/util/TableReport")

#convert the TableReport into an object that can be used by R
resultVectors <- J("net/maizegenetics/plugindef/GenerateRCode", "tableReportToVectors", tr)
data <- resultVectors$dataVector
dfmethods = data.frame(data$get(0L),data$get(1L),data$get(2L),data$get(3L),data$get(4L))
names(dfmethods) = resultVectors$columnNames

#maybe to use kotlin but does not work so far
.jaddClassPath("/Users/peterbradbury/git/tassel-5-standalone/lib/kotlin-stdlib-1.3.10.jar")

#build a graph
hgbp <- new(J("net/maizegenetics/pangenome/api/HaplotypeGraphBuilderPlugin"), .jnull("java/awt/Frame"), FALSE)
names(hgbp)
hgbp$getButtonName()
hgbp$configFile(config)
hgbp$configFile()
hgbp$methods("GATK_PIPELINE")
hgbp$methods()
phg <- hgbp$build()

#test phg RMethods
library("rJava")
library("rTASSEL")
.jaddClassPath("/Users/peterbradbury/git/practicalhaplotypegraph/dist/phg.jar")
rTASSEL::startLogger(fullPath = NULL, fileName = NULL)

#test data frame formatted import
df.source <- J("net.maizegenetics.pangenome.api/RMethods", "testDataFrame")
df.source$columnNames
df.source$rowNames
df <- data.frame(df.source$dataVectors$get(as.integer(0)), df.source$dataVectors$get(as.integer(1)), df.source$dataVectors$get(as.integer(2)))
if (!is.jnull(df.source$columnNames)) names(df) <- df.source$columnNames
if (!is.jnull(df.source$rowNames)) row.names(df) <- df.source$rowNames

#test matrix formatted import
matrix.source <- J("net.maizegenetics.pangenome.api/RMethods", "testMatrix")
matrix.test = matrix.source$matrix
if(!is.jnull(matrix.source$rowNames)) rownames(matrix.test) = matrix.source$rowNames
if(!is.jnull(matrix.source$columnNames)) colnames(matrix.test) = matrix.source$columnNames

#test phg methods
config = "/Users/peterbradbury/temp/phgSmallSeq/data/configSQLite.txt"
#create a graph builder object and set the required parameters
hgbp <- new(J("net/maizegenetics/pangenome/api/HaplotypeGraphBuilderPlugin"), .jnull("java/awt/Frame"), FALSE)
hgbp$configFile(config)
hgbp$methods("CONSENSUS")

#build the graph
phg <- hgbp$build()

#haplotype id matrix
hapids <- J("net.maizegenetics.pangenome.api/RMethods", "hapidTableAsMatrix", phg)
hapid.matrix = hapids$matrix
if(!is.jnull(hapids$rowNames)) rownames(hapid.matrix) = hapids$rowNames
if(!is.jnull(hapids$columnNames)) colnames(hapid.matrix) = hapids$columnNames

#reference range table
rr <- J("net.maizegenetics.pangenome.api/RMethods", "referenceRanges", phg)
refranges <- data.frame(rr$dataVectors$get(0L), rr$dataVectors$get(1L),rr$dataVectors$get(2L),rr$dataVectors$get(3L))
names(refranges) = rr$columnNames

#path hapids
#get some files
path.filenames <- list.files('temp/phgSmallSeq/hap_count_best_path/')
directory.name = '/Users/peterbradbury/temp/phgSmallSeq/hap_count_best_path/'
path.files = paste(directory.name, path.filenames, sep="")
jfiles <- .jarray(path.files)
path.hapids <- J("net.maizegenetics.pangenome.api/RMethods", "pathHapids", config, jfiles)
path.matrix <- path.hapids$matrix
if(!is.jnull(path.hapids$rowNames)) rownames(path.matrix) <- path.hapids$rowNames
if(!is.jnull(path.hapids$columnNames)) colnames(path.matrix) = path.hapids$columnNames

