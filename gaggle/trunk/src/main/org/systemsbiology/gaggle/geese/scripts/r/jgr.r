# jgr-init.r:  gaggle-related functions specific to the JGR goose -- that is, the gaggled
# version of the 'Java Gui for R' (JGR, or jaguar), whose specific gaggle-related functions
#--------------------------------------------------------------------------------------------------
print (.jcall("java/lang/System", "S", "getProperty", "os.name"))
goose = .jcall ("org/systemsbiology/gaggle/geese/jaguar/JaguarGoose", 
                "Lorg/systemsbiology/gaggle/geese/jaguar/JaguarGoose;", "getCurrentConsole")

print (.jcall (goose, "S", "getName"))
url = 'http://gaggle.systemsbiology.net/utilities.r'
cat (  'about to source ', url, '\n')
source (url)
#--------------------------------------------------------------------------------
scriptVersion <- function ()
{
  return ("jgr-init.r   $Revision: 148 $   $Date: 2005-08-15 17:11:29 -0400 (Mon, 15 Aug 2005) $");
}
#--------------------------------------------------------------------------------------------------
broadcast <- function (x, name=NULL) 
{
  if (is.matrix (x)) {
    variableName <- "matrixFromR"
    newMatrixName <- name
    assign (variableName, x, env = .GlobalEnv)
    rowNamesVar = paste (variableName, ".tmp.rownames", sep="")
    matrixNameVar = paste (variableName, ".newName", sep="")
    assign (rowNamesVar, rownames (x), env = .GlobalEnv)
    if (is.null (newMatrixName))
      name <- "from R"
    assign (matrixNameVar, newMatrixName,  env = .GlobalEnv)
    .jcall (goose, "V", "broadcast", variableName)
    }
  else if (is.vector (x)) {
    variableName <- "listFromR"
    assign (variableName, x, env = .GlobalEnv)
    .jcall (goose, "V", "broadcast", variableName)
    }
  else {
    cat ("no support yet for broadcasting variables of type ", typeof (x), "\n")
    }
  invisible (NULL)
  
}
#--------------------------------------------------------------------------------
newbroadcast <- function (x, y=NULL, name=NULL) 
{
  if (is.matrix (x)) {
    variableName <- "matrixFromR"
    newMatrixName <- name
    assign (variableName, x, env = .GlobalEnv)
    rowNamesVar = paste (variableName, ".tmp.rownames", sep="")
    matrixNameVar = paste (variableName, ".newName", sep="")
    assign (rowNamesVar, rownames (x), env = .GlobalEnv)
    if (is.null (newMatrixName))
      name <- "from R"
    assign (matrixNameVar, newMatrixName,  env = .GlobalEnv)
    .jcall (goose, "V", "broadcast", variableName)
    }
  else if (is.vector (x) && is.null (y)) {
    variableName <- "listFromR"
    assign (variableName, x, env = .GlobalEnv)
    .jcall (goose, "V", "broadcast", variableName)
    }
  else if (is.vector (x) && is.vector (y)) {
    assign ("rowNames", x, env = .GlobalEnv)
    assign ("colNames", y, env = .GlobalEnv)
    .jcall (goose, "V", "broadcast", "rowNames", "colNames")
    }
  else {
    cat ("no support yet for broadcasting variables of type ", typeof (x), "\n")
    }
  invisible (NULL)
  
}
#--------------------------------------------------------------------------------
