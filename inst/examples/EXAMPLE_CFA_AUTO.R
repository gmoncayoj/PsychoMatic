###############################################################################################
##                                                                                           ##
##  PsychoMatic can automatically detect which model you are about to run.                   ##
##  If it detects a second-order model, it will report the higher-order omega reliability.   ##
##  If it detects a two-factor model, it will also calculate the hierarchical omega,         ##
##  ECV, and PUC indices, and provide you with a final interpretation of the strength        ##
##  of the general factor!                                                                   ##
##                                                                                           ##
###############################################################################################
###############################################################################################


# Load a CSV file from a local folder
data <- read.csv("C:/Users/YourUsername/Documents/project/data/database.csv")

# Load an Excel file
library(readxl)
data <- read_excel("C:/Users/YourUsername/Documents/project/data/database.xlsx")


library(PsychoMatic)

# Define the model in Lavaan syntax
# Unidimensional
model <- 'F =~ i1 + i2 + i3 + i4 + i5'

# Correlated factors
model <- '
  F1 =~ i1 + i2 + i3
  F2 =~ i4 + i5 + i6
'

# Second order
model <- '
  F1 =~ i1 + i2 + i3
  F2 =~ i4 + i5 + i6
  SO  =~ F1 + F2
'

# Bifactor
model <- '
  GEN  =~ i1 + i2 + i3 + i4 + i5 + i6
  F1 =~ i1 + i2 + i3
  F2 =~ i4 + i5 + i6
'

# Run the automated AFC
result <- cfa_auto(data = data, model = model)

# View the complete psychometric report
print(result)