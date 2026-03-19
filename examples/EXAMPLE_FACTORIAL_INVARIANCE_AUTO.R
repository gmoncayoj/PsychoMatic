################################################################################
################################################################################
##                                                                            ##
##   When estimator is set to “auto”, the function automatically selects      ##
##   between WLSMV, ML, and MLR based on the number of response categories    ##
##   and multivariate normality (Mardia's test).                              ##
##                                                                            ##
################################################################################
################################################################################

# Load a CSV file from a local folder
data <- read.csv("C:/Users/YourUsername/Documents/project/data/database.csv")

# Load an Excel file
library(readxl)
data <- read_excel("C:/Users/YourUsername/Documents/project/data/database.xlsx")


library(PsychoMatic)

# Define your CFA model in lavaan syntax
model <- "
  F1 =~ item1 + item2 + item3
  F2 =~ item4 + item5 + item6
"

# Run — automatic estimator selection, output in Spanish (default)
result <- factorial_invariance_auto(
  data  = data,
  group = "group_var",
  model = model
)

# Common variants

# English/spanish output + save Excel report to working directory
result <- factorial_invariance_auto(
  data     = data,
  group    = "group_var",
  model    = model,
  language = "eng", #If you need the Spanish language, change “esp”
  report   = TRUE
)

# Force a specific estimator
result <- factorial_invariance_auto(
  data      = data,
  group     = "group_var",
  model     = model,
  estimator = "MLR"
)

# Save report to a custom path
result <- factorial_invariance_auto(
  data   = data,
  group  = "group_var",
  model  = model,
  report = "C:/results/invariance_report.xlsx"
)
