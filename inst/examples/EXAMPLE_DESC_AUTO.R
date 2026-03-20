##########################################################################
##                                                                      ##
## Using desc_auto(), you can calculate the mean, standard deviation,   ##
## skewness, kurtosis, and response rate per item.                      ##
##                                                                      ##
##########################################################################

# Load a CSV file from a local folder
data <- read.csv("C:/Users/YourUsername/Documents/project/data/datadata.csv")

# Load an Excel file
library(readxl)
data <- read_excel("C:/Users/YourUsername/Documents/project/data/datadata.xlsx")

library(PsychoMatic)

### Basic use (Spanish, no report)
desc_auto(data)

### With English labels
desc_auto(data, language = "eng") # If you want Spanish, enter “esp”

### Adjust decimal places
desc_auto(data, digits = 2)

### Export report to Excel
desc_auto(data, report = TRUE, language = "eng")