# Load a CSV file from a local folder
data <- read.csv("C:/Users/YourUsername/Documents/project/data/database.csv")

# Load an Excel file
library(readxl)
data <- read_excel("C:/Users/YourUsername/Documents/project/data/database.xlsx")

library(PsychoMatic)

result <- efa_auto(data,
                   rotation = "oblique", #You can change the rotation to “orthogonal” if needed
                   language = "eng") #For Spanish, enter “esp”
result

# Export results

# Export to Excel
export_efa(result, format = "excel", file_name = "results_afe")
# Export to Word
export_efa(result, format = "word", file_name = "results_afe")