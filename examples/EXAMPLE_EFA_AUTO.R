# Load a CSV file from a local folder
data <- read.csv("C:/Users/YourUsername/Documents/project/data/database.csv")

# Load an Excel file
library(readxl)
data <- read_excel("C:/Users/YourUsername/Documents/project/data/database.xlsx")

library(PsychoMatic)

result <- efa_auto(data,
                   rotacion = "oblicua", #You can change the rotation to “orthogonal” if needed
                   language = "eng") #For Spanish, enter “esp”
result

# Export results

# Export to Excel
exportar_efa(result, formato = "excel", archivo = "results_afe")
# Export to Word
exportar_efa(result, formato = "word", archivo = "results_afe")