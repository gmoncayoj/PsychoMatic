# Load a CSV file from a local folder
data <- read.csv("C:/Users/YourUsername/Documents/project/data/database.csv")

# Load an Excel file
library(readxl)
data <- read_excel("C:/Users/YourUsername/Documents/project/data/database.xlsx")


library(PsychoMatic)

## Option 1 — From raw item data and a grouping variable

# items_df: data frame with numeric item columns only
# group_var: vector indicating group membership (one value per row)

res <- inv_align_auto(
  data     = data,
  group    = group_var,
  language = "eng"
)

print(res)       # Full alignment report
summary(res)     # Summary table (R2, % noninvariant, invariance verdict)

## Option 2 — From pre-estimated λ and ν matrices

# lambda and nu: groups × items numeric matrices
# (rows = groups, columns = items)

res <- inv_align_auto(
  lambda   = lambda_mat,
  nu       = nu_mat,
  language = "eng"
)

print(res)