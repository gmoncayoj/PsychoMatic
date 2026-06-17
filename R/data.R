#' Simulated continuous item responses
#'
#' A small continuous dataset for examples, tests, and vignettes. The variables
#' approximate six positively related psychometric items generated from one
#' common factor with residual noise.
#'
#' @format A data frame with 120 rows and 6 variables:
#' \describe{
#'   \item{item1,item2,item3,item4,item5,item6}{Continuous item responses.}
#' }
"psychomatic_continuous"

#' Simulated ordinal item responses
#'
#' A small ordinal Likert-type dataset for reproducible examples. The variables
#' are five-category transformations of the continuous example items.
#'
#' @format A data frame with 120 rows and 6 variables:
#' \describe{
#'   \item{ord1,ord2,ord3,ord4,ord5,ord6}{Ordinal item responses from 1 to 5.}
#' }
"psychomatic_ordinal"

#' Simulated multigroup item responses
#'
#' A small multigroup dataset for measurement-invariance examples. It contains
#' two groups and six continuous item responses generated from a common latent
#' variable with a small group mean difference.
#'
#' @format A data frame with 140 rows and 7 variables:
#' \describe{
#'   \item{mg1,mg2,mg3,mg4,mg5,mg6}{Continuous item responses.}
#'   \item{group}{Grouping factor with levels `A` and `B`.}
#' }
"psychomatic_multigroup"
