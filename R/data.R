#' US Health Insurance Data
#'
#' A dataset with 1338 observations on health insurance charges and demographic factors.
#'
#' @format A data frame with 1338 rows and 7 variables:
#' \describe{
#' \item{age}{age of primary beneficiary in years}
#' \item{sex}{insurance contractor sex, male/female}
#' \item{bmi}{body mass index}
#' \item{children}{number of dependents}
#' \item{smoker}{smoker/non-smoker}
#' \item{region}{beneficiary's region of US}
#' \item{charges}{individual medical costs billed by health insurance}
#' }
#'
#' @source \url{https://www.kaggle.com/datasets/teertha/ushealthinsurancedataset}
"insurance"


#' Depression Treatment Data
#'
#' A longitudinal dataset comparing two treatments for depression.
#'
#' @format A data frame with 1020 rows and 5 variables:
#' \describe{
#' \item{diagnose}{diagnosed depression severity}
#' \item{drug}{treatment; standard or new}
#' \item{id}{patient id}
#' \item{time}{time point of treatment}
#' \item{depression}{depression response at time of treatment. 1 = Normal, 0 = Abnormal}
#' }
#'
#' @source \url{http://static.lib.virginia.edu/statlab/materials/data/depression.csv}
#' @references Agresti, A. (2002). Categorical Data Analysis. Wiley, 2nd Edition.
"depression"
