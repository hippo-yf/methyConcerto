

library(readr)
library(matrixStats)
library(argparser, quietly=TRUE)
library(assertthat)

## parse arguments

arg <- arg_parser("Bayesian Methylation Caller (BMC)")

arg <- add_argument(arg, "--CG", help="CG site or nonCG site", type="logical", default = TRUE)
arg <- add_argument(arg, "--input", help="input .CGmap file", type="character")
arg <- add_argument(arg, "--output", help="output .prob (probability) file", type="character")
arg <- add_argument(arg, "--summary", help="summary file", type="character")
arg <- add_argument(arg, "--prior", help="prior probabilities overriding the CG flag", 
                    type="numeric", nargs = 3, default = c(0.59, 0.4, 0.01)
)
arg <- add_argument(arg, "--error_rate", help="prior rate of C/T mis-detction", 
                    type="numeric", default = 0.01
)
arg <- add_argument(arg, "--Poisson_approx", help="threshold of depth the Poisson approximation applied", 
                    type="numeric", default = 20
)
arg <- add_argument(arg, "--rounds", help="rounds of iteration", 
                    type="integer", default = 2
)

argv <- parse_args(arg)


assert_that(!anyNA(argv$output))
assert_that(!anyNA(argv$CG))
assert_that(!anyNA(argv$prior))


if (length(argv$prior) == 3 & (!anyNA(argv$prior))) {
  r = argv$prior
} else if(CG){
  r = c(0.59, 0.4, 0.01)
}else{
  r = c(0.01, 0.98, 0.01)
}

# input: stdio or file

if(anyNA(argv$input)){
  # read counts from stdio
  count = base::as.matrix(readr::read_delim(stdin(), 
                                            delim = '\t', 
                                            col_types = cols(col_integer(), col_integer()),
                                            col_names = FALSE)
                          )
}else{
  count = base::as.matrix(readr::read_delim(argv$input, 
                                            delim = '\t', 
                                            col_types = cols(col_integer(), col_integer()),
                                            col_names = FALSE)
  )
}


# count = as.matrix(readr::read_delim('AF1_13.CGmap.chr.CG', delim = '\t', col_names = FALSE))
# count = readr::read_delim('AF1_13.CGmap.chr.CG',
#                           col_types = cols(
#                             X1 = col_character(),
#                             X2 = col_character(),
#                             X3 = col_double(),
#                             X4 = col_character(),
#                             X5 = col_character(),
#                             X6 = col_double(),
#                             X7 = col_double(),
#                             X8 = col_double()
#                           ),
#                           delim = '\t', col_names = FALSE)

# 
# count = prob.1.51.chg


## sqrt-transform of large DP 

DP_max = 60
k = DP_max - sqrt(DP_max)

count[count > DP_max] = round(sqrt(count[count > DP_max]) + DP_max)



# prior

# r = c(0.01, 0.98, 0.01)
# r = c(0.59, 0.4, 0.01)

# error rate

p = argv$error_rate

STATUS = c(1, 0, 0.5)

# Poisson approx

i = count[,2] > argv$Poisson_approx

probs = matrix(0, ncol = 3, nrow = nrow(count))

update.probs <- function(probs){
  
  r.m = matrix(0, nrow = 3, ncol = 3) 
  diag(r.m) = r
  
  probs[i,] = cbind(dpois(count[i,2] - count[i,1], count[i,2] * p), 
                    dpois(count[i,1], count[i,2] * p), 
                    dpois(count[i,1], count[i,2] * 0.5)) %*% r.m
  
  probs[!i,] = cbind(dbinom(count[!i,2] - count[!i,1], count[!i,2], p), 
                     dbinom(count[!i,1], count[!i,2], p),
                     dbinom(count[!i,1], count[!i,2], 0.5)) %*% r.m
  
  probs = probs / rowSums(probs)
  
  return(probs)
}


update.p <- function(p){
  
  j = (count[,2] > 2) & (!(status == 0.5 & max.prob > 0.8))
  
  p = sum(cbind(count[j,2] - count[j,1], count[j,1]) * 
            (probs[j, 1:2] / rowSums(probs[j, 1:2])))/
    sum(count[j,2])
  
  return(p)
}

update.r <- function(r){
  
  # return(colMeans(probs))
  
  a = c(sum(status==1), sum(status==0))/length(status)
  return(c(a, 1 - sum(a)))
}

## iterations

for(K in 1:argv$rounds){
  
  
  # posterior probs
  
  probs = update.probs(probs)
  
  # status
  
  status = STATUS[apply(probs, 1, which.max)]
  
  # maximal prob/confidence
  
  max.prob = signif(matrixStats::rowMaxs(probs), digits = 6)
  
  # update r, compositions
  
  r = update.r(r)
  
  # update p, error rate
  
  p = update.p(p)
  
}


## write p and r to summary file 

# write(sprintf('Rows: %d Columns: %d', dim(count)[1], dim(count)[2]), 
#       argv$summary
#       )
write(sprintf("Rows:\t%d", nrow(count)), argv$summary)


# write(format(Sys.time(), "%X %b %d %Y"), 
#       argv$summary, append = TRUE
#       )

# write(sprintf(
#   "methylation composition
# methylated(1)\tunmethylated(0)\tepi-heterozygous(0.5)
# %f\t%f\t%f
# error rate of measurement (0->1 or 1->0): %f\n", 
#   r[1], r[2], r[3], p), 
#   argv$summary, append = TRUE
#   )

write(sprintf(
  "methylated(1):\t%f
unmethylated(0):\t%f
epi-heterozygous(0.5):\t%f
error rate of measurement (0->1 or 1->0):\t%f", 
r[1], r[2], r[3], p), 
argv$summary, append = TRUE
)


## write status and probability
# col 1: predicted status
# col 2: corresponding prob

write.table(
  cbind(status, max.prob),
  file = argv$output,
  sep = '\t',
  col.names = FALSE,
  row.names = FALSE
)



