

library(readr)
library(matrixStats)

# d = read.delim(stdin(), header = F)
# d = read.delim('T29_WGBS-000', header = F)

d = readr::read_delim(stdin(), delim = '\t', col_names = FALSE,
                      na = c("", "NA", "na"),
                      col_types = cols(
                        X1 = col_character(),
                        X2 = col_character(),
                        X3 = col_double(),
                        X4 = col_character(),
                        X5 = col_character(),
                        X6 = col_double(),
                        X7 = col_double(),
                        X8 = col_double(),
                        X9 = col_double(),
                        X10 = col_double(),
                        X11 = col_double(),
                        X12 = col_double(),
                        X13 = col_double(),
                        X14 = col_double(),
                        X15 = col_double(),
                        X16 = col_double()
                      ))

Watson = as.matrix(d[,6:9])
Crick  = as.matrix(d[,11:14])[,c(2,1,4,3)]

# take sqrt transform of the counts
# considering the counts are actually positively correlated

t = 60
k = t - sqrt(t)

Watson[Watson > t] = round(sqrt(Watson[Watson > t]) + k)
Crick[Crick > t] = round(sqrt(Crick[Crick > t]) + k)

# mutation rate
pm = 1/1000/3

# error rate
# in oocyte samples, error rate is set triple
# pe = 1/100/3
pe = 3/100/3

# total mis rate
p = pm + pe

# methylation rate/proportion 
pr.cg = 0.6        # CG context
pr.ncg = 1/100     # non-CG context


# transition prob of haploidy

PA = function(pr){c(1-3*pm-3*pe, 2*pm-pm*pr+pe, pm*pr+pe, pm+pe)}
PT = function(pr){c(pm+pe, 1-2*pm-pm*pr-3*pe, pm*pr+pe, pm+pe)}
PC = function(pr){c(pm+pe, pm+pe+(1-3*pm-3*pe)*(1-pr), (1-3*pm-3*pe)*pr, pm+pe)}
PG = function(pr){c(pm+pe, 2*pm-pm*pr+pe, pm*pr+pe, 1-3*pm-3*pe)}


# STATUS
HOMO = c('A', 'T', 'C', 'G')
HETER = c('AC', 'AG', 'AT', 'CG', 'CT', 'GT')
STATUS = c(HOMO, HETER)

# prior

ps = c((1-3*p)^2, p^2, 2*p*(1-3*p))

pri.A = ps[c(1,2,2,2,3,3,3,2,2,2)]
pri.T = ps[c(2,1,2,2,2,2,3,2,3,3)]
pri.C = ps[c(2,2,1,2,3,2,2,3,3,2)]
pri.G = ps[c(2,2,2,1,2,3,2,3,2,3)]

pris = list(pri.A, pri.T, pri.C, pri.G)


postp <- function(ref = 'A', cg = TRUE, Watson = 1:4, Crick = 1:4) {
  
  if(ref == 'N') return(NA)
  
  if(cg) pr = pr.cg
  else pr = pr.ncg
  
  # prior
  
  theta = pris[[which.max(c('A', 'T', 'C', 'G') == ref)]]
  
  # conditional prob
  
  
  PA = PA(pr)
  PT = PT(pr)
  PC = PC(pr)
  PG = PG(pr)
  
  p.cond = c(dmultinom(c(Watson, Crick), prob = c(PA, PT)/2), # A
             dmultinom(c(Watson, Crick), prob = c(PT, PA)/2), # T
             dmultinom(c(Watson, Crick), prob = c(PC, PG)/2), # C
             dmultinom(c(Watson, Crick), prob = c(PG, PC)/2), # G
             dmultinom(c(Watson, Crick), prob = c(PA+PC, PT+PG)/4), # AC
             dmultinom(c(Watson, Crick), prob = c(PA+PG, PT+PC)/4), # AG
             dmultinom(c(Watson, Crick), prob = c(PA+PT, PT+PA)/4), # AT
             dmultinom(c(Watson, Crick), prob = c(PC+PG, PG+PC)/4), # CG
             dmultinom(c(Watson, Crick), prob = c(PC+PT, PG+PA)/4), # CT
             dmultinom(c(Watson, Crick), prob = c(PG+PT, PC+PA)/4)  # GT
  )
  
  # posterior prob
  p.post.unnorm = p.cond*theta
  p.post = p.post.unnorm/sum(p.post.unnorm)
  
  # prob not mutation (same with ref)
  # regarded as p.value
  
  p.value = p.post[1:4][which.max(c('A', 'T', 'C', 'G') == ref)]
  return(c(p.post, p.value, sum(Watson), sum(Crick), sum(p.post[1:4])))
}

# test
prob.post = matrix(0, nrow = nrow(d), ncol = length(STATUS) + 4)
# status.pred = rep('N', nrow(d))


for (i in 1:nrow(d)) {
  pp = postp(d[[i,2]], d[i,4] == 'CG', Watson[i,], Crick[i,])
  # pp = postp(ref = d[i,2], cg = d[i,4] == 'CG', c(0,0,2,2), c(0,0,2,2))
  
  # names(pp) = STATUS
  # barplot(pp, ylim = c(0,1))
  # 
  prob.post[i,] = pp
  # status.pred[i] = STATUS[which.max(pp)]
  
}


# max.prob = rowMaxs(prob.post[,1:10])
# sum(max.prob<0.95)


## allele frequencies

allele.weights = t(matrix(
  c(1, 0, 0, 0, 0.5, 0.5, 0.5, 0  , 0  , 0  ,
    0, 1, 0, 0, 0,   0  , 0.5, 0  , 0.5, 0.5,
    0, 0, 1, 0, 0.5, 0  , 0  , 0.5, 0.5, 0  ,
    0, 0, 0, 1, 0  , 0.5, 0  , 0.5, 0  , 0.5
  ),
  nrow = 4, byrow = T
))

allele.freq = prob.post[,1:10] %*% allele.weights

# 
# # validation
# 
# cg.bayes = read.delim('bayes.20m.snv')
# cg.binom =  read.delim('binom.20m.snv')
# 

# i.snv = prob.post[,11] < 0.05
# snv = d$V3[i.snv]
# 


## write p and r to stderr

# write(format(Sys.time(), "%X %b %d %Y"), stderr())
# write(sprintf(
#   "methylation composition
# methylated(1)\tunmethylated(0)\tepi-heterozygous(0.5)
# %f\t%f\t%f
# error rate of measurement (0->1 or 1->0): %f\n", 
#   r[1], r[2], r[3], p), stderr())


## write status and probability
# col 1:4, site info
# col 6, p-value testing SNV
# col 7:10, estimated frequencies of ATCG
# col 11:12, coverage depths of Wastson and Crick strands
# col 13, propbabilty of homozygous

res = cbind(d[,1:5], 
            sprintf('%.5e', prob.post[,11]),  # p-value is logged in raw format
            apply(allele.freq, 2, function(x) sprintf(fmt = '%.5e', x)),
            prob.post[,12:13],
            sprintf('%.5f', prob.post[,14])
)


readr::write_delim(
  res,
  file = stdout(),
  # 't',
  delim = '\t',
  col_names = FALSE
)
