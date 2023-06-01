


library(matrixStats)
library(tidyverse)
library(magrittr)

# library(MASS)
# library(e1071)
# library(quantreg)
# library(splines)
# library(Matrix)
# library(sparseMatrixStats)
# library(lattice)

#############################################################################################

snv = read_delim('snv-10K-merge.simple', 
                 delim = '\t', 
                 col_names = F, 
                 col_types = cols(X1 = col_character()),
                 num_threads = 4
)
names(snv)= c('chr', 'nuc', 'pos', 'sample')


# methpipe allelicmeth
asm.allelicmeth = read_delim('methpipe-10K-merge.asm', 
                             delim = '\t', 
                             col_names = F,
                             num_threads = 4
                             )
colnames(asm.allelicmeth) = 
  c('chr', 'pos', 'strand', 'CpG', 'p_value', 'dp', 'MM', 'MU', 'UM', 'UU', 'sample')

# asm.allelicmeth = asm.allelicmeth[sample.int(nrow(asm.allelicmeth), 5e4), ]

#############################################################################################
## effective coverage
## both covered at least 5 reads

sDP = 5
# asm.allelicmeth %<>% subset((MM>=sDP & UU>=sDP) | (MU>=sDP & UM>=sDP))

asm.allelicmeth %<>% 
  mutate(is_asm = (p_value < 1e-3) & ((MM>=sDP & UU>=sDP) | (MU>=sDP & UM>=sDP)),
         pos = pos + 1,
         # nuc = 'C',
         CpG_id = 1:nrow(.),
         strand = NULL,
         CpG = NULL
  )

# summary of each sample
asm.count = asm.allelicmeth %>% 
  group_by(sample) %>%
  summarise(count.CG = n(), count.asm = sum(is_asm), pi = count.asm/count.CG)

#############################################################################################
## C or G mutation in a CpG disables the CpG

asm.allelicmeth = asm.allelicmeth %>% 
  # subset(is_asm) %>% 
  bind_rows(mutate(., pos = pos + 1, 
                   # nuc = 'G'
  ))

# asm.allelicmeth %<>% arrange(sample, chr, pos)

snv.cg = snv %>% 
  mutate(is_snv = T) %>% 
  subset(nuc %in% c('C', 'G')) %>%
  dplyr::select(c(chr, pos, sample, is_snv))

asm.allelicmeth = asm.allelicmeth %>%
  group_by(sample) %>%
  group_modify(
    ~ merge(.x, snv.cg %>% subset(.$sample == .y$sample, select = -sample),
            by = c('chr', 'pos'),
            all.x = T,
            sort = F)) %>%
  ungroup


asm.allelicmeth %<>% replace_na(list(is_snv = F))

is_snv = asm.allelicmeth %>% group_by(CpG_id) %>%
  summarise(is_snv = any(is_snv)) %>%
  arrange(CpG_id)
  
# asm.allelicmeth = merge(asm.allelicmeth %>% subset(duplicated(.$CpG_id)),
#                         is_snv)

asm1 = asm.allelicmeth %>% subset(duplicated(.$CpG_id)) %>%
  arrange(CpG_id) %>%
  mutate(is_snv = NULL)

asm.allelicmeth = bind_cols(asm1, is_snv[, "is_snv"])

#############################################################################################

prop = asm.count$pi %>% set_names(asm.count$sample)


# mark pi of Ber(pi) of each CpG site covered

asm.allelicmeth %<>% mutate(pi = prop[sample])

## CpG table with asm marked

CpG_with_asm = asm.allelicmeth %>%
  group_by(chr, pos) %>%
  summarise(chi.sq = sum(-log(p_value)*(!is_snv)),
            count.cov = n(),
            count.asm = sum(is_asm),
            count.snv = sum(is_snv),
            lambda = sum(pi*(!is_snv))
            ) %>%
  ungroup

#############################################################################################

## test each CpG is a consistent ASM site 

# using Chi square test of the sum of p.values
CpG_with_asm %<>% 
  mutate(p.chisq = pchisq(chi.sq, df = count.cov - count.snv, lower.tail = F))


## using a Pois dist of the asm count
CpG_with_asm %<>% 
  mutate(p.pois = ppois(count.asm - count.snv, lambda = lambda, lower.tail = F))

## reorder columns

CpG_with_asm %<>% transmute(chr = chr, 
                            pos = pos,
                            count.cov = count.cov,
                            count.asm = count.asm,
                            count.snv = count.snv,
                            chi.sq = chi.sq,
                            p.chisq = p.chisq,
                            lambda = lambda,
                            p.pois = p.pois
                            )

#############################################################################################

## conditions of BCA ASM CpGs

asm.cons = CpG_with_asm %>% subset(p.pois < 0.001 & 
                                     lambda > 0 & 
                                     count.cov >= 10 &
                                     (count.asm-count.snv) >= 5 &
                                     (count.asm-count.snv)/(count.cov-count.snv) >= 0.7
                                   )

write_tsv(asm.cons, file = './output/asm.consensus.tsv')



