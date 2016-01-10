source('NGG_lp.r')

E_Set = readRDS('es.rData')
result = NGG_lp(E_Set, verbose = 1, is_sim = 1, M = 100)