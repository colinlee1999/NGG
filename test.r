source('NGG_lp.r')
source('gen_data_by_NGG_lp.r')

start_time = proc.time()

G = 1000
n = 30

psi = c(
	1,0.5,1,0.5,
	1,0.5,1,0.5,
	1,0.5)

t_pi = c(0.05, 0.05, 0.90)

E_Set = gen_data_by_NGG_lp(G, n, psi, t_pi)
saveRDS(E_Set, 'E_Set.Rdata')
# E_Set = readRDS('E_Set.Rdata')
result = NGG_lp(E_Set, verbose = 1, is_sim = 1, M = 100, limma_prior = 0)

end_time = proc.time()
print(end_time - start_time)
