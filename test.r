source('NGG_lp.r')
source('gen_data_by_NGG_lp.r')

start_time = proc.time()

G = 10000
n = 30

psi = c(
	2,1,2,1,
	2,1,2,1,
	2,1)

# psi = c(
# 	1,0.5,1,0.5,
# 	1,0.5,1,0.5,
# 	1,0.5)

# psi = c(
# 	1,0.5,-1,3,
# 	1,0.5,-1,3,
# 	1,0.5)

# psi = c(
# 	0,-2,1,0.5,
# 	0,-2,1,0.5,
# 	0,-2)

t_pi = c(0.05, 0.05, 0.90)

generate = 0

if (generate)
{
	E_Set = gen_data_by_NGG_lp(G, n, psi, t_pi)
	saveRDS(E_Set, 'E_Set.Rdata')	
}else
{
	E_Set = readRDS('E_Set.Rdata')
}
result = NGG_lp(E_Set, verbose = 1, is_sim = 1, M = 100, limma_prior = 0, cores = 4)

end_time = proc.time()
result$runnint_time = end_time - start_time
print(end_time - start_time)
