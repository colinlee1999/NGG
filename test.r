source('NGG_lp.r')
source('gen_data_by_NGG_lp.r')

G = 10000
n = 20

psi = c(
	2,1,2,1,
	2,1,2,1,
	2,1)

t_pi = c(0.05, 0.05, 0.90)

E_Set = gen_data_by_NGG_lp(G, n, psi, t_pi)
result = NGG_lp(E_Set, verbose = 1, is_sim = 1, M = 100)
