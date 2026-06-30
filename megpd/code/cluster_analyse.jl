include("markov_megpd_cluster.jl")

using Plots

clust_res = clust_sim(5.0, 10000, 2.0, 1.0, 0.2);
# boot = bootstrap_clust(clust_res.raw[1], clust_res.raw[2])

p1 = bar(
    1:length(clust_res.pi_dist),
    clust_res.pi_dist,
    xlabel = "Cluster size (i)",
    ylabel = "π(i)",
    legend = false
)

p2 = bar(
    1:length(clust_res.piC_dist),
    clust_res.piC_dist,
    xlabel = "Run length (i)",
    ylabel = "πC(i)",
    legend = false
)

plot(p1, p2, layout = (1, 2))