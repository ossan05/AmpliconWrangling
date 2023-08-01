using Plots

x10 = [i for i in 300:300:3000]
x20 = [i for i in 150:150:3000]
Reg = [1.613, 5.567, 13.056, 22.244, 35.604, 51.937, 73.169, 94.903, 120.708, 147.348]
Tri = [3.175, 12.005, 25.561, 45.292, 75.864, 109.775, 154.902, 193.013, 246.732, 319.301]
Affine = [2.577, 9.271, 20.648, 41.571, 60.159, 91.789, 122.633, 159.275, 213.057, 258.521]
ChainMatching = [0.3841, 1.278, 2.108, 3.777, 4.935, 5.530, 5.946, 7.253, 8.819, 9.298]
RegKmer = [0.3701, 1.280, 2.270, 3.786, 4.435, 4.578, 5.124, 6.838, 6.356, 9.142]

Old_nw = [5.003 , 23.141 , 53.833 , 95.217 , 190.597 , 267.570 , 441.479 , 387.069 , 480.260 , 604.483 ]
Old_aff_nw = [14.946 , 68.349 , 188.608 , 399.125 , 673.922 , 935.105 , 1309 , 1144 , 1566 , 1874 ]
Old_kmeralign = [0.394400 , 0.9072 , 4.847 , 8.476 , 18.812 , 25.321 , 26.667 , 30.455 , 96.079 , 39.447 , 97.425 , 17.375 , 67.941 , 220.274 , 134.627 , 122.520 , 226.438 , 125.014 , 312.716 , 305.850 ]
Old_triplet_nw =[1.368 , 5.595 , 12.813 , 25.096 , 40.734 , 58.606 , 86.346 , 120.272 , 181.555 , 248.476 , 290.408 , 360.956 , 454.040 , 444.792 , 607.247 , 654.854 , 475.524 , 532.713 , 612.120 , 670.440 ]

var_sum = 0.0005 + 0.003 + 0.00005 + 0.0001

xvar = [3i * var_sum for i in 1:10]
println(xvar)
Var = [1.617, 1.815, 1.760, 1.668, 1.992, 2.009, 2.003, 2.165, 1.849, 2.230]

# p = plot(x10, [Reg RegKmer],
#         xlims = (0, 3000),
#         xlabel="Length (nt bases)",
#         ylabel="Time (ms)",
#         label=["Reg NW" "Kmer1"],
#         lw = [1 2])

# plot!(x10, [Tri, Affine, ChainMatching],
#         label=["Triplet NW" "Affine NW" "ChainMatching"],
#         lw=[1 1 2])

# plot!(x10, [Old_nw Old_aff_nw],
#         label=["Old_NW" "Old_NW with affine gap" ],
#         lw=[1 1])
# plot!(x20, [Old_kmeralign Oldriplet_nw],
#         label=["Old kmer" "Old triplet NW"],
#         lw=[2 1])

# plot(x10, RegKmer,
#         title="K-mer aligner",
#         xlims = (0, 3000),
#         xlabel="Length (nt bases)",
#         ylabel="Time (ms)",
#         label="Triplet",
#         lw = 4,
#         lc = :red)


# plot!(x20, Old_kmeralign,
#         label="Old Triplet",
#         lw=4,
#         lc=:blue)

# savefig("KmerAligners.png")

plot(xvar, Var,
        title="Kmer matching vs variance",
        xlabel="Variance",
        ylabel="Time(ms)")

savefig("KmerVariance")