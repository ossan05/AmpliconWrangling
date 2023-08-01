using Plots

x = [i for i in 300:300:3000]
Reg = [1.613, 5.567, 13.056, 22.244, 35.604, 51.937, 73.169, 94.903, 120.708, 147.348]
Tri = [3.175, 12.005, 25.561, 45.292, 75.864, 109.775, 154.902, 193.013, 246.732, 319.301]
Affine = [2.577, 9.271, 20.648, 41.571, 60.159, 91.789, 122.633, 159.275, 213.057, 258.521]
ChainMatching = [0.3841, 1.278, 2.108, 3.777, 4.935, 5.530, 5.946, 7.253, 8.819, 9.298]
RegKmer = [0.3701, 1.280, 2.270, 3.786, 4.435, 4.578, 5.124, 6.838, 6.356, 9.142]

p = plot(x, [Reg RegKmer],
        xlims = (0, 3000),
        yscale=:log10,
        xlabel="Length (nt bases)",
        ylabel="Time (milliseconds)",
        label=["Reg NW" "Kmer1"],
        lw = [1 2])

plot!(x, [Tri, Affine, ChainMatching],
        label=["Triplet NW" "Affine NW" "ChainMatching"],
        lw=[1 1 2])