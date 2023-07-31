using Plots

x = [1000, 2000, 3000]
y1 = [2^10, 1^10, 5^10]
y2 = [3^10, 3^10, 5^10]

p = plot(x, [y1 y2],
        xlims = (0, 3000),
        yscale=:log10,
        xlabel="Length",
        ylabel="Time",
        label=["NW" "Kmer"],
        lw = [3 4])

