using Plots

data1 = [1000, 2000, 3000]
data2 = [2^10, 1^10, 5^10]

p = plot(data1, data2, xlims = (0, 3000), yscale=:log10)