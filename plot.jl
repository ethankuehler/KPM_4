using Plots
using DelimitedFiles

#for doing vectors in the ploting phase
function T_nv(x, n)
    return cos.(n .* acos.(x))
end


#jackson kernal
function Jackson(N, n)
    a1 = (N-n+1)*cos((pi*n)/(N+1))
    a2 = sin((pi*n)/(N+1))
    a3 = cot(pi/(N+1))
    return (1/(N + 1))*(a1 + a2 * a3)
end


Chemical_potential = vec(readdlm("data/mu_mag.csv", ','))

order = size(Chemical_potential)[1]
Length_of_output = 2*order
jack(n) = Jackson(order, n)

plot
#plot dos
x = collect(range(-1,stop=1,length=Length_of_output))
y = ((pi*sqrt.(-x.^2 .+ 1)).\1) .* (Chemical_potential[1]*jack(0) .+ sum(2*Chemical_potential[i]*T_nv(x, i-1)*jack(i) for i in (2:order)))
plot()
xlabel!("energy")
ylabel!("density")
title!("DOS weyl KPM")
display(plot!(x, y, legend=false))
savefig("dos.svg")
println("done")

writedlm( "data/mag.csv",  y, ',')
writedlm( "points.csv",  x, ',')

