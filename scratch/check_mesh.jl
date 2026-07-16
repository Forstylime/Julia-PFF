using FerriteGmsh
using Ferrite
grid = FerriteGmsh.togrid("data/mesh/l_shape.msh")
nodes = getnodes(grid)
x_coords = [n.x[1] for n in nodes]
y_coords = [n.x[2] for n in nodes]
println("X range: ", minimum(x_coords), " to ", maximum(x_coords))
println("Y range: ", minimum(y_coords), " to ", maximum(y_coords))
