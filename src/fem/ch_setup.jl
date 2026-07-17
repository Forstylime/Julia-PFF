# src/fem/ch_setup.jl

using Serialization, Ferrite, FerriteGmsh

"""
    CHSetup

Holds the finite element initialization results for the Cahn-Hilliard equation solver.
"""
Base.@kwdef struct CHSetup{G,DH,CH}
    grid::G
    dh::DH
    ch::CH
end

"""
    create_ch_grid(msh_file = "data/mesh/l_shape.msh")
    
Generates a grid from a Gmsh .msh file.
"""
function create_ch_grid(msh_file = "data/mesh/l_shape.msh")
    cache_file = msh_file * ".ch.jls"
    if isfile(cache_file) && mtime(cache_file) >= mtime(msh_file)
        println("从缓存加载网格: ", cache_file)
        return deserialize(cache_file)
    end
    println("解析 .msh 文件: ", msh_file)
    grid = FerriteGmsh.togrid(msh_file)
    serialize(cache_file, grid)
    println("网格已缓存至: ", cache_file)
    return grid
end

"""
    create_ch_dofhandler(grid)

Creates a mixed DofHandler for the Cahn-Hilliard equation containing two scalar fields: ϕ and μ.
"""
function create_ch_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    
    # We assume a 2D quadrilateral mesh, consistent with the rest of the project.
    ip = Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()
    
    Ferrite.add!(dh, :phi, ip)
    Ferrite.add!(dh, :mu, ip)
    
    Ferrite.close!(dh)

    println("Cahn-Hilliard 混合场总自由度数量: ", Ferrite.ndofs(dh))

    return dh
end

"""
    create_ch_constraints(dh)

Creates constraint handler for the Cahn-Hilliard equation.
Assuming homogeneous Neumann boundary conditions, no Dirichlet BCs are needed.
"""
function create_ch_constraints(dh)
    ch = Ferrite.ConstraintHandler(dh)
    Ferrite.close!(ch)
    Ferrite.update!(ch, 0.0)
    return ch
end

"""
    setup_ch(; msh_file = "data/mesh/l_shape.msh")

One-stop setup for the Cahn-Hilliard problem.
"""
function setup_ch(; msh_file = "data/mesh/l_shape.msh")
    grid = create_ch_grid(msh_file)
    dh = create_ch_dofhandler(grid)
    ch = create_ch_constraints(dh)
    return CHSetup(grid, dh, ch)
end
