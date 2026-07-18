using Serialization, Ferrite, FerriteGmsh

"""
    TensionSetup

保存方形拉伸相场断裂算例的有限元初始化结果。也可以用于其他几何相似的算例（如 L 形拉伸），只要保证网格和边界条件设置与求解器要求一致。

该结构体把后续求解器需要反复使用的对象集中在一起：计算网格、位移场和
相场的自由度处理器、两类场的约束处理器，以及用于设置预制裂纹的节点编号。
"""
Base.@kwdef struct TensionSetup{G,DHU,DHD,CHU,CHD,T<:AbstractFloat}
    dir::Int
    grid::G
    dh_u::DHU
    dh_d::DHD
    ch_u::CHU
    ch_d::CHD
    final_displacement::T
end

"""
    create_grid(msh_file = "data/mesh/l_shape.msh")
    生成 L 形算例使用的二维四边形结构网格。
    网格一般已在Gmsh中生成好，直接从 .msh 文件读取。
"""
function create_grid(msh_file = "data/mesh/l_shape.msh")
    cache_file = msh_file * ".jls"
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
    create_dofhandlers(grid)

为交错求解格式创建位移场和相场各自的自由度处理器。
"""
function create_dofhandlers(grid)
    # 位移自由度处理器：每个节点有两个分量，对应 ux 和 uy。
    dh_u = Ferrite.DofHandler(grid)
    ip_u = Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()^2
    Ferrite.add!(dh_u, :u, ip_u)
    # `close!` 会冻结自由度布局并建立单元到全局自由度的映射。
    Ferrite.close!(dh_u)

    # 相场自由度处理器：每个节点一个标量损伤变量 d。
    dh_d = Ferrite.DofHandler(grid)
    ip_d = Ferrite.Lagrange{Ferrite.RefQuadrilateral, 1}()
    Ferrite.add!(dh_d, :d, ip_d)
    # 关闭后才能查询自由度数量、组装矩阵或创建约束。
    Ferrite.close!(dh_d)

    println("位移自由度数量: ", Ferrite.ndofs(dh_u))
    println("相场自由度数量: ", Ferrite.ndofs(dh_d))

    return dh_u, dh_d
end

"""
    create_displacement_constraints(dh_u, grid; final_displacement = 0.0)

创建位移场的 Dirichlet 边界条件。
"""
function create_displacement_constraints(dh_u, grid, fixed_face = "top", final_displacement = 0.0, dir = 2)
    dir in (1, 2) || throw(ArgumentError("dir must be 1 (x) or 2 (y)"))
    ch_u = Ferrite.ConstraintHandler(dh_u)

    # 读取网格生成阶段创建的边界 facet set。
    fixed = Ferrite.getfacetset(grid, fixed_face)
    right = Ferrite.getfacetset(grid, "right")

    # ux、uy 全固定，作为拉伸试样的支承边界。
    Ferrite.add!(
        ch_u,
        Ferrite.Dirichlet(:u, fixed, (x, t) -> zeros(2), [1, 2]),
    )
    # 右边界只约束位移分量；加载幅值通过时间/载荷参数 t 缩放。
    Ferrite.add!(
        ch_u,
        Ferrite.Dirichlet(:u, right, (x, t) -> t * final_displacement, dir),
    )

    # 关闭约束处理器并在 t = 0 时初始化约束值。
    Ferrite.close!(ch_u)
    Ferrite.update!(ch_u, 0.0)
    return ch_u
end

"""
    create_phase_field_constraints(dh_d)

创建相场变量的约束处理器。
"""
function create_phase_field_constraints(dh_d)
    ch_d = Ferrite.ConstraintHandler(dh_d)
    Ferrite.close!(ch_d)
    Ferrite.update!(ch_d, 0.0)
    return ch_d
end

"""
    setup_tension(; msh_file = "data/mesh/l_shape.msh", final_displacement, fixed_face, dir)

一站式构建方形拉伸相场断裂算例的有限元初始化对象。
"""
function setup_tension(;
    msh_file = "data/mesh/l_shape.msh",
    final_displacement = 0.0,
    fixed_face = "top",
    dir = 2
)
    # 1. 创建计算网格，后续所有自由度和边界集合都基于同一个 grid。
    grid = create_grid(msh_file)
    # 2. 分别为位移场 u 和相场 d 建立自由度编号，适配交错求解流程。
    dh_u, dh_d = create_dofhandlers(grid)
    # 3. 设置力学边界条件：上边(top)固定，右边(right)施加竖直向下位移加载。
    ch_u = create_displacement_constraints(dh_u, grid, fixed_face, final_displacement, dir)
    # 4. 创建相场约束处理器。
    ch_d = create_phase_field_constraints(dh_d)

    # 统一封装初始化结果，减少求解脚本需要手动传递的对象数量。
    return TensionSetup(
        dir,
        grid,
        dh_u,
        dh_d,
        ch_u,
        ch_d,
        float(final_displacement),
    )
end
