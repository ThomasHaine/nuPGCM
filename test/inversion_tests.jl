using Test
using nuPGCM
using Gridap, GridapGmsh
using IncompleteLU, Krylov, LinearOperators, CuthillMcKee
using CUDA, CUDA.CUSPARSE, CUDA.CUSOLVER
using SparseArrays, LinearAlgebra, Statistics
using Printf
using PyPlot

pygui(false)
plt.style.use("plots.mplstyle")
plt.close("all")

function coarse_inversion(dim::AbstractDimension, arch::AbstractArchitecture)
    # tolerance and max iterations for iterative solvers
    tol = 1e-6
    itmax = 0

    # Vector type 
    VT = typeof(arch) == CPU ? Vector{Float64} : CuVector{Float64}

    # coarse model
    h = 0.05
    model = GmshDiscreteModel(@sprintf("meshes/bowl%s_%0.2f.msh", dim, h))

    # FE spaces
    X, Y, B, D = setup_FESpaces(model)
    Ux, Uy, Uz, P = unpack_spaces(X)
    nx = Ux.space.nfree
    ny = Uy.space.nfree
    nz = Uz.space.nfree
    nu = nx + ny + nz
    np = P.space.space.nfree
    N = nu + np - 1
    if typeof(dim) == TwoD
        @test N == 7451
    else
        @test N == 227166
    end

    # triangulation and integration measure
    Ω = Triangulation(model)
    dΩ = Measure(Ω, 4)

    # depth
    H(x) = 1 - x[1]^2 - x[2]^2

    # forcing
    ν(x) = 1
    κ(x) = 1e-2 + exp(-(x[3] + H(x))/0.1)

    # params
    ε² = 1e-2
    γ = 1/4
    f₀ = 1
    β = 0.5
    f(x) = f₀ + β*x[2]

    # assemble LHS inversion and test against saved matrix
    LHS_inversion_fname = @sprintf("test/data/LHS_inversion_%s_%e_%e_%e_%e_%e.h5", dim, h, ε², γ, f₀, β)
    # LHS_inversion, perm_inversion, inv_perm_inversion = assemble_LHS_inversion(arch, γ, ε², ν, f, X, Y, dΩ; fname=LHS_inversion_fname) # FIRST TIME ONLY: save LHS matrix
    LHS_inversion, perm_inversion, inv_perm_inversion = assemble_LHS_inversion(arch, γ, ε², ν, f, X, Y, dΩ; fname="LHS_inv_temp.h5")
    @test LHS_inversion ≈ read_sparse_matrix(LHS_inversion_fname)[1]

    # inversion RHS
    RHS_inversion = assemble_RHS_inversion(perm_inversion, B, Y, dΩ)

    # preconditioner
    if typeof(dim) == TwoD 
        if typeof(arch) == CPU
            P_inversion = lu(LHS_inversion)
            ldiv_P_inversion = true
        else
            P_inversion = Diagonal(on_architecture(arch, 1/h^2*ones(N)))
            ldiv_P_inversion = false
        end
    else
        P_inversion = I
        ldiv_P_inversion = false
    end

    # put on GPU, if needed
    LHS_inversion = on_architecture(arch, LHS_inversion)
    RHS_inversion = on_architecture(arch, RHS_inversion)

    # Krylov solver for inversion
    solver_inversion = GmresSolver(N, N, 20, VT)
    solver_inversion.x .= 0.

    # inversion functions
    function invert!(arch::AbstractArchitecture, solver, b)
        b_arch = on_architecture(arch, b.free_values)
        if typeof(arch) == GPU
            RHS = [CUDA.zeros(nx); CUDA.zeros(ny); RHS_inversion*b_arch; CUDA.zeros(np-1)]
        else
            RHS = [zeros(nx); zeros(ny); RHS_inversion*b_arch; zeros(np-1)]
        end
        Krylov.solve!(solver, LHS_inversion, RHS, solver.x, M=P_inversion, ldiv=ldiv_P_inversion,
                    atol=tol, rtol=tol, verbose=0, itmax=itmax, restart=true,
                    history=true)
        return solver
    end
    function update_u_p!(ux, uy, uz, p, solver)
        sol = on_architecture(CPU(), solver.x[inv_perm_inversion])
        ux.free_values .= sol[1:nx]
        uy.free_values .= sol[nx+1:nx+ny]
        uz.free_values .= sol[nx+ny+1:nx+ny+nz]
        p = FEFunction(P, sol[nx+ny+nz+1:end])
        return ux, uy, uz, p
    end

    # background state ∂z(b) = N^2
    N² = 1.

    # simple test buoyancy field: b = δ exp(-(z + H)/δ)
    b = interpolate_everywhere(x -> 0.1*exp(-(x[3] + H(x))/0.1), B)
    ux = interpolate_everywhere(0, Ux)
    uy = interpolate_everywhere(0, Uy)
    uz = interpolate_everywhere(0, Uz)
    p  = interpolate_everywhere(0, P) 

    # invert
    solver_inversion = invert!(arch, solver_inversion, b)
    ux, uy, uz, p = update_u_p!(ux, uy, uz, p, solver_inversion)

    # # plot for sanity check
    # sim_plots(dim, ux, uy, uz, b, N², H, 0, 0, "test")

    # # FIRST TIME ONLY: save state
    # save_state(ux, uy, uz, p, b, 0; fname=@sprintf("test/data/inversion_%s.h5", dim))

    # compare state with data
    ux_data, uy_data, uz_data, p_data, b_data, t_data = load_state(@sprintf("test/data/inversion_%s.h5", dim))
    @test isapprox(ux.free_values, ux_data, rtol=1e-2)
    @test isapprox(uy.free_values, uy_data, rtol=1e-2)
    @test isapprox(uz.free_values, uz_data, rtol=1e-2)
    @test isapprox(p.free_values,  p_data,  rtol=1e-2)
    @test isapprox(b.free_values,  b_data,  rtol=1e-2)
    println(norm(ux.free_values - ux_data)/norm(ux_data))
    println(norm(uy.free_values - uy_data)/norm(uy_data))
    println(norm(uz.free_values - uz_data)/norm(uz_data))
    println(norm(p.free_values - p_data)/norm(p_data))
    println(norm(b.free_values - b_data)/norm(b_data))

    # remove temporary files
    rm("LHS_inv_temp.h5", force=true)
end

@testset "Inversion Tests" begin
    @testset "2D CPU" begin
        coarse_inversion(TwoD(), CPU())
    end
    # @testset "2D GPU" begin
    #     coarse_inversion(TwoD(), GPU())
    # end
    # @testset "3D GPU" begin
    #     coarse_inversion(ThreeD(), GPU())
    # end
end