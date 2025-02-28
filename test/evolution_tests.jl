using Test
using nuPGCM
using Gridap
using Krylov
using CUDA, CUDA.CUSPARSE, CUDA.CUSOLVER
using SparseArrays
using LinearAlgebra
using ProgressMeter
using Printf
using PyPlot

pygui(false)
plt.style.use("plots.mplstyle")
plt.close("all")

function coarse_evolution(dim, arch)
    # params/funcs
    ε² = 1e-2
    γ = 1/4
    f₀ = 1
    β = 0.5
    f(x) = f₀ + β*x[2]
    μϱ = 1e-4
    Δt = 1e-4*μϱ/ε²
    T = 5e-2*μϱ/ε²
    α = Δt/2*ε²/μϱ
    H(x) = 1 - x[1]^2 - x[2]^2
    ν(x) = 1
    κ(x) = 1e-2 + exp(-(x[3] + H(x))/0.1)

    # coarse mesh
    h = 0.1
    mesh = Mesh(@sprintf("meshes/bowl%s_%0.2f.msh", dim, h))

    # dof
    nu, nv, nw, np, nb = get_n_dof(mesh)

    # dof perms
    p_u, p_v, p_w, p_p, p_b = compute_dof_perms(mesh)
    inv_p_b = invperm(p_b)
    p_inversion = [p_u; p_v .+ nu; p_w .+ nu .+ nv; p_p .+ nu .+ nv .+ nw]
    inv_p_inversion = invperm(p_inversion)

    # build LHS matrix for inversion and test against saved matrix
    A_inversion_fname = @sprintf("test/data/A_inversion_%s_%e_%e_%e_%e_%e.h5", dim, h, ε², γ, f₀, β)
    if !isfile(A_inversion_fname)
        @warn "A_inversion file not found, generating..."
        A_inversion = build_A_inversion(mesh, γ, ε², ν, f; fname=A_inversion_fname)
    else
        # just read inversion matrix instead of building and testing it; this is already tested in inversion_tests.jl
        A_inversion = load(A_inversion_fname, "A_inversion")
    end

    # re-order dof
    A_inversion = A_inversion[p_inversion, p_inversion]

    # build RHS matrix for inversion
    B_inversion = build_B_inversion(mesh)

    # re-order dof
    B_inversion = B_inversion[p_inversion, :]

    # preconditioner
    if typeof(arch) == CPU
        P_inversion = lu(A_inversion)
    else
        P_inversion = Diagonal(on_architecture(arch, 1/h^dim.n*ones(N)))
    end

    # put on GPU, if needed
    A_inversion = on_architecture(arch, A_inversion)
    B_inversion = on_architecture(arch, B_inversion)

    # setup inversion toolkit
    inversion_toolkit = InversionToolkit(A_inversion, P_inversion, B_inversion)

    function update_u_p!(ux, uy, uz, p, solver)
        sol = on_architecture(CPU(), solver.x[inv_p_inversion])
        ux.free_values .= sol[1:nx]
        uy.free_values .= sol[nx+1:nx+ny]
        uz.free_values .= sol[nx+ny+1:nx+ny+nz]
        p = FEFunction(P, sol[nx+ny+nz+1:end])
        return ux, uy, uz, p
    end

    # background state \partial_z b = N^2
    N² = 1.

    # initial condition: b = N^2 z, t = 0
    b = interpolate_everywhere(0, B)
    t = 0.
    ux = interpolate_everywhere(0, Ux)
    uy = interpolate_everywhere(0, Uy)
    uz = interpolate_everywhere(0, Uz)
    p  = interpolate_everywhere(0, P) 

    # assemble evolution matrices and test against saved matrices
    LHS_diff_fname = @sprintf("test/data/LHS_diff_%s_%e_%e_%e.h5", dim, h, α, γ)
    LHS_adv_fname = @sprintf("test/data/LHS_adv_%s_%e.h5", dim, h)
    if !isfile(LHS_diff_fname) || !isfile(LHS_adv_fname)
        @warn "LHS_diff or LHS_adv file not found, generating..."
        LHS_adv, LHS_diff, perm_b, inv_perm_b = assemble_LHS_adv_diff(CPU(), α, γ, κ, B, D, dΩ; fname_adv=LHS_adv_fname, fname_diff=LHS_diff_fname)
    else
        LHS_adv, LHS_diff, perm_b, inv_perm_b = assemble_LHS_adv_diff(CPU(), α, γ, κ, B, D, dΩ; fname_adv="LHS_adv_temp.h5", fname_diff="LHS_diff_temp.h5")
        @test LHS_adv ≈ read_sparse_matrix(LHS_adv_fname)[1]
        @test LHS_diff ≈ read_sparse_matrix(LHS_diff_fname)[1]
    end

    # diffusion RHS matrix and vector
    RHS_diff, rhs_diff = assemble_RHS_diff(perm_b, α, γ, κ, N², B, D, dΩ)

    # preconditioners
    if typeof(arch) == CPU
        P_diff = lu(LHS_diff)
        P_adv  = lu(LHS_adv)
    else
        P_diff = Diagonal(on_architecture(arch, Vector(1 ./ diag(LHS_diff))))
        P_adv  = Diagonal(on_architecture(arch, Vector(1 ./ diag(LHS_adv))))
    end

    # put on GPU, if needed
    LHS_diff = on_architecture(arch, LHS_diff)
    RHS_diff = on_architecture(arch, RHS_diff)
    rhs_diff = on_architecture(arch, rhs_diff)
    LHS_adv = on_architecture(arch, LHS_adv)

    # Krylov solver for evolution
    VT = typeof(arch) == CPU ? Vector{Float64} : CuVector{Float64}
    tol = 1e-6
    itmax = 0
    solver_evolution = CgSolver(nb, nb, VT)
    solver_evolution.x .= on_architecture(arch, copy(b.free_values)[perm_b])

    # evolution functions
    b_half = interpolate_everywhere(0, B)
    function evolve_adv!(inversion_toolkit, solver_evolution, ux, uy, uz, p, b)
        # determine architecture
        arch = architecture(solver_evolution.x)

        # half step
        l_half(d) = ∫( b*d - Δt/2*(ux*∂x(b) + uy*∂y(b) + uz*(N² + ∂z(b)))*d )dΩ
        RHS = on_architecture(arch, assemble_vector(l_half, D)[perm_b])
        Krylov.solve!(solver_evolution, LHS_adv, RHS, solver_evolution.x, M=P_adv, atol=tol, rtol=tol, verbose=0, itmax=itmax)

        # u, v, w, p, b at half step
        update_b!(b_half, solver_evolution)
        invert!(inversion_toolkit, b_half)
        ux, uy, uz, p = update_u_p!(ux, uy, uz, p, inversion_toolkit.solver)

        # full step
        l_full(d) = ∫( b*d - Δt*(ux*∂x(b_half) + uy*∂y(b_half) + uz*(N² + ∂z(b_half)))*d )dΩ
        RHS = on_architecture(arch, assemble_vector(l_full, D)[perm_b])
        Krylov.solve!(solver_evolution, LHS_adv, RHS, solver_evolution.x, M=P_adv, atol=tol, rtol=tol, verbose=0, itmax=itmax)

        return inversion_toolkit, solver_evolution
    end
    function evolve_diff!(solver, b)
        arch = architecture(solver.x)
        b_arch= on_architecture(arch, b.free_values)
        RHS = RHS_diff*b_arch + rhs_diff
        Krylov.solve!(solver, LHS_diff, RHS, solver.x, M=P_diff, atol=tol, rtol=tol, verbose=0, itmax=itmax)
        return solver
    end
    function update_b!(b, solver)
        b.free_values .= on_architecture(CPU(), solver.x[inv_p_b])
        return b
    end

    # solve function
    function solve!(ux, uy, uz, p, b, t, inversion_toolkit, solver_evolution, i_step, n_steps)
        @showprogress for i ∈ i_step:n_steps
            # advection step
            evolve_adv!(inversion_toolkit, solver_evolution, ux, uy, uz, p, b)
            update_b!(b, solver_evolution)

            # diffusion step
            evolve_diff!(solver_evolution, b)
            update_b!(b, solver_evolution)

            # invert
            invert!(inversion_toolkit, b)
            ux, uy, uz, p = update_u_p!(ux, uy, uz, p, inversion_toolkit.solver)

            # time
            t += Δt
        end
        return ux, uy, uz, p, b
    end

    # run
    i_step = Int64(round(t/Δt)) + 1
    n_steps = Int64(round(T/Δt))
    ux, uy, uz, p, b = solve!(ux, uy, uz, p, b, t, inversion_toolkit, solver_evolution, i_step, n_steps)

    # # plot for sanity check
    # sim_plots(dim, ux, uy, uz, b, N², H, t, 0, "test")

    # compare state with data
    datafile = @sprintf("test/data/evolution_%s.h5", dim)
    if !isfile(datafile)
        @warn "Data file not found, saving state..."
        save_state(ux, uy, uz, p, b, t; fname=@sprintf("test/data/evolution_%s.h5", dim))
    else
        ux_data, uy_data, uz_data, p_data, b_data, t_data = load_state(@sprintf("test/data/evolution_%s.h5", dim))
        @test isapprox(ux.free_values, ux_data, rtol=1e-2)
        @test isapprox(uy.free_values, uy_data, rtol=1e-2)
        @test isapprox(uz.free_values, uz_data, rtol=1e-2)
        @test isapprox(p.free_values,  p_data,  rtol=1e-2)
        @test isapprox(b.free_values,  b_data,  rtol=1e-2)
    end

    # remove temporary files
    rm("LHS_adv_temp.h5", force=true)
    rm("LHS_diff_temp.h5", force=true)
end

@testset "Evolution Tests" begin
    @testset "2D CPU" begin
        coarse_evolution(TwoD(), CPU())
    end
    @testset "2D GPU" begin
        coarse_evolution(TwoD(), GPU())
    end
    @testset "3D CPU" begin
        coarse_evolution(ThreeD(), CPU())
    end
    @testset "3D GPU" begin
        coarse_evolution(ThreeD(), GPU())
    end
end