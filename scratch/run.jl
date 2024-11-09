using nuPGCM
using Gridap, GridapGmsh
using IncompleteLU, Krylov, LinearOperators, CuthillMcKee
using CUDA, CUDA.CUSPARSE, CUDA.CUSOLVER
using SparseArrays, LinearAlgebra, Statistics
using Printf
using PyPlot

pygui(false)
plt.style.use("../plots.mplstyle")
plt.close("all")

# out_folder = "../sims/sim041"
out_folder = "../out"

if !isdir(out_folder)
    println("creating folder: ", out_folder)
    mkdir(out_folder)
end
if !isdir("$out_folder/images")
    println("creating subfolder: ", out_folder, "/images")
    mkdir("$out_folder/images")
end
if !isdir("$out_folder/data")
    println("creating subfolder: ", out_folder, "/data")
    mkdir("$out_folder/data")
end
flush(stdout)
flush(stderr)

# choose dimensions
dim = TwoD()
# dim = ThreeD()

# choose architecture
# arch = CPU()
arch = GPU()

# tolerance and max iterations for iterative solvers
tol = 1e-6
@printf("tol = %.1e\n", tol)
itmax = 0
@printf("itmax = %d\n\n", itmax)

# Vector type 
VT = typeof(arch) == CPU ? Vector{Float64} : CuVector{Float64}

# model
hres = 0.01
# model = GmshDiscreteModel(@sprintf("../meshes/bowl%s_%0.2f.msh", dim, hres))
model = GmshDiscreteModel(@sprintf("../meshes/bowl%s_%0.2f_dm.msh", dim, hres))

# full grid
m = Mesh(model)

# surface grid
m_sfc = Mesh(model, "sfc")

# mesh res
hs = [norm(m.p[m.t[i, j], :] - m.p[m.t[i, mod1(j+1, dim.n+1)], :]) for i ∈ axes(m.t, 1), j ∈ 1:dim.n+1]
hmin = minimum(hs)
hmax = maximum(hs)
h = mean(hs)
@printf("\n%.1e < h < %.1e (mean = %.1e)\n", hmin, hmax, h)

# FE spaces
X, Y, B, D = setup_FESpaces(model)
Ux, Uy, Uz, P = unpack_spaces(X)
nx = Ux.space.nfree
ny = Uy.space.nfree
nz = Uz.space.nfree
nu = nx + ny + nz
np = P.space.space.nfree
nb = B.space.nfree
N = nu + np - 1
@printf("\nN = %d (%d + %d) ∼ 10^%d DOF\n", N, nu, np-1, floor(log10(N)))

# triangulation and integration measure
Ω = Triangulation(model)
dΩ = Measure(Ω, 4)

# depth
H(x) = 1 - x[1]^2 - x[2]^2

# forcing
ν(x) = 1
κ(x) = 1e-2 + exp(-(x[3] + H(x))/0.1)

# params
ε² = 1e-4
γ = 1/4
f₀ = 1
# f₀ = 0
β = 0
f(x) = f₀ + β*x[2]
μϱ = 1e0
Δt = 1e-2
# Δt = 1e-4*μϱ/ε²
T = 5e-2*μϱ/ε²
α = Δt/2*ε²/μϱ # for timestep
println("\n---")
println("Parameters:\n")
@printf("ε² = %.1e (δ = %.1e, %.1e ≤ h ≤ %.1e)\n", ε², √(2ε²), hmin, hmax)
@printf("f₀ = %.1e\n", f₀)
@printf(" β = %.1e\n", β)
@printf(" γ = %.1e\n", γ)
@printf("μϱ = %.1e\n", μϱ)
@printf("Δt = %.1e\n", Δt)
@printf(" T = %.1e\n", T)
println("---\n")

# filenames for LHS matrices
LHS_inversion_fname = @sprintf("../matrices/LHS_inversion_%s_dm_%e_%e_%e_%e_%e.h5", dim, hres, ε², γ, f₀, β)
LHS_diff_fname = @sprintf("../matrices/LHS_diff_%s_dm_%e_%e_%e.h5", dim, hres, α, γ)
LHS_adv_fname = @sprintf("../matrices/LHS_adv_%s_dm_%e.h5", dim, hres)

# inversion LHS
if isfile(LHS_inversion_fname)
    LHS_inversion, perm_inversion, inv_perm_inversion = read_sparse_matrix(LHS_inversion_fname)
else
    LHS_inversion, perm_inversion, inv_perm_inversion = assemble_LHS_inversion(arch, dim, γ, ε², ν, f, X, Y, dΩ; fname=LHS_inversion_fname)
end

# inversion RHS
RHS_inversion = assemble_RHS_inversion(perm_inversion, B, Y, dΩ)

# preconditioner
if typeof(dim) == TwoD
    # if typeof(arch) == GPU 
    #     @time "lu(LHS_inversion)" P_inversion = CUSOLVER.lu(LHS_inversion) #FIXME: this gives an LU factor on the CPU, can't be used by Krylov
    # else
    #     @time "lu(LHS_inversion)" P_inversion = lu(LHS_inversion)
    # end
    # ldiv_P_inversion = true
    P_inversion = Diagonal(on_architecture(arch, 1/h^2*ones(N)))
    ldiv_P_inversion = false
else
    P_inversion = Diagonal(on_architecture(arch, 1/h^3*ones(N)))
    ldiv_P_inversion = false
end

# put on GPU, if needed
LHS_inversion = on_architecture(arch, LHS_inversion)
RHS_inversion = on_architecture(arch, RHS_inversion)
if typeof(arch) == GPU
    CUDA.memory_status()
    println()
end

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
                  atol=tol, rtol=tol, verbose=1, itmax=itmax, restart=true,
                  history=true)
    @printf("inversion GMRES solve: solved=%s, niter=%d, time=%f\n", solver.stats.solved, solver.stats.niter, solver.stats.timer)
    # solver.x = P_inversion \ RHS
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

flush(stdout)
flush(stderr)

# background state \partial_z b = N^2
N² = 1.
# N² = 0.

# # initial condition: b = N^2 z, t = 0
# i_save = 0
# b = interpolate_everywhere(0, B)
# t = 0.
# ux = interpolate_everywhere(0, Ux)
# uy = interpolate_everywhere(0, Uy)
# uz = interpolate_everywhere(0, Uz)
# p  = interpolate_everywhere(0, P) 
# save_state(ux, uy, uz, p, b, t; fname=@sprintf("%s/data/state%03d.h5", out_folder, i_save))

# # initial condition: b = gaussian, t = 0
# i_save = 0
# b = interpolate_everywhere(x -> -exp(-100*(x[1]^2 + x[2]^2 + (x[3] + 0.5)^2)), B)
# t = 0.
# ux = interpolate_everywhere(0, Ux)
# # ux = interpolate_everywhere(0.1, Ux)
# uy = interpolate_everywhere(0, Uy)
# uz = interpolate_everywhere(0, Uz)
# p  = interpolate_everywhere(0, P)
# solver_inversion = invert!(arch, solver_inversion, b)
# ux, uy, uz, p = update_u_p!(ux, uy, uz, p, solver_inversion)
# save_state(ux, uy, uz, p, b, t; fname=@sprintf("%s/data/state%03d.h5", out_folder, i_save))

# initial condition: load from file
i_save = 1
# statefile = @sprintf("%s/data/state%03d.h5", out_folder, i_save)
statefile = @sprintf("../sims/sim040/data/state001.h5")
ux, uy, uz, p, b, t = load_state(statefile)
solver_inversion.x .= on_architecture(arch, [ux; uy; uz; p][perm_inversion])
ux = FEFunction(Ux, ux)
uy = FEFunction(Uy, uy)
uz = FEFunction(Uz, uz)
p  = FEFunction(P, p)
# ux = interpolate_everywhere(0, Ux)
# uy = interpolate_everywhere(0, Uy)
# uz = interpolate_everywhere(0, Uz)
# p  = interpolate_everywhere(0, P)
b  = FEFunction(B, b)

# # intialize solver_inversion
# sol = [ux.free_values; uy.free_values; uz.free_values; p.free_values]
# solver_inversion.x .= on_architecture(arch, sol[perm_inversion])

# sif = load("$out_folder/data/solver_inversion_f.jld2", "solver_inversion")
# solver_inversion.x .= on_architecture(arch, sif.x)
solver_inversion = invert!(arch, solver_inversion, b)
ux, uy, uz, p = update_u_p!(ux, uy, uz, p, solver_inversion)

# plot initial condition
plots_cache = sim_plots(dim, ux, uy, uz, b, N², H, t, i_save, out_folder)
i_save += 1

# plot gmres error vs iterations
fig, ax = subplots(1)
solver_inversion_f = load("$out_folder/data/solver_inversion_f.jld2", "solver_inversion")
solver_inversion_f_init = load("$out_folder/data/solver_inversion_f_init.jld2", "solver_inversion")
solver_inversion = load("$out_folder/data/solver_inversion.jld2", "solver_inversion")
ax.plot(solver_inversion_f.stats.residuals, "-", label=L"$f = 1$")
ax.plot(solver_inversion_f_init.stats.residuals, "-", label=L"$f = 1$, true init")
ax.plot(solver_inversion.stats.residuals, "-", label=L"$f = 0$")
ax.legend()
ax.set_yscale("log")
ax.set_xlabel("Iteration")
ax.set_ylabel("Residual")
savefig("$out_folder/images/gmres_convergence.png")
println("$out_folder/images/gmres_convergence.png")
plt.close()

error()

# evolution LHSs
if isfile(LHS_adv_fname) && isfile(LHS_diff_fname)
    LHS_adv,  perm_b, inv_perm_b = read_sparse_matrix(LHS_adv_fname)
    LHS_diff, perm_b, inv_perm_b = read_sparse_matrix(LHS_diff_fname)
else
    LHS_adv, LHS_diff, perm_b, inv_perm_b = assemble_LHS_adv_diff(arch, α, γ, κ, B, D, dΩ; fname_adv=LHS_adv_fname, fname_diff=LHS_diff_fname)
end

# diffusion RHS matrix and vector
RHS_diff, rhs_diff = assemble_RHS_diff(perm_b, α, γ, κ, N², B, D, dΩ)

# preconditioners
if typeof(dim) == TwoD
    # if typeof(arch) == GPU 
    #     @time "lu(LHS_diff)" P_diff = CUSOLVER.lu(LHS_diff)
    #     @time "lu(LHS_adv)"  P_adv  = CUSOLVER.lu(LHS_adv)
    # else
    #     @time "lu(LHS_diff)" P_diff = lu(LHS_diff)
    #     @time "lu(LHS_adv)"  P_adv  = lu(LHS_adv)
    # end
    # ldiv_P_diff = true
    # ldiv_P_adv  = true
    P_diff = Diagonal(on_architecture(arch, Vector(1 ./ diag(LHS_diff))))
    P_adv  = Diagonal(on_architecture(arch, Vector(1 ./ diag(LHS_adv))))
    ldiv_P_diff = false
    ldiv_P_adv  = false
else
    P_diff = Diagonal(on_architecture(arch, Vector(1 ./ diag(LHS_diff))))
    P_adv  = Diagonal(on_architecture(arch, Vector(1 ./ diag(LHS_adv))))
    ldiv_P_diff = false
    ldiv_P_adv  = false
end

# put on GPU, if needed
LHS_diff = on_architecture(arch, LHS_diff)
RHS_diff = on_architecture(arch, RHS_diff)
rhs_diff = on_architecture(arch, rhs_diff)
LHS_adv = on_architecture(arch, LHS_adv)
if typeof(arch) == GPU
    CUDA.memory_status()
    println()
end

# Krylov solver for evolution
solver_evolution = CgSolver(nb, nb, VT)
solver_evolution.x .= on_architecture(arch, copy(b.free_values)[perm_b])

# evolution functions
∂x(u) = VectorValue(1.0, 0.0, 0.0)⋅∇(u)
∂y(u) = VectorValue(0.0, 1.0, 0.0)⋅∇(u)
∂z(u) = VectorValue(0.0, 0.0, 1.0)⋅∇(u)
# assembler = SparseMatrixAssembler(D, D)
# RHS_evolution = zeros(nb)
b_half = interpolate_everywhere(0, B)
function evolve_adv!(arch::AbstractArchitecture, solver_inversion, solver_evolution, ux, uy, uz, p, b)
    # half step
    l_half(d) = ∫( b*d - Δt/2*(ux*∂x(b) + uy*∂y(b) + uz*(N² + ∂z(b)))*d )dΩ
    @time "build RHS_evolution 1" RHS = on_architecture(arch, assemble_vector(l_half, D)[perm_b])
    # @time "build RHS_evolution" Gridap.FESpaces.assemble_vector!(l, RHS_evolution, assembler, D)
    RHS = on_architecture(arch, RHS_evolution[perm_b])
    Krylov.solve!(solver_evolution, LHS_adv, RHS, solver_evolution.x, M=P_adv, ldiv=ldiv_P_adv, atol=tol, rtol=tol, verbose=0, itmax=itmax)
    @printf("advection CG solve 1: solved=%s, niter=%d, time=%f\n", solver_evolution.stats.solved, solver_evolution.stats.niter, solver_evolution.stats.timer)
    # solver_evolution.x .= P_adv \ RHS

    # u, v, w, p, b at half step
    update_b!(b_half, solver_evolution)
    solver_inversion = invert!(arch, solver_inversion, b_half)
    ux, uy, uz, p = update_u_p!(ux, uy, uz, p, solver_inversion)

    # full step
    l_full(d) = ∫( b*d - Δt*(ux*∂x(b_half) + uy*∂y(b_half) + uz*(N² + ∂z(b_half)))*d )dΩ
    @time "build RHS_evolution 2" RHS = on_architecture(arch, assemble_vector(l_full, D)[perm_b])
    Krylov.solve!(solver_evolution, LHS_adv, RHS, solver_evolution.x, M=P_adv, ldiv=ldiv_P_adv, atol=tol, rtol=tol, verbose=0, itmax=itmax)
    @printf("advection CG solve 2: solved=%s, niter=%d, time=%f\n", solver_evolution.stats.solved, solver_evolution.stats.niter, solver_evolution.stats.timer)
    # solver_evolution.x .= P_adv \ RHS

    return solver_inversion, solver_evolution
end
function evolve_diff!(arch::AbstractArchitecture, solver, b)
    b_arch= on_architecture(arch, b.free_values)
    RHS = RHS_diff*b_arch + rhs_diff
    Krylov.solve!(solver, LHS_diff, RHS, solver.x, M=P_diff, ldiv=ldiv_P_diff, atol=tol, rtol=tol, verbose=0, itmax=itmax)
    @printf("diffusion CG solve: solved=%s, niter=%d, time=%f\n", solver.stats.solved, solver.stats.niter, solver.stats.timer)
    # solver.x .= P_diff \ RHS
    return solver
end
function update_b!(b, solver)
    b.free_values .= on_architecture(CPU(), solver.x[inv_perm_b])
    return b
end

# solve function
function solve!(arch::AbstractArchitecture, ux, uy, uz, p, b, t, solver_inversion, solver_evolution, i_save, i_step, n_steps)
    t0 = time()
    for i ∈ i_step:n_steps
        flush(stdout)
        flush(stderr)

        # advection step
        solver_inversion, solver_evolution = evolve_adv!(arch, solver_inversion, solver_evolution, ux, uy, uz, p, b)
        b = update_b!(b, solver_evolution)

        # diffusion step
        solver_evolution = evolve_diff!(arch, solver_evolution, b)
        b = update_b!(b, solver_evolution)

        # # invert
        # solver_inversion = invert!(arch, solver_inversion, b)
        # ux, uy, uz, p = update_u_p!(ux, uy, uz, p, solver_inversion)

        # blow up
        if any(isnan.(solver_inversion.x)) || any(isnan.(solver_evolution.x))
            # save and kill
            save_state(ux, uy, uz, p, b, t; fname=@sprintf("%s/data/state%03d.h5", out_folder, i_save))
            sim_plots(plots_cache, ux, uy, uz, b, t, i_save, out_folder)
            error("Solution diverged: NaN(s) found.")
        end

        # time
        t += Δt

        # info
        ux_max = maximum(abs.(ux.free_values))
        uy_max = maximum(abs.(uy.free_values))
        uz_max = maximum(abs.(uz.free_values))
        t1 = time()
        println("\n---")
        @printf("t = %f (i = %d/%d, Δt = %f)\n\n", t, i, n_steps, Δt)
        @printf("time elapsed: %02d:%02d:%02d\n", hrs_mins_secs(t1-t0)...)
        @printf("estimated time remaining: %02d:%02d:%02d\n", hrs_mins_secs((t1-t0)*(n_steps-i)/(i-i_step+1))...)
        @printf("|u|ₘₐₓ = %.1e, %.1e ≤ b′ ≤ %.1e\n", max(ux_max, uy_max, uz_max), minimum([b.free_values; 0]), maximum([b.free_values; 0]))
        @printf("CFL ≈ %f\n", min(hmin/ux_max, hmin/uy_max, hmin/uz_max))
        println("---\n")

        # blow up
        if max(ux_max, uy_max, uz_max) > 10
            # save and kill
            save_state(ux, uy, uz, p, b, t; fname=@sprintf("%s/data/state%03d.h5", out_folder, i_save))
            sim_plots(plots_cache, ux, uy, uz, b, t, i_save, out_folder)
            error("Solution diverged: |u|ₘₐₓ > 10.")
        end

        # save/plot
        # if mod(i, n_steps ÷ 50) == 0
        #     save_state(ux, uy, uz, p, b, t; fname=@sprintf("%s/data/state%03d.h5", out_folder, i_save))
        #     sim_plots(plots_cache, ux, uy, uz, b, t, i_save, out_folder)
        #     i_save += 1
        # end
    end
    return ux, uy, uz, p, b
end

# # run
# i_step = Int64(round(t/Δt)) + 1
# n_steps = Int64(round(T/Δt))
# ux, uy, uz, p, b = solve!(arch, ux, uy, uz, p, b, t, solver_inversion, solver_evolution, i_save, i_step, n_steps)

# # compute error
# b_final = interpolate_everywhere(x -> -exp(-100*((x[1] - 0.1*T)^2 + x[2]^2 + (x[3] + 0.5)^2)), B)
# err_L2 = sqrt(sum(∫( (b - b_final)*(b - b_final) )dΩ))
# err_Linf = maximum(abs.(b.free_values - b_final.free_values))
# @printf("L2 error = %.1e\n", err_L2)
# @printf("Max error = %.1e\n", err_Linf)
# plot_slice(abs(b - b_final), b, N²; y=0, t=T, cb_label=L"Error $|b - b_a|$", fname=@sprintf("%s/images/error.png", out_folder))

# dts = [0.04, 0.02, 0.01]
# errors_L2   = [5.7e-5, 2.8e-5, 2.5e-5]
# errors_Linf = [5.7e-4, 3.0e-4, 3.0e-4]
# # dts = [5e-1, 2.5e-1, 1e-1, 5e-2, 2.5e-2, 1.25e-2]
# # errors = [7.9e-2, 2.1e-2, 3.6e-3, 1.0e-3, 4.4e-4, 3.4e-4]
# fig, ax = subplots(1, figsize=(3.2, 3.2))
# ax.spines["top"].set_visible(true)
# ax.spines["right"].set_visible(true)
# ax.loglog(dts, errors_L2, "o-")
# ax.loglog(dts, errors_Linf, "o-")
# ax.loglog(dts, errors_L2[3]/dts[3]^2*dts.^2, "k--", label=L"O(\Delta t^2)")
# ax.legend()
# ax.set_xlabel(L"Timestep $\Delta t$")
# # ax.set_xlim(1e-2, 1e0)
# # ax.set_ylim(1e-4, 5e-1)
# ax.set_ylabel(L"Error $||b - b_a||_\infty$")
# ax.grid(true, which="both", color="k", alpha=0.5, linestyle=":", linewidth=0.25)
# ax.set_axisbelow(true)
# ax.set_title(L"Global Trucation Error at $T = 0.4$")
# savefig("../out/images/convergence.png")
# println("../out/images/convergence.png")
# plt.close()

# dt_min = 0.005
# statefile = @sprintf("%s/data/state_%.3f.h5", out_folder, dt_min)
# ux0, uy0, uz0, p0, b0, t = load_state(statefile)
# dts = [0.1, 0.05, 0.04, 0.02, 0.01]
# errors = zeros(size(dts))
# for i in eachindex(dts)
#     statefile = @sprintf("%s/data/state_%.3f.h5", out_folder, dts[i])
#     ux, uy, uz, p, b, t = load_state(statefile)
#     errors[i] = maximum(abs.(b - b0))
# end
# fig, ax = subplots(1, figsize=(3.2, 3.2))
# ax.spines["top"].set_visible(true)
# ax.spines["right"].set_visible(true)
# ax.loglog(dts, errors[4]/dts[4]^2*dts.^2, "k--", label=L"O(\Delta t^2)", lw=1.0)
# ax.loglog(dts, errors[3]/dts[3]^1*dts.^1, "k--", alpha=0.5, label=L"O(\Delta t)", lw=1.0)
# ax.loglog(dts, errors, "o-")
# ax.legend()
# ax.set_xlabel(L"Timestep $\Delta t$")
# ax.set_xlim(9e-3, 2e-1)
# ax.set_ylim(1e-6, 1e-3)
# ax.set_ylabel(L"Error $||b - b_a||_\infty$")
# ax.grid(true, which="both", color="k", alpha=0.5, linestyle=":", linewidth=0.25)
# ax.set_axisbelow(true)
# ax.set_title(L"Global Trucation Error at $T = 5$")
# savefig("../out/images/convergence_spinup.png")
# println("../out/images/convergence_spinup.png")
# plt.close()

println("Done.")