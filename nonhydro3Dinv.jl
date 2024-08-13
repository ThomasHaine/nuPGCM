using NonhydroPG
using Gridap, GridapGmsh
using CUDA, CUDA.CUSPARSE, CUDA.CUSOLVER
using Krylov, SparseArrays, LinearAlgebra
using HDF5, Printf
using PyPlot

pygui(false)
plt.style.use("plots.mplstyle")
plt.close("all")

out_folder = "out"

if !isdir(out_folder)
    mkdir(out_folder)
    mkdir("$out_folder/images")
    mkdir("$out_folder/data")
end

function save(ux, uy, uz, p, b)
    # fname = "$out_folder/data/nonhydro3D.vtu"
    # writevtk(Ω, fname, cellfields=["u"=>ux, "v"=>uy, "w"=>uz, "p"=>p, "b"=>b])
    # println(fname)
    fname = "$out_folder/data/nonhydro3D.h5"
    h5open(fname, "w") do file
        write(file, "ux", ux.free_values)
        write(file, "uy", uy.free_values)
        write(file, "uz", uz.free_values)
        write(file, "p", Vector(p.free_values))
    end
    println(fname)
end

# model
hres = 0.02
model = GmshDiscreteModel(@sprintf("meshes/bowl3D_%0.2f.msh", hres))

# mesh res
pts, conns = get_p_t(model)
h1 = [norm(pts[conns[i, 1], :] - pts[conns[i, 2], :]) for i ∈ axes(conns, 1)]
h2 = [norm(pts[conns[i, 2], :] - pts[conns[i, 3], :]) for i ∈ axes(conns, 1)]
h3 = [norm(pts[conns[i, 3], :] - pts[conns[i, 4], :]) for i ∈ axes(conns, 1)]
h4 = [norm(pts[conns[i, 4], :] - pts[conns[i, 1], :]) for i ∈ axes(conns, 1)]
hmin = minimum([h1; h2; h3; h4])
hmax = maximum([h1; h2; h3; h4])

# reference FE 
reffe_ux = ReferenceFE(lagrangian, Float64, 2)
reffe_uy = ReferenceFE(lagrangian, Float64, 2)
reffe_uz = ReferenceFE(lagrangian, Float64, 2)
reffe_p  = ReferenceFE(lagrangian, Float64, 1)

# test FESpaces
Vx = TestFESpace(model, reffe_ux, conformity=:H1, dirichlet_tags=["bot"])
Vy = TestFESpace(model, reffe_uy, conformity=:H1, dirichlet_tags=["bot"])
Vz = TestFESpace(model, reffe_uz, conformity=:H1, dirichlet_tags=["bot", "sfc"])
Q  = TestFESpace(model, reffe_p,  conformity=:H1, constraint=:zeromean)
Y = MultiFieldFESpace([Vx, Vy, Vz, Q])

# trial FESpaces with Dirichlet values
Ux = TrialFESpace(Vx, [0])
Uy = TrialFESpace(Vy, [0])
Uz = TrialFESpace(Vz, [0, 0])
P  = TrialFESpace(Q)
X  = MultiFieldFESpace([Ux, Uy, Uz, P])
nx = Ux.space.nfree
ny = Uy.space.nfree
nz = Uz.space.nfree
nu = nx + ny + nz
np = P.space.space.nfree
N = nu + np - 1
@printf("\nN = %d (%d + %d) ∼ 10^%d DOF\n", N, nu, np-1, floor(log10(N)))

# initialize vectors
ux = interpolate_everywhere(0, Ux)
uy = interpolate_everywhere(0, Uy)
uz = interpolate_everywhere(0, Uz)
p  = interpolate_everywhere(0, P)

# triangulation and integration measure
Ω = Triangulation(model)
dΩ = Measure(Ω, 4)

# gradients 
∂x(u) = VectorValue(1.0, 0.0, 0.0)⋅∇(u)
∂y(u) = VectorValue(0.0, 1.0, 0.0)⋅∇(u)
∂z(u) = VectorValue(0.0, 0.0, 1.0)⋅∇(u)

# depth
H(x) = 1 - x[1]^2 - x[2]^2

# forcing
ν(x) = 1

# params
ε² = 1e-2
γ = 1
f₀ = 1
β = 0
f(x) = f₀ + β*x[2]
println("\n---")
println("Parameters:\n")
@printf("ε² = %.1e (δ = %.1e, %.1e ≤ h ≤ %.1e)\n", ε², √(2ε²), hmin, hmax)
@printf("f₀ = %.1e\n", f₀)
@printf(" β = %.1e\n", β)
@printf(" γ = %.1e\n", γ)
println("---\n")

# filenames for LHS matrices
LHS_inversion_fname = @sprintf("matrices/LHS_inversion_%e_%e_%e_%e_%e.h5", hres, ε², γ, f₀, β)

# inversion LHS
function assemble_LHS_inversion()
    a_inversion((ux, uy, uz, p), (vx, vy, vz, q)) = 
        ∫( γ*ε²*∂x(ux)*∂x(vx)*ν +   γ*ε²*∂y(ux)*∂y(vx)*ν +   ε²*∂z(ux)*∂z(vx)*ν - uy*vx*f + ∂x(p)*vx +
           γ*ε²*∂x(uy)*∂x(vy)*ν +   γ*ε²*∂y(uy)*∂y(vy)*ν +   ε²*∂z(uy)*∂z(vy)*ν + ux*vy*f + ∂y(p)*vy +
         γ^2*ε²*∂x(uz)*∂x(vz)*ν + γ^2*ε²*∂y(uz)*∂y(vz)*ν + γ*ε²*∂z(uz)*∂z(vz)*ν +           ∂z(p)*vz +
                                                                      ∂x(ux)*q + ∂y(uy)*q + ∂z(uz)*q )dΩ
    @time "assemble LHS_inversion" LHS_inversion = assemble_matrix(a_inversion, X, Y)
    write_sparse_matrix(LHS_inversion_fname, LHS_inversion)
    return LHS_inversion
end

if isfile(LHS_inversion_fname)
    LHS_inversion = read_sparse_matrix(LHS_inversion_fname)
else
    LHS_inversion = assemble_LHS_inversion()
end

# Cuthill-McKee DOF reordering
@time "RCM perm" begin 
a_m(u, v) = ∫( u*v )dΩ
M_ux = assemble_matrix(a_m, Ux, Vx)
M_uy = assemble_matrix(a_m, Uy, Vy)
M_uz = assemble_matrix(a_m, Uz, Vz)
M_p  = assemble_matrix(a_m, P, Q)
perm_ux = CUSOLVER.symrcm(M_ux) .+ 1
perm_uy = CUSOLVER.symrcm(M_uy) .+ 1
perm_uz = CUSOLVER.symrcm(M_uz) .+ 1
perm_p  = CUSOLVER.symrcm(M_p)  .+ 1
perm_inversion = [perm_ux; 
                  perm_uy .+ nx; 
                  perm_uz .+ nx .+ ny; 
                  perm_p  .+ nx .+ ny .+ nz]
inv_perm_inversion = invperm(perm_inversion)
# plot_sparsity_pattern(LHS_inversion, fname="$out_folder/images/LHS_inversion.png")
LHS_inversion = LHS_inversion[perm_inversion, perm_inversion]
# plot_sparsity_pattern(LHS_inversion, fname="$out_folder/images/LHS_inversion_symrcm.png")
end

# put on GPU
LHS_inversion = CuSparseMatrixCSR(LHS_inversion)
CUDA.memory_status()

# Krylov solver
# solver_inversion = DqgmresSolver(N, N, 20, CuVector{Float64})
solver_inversion = GmresSolver(N, N, 20, CuVector{Float64})
solver_inversion.x .= 0
tol = 1e-8
itmax = 0
# tol = 1e-7
# itmax = 2000
@printf("tol = %.1e, itmax = %d\n", tol, itmax)

# inversion functions
function invert!(solver_inversion, b)
    l_inversion((vx, vy, vz, q)) = ∫( b*vz )dΩ
    @time "build RHS_inversion" RHS_inversion = CuVector(assemble_vector(l_inversion, Y)[perm_inversion])
    @time "invert!" Krylov.solve!(solver_inversion, LHS_inversion, RHS_inversion, solver_inversion.x, 
                                  verbose=1, atol=tol, rtol=tol, itmax=itmax, restart=true)
    return solver_inversion
end
function update_u_p!(ux, uy, uz, p, solver_inversion)
    sol = Vector(solver_inversion.x[inv_perm_inversion])
    ux.free_values .= sol[1:nx]
    uy.free_values .= sol[nx+1:nx+ny]
    uz.free_values .= sol[nx+ny+1:nx+ny+nz]
    p = FEFunction(P, sol[nx+ny+nz+1:end])
    return ux, uy, uz, p
end

b(x) = x[3] + 0.1*exp(-(x[3] + H(x))/0.1)
solver_inversion = invert!(solver_inversion, b)
ux, uy, uz, p = update_u_p!(ux, uy, uz, p, solver_inversion)
@time "profiles" plot_profiles(ux, uy, uz, b, 0.5, 0.0, H; t=0, fname="$out_folder/images/profiles.png")
@time "save" save(ux, uy, uz, p, b)