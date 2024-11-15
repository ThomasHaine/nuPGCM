# Experiment Set 1
- 2D
- $h = 0.01$
- $b$ from spinup simulation
- Preconditioners of the form $P = \left(\begin{array}{c c} \tilde A^{-1} & 0\\0 & \varepsilon^{2} \tilde M_p^{-1} \end{array}\right)$
    - $\tilde M_p^{-1} = $ 4 `cg` iterations preconditioned by `diag(M_p)` unless otherwise stated
    - Using GMRES with `atol = 1e-6`, `rtol = 1.5e-8`, `k=20` unless otherwise stated

### Classic Stokes ($\varepsilon^2 = 1$, $\gamma = 1$, $f = 0$)

__NOTE__ : $A$ is SPD

See `*1.png` for images of solution.

- $\tilde A^{-1} =$ `lu(A)`: converges in __15__ iterations 
- $\tilde A^{-1} =$ `ilu(A, τ=1e-5)`: converges in __17__ iterations
- $\tilde A^{-1} =$ `ilu(A, τ=1e-4)`: converges in __26__ iterations
- $\tilde A^{-1} =$ `ilu(A, τ=1e-3)`: converges in __72__ iterations
- no blocks, just `1/h^2` normalization: __53,581__ iterations (steady slope)

### Aspect Ratio Stokes ($\varepsilon^2 = 1$, $\gamma = 1/4$, $f = 0$)

__NOTE__ : $A$ is SPD

See `*2.png` for images of solution.

- $\tilde A^{-1} =$ `lu(A)`: converges in __46__ iterations 
- $\tilde A^{-1} =$ `ilu(A, τ=1e-3)`: converges in __218__ iterations

### Less Diff Stokes ($\varepsilon^2 = 10^{-4}$, $\gamma = 1$, $f = 0$)

__NOTE__ : $A$ is _neither_ symmetric _nor_ positive definite!

See `*3.png` for images of solution.

- $\tilde A^{-1} =$ `lu(A)`: converges in __15__ iterations 
- $\tilde A^{-1} =$ `ilu(A, τ=1e-10)`: converges in __26__ iterations
- $\tilde A^{-1} =$ `ilu(A, τ=1e-9)`: converges in __30__ iterations
- $\tilde A^{-1} =$ `ilu(A, τ=1e-8)`: converges in __55__ iterations
- $\tilde A^{-1} =$ `ilu(A, τ=1e-7)`: converges in __185__ iterations

### PG Thick BL ($\varepsilon^2 = 1$, $\gamma = 1$, $f = 1$)

__NOTE__ : $A$ is _neither_ symmetric _nor_ positive definite!

See `*4.png` for images of solution.

- $\tilde A^{-1} =$ `lu(A)`: converges in __15__ iterations 
- $\tilde A^{-1} =$ `ilu(A, τ=1e-3)`: converges in __72__ iterations

### PG Thin BL ($\varepsilon^2 = 10^{-4}$, $\gamma = 1$, $f = 1$)

__NOTE__ : $A$ is _neither_ symmetric _nor_ positive definite!

See `*5.png` for images of solution.

- $\tilde A^{-1} =$ `lu(A)`: "converges" in __1,568__ iterations 
    - rapid convergence in first 100 steps followed by very slow progress
    - not actually converged to true solution (see `*5a.png`)
- no blocks, just `1/h^2` normalization: __308,990__ iterations (steady slope)

### PG Thin BL, Aspect Ratio ($\varepsilon^2 = 10^{-4}$, $\gamma = 1/4$, $f = 1$)

---

Now try $P = \left(\begin{array}{c c} \tilde A^{-1} & B^T\\0 & -\varepsilon^{2} \tilde M_p^{-1} \end{array}\right)$

- Classic Stokes: __7__ iterations 
- Aspect Ratio Stokes: __11__ iterations 
- Less Diff Stokes: __10__ iterations 
-  PG Thick BL: __7__ iterations 
- PG Thin BL: __1,281__ iterations (same failure)