# Benchmark Plan — DPIM Higher Orders (`beam` demo)

## Goal

Instrument the higher-order loop in `source/dpim.jl` (lines 107–126) to produce
per-monomial timing and allocation tables, readable in Python.

## Configuration

| Setting | Value |
|---------|-------|
| mesh file | `input/beam.mphtxt` |
| info file | `input/beam.jl` (single master mode, `info.max_order = 7`) |
| FEM DOFs (`nK`) | determined at runtime from mesh |
| master modes (`nm`) | 1 |
| reduced vars (`nz`) | 2 |
| forced vars (`nzforce`) | 2 |

Set these in the launch script before the DPIM call:

```julia
mesh_file = "beam.mphtxt"
info_file = "beam.jl"
```

## What to measure

### Per-order operations (one call per order `p`)

| Symbol | Function | Notes |
|--------|----------|-------|
| `fillrhsG` | `fillrhsG!(mesh,U,Cp,p)` | Quadratic geometric nonlinearity RHS |
| `fillrhsH` | `fillrhsH!(mesh,U,Cp,p)` | Cubic geometric nonlinearity RHS |
| `fillWf` | `fillWf!(Cp,p,info)` | Autonomous forcing (Wf) accumulation |

### Per-monomial operations (inner loop over `i in 1:Cp[p].nc`)

| Symbol | Condition | Function | Notes |
|--------|-----------|----------|-------|
| `fillWfnonaut` | `corresp > 0` | `fillWfnonaut!(Cp,p,Avector,i,info)` | Non-autonomous Wf term |
| `homological` | `corresp > 0` | `homological_FULL!(...)` | Bordered linear system solve |
| `conj_copy` | `corresp < 0` | in-place conjugate assignment | Cheap symmetry shortcut |

## Instrumentation

Wrap each targeted call with Julia's `@timed` (returns `(value, time_s, bytes, gctime_s, gcstats)`).

### Modified higher-order block in `dpim.jl`

Replace lines 107–126 with:

```julia
println("Higher orders")
# --- open benchmark CSV files ---
f_order = open(info.output_dir*"/benchmark_per_order.csv","w")
f_mono  = open(info.output_dir*"/benchmark_per_monomial.csv","w")

write(f_order, "order,n_monomials,n_solve,n_conj,"*
               "fillrhsG_time_s,fillrhsG_alloc_bytes,"*
               "fillrhsH_time_s,fillrhsH_alloc_bytes,"*
               "fillWf_time_s,fillWf_alloc_bytes\n")

write(f_mono, "order,monomial_idx,alpha_vector,corresp,"*
              "fillWfnonaut_time_s,fillWfnonaut_alloc_bytes,"*
              "homological_time_s,homological_alloc_bytes,"*
              "conj_copy_time_s,conj_copy_alloc_bytes,"*
              "monomial_total_time_s,monomial_total_alloc_bytes\n")

for p = 2:info.max_order
    println("Order $p")

    # --- per-order timed calls ---
    (_, tG, bG) = @timed fillrhsG!(mesh,U,Cp,p)
    (_, tH, bH) = @timed fillrhsH!(mesh,U,Cp,p)
    (_, tWf,bWf)= @timed fillWf!(Cp,p,info)

    n_solve = count(x->x>0, Cp[p].corresp)
    n_conj  = count(x->x<0, Cp[p].corresp)
    write(f_order, "$p,$(Cp[p].nc),$n_solve,$n_conj,"*
                   "$tG,$bG,$tH,$bH,$tWf,$bWf\n")

    for i in 1:Cp[p].nc
        corresp = Cp[p].corresp[i]
        alpha   = string(Cp[p].Avector[i])
        t_fwnaut = NaN; b_fwnaut = 0
        t_hom    = NaN; b_hom    = 0
        t_conj   = NaN; b_conj   = 0

        if corresp > 0
            (_, t_fwnaut, b_fwnaut) = @timed fillWfnonaut!(Cp,p,Cp[p].Avector[i],i,info)
            (_, t_hom,    b_hom)    = @timed homological_FULL!(Sol_FULL,Rhs_FULL,Mat_FULL,
                                                                Λ,Cp[p],i,M,K,C,AY,XTA,info)
        elseif corresp < 0
            (_, t_conj, b_conj) = @timed begin
                Cp[p].W[:,i]                     = conj(Cp[p].W[:,-corresp])
                Cp[p].f[1:info.nm,i]             = conj(Cp[p].f[info.nm+1:2*info.nm,-corresp])
                Cp[p].f[info.nm+1:2*info.nm,i]  = conj(Cp[p].f[1:info.nm,-corresp])
            end
        end

        t_tot = sum(x -> isnan(x) ? 0.0 : x, [t_fwnaut, t_hom, t_conj])
        b_tot = b_fwnaut + b_hom + b_conj
        write(f_mono, "$p,$i,\"$alpha\",$corresp,"*
                      "$t_fwnaut,$b_fwnaut,"*
                      "$t_hom,$b_hom,"*
                      "$t_conj,$b_conj,"*
                      "$t_tot,$b_tot\n")
    end
end

close(f_order)
close(f_mono)
```

> **Note:** `@timed` in Julia ≥ 1.8 returns a `NamedTuple`; destructuring
> as `(_, time, bytes) = @timed expr` still works because the first three
> positional fields are `value`, `time`, `bytes`.

## Output files

Both files are written to `info.output_dir` (created by `make_output_dir`).

### `benchmark_per_order.csv`

One row per polynomial order.

| Column | Type | Description |
|--------|------|-------------|
| `order` | int | Polynomial order `p` |
| `n_monomials` | int | Total monomials `Cp[p].nc` |
| `n_solve` | int | Monomials solved via homological equation |
| `n_conj` | int | Monomials filled by conjugate symmetry |
| `fillrhsG_time_s` | float | Wall time for `fillrhsG!` (s) |
| `fillrhsG_alloc_bytes` | int | Heap allocation for `fillrhsG!` (bytes) |
| `fillrhsH_time_s` | float | Wall time for `fillrhsH!` (s) |
| `fillrhsH_alloc_bytes` | int | Heap allocation for `fillrhsH!` (bytes) |
| `fillWf_time_s` | float | Wall time for `fillWf!` (s) |
| `fillWf_alloc_bytes` | int | Heap allocation for `fillWf!` (bytes) |

### `benchmark_per_monomial.csv`

One row per monomial (all orders combined).

| Column | Type | Description |
|--------|------|-------------|
| `order` | int | Polynomial order `p` |
| `monomial_idx` | int | 1-based index within order |
| `alpha_vector` | str | Multi-index string, e.g. `[2, 1, 0, 0]` |
| `corresp` | int | Correspondence flag (>0 solve, <0 conj, 0 skip) |
| `fillWfnonaut_time_s` | float | Wall time (NaN if not applicable) |
| `fillWfnonaut_alloc_bytes` | int | Heap allocation (0 if not applicable) |
| `homological_time_s` | float | Wall time for bordered system solve (NaN if not applicable) |
| `homological_alloc_bytes` | int | Heap allocation for bordered system solve |
| `conj_copy_time_s` | float | Wall time for conjugate copy (NaN if not applicable) |
| `conj_copy_alloc_bytes` | int | Heap allocation for conjugate copy |
| `monomial_total_time_s` | float | Sum of applicable wall times |
| `monomial_total_alloc_bytes` | int | Sum of applicable allocations |

## Python reading snippet

```python
import pandas as pd
import numpy as np

df_order = pd.read_csv("benchmark_per_order.csv")
df_mono  = pd.read_csv("benchmark_per_monomial.csv", na_values="NaN")

# total time per order (per-order ops + sum of monomial ops)
order_mono_totals = df_mono.groupby("order")["monomial_total_time_s"].sum()
df_order = df_order.set_index("order")
df_order["monomial_loop_time_s"] = order_mono_totals
df_order["order_total_time_s"] = (
    df_order["fillrhsG_time_s"]
    + df_order["fillrhsH_time_s"]
    + df_order["fillWf_time_s"]
    + df_order["monomial_loop_time_s"]
)

# breakdown: which function dominates?
df_solve = df_mono[df_mono["corresp"] > 0].copy()
print(df_solve.groupby("order")[["fillWfnonaut_time_s","homological_time_s"]].describe())
```

## Execution

Run from the `MORFE2.0/` directory (first `cd` there in the Julia session):

```bash
julia --project launch_script.jl
```

Or from within a REPL after `include("launch_script.jl")`.

The benchmark adds negligible overhead: `@timed` is a macro that wraps the
original call with `time_ns()` and `Base.gc_num()` snapshots — no extra
allocations beyond those columns.
