# Condensate exhaustion

Liquid and ice going to zero is a Lambert (closed-form) root of

```
Ā·t + K·(1 − e^{−t/τ}) − c = 0
```

See [`_mm_smallest_depletion_time`](@ref), [`get_t_out_of_q_liq`](@ref), and [`get_t_out_of_q_ice`](@ref) in `src/depletion.jl`.

Do **not** iterative-root `q_l → 0` / `q_i → 0`. Saturation and freeze, for [`MM2015`](@ref) only, are residual roots; condensate out stays Lambert.
