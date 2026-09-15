# MM2015FixedT — frozen T

[`MM2015FixedT`](@ref) is the event-driven Appendix C solver at a **frozen** thermodynamic state.

1. Derive active phases `(la, ia)` from condensate and supersaturation.
2. Invert **C5** for the time to liquid or ice saturation (`Δ = q_sl − q_si` is constant at frozen T, so that invert **is** the saturation event).
3. Lambert for condensate exhaustion.
4. Apply C6 over the earliest event (or the remaining `Δt`).

Entry point: `morrison_milbrandt_2015(::MM2015FixedT, …)` or the same name on [`MM2015Inputs`](@ref).
The caller passes unpacked scalars. *T* is held fixed for the whole call; [`MM2015`](@ref)
updates *T* through the residual path.
