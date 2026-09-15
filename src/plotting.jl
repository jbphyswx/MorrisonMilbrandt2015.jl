"""Plot mean liquid/ice tendencies across timesteps. Requires `using CairoMakie`."""
function plot_rates end

"""Plot one T-updating trajectory and its frozen-coefficient comparison. Requires `using CairoMakie`."""
function plot_evolution end

"""Plot liquid and ice moisture trajectories. Requires `using CairoMakie`."""
function plot_condensate end

"""Plot liquid- and ice-frame supersaturation trajectories. Requires `using CairoMakie`."""
function plot_supersaturation end

"""Plot vapor and phase-saturation moisture coordinates. Requires `using CairoMakie`."""
function plot_specific_humidities end

"""Plot parcel temperature through one solver step. Requires `using CairoMakie`."""
function plot_temperature end
