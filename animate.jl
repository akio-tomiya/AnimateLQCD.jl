using Gaugefields
using Wilsonloop
using LinearAlgebra
using Plots
using Interpolations
using LaTeXStrings
using Printf

# Constants
const LtoSec = 10/3 # to change unit from fm to sec, we need this.
const myxlabel = L"$x$ [fm]"
const myylabel = L"$y$ [fm]"
const r_0 = 0.48 # [fm] . Sommer scale.
const plaquette_min = 0.5
const plaquette_max = 1.0

# Calculates ln_a(beta) based on the formula provided in:
# https://arxiv.org/abs/hep-lat/9806005
function ln_a(beta::Float64)::Float64
    if beta < 5.7 || beta > 6.57
        throw(ArgumentError("Beta should be in the range [5.7, 6.57]"))
    end

    delta_beta = beta - 6
    return -1.6805 - 1.7139*delta_beta + 0.8155*delta_beta^2 - 0.6667*delta_beta^3
end

# Subroutine to calculate 'a' from beta
function calculate_a(beta::Float64)::Float64
    return r_0 * exp(ln_a(beta))
end

function create_animation(U1, scale_factor)
    NX, NY, NZ, NT, NC = size(U1[1])[3], size(U1[1])[4], size(U1[1])[5], size(U1[1])[6], size(U1[1])[1]
    min_value, max_value = 0.0, 1.0
    loop = [(1, +1), (2, +1), (1, -1), (2, -1)]
    
    # Prepare Wilson Line
    w = Wilsonline(loop)
    Uloop = similar(U1[1])
    temps = [similar(U1[1]) for _ in 1:10]
    Gaugefields.evaluate_gaugelinks!(Uloop, w, U1, temps)

    anim = @animate for t in 1:NT
        plaqs = zeros(Float64, NX, NY)
        
        for y in 1:NY, x in 1:NX, z in 1:NZ
            plaqs[x, y] += real(tr(Uloop[:, :, x, y, z, t]))
        end
        plaqs .*= 1 / (NZ * NC)

        # Smoothing
        itp = interpolate(plaqs', BSpline(Cubic(Flat(OnGrid()))))
        x_fine = 1:0.05:NX
        y_fine = 1:0.05:NY
        x_physical = x_fine .* scale_factor
        y_physical = y_fine .* scale_factor
        plaqs_fine = [itp(x, y) for x in x_fine, y in y_fine]
        t_showing_val = round(t * scale_factor * LtoSec, digits=2)
        my_title = L"Energy density at $t=" * @sprintf("%.2f", t_showing_val) * L"\times 10^{-24}~\mathrm{s}$"
        surface(x_physical, y_physical, plaqs_fine, title=my_title, zlim=(min_value, max_value), xlabel=myxlabel, ylabel=myylabel, clims=(plaquette_min, plaquette_max))
    end

    gif(anim, "latticeQCD_animation.gif", fps=10)
end

function main()
    NX, NY, NZ, NT, NC, Nwing = 16, 16, 4, 32, 3, 1
    beta = 6.1
    filename = "confs_Heatbath_L16160432_beta6.1_quenched/conf_00000100.ildg"

    a = calculate_a(beta)
    scale_factor = a

    U1 = Initialize_Gaugefields(NC, Nwing, NX, NY, NZ, NT, condition="cold")
    ildg = ILDG(filename)
    load_gaugefield!(U1, 1, ildg, [NX, NY, NZ, NT], NC)

    create_animation(U1, scale_factor)
end

main()
