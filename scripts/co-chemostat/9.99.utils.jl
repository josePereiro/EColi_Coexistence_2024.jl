## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
function U_lines(aux_id)
    bb =  B[r"fba.glc.q_aa.phy"]

    U_LINES = Dict()
    for b in eachblob(bb)
        b["context", "aux_id"] == aux_id || continue

        # data
        glc_U = b["context", "glc_U"]
        aa_U = b["context", "aa_U"]
        opm_sol = b["opm.sol", "."]
        biom = opm_sol["BIOMASS_Ec_iML1515_WT_75p37M_f"]

        # lines
        glc_ulines = get!(U_LINES, "glc_U", Dict())
        glc_lines = get!(glc_ulines, glc_U, Dict())
        aa_ulines = get!(U_LINES, "aa_U", Dict())
        aa_lines = get!(aa_ulines, aa_U, Dict())
        
        # glc_ulines
        xs = get!(glc_lines, "xs", Float64[])
        push!(xs, aa_U)
        
        ys = get!(glc_lines, "ys", Float64[])
        push!(ys, biom)
        
        # aa_lines
        xs = get!(aa_lines, "xs", Float64[])
        push!(xs, glc_U)

        ys = get!(aa_lines, "ys", Float64[])
        push!(ys, biom)
    end
    return U_LINES
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
function U_envelop(U_LINES, nonlim_arg)
    ENV = Dict()
    for (fixxed, line) in U_LINES[nonlim_arg]
        idxs = sortperm(line["xs"])
        U1 = line["xs"][idxs]
        B1 = line["ys"][idxs]
        for (u1, b1) in zip(U1, B1)
            get!(ENV, u1, 0)
            ENV[u1] = max(maximum(b1), ENV[u1])
        end
    end

    xs = collect(Float64, keys(ENV))
    sort!(xs)
    ys = [ENV[x] for x in xs]
    return (;xs, ys)
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
function extrapolate(xt, xs::Vector, ys::Vector)
    # find closest
    # find closest
    sidx = sortperm(xs;
        by = (x) -> (x - xt)^2
    )
    i0, i1 = first(sidx, 2)
    x0, x1 = xs[i0], xs[i1]
    y0, y1 = ys[i0], ys[i1]
    
    m = (y1 - y0) / (x1 - x0)
    ye = y0 + m * (xt - x0)
    return ye
end



## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
nothing