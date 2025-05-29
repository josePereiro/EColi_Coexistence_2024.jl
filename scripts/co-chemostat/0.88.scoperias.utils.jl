
## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: init_block
macro init_block(lb, ver)
    quote
        empty!(B)
        empty!(G)
        empty!(C)    
        Scoperias.@sc_label $(lb)
        blockver = $(ver)
        blockttag = Dates.format(now(), "yyyymmdd-HHMMSSS")
        blockuuid =  uuid4()
        bb = blobbatch!(B, $(lb))

        @sc_call "init_block"
        
    end |> esc
end


## -- .. - .-- .-. . .... -- -- -- .. ...
# MARK: blockbatch!
function blockbatch!(sc::Scope)
    bbid = string(
        sc["blockttag"][], "-",
        sc["blockuuid"][]
    )
    blobbatch!(B, bbid)
end

macro blockbatch!()
    quote
        blockbatch!(@sc_scope())
    end |> esc
end

macro blockbatch!()
    quote
        blockbatch!(@sc_scope())
    end |> esc
end

## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
ascontext(v::Number) = v
ascontext(v::Bool) = v
ascontext(v::AbstractString) = v
ascontext(v::DateTime) = v
ascontext(v::Date) = v
ascontext(v::BlobyRef) = v

ascontext(v::Any) = hash(v)

function litecontext(sc::Scope)
    new_sc = Dict{String, Any}()
    for scv in values(sc)
        new_sc[scv.key] = ascontext(scv.val)
    end
    return new_sc
end


## --- .- -.. -- .-. . . .-- - -. . . .- .- .-.- .. .
nothing