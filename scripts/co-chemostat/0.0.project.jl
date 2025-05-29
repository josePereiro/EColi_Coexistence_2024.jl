using EColi_Coexistence_2024
using ProjFlows
using MetX
using MetX: MetXGEMs, MetXEP
using Bloberias

# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.
#  Project
# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.

PROJ = Project0(EColi_Coexistence_2024)
globproj!(PROJ)

# TODO: implement version into project0
# TODO: implement a `versioninfo`-like method for Projects to summarize 
PROJVER = v"0.1.0"

# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.
# proc
_procdir(args...) = projpath(PROJ, ["data", string(PROJVER)], args...)

# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.
# db
_blobsdir(args...) = _procdir(["blob"], args...)

# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.
# Bloberia
B = Bloberia(_blobsdir())
G = blobbatch!(B, "_sim.globals") 
C = blobbatch!(B, "_sim.cache") 
mkpath(B)

# .-- .- -.-.-.--. ...---. . . . -- .--. -. -. -.
nothing