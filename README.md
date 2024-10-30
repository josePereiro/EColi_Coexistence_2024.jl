# EColi_Coexistence_2024

[![CI](https://github.com/josePereiro/EColi_Coexistence_2024.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/josePereiro/EColi_Coexistence_2024.jl/actions/workflows/CI.yml)
<!-- TODO: Make CODECOV work -->
<!-- [![Coverage](https://codecov.io/gh/josePereiro/EColi_Coexistence_2024.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/josePereiro/EColi_Coexistence_2024.jl) -->

# Dev install

```bash
julia -e 'import Pkg; Pkg.Registry.add(url = "https://github.com/MetabolicXploration/MetX_Registry_jl"); Pkg.Registry.add(url = "https://github.com/FF-UH/CSC_Registry.jl"); Pkg.Registry.add(url = "https://github.com/JuliaRegistries/General.gitjl"); Pkg.develop(url="https://github.com/josePereiro/EColi_Coexistence_2024.jl");Pkg.develop("EColi_Coexistence_2024"); Pkg.instantiate()'
```