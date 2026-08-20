using Documenter
using OpenSpinDynamics

makedocs(
    sitename="OpenSpinDynamics.jl",
    modules=[OpenSpinDynamics],

    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
    ),

    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",

        "Examples" => [
            "Closed XXZ dynamics" =>
                "examples/closed_xxz.md",

            "Amplitude damping" =>
                "examples/amplitude_damping.md",

            "Dissipative XXZ chain" =>
                "examples/dissipative_xxz.md",

            "Trajectories vs Lindblad" =>
                "examples/trajectories_vs_lindblad.md",
        ],

        "API" => "api.md",
    ],
)

deploydocs(
    repo="github.com/javahedi/OpenSpinDynamics.jl.git",
    devbranch="main",
)