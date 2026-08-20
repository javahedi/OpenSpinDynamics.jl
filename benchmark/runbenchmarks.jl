using Dates
using InteractiveUtils

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RESULTS = joinpath(@__DIR__, "results")

mkpath(RESULTS)

function git_output(args...)
    try
        cmd = Cmd(["git", "-C", ROOT, args...])
        return strip(read(cmd, String))
    catch
        return "unknown"
    end
end


function write_system_info()
    filename = joinpath(RESULTS, "system_info.txt")

    open(filename, "w") do io
        println(io, "OpenSpinDynamics.jl benchmark environment")
        println(io, "=========================================")
        println(io)
        println(io, "Timestamp:      ", now())
        println(io, "Git commit:     ", git_output("rev-parse", "HEAD"))
        println(io, "Git branch:     ", git_output("rev-parse", "--abbrev-ref", "HEAD"))
        println(io, "Julia version:  ", VERSION)
        println(io, "Julia threads:  ", Threads.nthreads())
        println(io, "OS:             ", Sys.KERNEL)
        println(io, "Architecture:   ", Sys.ARCH)
        println(io, "Word size:      ", Sys.WORD_SIZE)
        println(io, "CPU threads:    ", Sys.CPU_THREADS)

        if !isempty(Sys.cpu_info())
            println(io, "CPU model:      ", Sys.cpu_info()[1].model)
        end
    end

    return filename
end

function run_benchmark(script, output)
    script_path = joinpath(@__DIR__, script)
    output_path = joinpath(RESULTS, output)

    println()
    println("Running ", script, " ...")

    cmd = `$(Base.julia_cmd()) --project=$(@__DIR__) $script_path`

    open(output_path, "w") do io
        run(pipeline(cmd, stdout=io))
    end

    println("Saved → ", output_path)

    return output_path
end

println("OpenSpinDynamics.jl benchmarks")
println("==============================")

info = write_system_info()
println("System information → ", info)

run_benchmark(
    "closed_system.jl",
    "closed_system.csv",
)

run_benchmark(
    "lindblad.jl",
    "lindblad.csv",
)

run_benchmark(
    "trajectories.jl",
    "trajectories.csv",
)

println()
println("Benchmark suite complete.")
println("Results are available in:")
println(RESULTS)
