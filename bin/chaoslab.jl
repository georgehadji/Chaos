# Command-line entry point: argument parsing lives in Chaoslab.Experiments.cli (panel, run, pilot) and Chaoslab.Analysis.cli (analyze).
using Chaoslab
try
    @eval using CUDA
catch
end
try
    @eval using CairoMakie     # the Figures extension of the analysis; present under the verify environment (D60)
catch
end
exit(!isempty(ARGS) && ARGS[1] in ("analyze", "calibrate") ? Chaoslab.Analysis.cli(ARGS) : Chaoslab.Experiments.cli(ARGS))
