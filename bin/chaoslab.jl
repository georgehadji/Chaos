# Command-line entry point: argument parsing lives in Chaoslab.Experiments.cli (panel, run, pilot).
using Chaoslab
try
    @eval using CUDA
catch
end
exit(Chaoslab.Experiments.cli(ARGS))
