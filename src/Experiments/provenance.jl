# Provenance of a panel file: git state of the package and its version.

"""
    git_state(dir) -> (sha, dirty)

Head commit of the repository containing `dir` and whether a tracked file is modified (`git status --porcelain
--untracked-files=no` is not empty). `("unknown", true)` when git cannot be run.
"""
function git_state(dir::AbstractString)
    try
        sha = strip(read(`git -C $dir rev-parse HEAD`, String))
        dirty = !isempty(strip(read(`git -C $dir status --porcelain --untracked-files=no`, String)))
        return (sha = String(sha), dirty = dirty)
    catch
        return (sha = "unknown", dirty = true)
    end
end

"Version of the package as written in its `Project.toml`."
package_version() = VersionNumber(TOML.parsefile(joinpath(_root(), "Project.toml"))["version"])

# Directory of the package (the repository root when developed in place).
_root() = pkgdir(parentmodule(@__MODULE__))::String
