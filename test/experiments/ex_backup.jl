# Amendment D59: retry and best effort for the manifest backups. Spec: experiments-spec.md Amendment D59.

bk_busy(code) = Base.IOError("open: held by another process", code)

# a copier that fails `n` times with `code` before it copies (the Base.cp with force)
function bk_flaky(n, code)
    calls = Ref(0)
    return (a, b; force = false) -> (calls[] += 1; calls[] <= n ? throw(bk_busy(code)) : cp(a, b; force = force)), calls
end

# a path with the versions v1 .. vN written by the real writer (keep = 3), so that backups exist
function bk_versions(tag, n)
    path = ex_path("bk-$tag.json")
    foreach(f -> rm(f; force = true), [path; ["$path.$k" for k in 1:3]; path * ".partial"])
    for v in 1:n
        EX._write_json(path, Dict("v" => v); keep = 3)
    end
    return path
end

bk_v(f) = EX.read_json(f)["v"]

@testset "D59 retry of backup copies and of the temporary write" begin
    EX.RETRY_PAUSE[] = 0.0
    try
        for code in (Base.UV_EBUSY, Base.UV_EACCES, Base.UV_EPERM)
            path = bk_versions("flaky-$code", 3)
            copier, calls = bk_flaky(2, code)
            EX.BACKUP_COPIER[] = copier
            try
                EX._write_json(path, Dict("v" => 4); keep = 3)
            finally
                EX.BACKUP_COPIER[] = cp
            end
            @test calls[] > 2
            @test (bk_v(path), bk_v(path * ".1"), bk_v(path * ".2"), bk_v(path * ".3")) == (4, 3, 2, 1)
        end
        # the temporary write fails twice, then works
        path = ex_path("bk-writer.json")
        n = Ref(0)
        writer = (f, s) -> (n[] += 1; n[] <= 2 ? throw(bk_busy(Base.UV_EACCES)) : write(f, s))
        @test EX._write_json(path, Dict("v" => 1); writer = writer) == path && n[] == 3 && bk_v(path) == 1
        # a writer that never works ends the write with the error, and the previous file stays whole
        @test_throws Base.IOError EX._write_json(path, Dict("v" => 2); writer = (f, s) -> throw(bk_busy(Base.UV_EACCES)))
        @test bk_v(path) == 1
    finally
        EX.RETRY_PAUSE[] = 0.1
        EX.BACKUP_COPIER[] = cp
    end
end

@testset "D59 a backup that cannot be copied is skipped" begin
    EX.RETRY_PAUSE[] = 0.0
    try
        path = bk_versions("skip", 3)          # v3 current, .1 = v2, .2 = v1
        EX.BACKUP_COPIER[] = (a, b; force = false) -> throw(bk_busy(Base.UV_EACCES))
        failures = String[]
        @test EX._write_json(path, Dict("v" => 4); keep = 3, failures = failures) == path
        @test bk_v(path) == 4
        @test length(failures) == 3 && all(f -> occursin("bk-skip", f), failures)
        @test (bk_v(path * ".1"), bk_v(path * ".2")) == (2, 1)             # the older backups stay
        # without a place for the failures the write still goes on
        @test EX._write_json(path, Dict("v" => 5); keep = 3) == path && bk_v(path) == 5
        # another kind of error in a backup copy is raised
        EX.BACKUP_COPIER[] = (a, b; force = false) -> throw(bk_busy(Base.UV_ENOSPC))
        @test_throws Base.IOError EX._write_json(path, Dict("v" => 6); keep = 3)
    finally
        EX.RETRY_PAUSE[] = 0.1
        EX.BACKUP_COPIER[] = cp
    end
end

@testset "D59 the manifest of a campaign lists skipped backups" begin
    fx = ex_fix(:small)
    EX.RETRY_PAUSE[] = 0.0
    try
        EX.BACKUP_COPIER[] = (a, b; force = false) -> throw(bk_busy(Base.UV_EACCES))
        res = ex_campaign(fx; store = ex_store("bk-fail"))
        @test res.status == "complete"
        man = EX.read_json(res.manifest)
        @test haskey(man, "backup_failures") && !isempty(man["backup_failures"])
    finally
        EX.RETRY_PAUSE[] = 0.1
        EX.BACKUP_COPIER[] = cp
    end
    res0 = ex_campaign(fx; store = ex_store("bk-ok"))
    @test res0.status == "complete" && !haskey(EX.read_json(res0.manifest), "backup_failures")
end
