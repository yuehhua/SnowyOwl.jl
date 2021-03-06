@testset "object" begin
    r, c = (100, 500)
    data = rand([0, 1], r, c)
    obs = DataFrame(A=rand(c), B=rand(c))
    var = DataFrame(C=rand(r), D=rand(r))
    @test_throws AssertionError Profile(data, obs, var)

    prof = Profile(data, var, obs)
    @test obsnames(prof) == ["A", "B"]
    @test varnames(prof) == ["C", "D"]
    @test nrow(prof) == r
    @test ncol(prof) == c
    @test nvar(prof) == r
    @test nobs(prof) == c
    @test maximum(prof) == 1
    @test minimum(prof) == 0
    @test size(prof) == (r, c)
    @test axes(prof) == (Base.OneTo(r), Base.OneTo(c))
    @test_throws AssertionError prof.obs = var
    @test_throws AssertionError prof.var = obs

    prof2 = prof[1:50, :]
    @test prof2.data == prof.data[:, 1:50]

    idx = rand([false,true], c)
    prof3 = prof[idx, :]
    @test prof3.data == prof.data[:, idx]
    @test prof3.var == prof.var
    @test prof3.obs == prof.obs[idx, :]
end
