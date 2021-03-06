using Statistics

"""
ref: CORRESPONDENCE ANALYSIS in PRACTICE
"""
function chisquare(x::AbstractVector, y::AbstractVector)
    x̄ = x - sum(x)
    ȳ = y - sum(y)
    n = sum(x) + sum(y)
    d = [n*(x̄[i] - ȳ[i])^2 / (x̄[i] + ȳ[i]) for i = 1:length(x)]
    sqrt(sum(d))
end

"column: feature, row: sample"
function central_chisquare(X::AbstractMatrix; dims=1)
    if dims == 1  # calculate χ² distance between column profile and central profile
        return _central_chisquare(X)
    elseif dims == 2  # calculate χ² distance between row profile and central profile
        return _central_chisquare(X')'
    end
end


function _central_chisquare(X::AbstractMatrix)
    c = sum(X, dims=2) / sum(X)
    X_ = @. (X - c)^2 / c
    sqrt.(sum(X_, dims=1))
end

""
function pairwise_chisquare(X::AbstractMatrix)
    c = sum(X, dims=1) / sum(X)
    X_ = X ./ sqrt(c)
    sqrt.(X_' * X_)
end

"column: feature, row: sample"
function chisquare_correction(X::AbstractMatrix)
    c = sum(X, dims=1) / sum(X)
    X ./ sqrt(c)
end


using StatsBase

"""
column: feature, row: sample
ref: propr: An R-package for Identifying Proportionally Abundant Features
Using Compositional Data Analysis
"""
phit(x::AbstractVector, y::AbstractVector) = var(x - y) / var(x)

perb(x::AbstractVector, y::AbstractVector) = 1 - var(x - y) / (var(x) + var(y))

phis(x::AbstractVector, y::AbstractVector) = var(x - y) / var(x + y)


"centered log-ratio transformation"
function clr(X::AbstractMatrix)
    Y = similar(X)
    for i = 1:size(X, 1)
        Y[i,:] .= log.(X[i,:] ./ StatsBase.geomean(X[i,:]))
    end
end

"additive log-ratio transformation"
function alr(X::AbstractMatrix)
    body
end
