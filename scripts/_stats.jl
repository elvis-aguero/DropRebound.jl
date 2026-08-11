# The four summary statistics the analysis scripts need, without a dependency.
#
# `Statistics` stopped being a bundled standard library in Julia 1.12 -- it is a downloadable
# package now, and declaring it makes `Pkg.build` fail in a clean environment before a single
# test runs. Since this repository already hand-rolls its CSV reading and its Legendre
# polynomials, four one-line reductions are more in keeping than a dependency that has to be
# fetched. Definitions match Statistics: `median` averages the middle pair for even counts,
# `std` is the sample standard deviation with the n-1 denominator, and `quantile` interpolates
# linearly between order statistics.

mean(v) = sum(v) / length(v)

function median(v)
    s = sort(collect(v)); n = length(s)
    isodd(n) ? s[(n+1) ÷ 2] : 0.5*(s[n ÷ 2] + s[n ÷ 2 + 1])
end

function std(v)
    n = length(v); n < 2 && return 0.0
    m = mean(v)
    sqrt(sum((x - m)^2 for x in v) / (n - 1))
end

function quantile(v, p)
    s = sort(collect(v)); n = length(s)
    n == 1 && return float(s[1])
    h = (n - 1) * p + 1
    lo = clamp(floor(Int, h), 1, n); hi = clamp(lo + 1, 1, n)
    s[lo] + (h - lo) * (s[hi] - s[lo])
end
