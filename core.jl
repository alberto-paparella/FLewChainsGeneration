"""
    struct FiniteFLewChain{N}
        cayleytable::Vector{Int8}
    end

A structure representing a finite FLew-chain over a set of N linearly ordered elements.

Given a set A of N linearly ordered elements A = {0, ..., N-1}, 
C_N = {0, 1/(N-1), 2/(N-1), ..., (N-2)/(N-1), (N-1)/(N-1) = 1}.

We could represent the Cayley's table for (C_N, ⋅) with an N × N matrix, but:
 - since the magma is abelian, we just need to represent its upper triangular
   using an (N*(N+1))÷2 vector;
 - since 1 is an identity element and 0 an absorbing element, we don't need
   to represent the first and last rows/columns, so we'd just need a (N-2)^2 matrix.
Taking both points into consideration, we just need an ((N-2)(N-1))÷2 vector.

E.g., with N=5:

   ⋅  |  0  1/4 2/4 3/4  1
 ===========================
   0  |  0   0   0   0   0
  1/4 |  0   A   B   C  1/4
  2/4 |  0   B   D   E  2/4
  3/4 |  0   C   E   F  3/4
   1  |  0  1/4 2/4 3/4  1

We just need {A, B, C, D, E, F} (6 elements instead of 25!!)
"""
struct FiniteFLewChain{N}
    cayleytable::Vector{Int8}
    
    function FiniteFLewChain{N}(cayleytable::Vector{Int8}) where {N}
        if length(cayleytable) != ((N-2)*(N-1))÷2
            error("Wrong number of elements provided")
        end
        return new{N}(cayleytable)
    end

    function FiniteFLewChain{N}() where {N}
        if N < 1
            error("N must be greater than 0")
        elseif N == 1 || N == 2
            return FiniteFLewChain{N}(Vector{Int8}())
        else
            return N > 2 && error("Please provide elements for N > 2")
        end
    end
end

import Base.show

function Base.show(io::IO, fm::FiniteFLewChain{N}) where {N}
    ct = Matrix{Union{Nothing, Int}}(nothing, N, N)
    for i in 1:N
        ct[1,i] = ct[i,1] = 0
        ct[N,i] = ct[i,N] = i-1
    end
    for i in 2:N-1
        for j in i:N-1
            ct[i,j] = ct[j,i] = evaluate(fm, i-1, j-1)
        end
    end
    print(io, "\n")
    print(io, " ⋅ |")
    for i in 1:N print(io, " $(Int(i-1))") end
    print(io, "\n")
    for i in 1:N+3 print(io, "==") end
    print(io, "\n")
    for i in 1:N
        print(io, " $(Int(i-1)) |")
        for j in 1:N print(io, " $(Int(ct[i,j]))") end
        print(io, "\n")
    end
end

function evaluate(fm::FiniteFLewChain{N}, a::Int, b::Int) where {N}
    a == 0 && return 0
    b == 0 && return 0
    a == N-1 && return b
    b == N-1 && return a
    (a > N-1 || a < 0) && error("a is not in the domain")
    (b > N-1 || b < 0) && error("b is not in the domain")
    a < b ? fm.cayleytable[(N-2)*(a-1)+b-(a*(a-1))÷2] : fm.cayleytable[(N-2)*(b-1)+a-(b*(b-1))÷2]
end

function weaklyincreasing(q, seq::Vector{Int8}, min::Int, max::Int, l::Int)
    if l == 0
        push!(q, copy(seq))
    else
        for i in min:max
            push!(seq, i)
            weaklyincreasing(q, seq, i, max, l-1)
            pop!(seq)
        end
        return q
    end
end


weaklyincreasing(min::Int, max::Int, l::Int) = weaklyincreasing(Vector{Vector{Int8}}(), Vector{Int8}(), min, max, l)

function generateflewchain(q::Vector{Vector{Int8}}, o::Vector{Int8}, min::Int, max::Int, l::Int, n::Int, u)
    wi = weaklyincreasing(min, max, l) 
    if l == 1
        for i in wi
            if checkassociativity(vcat(o, i), n)
                Threads.lock(u) do
                    push!(q, vcat(o, i))
                end
            end
        end
    else
        wi = weaklyincreasing(min, max, l)
        Threads.@threads for i in wi
            iswi = true
            if !isempty(o)
                for j in 1:length(i)
                    if reverse(i)[j] < reverse(o)[j]
                        iswi = false
                        break
                    end
                end
            end
            iswi && generateflewchain(q, vcat(o, i), Int64(i[2]), max+1, l-1, n, u)
        end
    end
end

function generateflewchain(n)
    q = Vector{Vector{Int8}}()
    n < 3 && return q
    generateflewchain(q, Vector{Int8}(), 0, 1, n-2, n, Threads.SpinLock())
    return q
end

function evaluate(o::Vector{Int8}, a::Int, b::Int, n::Int)
    a == 0 && return 0
    b == 0 && return 0
    a == n-1 && return b
    b == n-1 && return a
    a < b ? o[(n-2)*(a-1)+b-(a*(a-1))÷2] : o[(n-2)*(b-1)+a-(b*(b-1))÷2]
end

function checkassociativity(o::Vector{Int8}, n::Int)
    for a in 1:n-2
        for b in 1:n-2
            for c in 1:n-2
                evaluate(o, Int64(evaluate(o, a, b, n)), c, n) != evaluate(o, a, Int64(evaluate(o, b, c, n)), n) && return false
            end
        end
    end
    return true
end

function generateflewchains(n::Int)
    if n < 3
        [FiniteFLewChain{n}()]
    else
        # FiniteFLewChain{n}.([filter(x -> checkassociativity(x, n)==1, generateflewchain(n))...])
        FiniteFLewChain{n}.([generateflewchain(n)...])
    end
end
