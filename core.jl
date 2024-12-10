"""
Given a set A of N linearly ordered elements, A = {0, ..., N-1}, and
C_N = {0, 1/(N-1), 2/(N-1), ..., (N-2)/(N-1), (N-1)/(N-1) = 1}.

We could represent the magma (C_N) with an N × N matrix, using N^2 bytes.
 - Since the magma is abelian, we just need to represent its upper triangular
   using (N*(N+1))÷2 bytes.
 - Since 1 is an identity element and 0 an absorbing element, we don't need
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
struct FLewMonoid{N}
    o::Vector{Int}
    
    function FLewMonoid{N}(o::Vector{Int}) where {N}
        length(o) != ((N-2)*(N-1))÷2 && error("Wrong number of elements provided")
        new{N}(o)
    end

    function FLewMonoid{N}() where {N}
        N < 1 && error("N must be greater than 0")
        (N == 1 || N == 2) && return FLewMonoid{N}(Vector{Int}())
        N > 2 && error("Please provide elements for N > 2")
    end
end

import Base.show

function Base.show(io::IO, fm::FLewMonoid{N}) where {N}
    ct = Matrix{Union{Nothing, Int}}(nothing, N, N)
    for i in 1:N
        ct[1,i] = ct[i,1] = 0
        ct[N,i] = ct[i,N] = i-1
    end
    for i in 2:N-1
        for j in i:N-1
            # ct[i,j] = ct[j,i] = fm.o[(N-2)*(i-2)+j-1-((i-1)*(i-2))÷2]
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

function evaluate(fm::FLewMonoid{N}, a::Int, b::Int) where {N}
    a == 0 && return 0
    b == 0 && return 0
    a == N-1 && return b
    b == N-1 && return a
    a > N-1 || a < 0 && error("a is not in the domain")
    b > N-1 || b < 0 && error("b is not in the domain")
    a < b ? fm.o[(N-2)*(a-1)+b-(a*(a-1))÷2] : fm.o[(N-2)*(b-1)+a-(b*(b-1))÷2]
end

function weaklyincreasing(q, seq::Vector{Int}, min::Int, n::Int, l::Int)
    if l == 1
        for i in min:n-1
            push!(seq, i)
            push!(q, copy(seq))
            pop!(seq)
        end
        return q
    else
        for i in min:n-1
            push!(seq, i)
            weaklyincreasing(q, seq, i, n, l-1)
            pop!(seq)
        end
        return q
    end
end


weaklyincreasing(min::Int, n::Int,l::Int) = weaklyincreasing([], Vector{Int}(), min, n, l)
weaklyincreasing(n::Int,l::Int) = weaklyincreasing([], Vector{Int}(), 0, n, l)

function generateflewchain(q::Vector{Vector{Int}}, o::Vector{Int}, min::Int, n::Int, l::Int)
    if l == 1
        wi = weaklyincreasing(min, n, l)
        for i in wi
            push!(q, vcat(o, i))
        end
        return q    
    else
        wi = weaklyincreasing(min, n, l)
        for i in wi
            iswi = true
            if l != n-2
                for j in 1:length(i)
                    if reverse(i)[j] < reverse(o)[j]
                        iswi = false
                        break
                    end
                end
            end
            iswi && generateflewchain(q, vcat(o, i), i[2], n, l-1)
        end
        return q
    end
end

function generateflewchain(n)
    n < 3 && return Vector{Vector{Int}}()
    generateflewchain(Vector{Vector{Int}}(), Vector{Int}(), 0, n, n-2)
end

function evaluate(o::Vector{Int}, a::Int, b::Int, n::Int)
    a == 0 && return 0
    b == 0 && return 0
    a == n-1 && return b
    b == n-1 && return a
    a < b ? o[(n-2)*(a-1)+b-(a*(a-1))÷2] : o[(n-2)*(b-1)+a-(b*(b-1))÷2]
end

function checkassociativity(o::Vector{Int}, n::Int)
    for α in 1:n-2
        Δ = Vector{Int}()
        ∇ = Vector{Int}()
        for x in 1:n-2
            for y in 1:n-2
                push!(Δ, evaluate(o, evaluate(o, x, α, n), y, n))
                push!(∇, evaluate(o, x, evaluate(o, α, y, n), n))
            end
        end
        Δ != ∇ && return false
    end
    return true
end
