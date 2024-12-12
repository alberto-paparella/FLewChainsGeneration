# Finite FLew-Chains Generation

## Getting started

Launch a Julia command-line interface in your terminal:

```
$ julia
```

If you want to make use of multi-thread execution, add the `-t` option specifying the number `n` of cores to be used. Here's and example for launching the session using 8 cores:

```
$ julia -t 8
```

Then include `core.jl` file

```
julia> include("core.jl")
```

To generate all possible finite FLew-chains with `n` values, use the `generateflewchains` function. Putting a semicolon `;` after the function will prevent it to output all generated finite FLew-chains. One could also make use of the `@time` macro to take note of the execution time as well as other useful information. Here's an example generating all possible finite FLew-chains with 9 elements and storing them in the variable `q` printing the execution time:

```
julia> @time q = generateflewchains(9);
```

To know how many finite FLew-chains have been generated, simply use the `length` function (`generateflewchains` returns a `Vector` of `FiniteFLewChain`s):

```
julia> length(q)
```

One clean way to systematically print all generate finite FLew-chains is the following:

```
julia> for i in q
           println(i)
           sleep(0.1)
           Base.run(`clear`)
       end
```

Note that if you launched the julia command-line interface with the `-t` option, the generated finite FLew-chains are not stored in order.

## Results

| Values | FLew-chains |
|--------|-------------|
| 1      | 1           |
| 2      | 1           |
| 3      | 2           |
| 4      | 6           |
| 5      | 22          |
| 6      | 94          |
| 7      | 451         |
| 8      | 2386        |
| 9      | 12775       |
| 10     | 86417       |

