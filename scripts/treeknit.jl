using Pkg
Pkg.activate(".")
using ArgParse
using TreeKnit
using TreeTools

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--treel"
            help = "L raw nwk input file"
        "--treerest"
            help = "rest raw nwk input file"
        "--outputrest"
            help = "tree output rest resolved nwk"
        "--outputl"
            help = "tree output l resolved nwk"
    end
    return parse_args(s)
end

function main()
    @show parsed_args = parse_commandline()

    l_tree = parsed_args["treel"]
    rest_tree = parsed_args["treerest"]

    t_without = read_tree(rest_tree, label="rest")
    t_L = read_tree(l_tree, label="L")
    Knit = run_treeknit!(t_L, t_without, OptArgs(;pre_resolve=true, resolve=true, strict=true))

    write_newick(parsed_args["outputrest"], t_without) 
    write_newick(parsed_args["outputl"], t_L)
end


main()