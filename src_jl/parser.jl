
function parse_input(dirname::String)
    f = open(dirname)
    lines = readlines()
    @assert lines[1] == "problemData"
    chainlength = split(lines[2], ',')
    @assert chainlength[1] == "maxChainLength"
    if chainlength[2] == "Infinity"
        max_chain_length = Inf64
    else
        max_chain_length = Int64(chainlength[2])
    end
    close(f)
    return max_chain_length
end
