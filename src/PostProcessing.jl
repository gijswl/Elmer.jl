using CSV
using DataFrames

function load_dat(filename)
    header_file = filename * ".names"

    data = CSV.File(filename; delim=' ', ignorerepeated=true, header=false)

    header = open(header_file, "r") do file
        parse_names(file)
    end

    return DataFrame(data, header)
end

function parse_names(file::IOStream)
    header = String[]
    names_start = false

    for line ∈ eachline(file)
        if (names_start)
            name = join(split(line, ": ")[2:end], ": ") # Remove column number from name
            push!(header, name)
        end

        if (startswith(line, "Variables in columns of matrix"))
            names_start = true
        end
    end

    return header
end