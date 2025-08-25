
function strip_final_slash(s::String)
    if (s[end] == '/')
        return chop(s)
    end
    return s
end