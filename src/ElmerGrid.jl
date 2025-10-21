
"""
    function elmergrid_gmsh(mesh_path, mesh_file; elmergrid_path, log_file, flags)

Runs ElmerGrid (on the `elmergrid_path`) on the gmsh mesh stored in `mesh_path/mesh_file`.
Additional flags (by default `-autoclean`) are provided via the optional keyword `flags`.
The `stdout` logging is piped to the `log_file`.
"""
function elmergrid_gmsh(mesh_path::String, mesh_file::String; elmergrid_path="ElmerGrid", log_file="elmergrid.log", flags="-autoclean")
    logfile = joinpath(mesh_path, log_file)
    cmd = Cmd(`$elmergrid_path 14 2 $mesh_file $flags`, dir=mesh_path)
    pipe = pipeline(cmd, stdout=logfile)
    run(pipe)
end

"""
    function partition_mesh(mesh_path, mesh_dir, partitions; elmergrid_path, log_file, flags)

Partitions the mesh in directory `mesh_path/mesh_dir` using ElmerGrids' `-metiskway` flag.
Additional flags are provided via the optional keyword `flags`.
The `stdout` logging is piped to the `log_file`.
"""
function partition_mesh(mesh_path::String, mesh_dir::String, partitions::Int; elmergrid_path="ElmerGrid", log_file="elmergrid.log", flags="")
    logfile = joinpath(mesh_path, log_file)
    cmd = Cmd(`$elmergrid_path 2 2 $mesh_dir -partdual -metiskway $partitions $flags`, dir=mesh_path)
    pipe = pipeline(cmd, stdout=logfile)
    run(pipe)
end