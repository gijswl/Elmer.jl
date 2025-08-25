# function run_elmer_grid(mesh_path::String, mesh_file::String; elmergrid_path="ElmerGrid", log_file="elmergrid.log")
#     logfile = joinpath(mesh_path, log_file)
#     cmd = Cmd(`$elmergrid_path 14 2 $mesh_file`, dir=mesh_path)
#     pipe = pipeline(cmd, stdout=logfile)
#     run(pipe)
# end

function elmergrid_gmsh(mesh_path::String, mesh_file::String; elmergrid_path="ElmerGrid", log_file="elmergrid.log", flags="-autoclean")
    logfile = joinpath(mesh_path, log_file)
    cmd = Cmd(`$elmergrid_path 14 2 $mesh_file`, dir=mesh_path)
    pipe = pipeline(cmd, stdout=logfile)
    run(pipe)
end

function partition_mesh(mesh_path::String, mesh_dir::String, partitions::Int; elmergrid_path="ElmerGrid", log_file="elmergrid.log", flags="")
    logfile = joinpath(mesh_path, log_file)
    cmd = Cmd(`$elmergrid_path 2 2 $mesh_dir -partdual -metiskway $partitions $flags`, dir=mesh_path)
    pipe = pipeline(cmd, stdout=logfile)
    run(pipe)
end