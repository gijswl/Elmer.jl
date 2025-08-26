using YAML

function load_solver!(sif::SolverInformationFile, name::String, exec::ExecSolver, database::String)
    data = YAML.load_file(database, dicttype=OrderedDict{String,Any})
    @assert haskey(data, name) "No solver $name exists in database $database"

    solver_data = data[name]
    return add_solver!(sif, name, exec; data = solver_data)
end

function load_material!(sif::SolverInformationFile, name::String, database::String)
    data = YAML.load_file(database, dicttype=OrderedDict{String,Any})
    @assert haskey(data, name) "No material $name exists in database $database"

    material_data = data[name]
    return add_material!(sif, name; data = material_data)
end