using Gmsh

gmsh.initialize(["-noenv"]) # -noenv needed because gmsh butchers the PATH otherwise
model = gmsh.model
occ = model.occ
mesh = model.mesh

model.add("bent_waveguide")

path = abspath(joinpath(@__DIR__, "bent_waveguide.step"))
v = occ.importShapes(path)
occ.synchronize()

# Set mesh sizes
nodes_all = model.getEntities(0);
mesh.setSize(nodes_all, 10e-2 / 10);

gmsh.option.setNumber("Mesh.Algorithm3D", 10)
mesh.generate(3)

# Physical domains
model.addPhysicalGroup(2, [14], 1, "Input")
model.addPhysicalGroup(2, [1], 2, "Output")
model.addPhysicalGroup(2, collect(2:13), 3, "PEC")

model.addPhysicalGroup(3, [v[1][2]], 1, "Waveguide")

gmsh.write(joinpath(@__DIR__, "simdata" , "bent_waveguide.msh"))
