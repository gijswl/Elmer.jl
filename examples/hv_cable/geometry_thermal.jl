using Gmsh

##
OD_cond = 38.2

t_sheath  = 1
OD_sheath = 79
ID_sheath = OD_sheath - 2 * t_sheath

OD_core = 95
dc_core = OD_core * 1.01

OD_domain = 5000

h_burial = 1200

##
function gmsh_circle(x, d, l)
    geo = gmsh.model.geo;
    r = d / 2;
    
    p1 = geo.addPoint(x[1] + r, x[2], 0, l)
    p2 = geo.addPoint(x[1] - r, x[2], 0, l)
    p3 = geo.addPoint(x[1], x[2], 0, l)
    
    l1 = geo.addCircleArc(p1, p3, p2)
    l2 = geo.addCircleArc(p2, p3, p1)
    
    cl = geo.addCurveLoop([l1, l2])
    
    return cl
end

function gmsh_core(x, params)
    cond   = gmsh_circle(x, params["OD_cond"], params["l_cond"])
    ins    = gmsh_circle(x, params["ID_sheath"], params["l_sheath"])
    sheath = gmsh_circle(x, params["OD_sheath"], params["l_sheath"])
    core   = gmsh_circle(x, params["OD_core"]  , params["l_core"])

    cond_surf   = geo.addPlaneSurface([cond]);
    ins_surf    = geo.addPlaneSurface([ins, cond[1]]);
    sheath_surf = geo.addPlaneSurface([sheath, ins]);
    core_surf   = geo.addPlaneSurface([core, sheath]);
    
    return core, [cond_surf, ins_surf, sheath_surf, core_surf];
end

##
gmsh.initialize(["-noenv"]) # -noenv needed because gmsh butchers the PATH otherwise

gmsh.model.add("geometry_thermal")
geo = gmsh.model.geo;
mesh = gmsh.model.mesh;

x0 = [0, 0];
x1 = [           0,  dc_core / sqrt(3)      ];
x2 = [ dc_core / 2, -dc_core / (2 * sqrt(3))];
x3 = [-dc_core / 2, -dc_core / (2 * sqrt(3))];

params = Dict(
    "OD_cond" => OD_cond,
    "ID_sheath" => ID_sheath,
    "OD_sheath" => OD_sheath,
    "OD_core" => OD_core,
    "l_cond" => 3,
    "l_sheath" => (OD_sheath - ID_sheath) / 1.5,
    "l_core" => (OD_core - OD_sheath) / 1.5
);

core1 = gmsh_core(x1, params)
core2 = gmsh_core(x2, params)
core3 = gmsh_core(x3, params)

lmsh_domain = h_burial / 5;
p1 = geo.addPoint(OD_domain / 2, h_burial, 0, lmsh_domain);
p2 = geo.addPoint(0, h_burial - OD_domain / 2, 0, lmsh_domain);
p3 = geo.addPoint(-OD_domain / 2, h_burial, 0, lmsh_domain);
p4 = geo.addPoint(0, h_burial, 0, 1);
l1 = geo.addCircleArc(p1, p4, p2)
l2 = geo.addCircleArc(p2, p4, p3)
l3 = geo.addLine(p3, p1)
domain_cl   = geo.addCurveLoop([l1, l2, l3])
domain_surf = geo.addPlaneSurface([domain_cl, core1[1], core2[1], core3[1]])

geo.synchronize()

geo.addPhysicalGroup(1, [l3], 1)
gmsh.model.setPhysicalName(1, 1, "Surface")

geo.addPhysicalGroup(1, [l1, l2], 2)
gmsh.model.setPhysicalName(1, 2, "Soil_infty")

geo.addPhysicalGroup(1, [1, 2, 9, 10, 17, 18], 3)
gmsh.model.setPhysicalName(1, 3, "Conductors")

geo.addPhysicalGroup(2, [core1[2][1]], 1)
gmsh.model.setPhysicalName(2, 1, "Conductor1")

geo.addPhysicalGroup(2, [core2[2][1]], 2)
gmsh.model.setPhysicalName(2, 2, "Conductor2")

geo.addPhysicalGroup(2, [core3[2][1]], 3)
gmsh.model.setPhysicalName(2, 3, "Conductor3")

geo.addPhysicalGroup(2, [core1[2][2], core2[2][2], core3[2][2]], 4)
gmsh.model.setPhysicalName(2, 4, "Insulation")

geo.addPhysicalGroup(2, [core1[2][3], core2[2][3], core3[2][3]], 5)
gmsh.model.setPhysicalName(2, 5, "Sheath")

geo.addPhysicalGroup(2, [core1[2][4], core2[2][4], core3[2][4]], 6)
gmsh.model.setPhysicalName(2, 6, "Core jacket")

geo.addPhysicalGroup(2, [domain_surf], 7)
gmsh.model.setPhysicalName(2, 7, "Air")

geo.synchronize()

bl = mesh.field.add("BoundaryLayer")
mesh.field.setNumbers(bl, "CurvesList", [1, 2, 9, 10, 17, 18])
mesh.field.setNumbers(bl, "ExcludedSurfacesList", [2, 6, 10])
mesh.field.setNumber(bl, "SizeFar", 1)
mesh.field.setNumber(bl, "Size", 0.5)
mesh.field.setNumber(bl, "Ratio", 1.1)
mesh.field.setNumber(bl, "Thickness", 15)
mesh.field.setNumber(bl, "Quads", 1)
mesh.field.setAsBoundaryLayer(bl)

mesh.generate(2)

gmsh.fltk.run()
gmsh.write(joinpath(@__DIR__, "simdata", "thermal.msh"))

gmsh.finalize()