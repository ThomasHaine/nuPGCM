# Build a cylinder with a cone at the bottom to represent an Arctic basin.
# Add WarmingSalinifyingPatch.
# Script partly built with GPT-4, with some manual adjustments.
# twnh Aug '25

using Gmsh
using GridapGmsh
using Gridap

function generate_arctic_basin3D(h, α; savefile="TestArcticBasin", show_gui=false, DEBUG=false)
    gmsh.initialize()
    model = gmsh.model
    occ = model.occ
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Geometry.OCCUnionUnify",0)     # This is important. 
    model.add("ArcticBasin")
    
    # Parameters
    R  = 1.0
    H  = -α*R
    dH = 0.2 * H       # Can't have dH equal to zero, but it can be very small.
    # dH = 0.2 * H
    r  = R * 0.5
    WSP_theta = pi/6 ;   # Angle for the WarmingSalinifying patch
    slope_gradient = (H - dH)/(R - r)
    @info "Slope gradient (must be negative): " slope_gradient

    # Compute the coordinates for the WarmingSalinifying patch
    WSP_r = R*sqrt(2)*sqrt(1 - cos(WSP_theta/2));  # Radius of the wall patch
    WSP_x = R*cos(WSP_theta/2);                    # X coordinate of the wall patch
    WSP_y = R*sin(WSP_theta/2);                    # Y coordinate of the wall patch

    # MAIN CYLINDER
    cyl1 = occ.addCylinder(0.0, 0.0,  0.0, 0.0, 0.0, dH, R)
    
    # TOP CYLINDER (cutting for patch)
    cyl2 = occ.addCylinder(0.0, 0.0, -dH/2, 0.0, 0.0,dH, r)
    vol2, _ = occ.intersect([(3, cyl1)], [(3, cyl2)], -1, false, true)
    vol3, _ = occ.fuse(vol2,[(3,1)],-1,true,true)

    # Side cylinder for the WarmingSalinifying patch
    cyl3 = occ.addCylinder(WSP_x,WSP_y,-dH/2,0,0,2*dH,WSP_r)
    vol4, _ = occ.intersect(vol3, [(3, cyl3)], -1, false, true)
    vol5, _ = occ.fuse(vol4,vol3,-1,true,true)
 
    # DEEP BASIN CONE
    cone1 = occ.addCone(0.0, 0.0, dH, 0.0, 0.0, H-dH, R, r)
    vol6, _ = occ.fuse([(3,cone1)], vol5, -1, true, true)
    occ.synchronize()

    # Mesh sizing for all points in volume
    pts = model.getEntities(0)
    for pt in pts
        model.mesh.setSize([pt], h)
    end

    # Debug: Shows geometry entities before they get associated with PhysicalGroups and PhysicalNames.
    if DEBUG
        gmsh.fltk.run()
    end

    # PHYSICAL GROUPS
    arctic_vols = [v[2] for v in vol6]
    model.addPhysicalGroup(3, arctic_vols, 1)
    model.setPhysicalName(3, 1, "ArcticBasin")
    surf_dict = Dict(
    # "Slope" => [1],
    # "CoolingPatch" => [8],
    # "WarmingSalinifyingPatch" => [3],
    # "FresheningPatch" => [6,7],
    # "Wall" => [2, 4],
    # "Abyss" => [5],
    "Surface" => [6,7,8],
    "Bottom" => [1,2,3,4,5])

    # Actually add physical groups (skip empty ones)
    for (name, tags) in surf_dict
        if !isempty(tags)
            num = Dict(
            # "Slope" => 1, 
            # "CoolingPatch" => 2, 
            # "FresheningPatch" => 3, 
            # "Wall" => 4, 
            # "Abyss" => 5, 
            # "WarmingSalinifyingPatch" => 6, 
            "Bottom" => 7, 
            "Surface" => 8)[name]
            model.addPhysicalGroup(2, tags, num)
            model.setPhysicalName(2, num, name)
        end
    end

    if show_gui
        gmsh.fltk.run()
    end

    model.mesh.generate(3)
    gmsh.write("meshes/" * savefile * ".msh")
    gmsh.finalize()

    return H, dH, R, r
end

function domain_depth( x, H, dH, R, r)          # This must be consistent with the geometry in generate_arctic_basin3D()!
    rad = sqrt(x[1]^2 + x[2]^2)
    if rad < r
        height = H
    elseif rad >= r && rad < R
        m = (dH - H)/(R - r)
        c = H - m*r
        height = m*rad + c
    elseif rad >= R
        height = 0
    end
    return -height   # Height is negative, depth is positive.
end

depth(x) = domain_depth( x, H, dH, R, r)

function test_meshArcticBasin3D()
# Usage example:
    filename = "test_ArcticBasin3D.msh"
    generate_arctic_basin3D(0.1,1.0; show_gui=false, savefile=filename, DEBUG=false)

# Test read to make sure that the Physical labels are correct. View the .msh file in Gmsh too.
    model = GmshDiscreteModel(filename)
    label_names = model.face_labeling.tag_to_name
    @info "GridapGmsh reads labels:"
    @info label_names

end
