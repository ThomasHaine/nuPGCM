# Build a 2D version of a cylinder with a truncated cone at the bottom to represent an Arctic basin.
# Based on mesh_bowl2D.jl by Henry.
# twnh Aug '25

using Gmsh: gmsh
using GridapGmsh
using Gridap

function generate_arctic_basin2D(h, α; savefile="TestArcticBasin2D.msh", show_gui=false, DEBUG=false)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("ArcticBasin2D")

    # Parameters
    R  = 1.0
    H  = -α*R
    dH = 0.2 * H       # Can't have dH equal to zero, but it can be very small.
    # r  = R * 0.4
    r  = R * 0.95   # Make r very close to R so that the slope is nearly vertical.
    slope_gradient = (H - dH)/(R - r)
    @info "Slope gradient (must be negative): " slope_gradient

    # setup domain with H = α*(1 - x^2)
    gmsh.model.geo.addPoint(-R, 0, 0, h)
    gmsh.model.geo.addPoint( R, 0, 0, h)
    gmsh.model.geo.addPoint(-R, 0,dH, h)
    gmsh.model.geo.addPoint( R, 0,dH, h)
    gmsh.model.geo.addPoint(-r, 0, H, h)
    gmsh.model.geo.addPoint( r, 0, H, h)

    
    gmsh.model.geo.addLine(1, 3)
    gmsh.model.geo.addLine(3, 5)
    gmsh.model.geo.addLine(5, 6)
    gmsh.model.geo.addLine(6, 4)
    gmsh.model.geo.addLine(4, 2)
    gmsh.model.geo.addLine(2, 1)
    gmsh.model.geo.addCurveLoop(1:6, 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.synchronize()
    # Debug: Shows geometry entities before they get associated with PhysicalGroups and PhysicalNames.
    if DEBUG
        gmsh.fltk.run()
    end

    gmsh.model.addPhysicalGroup(0, [1, 2], 1, "Bottom")          # Needed?
    gmsh.model.addPhysicalGroup(1, [1, 2, 3, 4, 5], 1, "Bottom")
    gmsh.model.addPhysicalGroup(2, [1], 3, "int")
    gmsh.model.addPhysicalGroup(1, [6], 2, "Surface")
    
    if show_gui
        gmsh.fltk.run()
    end

    # generate and save
    gmsh.model.mesh.generate(2)
    gmsh.write("meshes/" * savefile * ".msh")
    gmsh.finalize()

    return H, dH, R, r
end

function domain_depth( x, H, dH, R, r)          # This must be consistent with the geometry in generate_arctic_basin2D()!
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

function test_meshArcticBasin2D()
# Usage example:
    filename = "test_ArcticBasin2D.msh"
    generate_arctic_basin2D(0.1,1.0; show_gui=true, savefile=filename, DEBUG=false)

# Test read to make sure that the Physical labels are correct. View the .msh file in Gmsh too.
    model = GmshDiscreteModel(filename)
    label_names = model.face_labeling.tag_to_name
    @info "GridapGmsh reads labels:"
    @info label_names

end