using Gmsh: gmsh
using Printf
using GridapGmsh
using Gridap

function generate_exp_bowl_mesh_2D(h, depth; savefile="Test_exp_bowl2D.msh", DEBUG=false, show_gui=false )
    # This function generates a 2D mesh for a bowl-shaped domain with depth profile given by the function `depth`. 
    # It's a template function: you provide the depth function as an argument, although it must take care of the edges at x = ±1.
    
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("exp_bowl2D")

    # Range of r and corresponding depth values for mesh refinement.
    # If depth goes to zero at the edges, like for Henry's parabolic bowl, then the mesh will be very fine there, and Gmsh may error.
    # So we avoid the endpoints -1 and 1 by a small amount if the depth goes to zero there.
        # xs = range(-1.0+h/2.0, 1.0-h/2.0, length=round(Int,2.0/h))     # radial points for mesh refinement
    xs = range(-1.0, 1.0, length=round(Int,2.0/h))                       # radial points for mesh refinement
    pts = [[x, 0.0] for x in xs]                                         # E.g. along x-axis
    Hs = depth.(pts)                                                     # corresponding depths

    # Setup domain with points at (-1,0,0) and (1,0,0) on the surface.
    point1 = gmsh.model.geo.addPoint(-1, 0, 0, h)
    point2 = gmsh.model.geo.addPoint( 1, 0, 0, h)
    
    # Add points along the depth curve
    point_tags = Int[]

    for (x, H) in zip(xs, Hs)
        push!(point_tags, gmsh.model.geo.add_point(x, 0, -H, h))
    end
    if(DEBUG)
        println("Point tags along depth curve: ", point_tags)
        println("Corresponding radial positions: ", xs)
        println("Corresponding depths: ", Hs)
    end

    # Create lines to connect the points, and a spline
    curve_tags = Int[]

    # 1. Line from point1 to first depth curve point
    push!(curve_tags, gmsh.model.geo.add_line(point1, point_tags[1]))

    # 2. Spline along depth curve points
    push!(curve_tags, gmsh.model.geo.add_spline(point_tags))

    # 3. Line from last depth curve point to point2
    push!(curve_tags, gmsh.model.geo.add_line(point_tags[end], point2))

    # 4. Line from point2 back to point1 (to close the loop)
    push!(curve_tags, gmsh.model.geo.add_line(point2, point1))

    if DEBUG
        gmsh.model.geo.synchronize()
        gmsh.fltk.run()
    end

    # Now add the loop using those curves
    gmsh.model.geo.addCurveLoop(curve_tags, 1)

    # Debug: Shows geometry entities before they get associated with PhysicalGroups and PhysicalNames.
    if DEBUG
        gmsh.model.geo.synchronize()
        gmsh.fltk.run()
    end

    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(1, [1, 2, 3], 1, "Bottom")
    gmsh.model.addPhysicalGroup(1, [4], 2, "Surface")
    gmsh.model.addPhysicalGroup(2, [1], 4, "int")

    if show_gui
        gmsh.fltk.run()
    end

    # generate and save
    gmsh.model.mesh.generate(2)
    gmsh.write("meshes/" * savefile * ".msh")
    gmsh.finalize()
end

# Based on Nøst & Isachsen (2003) Appendix B.
# Their topography is H(x,y) = -H₀ exp(-(x² + y²)/L²) - H₁
# with H₀ = 3700 m, L = 600 km, H₁ = 300 m. 
# Note that H₀ here is different to the H₀ used elsewhere in the model.
# Here we non-dimensionalize with distance L, aspect ratio α = H₀/L and β = H₀/H₁.
# So in non-dimensional form, the topography is
# H(x,y) = -α/(β + 1) (β exp(-(x² + y²)) - 1)
# The maximum depth is α * (β - 1)/(β + 1) at the center (x,y) = (0,0).

# Parameters.
β  = 3700.0 / 300.0  # H₀ / H₁ from Nøst & Isachsen (2003)
depth(x) = α/(β + 1) * (β * exp(-(x[1]^2 + x[2]^2)) - 1)
# depth(x) = α*(1 - x[1]^2 - x[2]^2)          # Parabolic bowl for testing. Note that this goes to zero at the edges x = ±1.

function test_mesh_exp_bowl2D()
# Usage example:
    hh = 0.1
    filename = "test_exp_bowl2D.msh"
    generate_exp_bowl_mesh_2D(hh, depth; savefile=filename,DEBUG=false, show_gui=false)

# Test read to make sure that the Physical labels are correct. View the .msh file in Gmsh too.
    model = GmshDiscreteModel(filename)
    label_names = model.face_labeling.tag_to_name
    @info "GridapGmsh reads labels:"
    @info label_names
end