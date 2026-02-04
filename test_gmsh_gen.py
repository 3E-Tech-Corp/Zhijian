#!/usr/bin/env python3
"""Generate test Gmsh MSH files in various formats and versions to test the reader."""

import gmsh
import os

output_dir = "/home/feng/repos/zhijian-fix/test_meshes"
os.makedirs(output_dir, exist_ok=True)

def gen_simple_quad_v41():
    """Generate a simple 2D quad mesh in MSH 4.1 ASCII format."""
    gmsh.initialize()
    gmsh.model.add("quad_v41")
    
    # Simple square
    gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
    gmsh.model.occ.synchronize()
    
    # Physical groups
    gmsh.model.addPhysicalGroup(2, [1], tag=1, name="domain")
    # Get boundary curves
    curves = gmsh.model.getBoundary([(2,1)], oriented=False)
    curve_tags = [c[1] for c in curves]
    gmsh.model.addPhysicalGroup(1, curve_tags, tag=2, name="wall")
    
    # Mesh parameters
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
    gmsh.option.setNumber("Mesh.Algorithm", 8)  # Frontal-Delaunay for quads
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    
    gmsh.model.mesh.generate(2)
    
    # Save as MSH 4.1 ASCII
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    path = os.path.join(output_dir, "quad_v41.msh")
    gmsh.write(path)
    print(f"Written: {path}")
    gmsh.finalize()
    return path

def gen_simple_tri_v41():
    """Generate a simple 2D triangular mesh in MSH 4.1 ASCII format."""
    gmsh.initialize()
    gmsh.model.add("tri_v41")
    
    # Simple square
    gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
    gmsh.model.occ.synchronize()
    
    # Physical groups
    gmsh.model.addPhysicalGroup(2, [1], tag=1, name="domain")
    curves = gmsh.model.getBoundary([(2,1)], oriented=False)
    curve_tags = [c[1] for c in curves]
    gmsh.model.addPhysicalGroup(1, curve_tags, tag=2, name="wall")
    
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
    
    gmsh.model.mesh.generate(2)
    
    # Save as MSH 4.1 ASCII
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    path = os.path.join(output_dir, "tri_v41.msh")
    gmsh.write(path)
    print(f"Written: {path}")
    gmsh.finalize()
    return path

def gen_simple_tri_v22():
    """Generate a simple 2D triangular mesh in MSH 2.2 ASCII format."""
    gmsh.initialize()
    gmsh.model.add("tri_v22")
    
    gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
    gmsh.model.occ.synchronize()
    
    gmsh.model.addPhysicalGroup(2, [1], tag=1, name="domain")
    curves = gmsh.model.getBoundary([(2,1)], oriented=False)
    curve_tags = [c[1] for c in curves]
    gmsh.model.addPhysicalGroup(1, curve_tags, tag=2, name="wall")
    
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
    
    gmsh.model.mesh.generate(2)
    
    # Save as MSH 2.2 ASCII
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    path = os.path.join(output_dir, "tri_v22.msh")
    gmsh.write(path)
    print(f"Written: {path}")
    gmsh.finalize()
    return path

def gen_order2_tri_v41():
    """Generate a second-order triangular mesh in MSH 4.1."""
    gmsh.initialize()
    gmsh.model.add("tri_order2_v41")
    
    gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
    gmsh.model.occ.synchronize()
    
    gmsh.model.addPhysicalGroup(2, [1], tag=1, name="domain")
    curves = gmsh.model.getBoundary([(2,1)], oriented=False)
    curve_tags = [c[1] for c in curves]
    gmsh.model.addPhysicalGroup(1, curve_tags, tag=2, name="wall")
    
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    
    gmsh.model.mesh.generate(2)
    
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    path = os.path.join(output_dir, "tri_order2_v41.msh")
    gmsh.write(path)
    print(f"Written: {path}")
    gmsh.finalize()
    return path

def gen_serendipity_quad_v41():
    """Generate a second-order serendipity quad mesh in MSH 4.1."""
    gmsh.initialize()
    gmsh.model.add("serendipity_v41")
    
    gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
    gmsh.model.occ.synchronize()
    
    gmsh.model.addPhysicalGroup(2, [1], tag=1, name="domain")
    curves = gmsh.model.getBoundary([(2,1)], oriented=False)
    curve_tags = [c[1] for c in curves]
    gmsh.model.addPhysicalGroup(1, curve_tags, tag=2, name="wall")
    
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)  # serendipity
    
    gmsh.model.mesh.generate(2)
    
    # Check element types
    types, _, _ = gmsh.model.mesh.getElements(2)
    print(f"  2D element types in serendipity mesh: {types}")
    
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    path = os.path.join(output_dir, "serendipity_quad_v41.msh")
    gmsh.write(path)
    print(f"Written: {path}")
    gmsh.finalize()
    return path

def gen_no_physgroups_v41():
    """Generate mesh WITHOUT physical groups (common user mistake)."""
    gmsh.initialize()
    gmsh.model.add("no_phys_v41")
    
    gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
    gmsh.model.occ.synchronize()
    
    # No physical groups defined!
    
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
    gmsh.option.setNumber("Mesh.SaveAll", 1)  # Save all elements
    
    gmsh.model.mesh.generate(2)
    
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    path = os.path.join(output_dir, "no_phys_v41.msh")
    gmsh.write(path)
    print(f"Written: {path}")
    gmsh.finalize()
    return path

def gen_no_physgroups_nosaveall_v41():
    """Generate mesh WITHOUT physical groups and WITHOUT SaveAll."""
    gmsh.initialize()
    gmsh.model.add("no_phys_nosaveall_v41")
    
    gmsh.model.occ.addRectangle(0, 0, 0, 1, 1)
    gmsh.model.occ.synchronize()
    
    # No physical groups, no SaveAll
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
    
    gmsh.model.mesh.generate(2)
    
    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.option.setNumber("Mesh.Binary", 0)
    path = os.path.join(output_dir, "no_phys_nosaveall_v41.msh")
    gmsh.write(path)
    print(f"Written: {path}")
    gmsh.finalize()
    return path

# Generate all test meshes
gen_simple_tri_v22()
gen_simple_tri_v41()
gen_simple_quad_v41()
gen_order2_tri_v41()
gen_serendipity_quad_v41()
gen_no_physgroups_v41()
gen_no_physgroups_nosaveall_v41()

print("\nAll test meshes generated!")
