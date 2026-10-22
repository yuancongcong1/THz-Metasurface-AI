"""Procedurally build and render the 4x4 THz metasurface in 图片1/图片2.

Run with Blender 5.2+:
    blender -b --factory-startup --python render_metasurface.py -- --output-dir PATH

The source illustration has no dimension labels. All dimensions below are
relative visualization units, not fabrication or simulation parameters.
"""

import argparse
import math
import os
import sys

import bpy
import bmesh
from mathutils import Vector


def args_after_double_dash():
    return sys.argv[sys.argv.index("--") + 1 :] if "--" in sys.argv else []


parser = argparse.ArgumentParser()
parser.add_argument("--output-dir", default=os.path.join(os.getcwd(), "outputs"))
parser.add_argument("--resolution-x", type=int, default=2200)
parser.add_argument("--resolution-y", type=int, default=1500)
args = parser.parse_args(args_after_double_dash())
output_dir = os.path.abspath(args.output_dir)
os.makedirs(output_dir, exist_ok=True)

blend_path = os.path.join(output_dir, "THz_metasurface_4x4.blend")
glb_path = os.path.join(output_dir, "THz_metasurface_4x4.glb")
image_path = os.path.join(output_dir, "THz_metasurface_4x4_render.png")
detail_path = os.path.join(output_dir, "THz_metasurface_unit_detail.png")

# Relative proportions reconstructed from the illustration. The substrate is
# four pitches wide, eight pitches long; the 4x4 motifs occupy its front half.
PITCH = 1.0
WIDTH = 4.0 * PITCH
LENGTH = 8.0 * PITCH
SUBSTRATE_THICKNESS = 0.60 * PITCH
METAL_THICKNESS = 0.075 * PITCH
OUTER_RADIUS = 0.39 * PITCH
RING_WIDTH = 0.10 * PITCH       # w in the unit-cell diagram
INNER_RADIUS = OUTER_RADIUS - RING_WIDTH  # r in the diagram
CENTER_BAND = 0.10 * PITCH      # d: continuous gold strip between halves
INNER_CAP_CHORD = 0.14 * PITCH  # location of each inner semicircle's flat edge
THETA_DEG = 22.0               # upper radial split, measured toward +y
PHI_DEG = 22.0                 # lower radial split, measured toward +y
SPLIT_SPAN_DEG = 15.0          # visual width of each radial split
MOTIF_X = (-1.5, -0.5, 0.5, 1.5)
MOTIF_Y = (-3.5, -2.5, -1.5, -0.5)
SEGMENTS = 96

bpy.ops.object.select_all(action="SELECT")
bpy.ops.object.delete(use_global=False)
for collection in list(bpy.data.collections):
    if collection.name != "Collection" and collection.users == 0:
        bpy.data.collections.remove(collection)


def make_collection(name):
    collection = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(collection)
    return collection


model_collection = make_collection("Metasurface model")
studio_collection = make_collection("Studio")


def move_to_collection(obj, collection):
    for old in list(obj.users_collection):
        old.objects.unlink(obj)
    collection.objects.link(obj)
    return obj


def material(name, color, metallic=0.0, roughness=0.35):
    mat = bpy.data.materials.new(name)
    mat.diffuse_color = (*color, 1.0)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    bsdf.inputs["Base Color"].default_value = (*color, 1.0)
    bsdf.inputs["Metallic"].default_value = metallic
    bsdf.inputs["Roughness"].default_value = roughness
    return mat


gold = material("Satin gold | patterned metal", (0.92, 0.67, 0.13), 0.72, 0.27)
cyan = material("Cyan dielectric substrate", (0.045, 0.52, 0.72), 0.10, 0.28)
floor_mat = material("Deep slate studio floor", (0.018, 0.034, 0.057), 0.13, 0.48)


def cube(name, dimensions, location, mat, collection, bevel=0.0):
    bpy.ops.mesh.primitive_cube_add(size=1, location=location)
    obj = bpy.context.object
    obj.name = name
    obj.dimensions = dimensions
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)
    obj.data.materials.append(mat)
    move_to_collection(obj, collection)
    if bevel:
        mod = obj.modifiers.new("Subtle machined edges", "BEVEL")
        mod.width = bevel
        mod.segments = 3
        mod.affect = "EDGES"
        obj.modifiers.new("Weighted normals", "WEIGHTED_NORMAL")
    return obj


# The top of the substrate is z=0; the cutouts reveal this cyan top face.
substrate = cube(
    "Dielectric slab | 4 x 8 x 0.6",
    (WIDTH, LENGTH, SUBSTRATE_THICKNESS),
    (0, 0, -SUBSTRATE_THICKNESS / 2),
    cyan,
    model_collection,
    bevel=0.025,
)
metal = cube(
    "Patterned gold sheet | 4 x 8",
    (WIDTH, LENGTH, METAL_THICKNESS),
    (0, 0, METAL_THICKNESS / 2),
    gold,
    model_collection,
)


def extruded_polygon(poly, cx, cy, zmin, zmax, vertices, faces):
    signed_area = sum(a[0] * b[1] - b[0] * a[1]
                      for a, b in zip(poly, poly[1:] + poly[:1]))
    if signed_area < 0:
        poly = list(reversed(poly))
    start = len(vertices)
    n = len(poly)
    if n < 3:
        raise ValueError("A unit-cell region collapsed during clipping")
    for z in (zmin, zmax):
        vertices.extend((cx + px, cy + py, z) for px, py in poly)
    faces.append(tuple(start + i for i in reversed(range(n))))
    faces.append(tuple(start + n + i for i in range(n)))
    for i in range(n):
        j = (i + 1) % n
        faces.append((start + i, start + j, start + n + j, start + n + i))


def circular_arc(radius, start, end):
    steps = max(4, math.ceil(SEGMENTS * abs(end - start) / (2 * math.pi)))
    return [(radius * math.cos(start + (end - start) * i / steps),
             radius * math.sin(start + (end - start) * i / steps))
            for i in range(steps + 1)]


def split_semicircle_opening(upper, center_deg):
    """Trace one cyan opening as a simple contour with the gold notch removed.

    Subtracting this contour directly from the sheet leaves the inner gold cap
    and its radial connection as part of the SAME solid, without a mesh union.
    x is perpendicular to the center band; y runs along that band.
    """
    outer_end = math.acos((CENTER_BAND / 2) / OUTER_RADIUS)
    inner_end = math.acos(INNER_CAP_CHORD / INNER_RADIUS)
    left_split = math.radians(center_deg - SPLIT_SPAN_DEG / 2)
    right_split = math.radians(center_deg + SPLIT_SPAN_DEG / 2)
    assert 0 < left_split < right_split < inner_end < outer_end

    contour = circular_arc(OUTER_RADIUS, -outer_end, left_split)
    contour.append((INNER_RADIUS * math.cos(left_split),
                    INNER_RADIUS * math.sin(left_split)))
    contour.extend(circular_arc(INNER_RADIUS, left_split, -inner_end)[1:])
    # A straight chord is the lower edge of the inner gold cap.
    contour.append((INNER_RADIUS * math.cos(inner_end),
                    INNER_RADIUS * math.sin(inner_end)))
    contour.extend(circular_arc(INNER_RADIUS, inner_end, right_split)[1:])
    contour.append((OUTER_RADIUS * math.cos(right_split),
                    OUTER_RADIUS * math.sin(right_split)))
    contour.extend(circular_arc(OUTER_RADIUS, right_split, outer_end)[1:])
    # Closing the contour supplies the straight edge beside the center band.
    return contour if upper else [(-x, y) for x, y in contour]


# Each unit has two semicircular openings around one continuous gold strip.
# A single subtraction per opening leaves one connected gold sheet.
upper_opening = split_semicircle_opening(True, THETA_DEG)
lower_opening = split_semicircle_opening(False, PHI_DEG)

cut_vertices, cut_faces = [], []
for y in MOTIF_Y:
    for x in MOTIF_X:
        for opening in (upper_opening, lower_opening):
            extruded_polygon(
                opening, x, y, -0.04, METAL_THICKNESS + 0.04,
                cut_vertices, cut_faces,
            )
cut_mesh = bpy.data.meshes.new("32 semicircular opening cutters")
cut_mesh.from_pydata(cut_vertices, [], cut_faces)
cut_mesh.update()
cutter = bpy.data.objects.new("Temporary paired-opening cutter", cut_mesh)
model_collection.objects.link(cutter)
boolean = metal.modifiers.new("32 split-ring openings", "BOOLEAN")
boolean.operation = "DIFFERENCE"
boolean.solver = "EXACT"
boolean.object = cutter
bpy.context.view_layer.objects.active = metal
bpy.ops.object.modifier_apply(modifier=boolean.name)
bpy.data.objects.remove(cutter, do_unlink=True)

# Enforce connectivity before saving either model file.
if metal.data.validate(verbose=False):
    raise RuntimeError("Gold mesh required validation repairs")
topology = bmesh.new()
topology.from_mesh(metal.data)
if any(not item.is_manifold for item in topology.edges):
    raise RuntimeError("Gold mesh has an open or non-manifold edge")
unseen = set(topology.verts)
components = 0
while unseen:
    components += 1
    stack = [unseen.pop()]
    while stack:
        vertex = stack.pop()
        for item in vertex.link_edges:
            neighbor = item.other_vert(vertex)
            if neighbor in unseen:
                unseen.remove(neighbor)
                stack.append(neighbor)
topology.free()
if components != 1:
    raise RuntimeError(f"Gold mesh is split into {components} components")
print("GOLD_CONNECTED_COMPONENTS=1; GOLD_NONMANIFOLD_EDGES=0")

edge = metal.modifiers.new("Soft continuous metal edges", "BEVEL")
edge.width = 0.006
edge.segments = 3
edge.affect = "EDGES"
metal.modifiers.new("Weighted normals", "WEIGHTED_NORMAL")

# Studio environment is excluded from the portable GLB model.
floor = cube(
    "Studio floor", (200, 200, 0.05), (0, 0, -0.75),
    floor_mat, studio_collection,
)


def area_light(name, location, power, size, color):
    data = bpy.data.lights.new(name, "AREA")
    data.energy = power
    data.shape = "DISK"
    data.size = size
    data.color = color
    obj = bpy.data.objects.new(name, data)
    studio_collection.objects.link(obj)
    obj.location = location
    direction = Vector((0, 0, -0.10)) - obj.location
    obj.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()
    return obj


area_light("Warm key", (-4, -6, 9), 1700, 7, (1.0, 0.86, 0.67))
area_light("Cool fill", (6, 1, 7), 1150, 8, (0.68, 0.85, 1.0))
area_light("Back rim", (-1, 6, 6), 1050, 6, (0.83, 0.94, 1.0))

camera_data = bpy.data.cameras.new("Orthographic isometric camera")
camera = bpy.data.objects.new("Orthographic isometric camera", camera_data)
studio_collection.objects.link(camera)
camera.location = (12, -12, 12)
target = Vector((0, 0, -0.12))
camera.rotation_euler = (target - camera.location).to_track_quat("-Z", "Y").to_euler()
camera_data.type = "ORTHO"
camera_data.ortho_scale = 10.5
bpy.context.scene.camera = camera

scene = bpy.context.scene
scene.render.engine = "BLENDER_EEVEE"
scene.render.resolution_x = args.resolution_x
scene.render.resolution_y = args.resolution_y
scene.render.resolution_percentage = 100
scene.render.image_settings.file_format = "PNG"
scene.render.film_transparent = False
scene.render.filepath = image_path
scene.world.color = (0.015, 0.025, 0.045)
scene.view_settings.view_transform = "AgX"
scene.view_settings.look = "AgX - Medium High Contrast"
scene.render.image_settings.color_mode = "RGBA"

# Save an editable Blender scene and a portable model containing just geometry.
bpy.context.preferences.filepaths.save_version = 0
bpy.ops.wm.save_as_mainfile(filepath=blend_path)
bpy.ops.object.select_all(action="DESELECT")
for obj in (substrate, metal):
    obj.select_set(True)
bpy.context.view_layer.objects.active = metal
bpy.ops.export_scene.gltf(
    filepath=glb_path,
    export_format="GLB",
    use_selection=True,
    # Export the manifold base mesh. The Blender scene retains the visual
    # bevel modifiers, while GLB avoids vertex splits from sharp normals.
    export_apply=False,
    export_normals=False,
)
bpy.ops.object.select_all(action="DESELECT")
bpy.ops.render.render(write_still=True)

# A direct top view makes the corrected inner arcs and angular splits easy to
# compare with the user's annotated unit-cell diagram.
camera.location = (0.5, -2.5, 10)
# Rotate the top view so the central strip is horizontal, as in 图片2.
camera.rotation_euler = (0, 0, math.pi / 2)
camera_data.ortho_scale = 1.18
scene.render.resolution_x = 1200
scene.render.resolution_y = 1200
scene.render.filepath = detail_path
bpy.ops.render.render(write_still=True)
print(f"BLEND={blend_path}")
print(f"GLB={glb_path}")
print(f"RENDER={image_path}")
print(f"DETAIL={detail_path}")

