"""Procedurally build and render the 4x4 THz metasurface in 图片1.png.

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

# Relative proportions reconstructed from the illustration. The substrate is
# four pitches wide, eight pitches long; the 4x4 motifs occupy its front half.
PITCH = 1.0
WIDTH = 4.0 * PITCH
LENGTH = 8.0 * PITCH
SUBSTRATE_THICKNESS = 0.60 * PITCH
METAL_THICKNESS = 0.075 * PITCH
APERTURE_RADIUS = 0.39 * PITCH
BAR_WIDTH = 0.115 * PITCH
BAR_CENTERS = (-0.145 * PITCH, 0.145 * PITCH)
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


def cylinder_mesh_at(cx, cy, radius, zmin, zmax, vertices, faces):
    start = len(vertices)
    for z in (zmin, zmax):
        for i in range(SEGMENTS):
            a = 2 * math.pi * i / SEGMENTS
            vertices.append((cx + radius * math.cos(a), cy + radius * math.sin(a), z))
    n = SEGMENTS
    faces.append(tuple(start + i for i in reversed(range(n))))
    faces.append(tuple(start + n + i for i in range(n)))
    for i in range(n):
        j = (i + 1) % n
        faces.append((start + i, start + j, start + n + j, start + n + i))


# One cutter with 16 disconnected cylinders gives one clean Boolean operation.
cut_vertices, cut_faces = [], []
for y in MOTIF_Y:
    for x in MOTIF_X:
        cylinder_mesh_at(
            x, y, APERTURE_RADIUS,
            -0.04, METAL_THICKNESS + 0.04,
            cut_vertices, cut_faces,
        )
cut_mesh = bpy.data.meshes.new("16 circular aperture cutters")
cut_mesh.from_pydata(cut_vertices, [], cut_faces)
cut_mesh.update()
cutter = bpy.data.objects.new("Temporary aperture cutter", cut_mesh)
model_collection.objects.link(cutter)
boolean = metal.modifiers.new("16 circular apertures", "BOOLEAN")
boolean.operation = "DIFFERENCE"
boolean.solver = "EXACT"
boolean.object = cutter
bpy.context.view_layer.objects.active = metal
bpy.ops.object.modifier_apply(modifier=boolean.name)
bpy.data.objects.remove(cutter, do_unlink=True)
edge = metal.modifiers.new("Soft metal edges", "BEVEL")
edge.width = 0.009
edge.segments = 3
edge.affect = "EDGES"
metal.modifiers.new("Weighted normals", "WEIGHTED_NORMAL")


def clip_x(poly, boundary, keep_greater):
    """Clip a counterclockwise polygon to x >= or x <= boundary."""
    result = []
    for a, b in zip(poly, poly[1:] + poly[:1]):
        ain = a[0] >= boundary if keep_greater else a[0] <= boundary
        bin = b[0] >= boundary if keep_greater else b[0] <= boundary
        if ain:
            result.append(a)
        if ain != bin:
            t = (boundary - a[0]) / (b[0] - a[0])
            result.append((boundary, a[1] + t * (b[1] - a[1])))
    return result


def extruded_polygon(poly, cx, cy, zmin, zmax, vertices, faces):
    start = len(vertices)
    n = len(poly)
    if n < 3:
        raise ValueError("A metal bridge collapsed during circle clipping")
    for z in (zmin, zmax):
        vertices.extend((cx + px, cy + py, z) for px, py in poly)
    faces.append(tuple(start + i for i in reversed(range(n))))
    faces.append(tuple(start + n + i for i in range(n)))
    for i in range(n):
        j = (i + 1) % n
        faces.append((start + i, start + j, start + n + j, start + n + i))


# The reference has two parallel gold bridges running along the long axis of
# each circular opening. Clipping them to the circle preserves three cyan gaps.
circle = [
    (APERTURE_RADIUS * math.cos(2 * math.pi * i / SEGMENTS),
     APERTURE_RADIUS * math.sin(2 * math.pi * i / SEGMENTS))
    for i in range(SEGMENTS)
]
bridge_vertices, bridge_faces = [], []
for y in MOTIF_Y:
    for x in MOTIF_X:
        for center in BAR_CENTERS:
            bridge = clip_x(circle, center - BAR_WIDTH / 2, True)
            bridge = clip_x(bridge, center + BAR_WIDTH / 2, False)
            extruded_polygon(
                bridge, x, y, 0, METAL_THICKNESS,
                bridge_vertices, bridge_faces,
            )
bridge_mesh = bpy.data.meshes.new("32 clipped gold bridges")
bridge_mesh.from_pydata(bridge_vertices, [], bridge_faces)
bridge_mesh.update()
bridges = bpy.data.objects.new("Two gold bridges per aperture", bridge_mesh)
model_collection.objects.link(bridges)
bridges.data.materials.append(gold)
bridge_bevel = bridges.modifiers.new("Soft bridge edges", "BEVEL")
bridge_bevel.width = 0.006
bridge_bevel.segments = 2
bridge_bevel.affect = "EDGES"
bridges.modifiers.new("Weighted normals", "WEIGHTED_NORMAL")

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
bpy.ops.wm.save_as_mainfile(filepath=blend_path)
bpy.ops.object.select_all(action="DESELECT")
for obj in (substrate, metal, bridges):
    obj.select_set(True)
bpy.context.view_layer.objects.active = metal
bpy.ops.export_scene.gltf(
    filepath=glb_path,
    export_format="GLB",
    use_selection=True,
    export_apply=True,
)
bpy.ops.object.select_all(action="DESELECT")
bpy.ops.render.render(write_still=True)
print(f"BLEND={blend_path}")
print(f"GLB={glb_path}")
print(f"RENDER={image_path}")

