# 4×4 metasurface visualization

`render_metasurface.py` uses Blender 5.2.1's built-in Python API to reconstruct the structure shown in the supplied overview and annotated unit-cell diagrams. It builds a blue dielectric slab and thin gold sheet. Each of the 16 units has **two semicircular cyan openings**, a continuous gold center strip of width `d`, inner gold semicircular caps, and an angular gold split in each outer arc. The metal is cut from one connected manifold sheet: the inner caps join the surrounding metal through the radial splits. The array occupies one half of the plate; the other half remains plain.

The diagrams give symbols `a`, `d`, `w`, `r`, `theta`, and `phi`, but no numeric values. The script exposes corresponding relative parameters near the top. These are visualization estimates, not fabrication or CST simulation dimensions.

Run from the repository root:

```powershell
blender -b --factory-startup --python visualization/render_metasurface.py -- --output-dir outputs
```

This creates an isometric PNG, a top-view PNG of one unit, an editable `.blend` scene, and a portable `.glb` model in the chosen output folder. Generated renders and models are not committed to this code repository.
