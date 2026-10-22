# 4×4 metasurface visualization

`render_metasurface.py` uses Blender 5.2.1's built-in Python API to reconstruct the structure shown in the supplied illustration. It builds a blue dielectric slab, a thin gold sheet, and a 4×4 array of circular openings with two parallel metal bridges in each opening. The array occupies one half of the plate; the other half remains plain.

The illustration provides no physical dimensions. The script therefore uses relative units for visualization only; its model is not a CST or fabrication model.

Run from the repository root:

```powershell
blender -b --factory-startup --python visualization/render_metasurface.py -- --output-dir outputs
```

The command creates `THz_metasurface_4x4_render.png`, an editable `.blend` scene, and a portable `.glb` model in the chosen output folder. Dimensions, pattern spacing, aperture radius, and bridge widths are constants near the top of the script. The rendered image and model are generated files and are not committed to this code repository.
