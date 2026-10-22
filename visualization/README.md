# 4×4 metasurface visualization

`render_metasurface.py` uses Blender 5.2.1's built-in Python API to reconstruct the structure shown in the supplied overview and annotated unit-cell diagrams. It builds a blue dielectric slab and thin gold sheet. Each of the 16 units has two semicircular cyan openings, a continuous gold center strip of width `d`, inner gold semicircular caps, and an angular gold split in each outer arc. The metal is one connected manifold sheet.

The diagrams give symbols `a`, `d`, `w`, `r`, `theta`, and `phi`, but no numeric values. Parameters in the script are relative visualization estimates, not fabrication or CST simulation dimensions.

Run from the repository root:

```powershell
blender -b --factory-startup --python visualization/render_metasurface.py -- --output-dir outputs
```

The default is the dark `midnight` style with an elevated orthographic camera. For a light figure with a higher camera:

```powershell
blender -b --factory-startup --python visualization/render_metasurface.py -- --output-dir outputs --style paper --view plan
```

Use separate output folders when keeping both versions. `--style` accepts `midnight` or `paper`; `--view` accepts `elevated` or `plan`. Camera location, target, and orthographic scale are near the end of the script; material colors are in the style block.

Each run creates an overview PNG, an oblique array close-up PNG, a top-view unit PNG, an editable `.blend` scene, and a portable `.glb` model in the chosen output folder. The `.blend` saves the overview camera. Generated renders and models are not committed to this code repository.
