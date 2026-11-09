TELEMAC MESH WORKFLOW

Requirements
- MATLAB, Gmsh, and AutoCAD for the drawing step.
- Use the same projected XY coordinates in metres for every input.

AutoCAD
- Run data_autocad in MATLAB and select a three-column X Y Z file.
- Points open in AutoCAD with elevation-based layers and colours; original XYZ values are preserved.
- Draw the outline with PLINE. Type C and press Enter to close it.
- Load outline_to_txt.lsp with APPLOAD, then run OUTLINETXT.
- Select the closed polyline and save the ordered XY outline as TXT.

Preprocess
- Run telemac_preprocess.
- Select the outline, bathymetry XYZ, and optional island boundary file.
- Island TXT: use closed XY rings, separated by NaN NaN for multiple islands.
- Enter global boundary/interior mesh sizes and bathymetry settings.
- Optional refinement: select Rectangle, Circle, or Polygon; add multiple regions if needed.
- Rectangle: two opposite corners. Circle: centre and edge. Polygon: vertices, then Enter.
- Set each region's target size and transition coefficient; smaller coefficients above 1 give wider transitions.
- Mouse wheel: zoom. Right-drag or Pan: move the view. Esc: cancel the current selection.
- Intermediate grid spacing affects bathymetry interpolation, not triangle size; 0 uses direct interpolation.
- Select the Z multiplier: 1 for bed elevations, -1 for positive depths; set an offset if needed.
- Choose the output path. Gmsh builds the mesh; only the SLF is saved permanently.

Boundary
- Run telemac_boundary and select the generated SLF; no remeshing is required.
- Click Add segment, select two boundary nodes, and confirm the highlighted path.
- Enter the segment name, hydraulic boundary type, tracer type, and values.
- Unassigned nodes remain closed walls (222); add further segments as needed.
- Click Save CLI + BC2. Both files are saved; named CLI rows end with # segment_name.
- Recreate boundary files whenever the mesh changes.

Status
- Early version; tested with MATLAB R2024b and Gmsh 4.13.1 on macOS.
- The boundary editor supports 2D triangular SLF geometry.
- AutoCAD outline export supports closed planar 2D polylines.
- Verify coordinate systems, vertical datum, mesh quality and boundary settings for each application.

Provided case datasets
- examples/sea_of_marmara.zip: Sea of Marmara outline, bathymetry and islands.
- examples/gulf_of_saros.zip: Gulf of Saros outline, bathymetry and islands.
- examples/black_sea.zip: Black Sea outline, bathymetry and islands.
- Extract a case ZIP inside examples; select its outline, bathymetry and island files in telemac_preprocess.
- Each extracted case folder includes a short README mapping the input files.
