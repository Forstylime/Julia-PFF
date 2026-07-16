using WriteVTK
using LightXML

# VTK files produced by WriteVTK are standard XML. We can parse them using LightXML.
xdoc = parse_file("data/sims/rlm_amor/fracture_step_500.vtu")
xroot = root(xdoc)
# Since reading appended binary data is hard, let's just use Ferrite to read the VTU? 
# Ferrite doesn't have a reader easily available.
