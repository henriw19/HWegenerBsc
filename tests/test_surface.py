from hwbsc.codes import repetition_code
from hwbsc.stab_codes import hypergraph_prod
from hwbsc.surface_code import *

rep = repetition_code(3)
sur = hypergraph_prod(rep).H
SurfCode(sur).draw(show_logicals=True)

