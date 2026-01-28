import pyqpp as pq
from aux import overloader, tbl2clms, tbl2raws
from symm import symm_add, build_multab, symm_order, symm_invert, is_symm_group, \
    is_normal_subgroup, mul_groups, idx2grp, grp2idx, abelian_sub, find_generators

REAL = 'double'

vector_types = [pq.vector3f, pq.vector3d]

matrix_types=[pq.matrix3f,pq.matrix3d]

array_pg_types = [pq.array_point_group_f,pq.array_point_group_d]


cell_types = [pq.periodic_cell_f, pq.periodic_cell_d]

geom_types = [pq.geometry_f,pq.geometry_d]

#xgeom_types = [pq.xgeometry_f, pq.xgeometry_d]

sphere_shape_types = [pq.sphere_shape_f, pq.sphere_shape_d]

box_shape_types = [pq.box_shape_f, pq.box_shape_d]

btable_types = [pq.bonding_table_f, pq.bonding_table_d]

find_core_shells = pq.find_core_shells

def add_fd(s,r):
    fd=['_f','_d']
    return 'pq.'+s+fd[r]
# ------------------------------------------

def real_type(s):
    if s.lower() in ["f","float","f32","float32"]:
        return 0
    elif s.lower() in ["d","double","f64","float64"]:
        return 1
    else:
        raise TypeError("qpp: unknown type name for real number: " + s)

# ------------------------------------------

def vector3(*args, real = None):
    if real is None:
        r = real_type(REAL)
    else:
        r = real_type(real)
    return vector_types[r](*args)

def matrix3(*args, real = None):
    if real is None:
        r = real_type(REAL)
    else:
        r = real_type(real)
    return matrix_types[r](*args)

def array_point_group(*args, real = None):
    if real is None:
        r = real_type(REAL)
    else:
        r = real_type(real)
    return eval(add_fd('array_point_group',r))(*args)

def geometry(*args,**kwds):
    r = real_type(REAL)
    return eval(add_fd('geometry',r))(*args,**kwds)

def mm_calculator(*args,**kwds):
    r = real_type(REAL)
    return eval(add_fd('mm_calculator',r))(*args,**kwds)

def xgeometry(*args,**kwds):
    r = real_type(REAL)
    #print(args,kwds)
    xg = eval(add_fd('xgeometry',r))(*args,kwds)
    for f in range(xg.nfields()):
        fname = xg.fnames[f]
        setattr(xg,fname, xgeometry_field_attr(xg,f))
    return xg

def coulomb_point_charges(*args,**kwds):
    real1 = kwds.get('real',REAL)
    #print(kwds)
    #print('real=',real1)
    if real1 is not None:
    #    print('WTF real=', real1)
        r = real_type(real1)
    else:
    #    print('real is None')
        r = real_type(REAL)
    return eval(add_fd('coulomb_point_charges',r))(*args,**kwds)

globals = pq.globals
pp = pq.pp

#def periodic_cell(*args,**kwds):
#    r = real_type(REAL)
#    return periodic_cell[r](*args,**kwds)

def sphere_shape(*args,**kwds):
    r = real_type(REAL)
    return sphere_shape_types[r](*args,**kwds)

def bonding_table(*args,**kwds):
    r = real_type(REAL)
    return eval(add_fd('bonding_table',r))(*args,**kwds)

def ngbr_table(*args,**kwds):
    r = real_type(REAL)
    return eval(add_fd('ngbr_table',r))(*args,**kwds)

 #def index_empty(d):
#    if type(d) is int:
#        return index([0]*d)
#    else: raise TypeError()
class xgeometry_field_attr(object):
    
    def __init__(self,geom1,f1):
        self.geom = geom1
        self.f=f1

    def __getitem__(self,i):
        return self.geom.field[self.f, i]

    def __setitem__(self,i,val):
        self.geom.field[self.f,i] = val

#---------------------------------------------------------



#---------------------------------------------------------

def cell_empty(dim,real="d"):
    if type(real) is str:
        r = real_type(real)
        return cell_types[r][dim]()
    else:
        raise TypeError()

def cell_copy(x):
    if type(x) in cell_types:
        return type(x)(x)
    else:
        raise TypeError()

def cell_from_1vectors(a):
    if type(a) in vector_types:
        r = vector_types.index(type(a))
        return cell_types[r](a)
    else:
        raise TypeError()

def cell_dim(dim):
    if not isinstance(dim,int):
        raise TypeError()
    else:
        r = real_type(REAL)
        return cell_types[r](dim)

def cell_from_2vectors(a,b):
    if type(a) in vector_types and type(b) in vector_types:
        r = vector_types.index(type(a))
        return cell_types[r](a,b)
    else:
        raise TypeError()

def cell_from_3vectors(a,b,c):
    if type(a) in vector_types and type(b) in vector_types and type(c) in vector_types :
        r = vector_types.index(type(a))
        return cell_types[r](a,b,c)
    else:
        raise TypeError()

def periodic_cell(*args,**kwds):
    return overloader([cell_empty, cell_copy, cell_dim, cell_from_1vectors, cell_from_2vectors, cell_from_3vectors], \
                      args, kwds, \
                      "qpp: invalid arguments of qpp.cell call"+str(args)+str(kwds))

#---------------------------------------------------------
'''
class xgeometry_field_attr(object):
    
    def __init__(self,geom1,f1):
        self.geom = geom1
        self.f=f1

    def __getitem__(self,i):
        return self.geom.field[self.f, i]

    def __setitem__(self,i,val):
        self.geom.field[self.f,i] = val

def xgeometry_from_list(real,cell,fields,name=''):
    rt = real_type(real)
    ctypes = [periodic_cell_f, point_group_f, crystal_group_f, periodic_cell_d, point_group_d, crystal_group_d]
    gtypes=[xgeometry_f, xgeometry_pgf, xgeometry_cgf, xgeometry_d, xgeometry_pgd, xgeometry_cgd]
    rtypes = [0,0,0,1,1,1]
    
    t = ctypes.index(type(cell))
    if rtypes[t]!=rt:
        raise ValueError("Real type in xgeometry.__init__ does not much the supercell type")

    print(gtypes[t], fields,name)
    xgeom = gtypes[t](cell,fields,name)

    for f in range(xgeom.nfields()):
        fname = xgeom.field_name(f)
        if not fname in ['atom','x','y','z']:
            setattr(xgeom,fname, xgeometry_field_attr(xgeom,f))

    return xgeom

def xgeometry_from_dict(real,cell,**kwds):    
    name = kwds.pop('name','')    
    return xgeometry_from_list(real,cell,list(kwds.items()),name)
    
def xgeometry(*args,**kwds):
    return overloader([xgeometry_from_list,xgeometry_from_dict], \
                      args, kwds, \
                      "qpp: invalid arguments of qpp.xgeometry call"+str(args)+str(kwds))
'''
index_range = pq.index_range

io = pq.io

from pyqpp import fill
