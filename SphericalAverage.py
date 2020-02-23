#! /usr/bin/env python
import load_cube
#import NewPotentialModule as pot
import numpy as np
from ase.io import read, write
from ase.io.trajectory import Trajectory

def cube_potential(origin,travelled,cube,Grid,nx,ny,nz):
    """Populates the sampling cube with the potential required"""

# Recalc the origin as grid point coordinates                                                            
    n_origin = np.zeros(shape=(3))
    n_origin[0] = int(origin[0]*nx)
    n_origin[1] = int(origin[1]*ny)
    n_origin[2] = int(origin[2]*nz)
    potential_cube = np.zeros(shape=(cube[0],cube[1],cube[2]))
    for x in range(0,cube[0]):
        for y in range(0,cube[1]):
            for z in range(0,cube[2]):
# Assign the values of coordinates in the original grid                                                  
                xv = n_origin[0]+travelled[0]+x
                yv = n_origin[1]+travelled[1]+y
                zv = n_origin[2]+travelled[2]+z
# Minimum image convention                                                                               
                zv = zv - nz*round(zv/nz)
                yv = yv - ny*round(yv/ny)
                xv = xv - nx*round(xv/nx)
                potential_cube[x,y,z] = Grid[xv,yv,zv]
    print 'potential_cube'  
    print potential_cube 
    return potential_cube.mean(), np.var(potential_cube)

def fractional(coordinates,atoms):
    """Converts an array of xyz coordinates to fractional coordinates for a given Atoms object """
    #import pdb;pdb.set_trace()
    coord_list =[]
    for coordinate in coordinates:
        fcoord = np.dot(coordinate,atoms.get_reciprocal_cell().T)
        coord_list.append(fcoord)
        fractional_coordinates = np.vstack(coord_list)
    return fractional_coordinates

input_file = '/Users/hxin/blueridge/hxin/mofs/PCN-223/Co/PBE-D3/bandstru/PBE-D3/PCN-223-Co-v_hartree-1_0.cube'
traj_file = 'PCN-223-Co.traj'
cube_size = [2,2,2]    # This size is in units of mesh points                                                                               
## origin defines the bottom left point of the cube the "0,0,0" point in fractional coordinates                                             
atoms_for_origin = read('for_origin.xyz')
positions =  atoms_for_origin.get_scaled_positions()
#cube_origin = (positions[184] + positions[185] + positions[729])/3
cube_origin = (positions[184] + positions[183] + positions[458])/3
cube_origin = [0.5,np.sqrt(3.0)/6,0.5]
#cube_origin = [0, 12.42, 8.445]
#cube_origin = [10.8, 6.3, 8.52]
print "cube_origin is" + str(cube_origin)

# No need to alter anything after here                                                                                                      
#------------------------------------------------------------------                                                                         
# Get the potential                                                                                                                         
# This section should not be altered                                                                                                        
#------------------------------------------------------------------                                                                        # create an object and read in data from file 
cube = load_cube.CUBE(input_file)

#print cube.data
potential = cube.data
#print potential.shape
#print potential
# calculate total number of data points
#print "Number of voxels: %d"%(cube.NX*cube.NY*cube.NZ)
# calculate box volume
#print "Total volume, Angs^3: %d"%((cube.NX-1)*cube.X[0]*(cube.NY-1)*cube.Y[1]*(cube.NZ-1)*cube.Z[2])
# calculate total electron density by summing up squared values, contained in 3d array cube.data 
#print "Number of electrons (sort of): %.2f"%((cube.data**2).sum()*cube.X[0]*cube.Y[1]*cube.Z[2])

traj = Trajectory(traj_file)
atoms = traj[-1]
print atoms.get_cell()
print 'len(positions of atoms) are ' + str(len(atoms.get_positions()))
print atoms.get_scaled_positions()
#print fractional(atoms.get_positions(),atoms)
vector_a = np.linalg.norm(atoms.cell[1])
vector_b = np.linalg.norm(atoms.cell[1])
vector_c = np.linalg.norm(atoms.cell[2])
NGX = len(potential)
NGY = len(potential[0])
NGZ = len(potential[0][0])
resolution_x = vector_a/NGX
resolution_y = vector_b/NGY
resolution_z = vector_c/NGZ
print NGX,NGY,NGZ
#------------------------------------------------------------------                                                                         

##------------------------------------------------------------------                                                                        
# Getting the average potential in a single cube of arbitrary size                                                                          
##------------------------------------------------------------------                                                                        
## cube defines the size of the cube in units of mesh points (NGX/Y/Z)                                                                      
cube = cube_size
## origin defines the bottom left point of the cube the "0,0,0" point in fractional coordinates                                             
origin = cube_origin
## travelled; do not alter this variable                                                                                                    
travelled = [0,0,0]
## Uncomment the lines below to do the business                                                                                             
cube_potential, cube_var = cube_potential(origin,travelled,cube,potential,NGX,NGY,NGZ)
print "Potential            Variance"
print "--------------------------------"
print cube_potential,"   ", cube_var
##------------------------------------------------------------------                                                                        
##------------------------------------------------------------------                                                                        
