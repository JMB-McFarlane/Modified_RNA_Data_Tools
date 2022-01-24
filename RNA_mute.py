"""
James McFarlane 2020 – Wetmore Group – University of Lethbridge
mcfarlane.james.mb@gmail.com

Mutator program that takes PDBs saved from LEaP non-standard residue units (uncapped)
and replaces defined residues in a PDB. Currently working with DNA duplexes and requires manual
variable changes and all necessary PDBs in local directory of execution.

"""
import signal
import os
import numpy as np
import numpy.linalg as linalg
#import matplotlib.pyplot as plt
import argparse
#mod_nuc_dir = "/Users/james/Data/Modified_Nucleotides/Edited-Post-RED/"
mod_nuc_dir = "./residues/"

parser = argparse.ArgumentParser(description='A program that mutates residues')
parser.add_argument("-i",type=str,action='store',dest='unmod_complex',help="Original PDB to be mutated")
parser.add_argument("-r",type=str,action='store',dest='mut_sel',help="Residue for insertion")
parser.add_argument("-s",type=int,action='store',dest='mut_site',help="Mutation site")
parser.add_argument("-m",type=str,action='store',dest='mode',help="Insertion Mode stack or replace")
parser.add_argument("-o",type=str,action='store',dest='output',help="Output file name")
args = parser.parse_args()

#mut_site = 4
#mut_sel = "I3A"   #will be connected to a library of residue names#
unmod_complex = args.unmod_complex
output_filename = args.output
mut_site = args.mut_site
mut_sel = args.mut_sel

#global vars
prep_pdb = open(mod_nuc_dir + mut_sel+".pdb")
old_pdb = open(unmod_complex + ".pdb")
mod_coords = []
prev_mut_coords = []
after_mut_coords = []
prev_mut_N_coords = []
after_mut_N_coords = []
orig_mut_N_coords = []
opp =[]
res_n_list = []



for line in old_pdb:
	if len(line.split()) >= 5:


		if line.split()[4] == str(mut_site):
			print("halp")
			if line.split()[2] == "P":
				P_coords = [float(line.split()[5]),float(line.split()[6]),float(line.split()[7])]
				original_P = [float(line.split()[5]),float(line.split()[6]),float(line.split()[7])]
		if line.split()[4] == str(mut_site):
			if "O3'" in line.split()[2]:
				O_coords = [float(line.split()[5]),float(line.split()[6]),float(line.split()[7])]
		if line.split()[1] == "254":
			opp_ref = [float(line.split()[5]),float(line.split()[6]),float(line.split()[7])]
		if line.split()[4] == str(mut_site -1):
			split = line.split()
			print("error")
			line_coords = [float(split[5]),float(split[6]),float(split[7])]
			prev_mut_coords.append(line_coords)
		if line.split()[4] == str(mut_site +1):
			split = line.split()
			line_coords = [float(split[5]),float(split[6]),float(split[7])]
			after_mut_coords.append(line_coords)
		if line.split()[4] == str(8):
			split = line.split()
			line_coords = [float(split[5]),float(split[6]),float(split[7])]
			opp.append(line_coords)

		#### Nitrogen coords for N-centroid ####
		if line.split()[4] == str(mut_site +1):
			if "N" in line.split()[2]:
				split = line.split()
				line_coords = [float(split[5]),float(split[6]),float(split[7])]
				#print("N found")
				after_mut_N_coords.append(line_coords)

		if line.split()[4] == str(mut_site -1):
			if "N" in line.split()[2]:
				split = line.split()
				line_coords = [float(split[5]),float(split[6]),float(split[7])]
				prev_mut_N_coords.append(line_coords)

		if line.split()[4] == str(mut_site):
                        if "N" in line.split()[2]:
                                split = line.split()
                                line_coords = [float(split[5]),float(split[6]),float(split[7])]
                                orig_mut_N_coords.append(line_coords)
		
		#print(orig_mut_N_coords)
for line in prep_pdb:
	if len(line.split()) >= 8:
		split = line.split()
		line_coords = [float(split[5]),float(split[6]),float(split[7])]
		mod_coords.append(line_coords)
		if line.split()[2][0]== "N":
			res_n_list.append(int(line.split()[1])-1)
prep_pdb.close()
#print(res_n_list)
#print("HERE IS THE ORIGINAL O COORDS" + str(mod_coords[-1]))

""" Translation of unit to 5'P coordinates"""
mod_coords = np.asarray(mod_coords)
delta_P = np.asarray(P_coords) - np.asarray(mod_coords[0])

for i in range(len(mod_coords)):
	mod_coords[i] = mod_coords[i] + delta_P 
print("Translating unit to 5'P...")



""" Reference vectors for PO rotation matrix calculation"""

P_coords = np.asarray(mod_coords[0])
O_res_coords = np.asarray(mod_coords[4])### EDIT THIS FOR PYRED PDB 
O_sys_coords = np.asarray(O_coords)
v_PO_res = P_coords - O_res_coords
v_PO_sys = P_coords - O_sys_coords

print("Rotating unit about 5'P to 3'O...")

""" Cross product method """

def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

rot_vector= rotation_matrix_from_vectors(v_PO_res,v_PO_sys)

def apply_rotation(p,rot_matrix):
	new_p = np.asarray(p)
	new_p = rot_matrix.dot(new_p)
	return new_p

for i in range(len(mod_coords)):
        mod_coords[i] = apply_rotation(mod_coords[i],rot_vector)
delta_P = np.asarray(original_P) - np.asarray(mod_coords[0])

for i in range(len(mod_coords)):
        mod_coords[i] = mod_coords[i] + delta_P 

""" Rotate about P-O axis to align other reference atoms """

def makeUnit(x):
    """Normalize entire input to norm 1. Not what you want for 2D arrays!"""
    return x / linalg.norm(x)

def xParV(x, v):
    """Project x onto v. Result will be parallel to v."""
    # (x' * v / norm(v)) * v / norm(v)
    # = (x' * v) * v / norm(v)^2
    # = (x' * v) * v / (v' * v)
    return np.dot(x, v) / np.dot(v, v) * v

def xPerpV(x, v):
    """Component of x orthogonal to v. Result is perpendicular to v."""
    return x - xParV(x, v)

def xProjectV(x, v):
    """Project x onto v, returning parallel and perpendicular components
    >> d = xProject(x, v)
    >> np.allclose(d['par'] + d['perp'], x)
    True
    """
    par = xParV(x, v)
    perp = x - par
    return {'par': par, 'perp': perp}

def rotateAbout(a, b, theta):
    """Rotate vector a about vector b by theta radians."""
    # Thanks user MNKY at http://math.stackexchange.com/a/1432182/81266
    proj = xProjectV(a, b)
    w = np.cross(b, proj['perp'])
    return (proj['par'] +
            proj['perp'] * np.cos(theta) +
            linalg.norm(proj['perp']) * makeUnit(w) * np.sin(theta))

def centroid(coordinate_set=[]):
	"""Calculates the centroid given a set of cartesian coordinates
	:param coordinate_set: An array of 3D coordinates
	:return centroid: Centroid of the coordinate set of input array
	"""
	x_set = []
	y_set = []
	z_set = []	
	for i in range(len(coordinate_set)):
		x_set.append(coordinate_set[i][0])
		y_set.append(coordinate_set[i][1])
		z_set.append(coordinate_set[i][2])	
	x_av = np.average(x_set)
	y_av = np.average(y_set)
	z_av = np.average(z_set)
	centroid = [x_av,y_av,z_av]
	return centroid

delta_P = np.asarray(original_P) - np.asarray(mod_coords[0])

for i in range(len(mod_coords)):
        mod_coords[i] = mod_coords[i] + delta_P

#print("Centroid test: " + str(centroid(after_mut_coords)))

def write_new_test_coords(pose_coords,name):
	"""Writes new PDB block of coordinate set relative to the original PDB
	:param pose_coords: Final translated and rotate coordinates of residue
	:param name: Name of file to be outputted
	:return none:
	"""
	print("Writing new PDB...")
	prep_pdb = open(mod_nuc_dir+ mut_sel+".pdb")
	new_pdb = open(name+".pdb",'w')
	i=0
	for line in prep_pdb:
		if len(line.split()) >= 5:
			line=line.split()
			new_pdb.write('{:>4}'.format(str(line[0])))
			new_pdb.write('{:>7}'.format(str(line[1])))
			new_pdb.write('{:<4}'.format(" " + str(line[2])))
			new_pdb.write('{:>5}'.format(str(line[3])))
			new_pdb.write('{:>6}'.format(str(mut_site)))

			new_pdb.write('{:>8}'.format(str(pose_coords[i][0]).split('.')[0]) + ".")
			new_pdb.write('{:<4.3}'.format(str(pose_coords[i][0]).split('.')[1]))

			new_pdb.write('{:>3}'.format(str(pose_coords[i][1]).split('.')[0]) + ".")
			new_pdb.write('{:<4.3}'.format(str(pose_coords[i][1]).split('.')[1]))

			new_pdb.write('{:>3}'.format(str(pose_coords[i][2]).split('.')[0]) + ".")
			new_pdb.write('{:<4.3}'.format(str(pose_coords[i][2]).split('.')[1]))

			new_pdb.write('{:>5}'.format(str(line[8])))
			new_pdb.write('{:>6}'.format(str(line[9])))
#		new_pdb.write('{:>2}'.format(str(line[10])))
			new_pdb.write("\n")
			i = i+1

def rot_trans(rads):
	test_rot_coords = mod_coords
	for i in range(len(test_rot_coords)):
		test_rot_coords[i] = rotateAbout(test_rot_coords[i],v_PO_sys,rads)
	delta_P = np.asarray(original_P) - np.asarray(test_rot_coords[0])
	for i in range(len(test_rot_coords)):
		test_rot_coords[i] = test_rot_coords[i] + delta_P
	return test_rot_coords

sweep_coords = []

for phi in range(0,1080,1):
#	print(phi)
	rads = phi*(np.pi/180)/3
#	print(rads)
	sweep_coords.append(rot_trans(rads).tolist())

def intercent_dist(frame):
	after_mut_cent = centroid(after_mut_coords)
	prev_mut_cent = centroid(prev_mut_coords)
	res_cent = centroid(frame)
	distance1 = ((res_cent[0] - after_mut_cent[0]) ** 2
	             + (res_cent[1] - after_mut_cent[1]) ** 2
	             + (res_cent[2] - after_mut_cent[2]) ** 2) ** 0.5
	distance2 = ((res_cent[0] - prev_mut_cent[0]) ** 2
	             + (res_cent[1] - prev_mut_cent[1]) ** 2
	             + (res_cent[2] - prev_mut_cent[2]) ** 2) ** 0.5
	dist = distance1 + distance2
	return dist

def intercent_dist_N(frame):
	after_mut_cent = centroid(after_mut_N_coords)
	prev_mut_cent = centroid(prev_mut_N_coords)
	n_list = [frame[index] for index in res_n_list]
	res_cent = centroid(n_list)
	distance1 = ((res_cent[0] - after_mut_cent[0]) ** 2
	             + (res_cent[1] - after_mut_cent[1]) ** 2
	             + (res_cent[2] - after_mut_cent[2]) ** 2) ** 0.5
	distance2 = ((res_cent[0] - prev_mut_cent[0]) ** 2
	             + (res_cent[1] - prev_mut_cent[1]) ** 2
	             + (res_cent[2] - prev_mut_cent[2]) ** 2)  ** 0.5
	dist = distance1 + distance2
	eq_dist = np.absolute(distance1 - distance2)
	return eq_dist + dist
#print(len(sweep_coords))

def overlap_orig_dist_N(frame):
	after_mut_cent = centroid(orig_mut_N_coords)
	prev_mut_cent = centroid(orig_mut_N_coords)
	n_list = [frame[index] for index in res_n_list]
	res_cent = centroid(n_list)
	distance1 = ((res_cent[0] - after_mut_cent[0]) ** 2
			+ (res_cent[1] - after_mut_cent[1]) ** 2
			+ (res_cent[2] - after_mut_cent[2]) ** 2) ** 0.5
	distance2 = ((res_cent[0] - prev_mut_cent[0]) ** 2
			+ (res_cent[1] - prev_mut_cent[1]) ** 2
			+ (res_cent[2] - prev_mut_cent[2]) ** 2)  ** 0.5
	#dist = distance1 + distance2
	#eq_dist = np.absolute(distance1 - distance2)
	#return eq_dist + dist
	return distance1	

def find_min():
	min_val = 100.0
	min_int = 0
	print("Finding inter-centroid distance minimum...")
	for i in range(len(sweep_coords)):
#		if intercent_dist_N(sweep_coords[i]) <= min_val:
		if overlap_orig_dist_N(sweep_coords[i]) <= min_val:
			min_val = overlap_orig_dist_N(sweep_coords[i])
			print(min_val)
			min_int = i
	return(min_int)


def plot_dist():
	x = []
	y = []
	for i in range(len(sweep_coords)):
		x.append(intercent_dist(sweep_coords[i]))
		y.append(i)
	plt.plot(x,y)
	plt.show()

#print(find_min())

def write_new_pdb():
#	newpdb = open(str(mut_sel)+"_"+str(mut_site)+"_"+unmod_complex+".pdb",'w')
	newpdb = open(output_filename,'w')
	old_pdb = open(unmod_complex + ".pdb")
	for line in old_pdb:
		if len(line.split()) <= 3:
			newpdb.write(line)
		if len(line.split()) >= 3:
			if int(line.split()[4]) < mut_site:
				newpdb.write(line)
			if int(line.split()[4]) == mut_site:
				break
	for line in open("minpose.pdb",'r'):
		newpdb.write(line)
	old_pdb = open(unmod_complex + ".pdb")
	for line in old_pdb:
		#if len(line.split()) <= 3:
		#	newpdb.write(line)
		if len(line.split()) >= 3:
			if int(line.split()[4]) > mut_site:
				newpdb.write(line)
	newpdb.write('END')
	newpdb.close()
	old_pdb.close()
### conditionals for what type of placement
write_new_test_coords(sweep_coords[find_min()],"minpose")
write_new_pdb()
prep_pdb.close()
old_pdb.close()
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
