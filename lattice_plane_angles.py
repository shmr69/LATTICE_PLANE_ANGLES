'''
Author: P.Dirk

This script can be used to help with indexing of diffraction patters. 
Input parameters are controlled in the area marked by the ######'s at the top of the script.

----------------------
Inputs:
 - lattice parameters of candidate structure (obtained e.g. from powder diffraction refinement)
 - type of Bravais lattice (currently only primitive, body-centred, and face-centred are implemented)
 - measured angle between two diffraction spots and the ratio between the lengths of the corresponding vectors
 - some other control parameters, which can be used to fine-tune the output (default values should be okay for most cases)
 
 Output:
 - miller indices of candidate lattice planes
 - angle and ratio between the two planes (note that d1 is the spacing for (h1,k1,l1) and vice versa)
 - calculated viewing orthogonal direction
 - candidate combinations of planes are ranked based on the squared deviation from the input angle and ratio, weighted by the respective tolerance (the 'score') 
----------------------

Notes:
 - this code currently only accounts for systematic absences due to Bravais lattice, so there are additional space group-dependent reflection conditions (see ITC vol A, table 3.2)
 - for electron diffraction you minght have to consider double diffraction, so there might be spots that appear to break the reflection conditions.
   If a spot can't be indexed with I- or F-centering, it might be wirth trying P-centering (since that shows all forbidden reflections).
 - the code can calculate the lattice spacings for a candidate lattice plane, but these values might differ significantly from measured values (e.g. due to measurement calibration or deviations in lattice parameters)
   so it's usually more reliable to just compare the ratios and angles.

'''


#############INPUT PARAMETERS###########################
LPs = np.array([14.995400,5.478150,5.245150]) #input lattice parameters (e.g. from powder diffraction)

centering = 'P' #P: primitive, I: body-centred, F: face-centred

target_angle = 87.3785 #interplanar angle in degrees

target_ratio = 1.00263 #lattice plane spacings ratio


#control parameters:

max_hkl = (4,4,4) #maximum |h,k,l|-values for calculating lattice planes

print_spacings = True #option to output the interplanar spacings

#these parameters control the search radius around target values (those default values should be fine, but can be increased if no matching pairs of planes are found)
angle_tol = 1.0 #absolute tolerance in the angle (in degrees)

tol_rel = 0.05 #relative tolerance for ratio, or
tol_abs = 0.00 #absolute tolarance (whichever value is larger will automatically be chosen)

score_cutoff = 0.3 #maximum value for the chi-squared score to show in output
######################################



import numpy as np
import math
from itertools import permutations
import warnings
warnings.filterwarnings('ignore')



# GCD of more than two (or array) numbers
# This function implements the Euclidean 
# algorithm to find H.C.F. of two number
def find_gcd(x, y):
    while(y):
        x, y = y, x % y
 
    return x

def is_colliniear(v1,v2):
    cross = np.cross(v1,v2)
    if cross[0]==0 and cross[1]==0 and cross[2]==0:
        return True
    else:
        return False

def calc_score(alpha, target_alpha, alpha_err, r, target_r, r_err):
    return ( (alpha - target_alpha)**2/(alpha_err**2) + (r - target_r)**2/(r_err**2) )/2

def calculate_angle(lat,p1,p2):
    v1=np.nan_to_num(lat/np.array(p1),nan=0,posinf=0,neginf=0)
    v2=np.nan_to_num(lat/np.array(p2),nan=0,posinf=0,neginf=0)
    arg = np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
    return math.degrees(math.acos(arg))


def calculate_d_hkl(a,b,c,h,k,l):
    return 1/(h**2/a**2+k**2/b**2+l**2/c**2)**(0.5)


def get_indices(cell,LPs,h_max,k_max,l_max,cutoff_d=1.7):
    a,b,c = LPs
    ind = {}

    if cell == 'P': #any combination of hkl allowed
        for h in range(-1*h_max,h_max+1):
            for k in range(-1*k_max,k_max+1):
                for l in range (-1*l_max,l_max):
                    comb = permutations([h,k,l],3)
                    for i in comb:
                        if i != (0,0,0):
                            #print(i)
                            d = calculate_d_hkl(a,b,c,i[0],i[1],i[2])
                            if d >= cutoff_d:
                                ind.update({i:d})

    if cell == 'I': #only h+k+l=even reflections allowed
        for h in range(h_max+1):
            for k in range(k_max+1):
                for l in range (l_max+1):
                    #print(h,k,l)
                    if (h+k+l)%2 == 0:
                        comb = permutations([h,k,l],3)
                        for i in comb:
                            if i != (0,0,0):
                                #print(i)
                                d = calculate_d_hkl(a,b,c,i[0],i[1],i[2])
                                if d >= cutoff_d:
                                    ind.update({i:d})
                    else: continue

    if cell == 'F': #only h,k,l all even or all odd reflections allowed
        for h in range(-1*h_max,h_max+1):
            for k in range(-1*k_max,k_max+1):
                for l in range (-1*l_max,l_max+1):
                    #print(h,k,l)
                    if ((h%2) == 0 and (k%2) == 0 and (l%2) == 0) or ((h%2) != 0 and (k%2) != 0 and (l%2) != 0):
                        comb = permutations([h,k,l],3)
                        for i in comb:
                            if i != (0,0,0):
                                #print(i)
                                d = calculate_d_hkl(a,b,c,i[0],i[1],i[2])
                                if d >= cutoff_d:
                                    ind.update({i:d})
                    else: continue


    return ind



r_tol = max(tol_rel*target_ratio,tol_abs)

print('absolute tolerance:')
print(f'ratio: {r_tol:.2f}')
print(f'angle: {angle_tol} deg')
print()
print('(h1,k1,l1)\t(h2,k2,l2)\tangle\td1/d2\tview\t\tscore')
print('-------------------------------------------------------------------------')

planes = get_indices(centering,LPs,*max_hkl) 
miller  = list(planes.keys())
results = []




for i in range(len(miller)):
    for j in range(i+1,len(miller)):

        if is_colliniear(miller[i],miller[j]) == True:
            continue
        else:
            #print(miller[i],miller[j])
            angle = calculate_angle(LPs,miller[i],miller[j])
            ratio = planes[miller[i]]/planes[miller[j]]

        if (abs(angle-target_angle) <= angle_tol):

            if (abs(ratio-target_ratio) <= r_tol):
                dir = np.cross(miller[i],miller[j]) #find perpendicular direction
                gcd = find_gcd(dir[0],dir[1]) #rescale vector to smallest integer components
                gcd = find_gcd(gcd,dir[2])
                score = calc_score(angle,target_angle,angle_tol,ratio,target_ratio,r_tol)
                results.append([miller[i],miller[j],f'{angle:.2f}',f'{ratio:.4f}',dir/gcd,f'{score:.6f}'])

            elif (abs(1/ratio-target_ratio) <= r_tol):
                dir = np.cross(miller[j],miller[i]) #find perpendicular direction
                gcd = find_gcd(dir[0],dir[1]) #find GCD to rescale vector to smallest integer components
                gcd = find_gcd(gcd,dir[2])
                score = calc_score(angle,target_angle,angle_tol,1/ratio,target_ratio,r_tol)
                results.append([miller[j],miller[i],f'{angle:.2f}',f'{1/ratio:.4f}',dir/gcd,f'{score:.6f}'])

            else: continue
        else: continue


#visualise results
sorted_results = sorted(results,key=lambda x: x[5])
for result in sorted_results:
    if float(result[5]) <= score_cutoff:
        print(*result,sep='\t')
        if print_spacings == True:
            print(f'd1 = {planes[result[0]]:.4f} A\td2 = {planes[result[1]]:.4f} A')
