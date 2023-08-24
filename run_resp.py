import psi4
import resp
import numpy as np

mol = psi4.geometry("""
O       1.0990     0.7800     0.7260
F      -2.3260     0.0910    -0.0950
O      -0.7320    -1.6160     1.3450
O       3.3640    -0.1730    -0.6070
C      -0.9310     0.1490    -0.3050
C      -0.2060    -1.1140     0.1170
C       1.2150    -0.6080     0.3410
C      -0.2670     1.2030     0.5660
C       2.0740    -0.6710    -0.9150
C      -0.3010     2.5930    -0.0450
H      -0.7120     0.3670    -1.3600
H      -0.2770    -1.9180    -0.6220
H       1.7180    -1.1240     1.1670
H      -0.7280     1.2390     1.5600
H       1.6620    -0.0500    -1.7160
H       2.1710    -1.7020    -1.2700
H       0.2240     2.6140    -1.0060
H       0.2140     3.3060     0.6070
H      -1.3290     2.9360    -0.1970
H      -0.1990    -2.3900     1.5930
H       3.8940    -0.2280    -1.4210
""")

# N.B. Atom numbering is 1-based!
backbone_indices = [9,15,16,7,13,1,6,12,5,11,2,]
backbone_charge = 0.2933

# I think most of these are the default values anyway, but i listed them
# just so it's explicit and clear which values are being used.
options = {'vdw_scale_factors' : [1.4, 1.6, 1.8, 2.0],
           'vdw_point_density' : 1.0,
           'resp_a'            : 0.0005,
           'resp_b'            : 0.1,
           'constraint_charge' : [[backbone_charge, backbone_indices]],
           }

# Call for first stage fit, which applies the standard parabolic restraints
# and the backbone charge constraints
unconstrained_charges1, constrained_charges1 = resp.resp([mol], options)

# Automagically determine aliphatic hydrogens and constrain them to be equivalent
resp.set_stage2_constraint(mol, constrained_charges1, options)
# The above call resets the options and replaces them with constraints that are designed
# to preserve 
options['constraint_charge'].append([backbone_charge, backbone_indices])

# Run the fit again with the new group / charge constraints introduced
constrained_charges, constrained_charges = resp.resp([mol], options)

print("Results")
print("=======")
print(f"\nconstrained charges:\n{constrained_charges}")
print(f"\nSum of constrained charges: {np.sum(constrained_charges)}")
print(f"\nAutomatically determined equivalence groups: {options['constraint_group']}")
for equiv in options['constraint_group']:
    if len(equiv) > 1:
        print(f"These should be the same: {constrained_charges[[i-1 for i in equiv]]}")
print(f"Sum of backbone atom charges: {np.sum(constrained_charges[[i-1 for i in backbone_indices]])}")
