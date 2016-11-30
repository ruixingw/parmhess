=========
Structure
=========

Prepare Input file
------------------

Simply starting from a model that contains the parameters you are going to parameterize, following processes must be done.

1. Geometry Optimization by Gaussian
2. Frequency Calculation to obtain Hessian by Gaussian
3. MK Charge Calculation to obtain potential fitting data by Gaussian
4. RESP Charge Calculation by AmberTools
5. Identify atom types by AmberTools


A program **Tsubasa** is made to automated these process. **Tsubasa** accepts a geometry input and a config file, then do these  steps in sequence.

Although, for large or complex system, it is recommended that you should do optimization and frequency calculation manually, and also charge calculation if necessary. These outputs can be provided for the following steps.

Provided a "H2O2.gau" and "H2O2.cfg", **Tsubasa** will produce several things:

1. mmH2O2.com     # MM input file. Includes the following information: geometry, atom type, charge, MM functions 
2. freqH2O2.fchk  # Frequency FCHK file. Hessian is read from this file.
3. freqH2O2.log   # Frequency LOG file. 
4. input.inp      # All-in-one input file.
5. Tsubasa        # Folder in which stores the intermediate files.

We are now ready for parameterization by Parmhess.

Parameterization
----------------

It is very common that several bonds/angles/dihedrals/impropers have the same type and share a same MM function. For example, two H-O bonds in H2O2 have the same bond type of "ho-oh". Parameterization process will determine force constants for both of them separately, and then average them to yield the force constant of "ho-oh" type.

For convinience, in the following discussion, term "internal coordinates" refers to all possible bonds/angles/dihedrals/impropers in the molecule. 

1. Create and switch to "hffiles" folder as the working space and copy all input files here.
2. From mmH2O2.com, read Geometry, atom types, charges and MM functions. From these informations, internal coordinates are identified and the corresponding objects will be created in the memory.
   In this example, the following MM functions are readed: ho-oh, oh-oh, ho-oh-oh, ho-oh-oh-ho. The following internal coordinates are identified and the objects are created: h1-o1; o1-o2; o2-h2; h1-o1-o2; o1-o2-h2; h1-o1-o2-h2.
3. Identify and count unknown force constants. If there are "improper" functions, the corresponding improper objects will be created. Here, all MM functions are unknown and the number of unknowns is 4.
4. Match all internal coordinates with MM functions. For example, o1-o2 will match with oh-oh MM function. 
5. Identify and Count all internal coordinates whose force constant is unknown.
6. Prepare and submit Gaussian jobs named as "hessXX.com". Each file corresponds to a internal coordinate. They are used to calculate the "lowercase h" in the formula.
7. Identify coupled terms, which are involved in 6- and 5-member member rings.
8. 1-4 pairs of dihedrals and impropers, 1-3 pairs of angles, 1-2 pairs of bonds are identified and stored.
9. All internal coordinates are classified into 5 groups based on their "pairs":

   a. one-four group: This group contains all "pure" dihedrals that are not coupled with others.
   b. one-tri-four group: This group contains dihedrals that are coupled with angles. This is the 5-member ring case, in which 1-4 atom pair of a dihedral is coupled with 1-3 atom pair of an angle.
   c. one-tri-coupled group: This group contains angles that are coupled with impropers. The 1-4 atom pair of impropers is always coupled with 1-3 atom pair of an angle.
   d. one-tri-uncoupled group: This group contains all "pure" angles that are not coupld with others.
   e. one-two group: This group contains all "pure" bonds.

10. Prepare the initial Hprime file, which yield the Hessian produced by known terms and nonbonded terms. These contribution will then be subtracted in the following processes.
11. Start solving equation systems for each group in sequence. After the end of each group, force constants of that group are solved. These solved parameters are then treated as "known terms" in the new Hprime file. Hence, the contribution of former groups will be excluded in the following processes.
12. Now, all force constants for each internal coordinate are solved. These results will then be averaged based on MM functions.



