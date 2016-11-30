================
How's it working
================

Prepare Input file
------------------

Simply starting from a model, following further processes must be done.

1. Geometry Optimization by Gaussian
2. Frequency Calculation by Gaussian to obtain Hessian 
3. MK Charge Calculation by Gaussian to obtain potential fitting data 
4. RESP Charge Calculation by AmberTools to obtain atom charge
5. Identify atom types by AmberTools


A program **Tsubasa** is made to automated these processes. **Tsubasa** accepts a geometry input and a config file, then do these  steps in sequence. Although, for large or complex system, it is recommended that you should do optimization and frequency calculation manually (also charge calculation if necessary). These outputs can be provided for the following steps.

Suppose we start from a "H2O2.gau" and "H2O2.cfg", **Tsubasa** will produce these files:

- :code:`mmH2O2.com`     : MM input file. Includes the following information: geometry, atom type, charge, types of bonded terms
- :code:`freqH2O2.fchk`  :Frequency FCHK file. Hessian is read from this file.
- :code:`freqH2O2.log`   : Frequency log file   
- :code:`input.inp`      : All-in-one input file.
- :code:`Tsubasa`        : Folder that stores intermediate files.

Now it is ready for parameterization by Parmhess.

Parameterization
----------------

It is very common that several bonds/angles/dihedrals/impropers have the same type. For example, two H-O bonds in H2O2 have the same bond type "ho-oh". Parmhess will solve force constant for each of them separately, and then average force constants of both to yield the force constant for "ho-oh" type.

For convinience, in the following text, term **internal coordinates** refers to all bonds/angles/dihedrals/impropers in the system. 

\1. Create and switch to :code:`hffiles` folder as the working space and then copy all input files here.


\2. From the Input files (:code:`mmH2O2.com`/:code:`freqH2O2.fchk` in the previous example), read geometry, atom types, charges , Hessian, and internal coordinate types.
   From these informations, internal coordinates are identified.
   In this example, the following coordinate types are readed: ho-oh, oh-oh, ho-oh-oh, ho-oh-oh-ho.
   The following internal coordinates are identified: h1-o1; o1-o2; o2-h2; h1-o1-o2; o1-o2-h2; h1-o1-o2-h2.


\3. Identify and count unknown coordinate types. If a **improper** type exists, the corresponding improper coordinate will be identified, too.
   In H2O2, all types are unknown and the number of unknowns is 4.


\4. Identify and match the type of internal coordinates.
   Here, The relation between o1-o2 and oh-oh type will be identified. 


\5. Identify and Count all internal coordinates whose force constant is unknown.
   Here, all internal coordinates are unknown and the number of unknowns is 6.


\6. Prepare and run Gaussian jobs named as :code:`hessXX.com`. Each job corresponds to a internal coordinate. They are used to calculate the **lowercase h**.


\7. 1-4 pair of dihedrals and impropers, 1-3 pair of angles, 1-2 pair of bonds are identified and stored.
   These pairs are corresponding to 3x3 partial Hessians.


\8. Identify cross-coupled pairs that are involved in 6- and 5-member rings.


\9. All pairs are classified into 5 groups based on the coupling situation:
   \a. one-four group:
      This group contains all 1-4 pairs that are not coupled with 1-3 pairs. This means, the contribution to the 3x3 partial-Hessian is only from dehedral terms and nonbonded terms.
      Please note if an 1-4 pair couples with a 1-2 pair, it will be a 4-member ring case, which is not supported by PHF.


   \b. one-tri-four group:
      This group contains 1-4 pairs that are coupled with 1-3 pairs. This is the 5-member ring case, in which the 1-4 pair of a dihedral is also the 1-3 pair of an angle. Hence, the 3x3 partial Hessian is contributed both from a dihedral and an angle.
      Since other dihedrals may also have contribution to this partial Hessian, this group is placed after group a) to subtract the contribution from group a).


   \c. one-tri-coupled group:
      This group contains angles that are coupled with impropers. The 1-4 atom pair of an improper is always also the 1-3 atom pair of an angle. Also, contributions from former groups should be subtracted.


   \d. one-tri-uncoupled group:
      This group contains all "pure" 1-3 pairs that are not coupld with others. After subtracting the former groups' contribution, the 3x3 partial Hessian is only contributed by a single angle.
      Please note if an 1-3 pair couples with a 1-3 or 1-2 pair, it will be a 4- or 3-member ring case, which is not supported by PHF.


   \e. one-two group:
      This group contains all "pure" bonds. After subtracting contribution from former groups, the 3x3 partial Hessian is purly contributed by a single bond.


\10. Prepare the initial *hprime* file to calculate **H'**, which is the Hessian contributed by *known terms* and nonbonded terms. Their contribution will then be subtracted in group a) process.


\11. Construct and solve equation systems for each group in sequence. After solving  force constants for each group, these solved parameters will then be treated as *known terms* and form a new :code:`hprime` file. Hence, the contribution of former groups will be excluded in the following processes.


\12. Now, all force constants for each internal coordinate are solved. These results will then be averaged based on coordinate types.



