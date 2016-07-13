%mem=12gb
#p amber=softonly geom=connectivity nosymm
iop(4/33=3,7/33=1)
freq=noraman

final

0 1
C-ca--0.117738   -0.196066032299     0.083954610270     0.000060359368
C-ca--0.117738    1.202825076336     0.083884471006     0.000506018868
C-ca--0.117738    1.902331654047     1.295325035437    -0.000047136416
C-ca--0.117738    1.202946707718     2.506835148041    -0.001045596240
C-ca--0.117738   -0.195944409385     2.506905316939    -0.001491091778
C-ca--0.117738   -0.895450997150     1.295464738220    -0.000938299116
H-ha-0.117738   -0.739751345816    -0.857625934164     0.000489987063
H-ha-0.117738    1.746416421211    -0.857750015636     0.001281451623
H-ha-0.117738    2.989606822804     1.295270625436     0.000299297106
H-ha-0.117738    1.746632277886     3.448415533193    -0.001475037473
H-ha-0.117738   -0.739536124154     3.448539598789    -0.002266379146
H-ha-0.117738   -1.982726171199     1.295519560979    -0.001284573853

1 2 1.5 6 1.5 7 1.0
2 3 1.5 8 1.0
3 4 1.5 9 1.0
4 5 1.5 10 1.0
5 6 1.5 11 1.0
6 12 1.0
7
8
9
10
11
12

AmbTrs   *   ca  ca  *    0  180  0   0    0.000 14.494  0.000  0.000   4.0
HrmBnd1   ca  ca  ha   46.700  120.00000
HrmBnd1   ca  ca  ca   91.704  120.00000
HrmStr1   ca  ha   378.918  1.08728
HrmStr1   ca  ca   377.250  1.39889
Nonbon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 -1.2
VDW    ca  1.9080  0.0860
VDW    ha  1.4590  0.0150


