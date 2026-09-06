
### Baseline (before any optimisation)

x86_64, 2.4.6 numpy

```
case      wall [s]  assembly  factorise   solve   other  note
steady        9.10      2.35       6.08    0.00    0.67  6 Newton it.
theta        73.16     43.44      22.40    0.00    7.32  20 steps
rk45         43.77     37.48       1.82    0.00    4.46  20 steps
rans         33.40      5.31      19.62    0.00    8.47  2 outer it.
```

### Phase 1: residual-only evaluations, fixed sparsity pattern (NumPy, SuperLU)

x86_64, 2.4.6 numpy

```
case      wall [s]  assembly  factorise   solve   other  note
steady        9.07      1.73       6.23    0.00    1.11  6 Newton it.
theta        32.38     10.57      20.24    0.00    1.57  20 steps
rk45          8.45      4.55       1.63    0.00    2.27  20 steps
rans         27.38      3.62      16.07    0.00    7.69  2 outer it.
```
