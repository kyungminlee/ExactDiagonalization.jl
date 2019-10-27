# ExactDiagonalization

Implements exact diagonalization.

Schematics for the structure of the package
```
                State
                  ↓
                Site
                  ↓
                HilbertSpace → HilbertSpaceSector    Operator
                  ↓              ↓                     ↓
                HilbertSpaceRepresentation         → OperatorRepresentation
                  ↓                                    ↓
SymmetryGroup → ReducedHilbertSpaceRepresentation  → ReducedOperatorRepresentation
```
