# Operator

## Operator Types

### NullOperator

A [`NullOperator`](@ref), as the name suggests, represents a null operator. It contains
no fields, and is thus a singleton.

### PureOperator

A [`PureOperator`](@ref) represents an operator of the following form:
```math
\hat{O} = \alpha \hat{P}_1 \otimes \hat{P}_2 \otimes \ldots \otimes \hat{P}_N
```
where $\alpha$ is a complex number, and $\hat{P}_i$ is either identity, or projection
``|rᵢ⟩⟨cᵢ|``.
It serves as a building block for all the operators used for the construction of
the representation of the operators.


Internally, `PureOperator` has fields `bitmask`, `bitrow`, `bitcol`, and `amplitude`.
The `bitmask` marks whether the $\hat{P}_i$ is identity or projection:
If the bitmask for site `i` is unset, then $\hat{P}_i$ is an identity operator;
if it is set, then $\hat{P}_{i}$ is a projection.
The fields `bitrow` and `bitcol` can contain information on rᵢ and cᵢ: they can
contain nonzero bit-field only at sites with nonzero `bitmask`.

### SumOperator

A [`SumOperator`](@ref) represents a sum of `PureOperator`. The scalar types of
the `PureOperator`s are required to be the same. While a `SumOperator` can be constructed
from the `PureOperators`, it can also be constructed using additions/subtractions
(See [Binary Operations](@ref binary-operations)).

## Mathematical Operations for Operators

### Unary Operations

Unary operations `+` and `-` are defined for the operators. These simply act on
the overall amplitude of the operators. `+` does not do anything and simply returns
the original operator, while `-` changes sign only. The type of the resulting operator,
is therefore the same as the original operator. There is one exception: when acting
`-` on an operator whose scalar type is `Bool`, the resulting type has scalar type `Int`.


In addition to `+` and `-`, functions `real` and `imag` are also defined for the
operators. Depending on the scalar type, the resulting operator has a different type:

|    .     |  N  | PR | PC | SR | SC |
|:--------:|:---:|:--:|:--:|:--:|:--:|
|  `real`  |  N  | PR | PR | SR | SR |
|  `imag`  |  N  | N  | PR | N  | SR |

- N: `NullOperator`
- PR: `PureOperator` with real scalar type
- PC: `PureOperator` with complex scalar type
- SR: `SumOperator` with real scalar type
- SC: `SumOperator` with complex scalar type

### [Binary Operations](@id binary-operations)

Binary operations are also defined for the operators. Since `PureOperator`s are
closed under multiplication, while product with `NullOperator` is always `NullOperator`,
we get the following multiplication table

| `*`    | N | P | S |
|:------:|:-:|:-:|:-:|
| N      | N | N | N |
| P      | N | P | S |
| S      | N | S | S |

- N: `NullOperator`, P: `PureOperator`, S: `SumOperator`

Additions or subtractions of two `PureOperator`s, on the other hand, produce `SumOperator`s

|`+`/`-` | N | P | S |
|:------:|:-:|:-:|:-:|
| N      | N | P | S |
| P      | P | S | S |
| S      | S | S | S |
