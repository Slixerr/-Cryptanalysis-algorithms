# -Cryptanalysis-algorithms

In this repository you can multiple algorithms developed for the Cryptology and Security Data class during my Master's Degree. Python is the programming language used, and Pandas is imported to generate the result tables.

You can find the full analysis of the algorithms in the PDF file [within this repository](https://github.com/Slixerr/-Cryptanalysis-algorithms/blob/main/TrabajoCriptoanalisis-SilviuManolescu.pdf) (in Spanish).

This assignment has been graded by the teacher, receiving a score of 10.

## Factorization algorithms

1. Fermat: method is based on the difference of squares. It aims to express an odd composite number as the difference of two squares.
2. Pollard-rho: probabilistic algorithm used for integer factorization. It employs Floyd's cycle-finding algorithm to find non-trivial factors.
3. Pollard P-1: another integer factorization method that focuses on finding factors by choosing random values and computing their powers modulo the number to be factored.
4. Lenstra: utilizes elliptic curves to find non-trivial factors of a composite number. It is particularly effective for numbers with small factors.

## Discrete Logarithm algorithms

1. Baby step Giant step: used for solving discrete logarithm problems. It employs a precomputation step to efficiently find the logarithm in a cyclic group.
2. Pollard-rho: algorithm can be adapted to find logarithms in cyclic groups, providing a probabilistic approach to solving the discrete logarithm problem
3. Pollard P-1: Similar to its application in factorization, Pollard's P-1 algorithm can be adapted for solving discrete logarithm problems, targeting cyclic groups with known order.
