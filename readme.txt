This is a Matlab library for detecting and approximation decision lines/faults in 2D and 3D. We describe the underlying algorithm in
"Detecting and approximating decision boundaries in low dimensional spaces", which is availabe at arXiv.org.

This library requires the Matlab library mesh2D by Darren Engwirda (https://github.com/dengwirda/mesh2d).
There is no explicit documentation yet, but all functions are documented in the source code.

Organisation of faultapprox-matlab

├── src  : the actual Matlab implementation of the algorithm. The "main" function is faultApprox.m.
├── tests: ┐
│          ├─ generic: a number of 2D and 3D test cases, unrelated to the aforementioned paper
│          └─ GrajewskiKleefeld_2023: scripts of (some of) the tests presented in the paper
│
└── utils: utility Matlab functions