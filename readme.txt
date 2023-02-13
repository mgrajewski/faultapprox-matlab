This is a Matlab library for detecting and approximation decision lines/faults in 2D and 3D.
We describe the underlying algorithm in "Detecting and approximating decision boundaries in low dimensional spaces", which is availabe at arXiv.org.

This library requires the Matlab library mesh2D by Darren Engwirda (https://github.com/dengwirda/mesh2d).
There is no explicit documentation yet, but all functions are documented in the source code.

Organisation of faultapprox-matlab

├── src  : the actual Matlab implementation of the algorithm. The "main" function is faultApprox.m.
├── tests: ┐
│          ├─ generic: a number of 2D and 3D test cases, unrelated to the aforementioned paper
│          └─ GrajewskiKleefeld_2023: scripts of (some of) the tests presented in the paper
│
└── utils: utility Matlab functions


Running the test cases
──────────────────────

- For running the tests, make sure that the faultapprox-matlab-folder with its subfolders is added to the Matlab path.
- For running the tests in tests/generic, just execute the files testCaseFaultApprox[n]D_[test index]. The resulting vtu-files are stored in the subfolder tests/generic/results.
- For running the test cases presented in "Detecting and approximating decision boundaries in low   dimensional spaces", section 2.1, 2.2 and 4.1, run the Matlab files testCase2D_[test index].m, for   the 3D test case of section 2.2, run testFunc3D.m in tests/GrajewskiKleefeld2023.
- For the reconstruction of scatters using Kirsch's factorization method (see "Detecting and   approximating decision boundaries in low dimensional spaces", section 4.3, first unzip the files   in farFiledData1026, which contain the farfield data. Then, run "reconstruct[scatterer].m". This script creates two txt-files in the subfolder raw_results, "normal[scatterer].m" and "points[scatterer].m" which contain the points and the outer normals of the surface at these points. For visualisation purposes, we created a mesh using open3d (http://www.open3d.org/). A template file
for this is available ("triangulate_points.py") in the directory tests/GrajewskiKleefeld_2023. 
