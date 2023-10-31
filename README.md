# virus-research
This is the repository for the C++ program used to gather initial data for my summer 2023 research.
It tries to find a matrix equation of the form $`TB_0 = B_1`$ that can be solved, which represents an icosahedral viral transition that preserves some or all of icosahedral symmetry.

It is not as comprehensive or exhaustive as the Python version used in my other repo [https://github.com/RochX/virus-research-python](https://github.com/RochX/virus-research-python).
That version should be used instead of this one, but I leave this one public for posterity.

## Running the Program
### Eigen Library
To run the cpp program, note the use of the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library for fast matrix operations.

This library is simply a C++ header library, and ***is not included in the repo***, so to install it just download version 3.4.0 and put it in a folder named `eigen-3.4.0` (it should already come named this after extracting the zip).

### CMake
After installing Eigen run `cmake --build /path/to/build-dir` when inside the project directory to produce the make file.
Then run `make` within the build directory to compile the program.

It produces 3 executables:
- `virus_research` is the main program;
- `print_output_results` is a helper program to read the data created by `virus_research`
- `entry_sampling` is another helper program used to gather information that could make the main program run faster.

All programs take one optional argument that specifies a directory which all files will be assumed to be in.
If not given it assumes current directory.
