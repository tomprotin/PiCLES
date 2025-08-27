# Particle-in-Cell for Efficient Swell - PiCLES
PiCLES is a fast and efficient wave model for Earth System Models, using Particle-in-Cell methods for better performance. This version was developped to introduce a sequential importance resampling method to swell dynamics in PiCLES.

PiCLES:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13799205.svg)](https://doi.org/10.5281/zenodo.13799205)

## Quick Start
A brief guide on how to use PiCLES.

## How to Install

### Step 1: Install Julia
You will need Julia version 1.10.0 installed. For this, we recommend using juliaup
1. Open a terminal or command prompt.
21. Install juliaup using:
   ```
   curl -fsSL https://install.julialang.org | sh
   ```
3. Restart your shell or add Juliaup to your PATH if needed.
4. Install Julia version 1.10.0:
   ```
   juliaup add 1.10.0
   ```

### Step 2: Download PiCLES Repository

(Once registered as a Package this will be simpler, sry for the delay)

1. Open a terminal or command prompt.
2. Clone the PiCLES repository using Git:
   ```
   git clone https://github.com/tomprotin/PiCLES.git
   ```
3. Navigate into the cloned repository directory:
   ```
   cd PiCLES
   ```

### Step 3: Install Dependencies and Activate Environment
1. Start Julia version, 1.10.0 by typing `julia +1.10.0` in your terminal within the PiCLES directory where the `Project.toml` file is located.
2. Activate the project environment:
   ```julia
   using Pkg
   Pkg.activate(".")
   ```
3. Install the required dependencies:
   ```julia
   Pkg.instantiate()
   ```
   This may take some time.

You are now ready to reproduce the test cases for your simulations.

## How to Run Test cases
To run the two test cases from the command line, follow these steps:

1. Open a terminal or command prompt.
2. Navigate to the PiCLES directory where the `Project.toml` file is located.
3. Start Julia by typing `julia +1.10.0` in your terminal.
4. In the Julia REPL, include the `test_case_1.jl` file:
    ```julia
    include("tests/test_case_1.jl")
    ```
5. The test will run and save the results in a new repertory at `PiCLES/plots/test_case_1/`.
6. In the Julia REPL, include the `test_case_2.jl` file:
    ```julia
    include("tests/test_case_2.jl")
    ```
7. The test will run and save the results in `PiCLES/plots/test_case_2/`.

## How to Cite
