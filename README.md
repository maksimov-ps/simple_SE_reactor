# Simple SE Reactor Model Documentation

## Overview
This repository contains MATLAB scripts for simulating a simple SE (Sorption Enhanced) reactor model. The code is designed to solve a system of partial differential equations (PDEs) and ordinary differential equations (ODEs) to model the behavior of a chemical reactor. The model includes calculations for reaction kinetics, heat transfer, adsorption, and mass transfer.

## Repository Structure
The repository consists of the following files:

- **calculation_script.m**: The main script to run the simulation. This script initializes parameters, calls the solver, and processes results.
- **ODEsolver_Graaf.m**: Contains the main ODE solver function that calculates the time derivatives of the system's state variables.
- **dss020.m**: Implements a discretization scheme for the first derivative (used for convection terms).
- **dss044.m**: Implements a discretization scheme for the second derivative (used for diffusion terms).
- **efficiency.m**: Calculates the efficiency factors for reactions and diffusion coefficients to account for the mass transfer limitations on the rate of the reactions. Calculation of efficiency factors is based on https://doi.org/10.1016/S0009-2509(00)00194-9
- **equilibrium_calculation.m**: Computes equilibrium constants for the reactions.
- **heat_transfer.m**: Calculates heat transfer coefficients and thermal conductivity. Calculation is based on https://doi.org/10.1016/0009-2509(77)80143-7
- **rate_calculation_Graaf.m**: Defines the reaction rate equations based on the Graaf et al. model. Rate equations are based on: https://doi.org/10.1016/0009-2509(88)85127-3
- **RKS.m**: Computes fugacity coefficients using the Redlich-Kwong-Soave equation of state.
- **viscosity.m**: Calculates the viscosity and thermal conductivity of the gas mixture.
- **PSA_script.m** (optional): Implements the Pressure-Swing-Absorption (PSA) logic for the SE reactor model. This script simulates the cyclic operation of the reactor, including adsorption and desorption phases.
- **SERP_PDE.m** (optional): Solves the system of partial differential equations (PDEs) for the SE reactor model, providing a detailed spatial and temporal resolution of the reactor's behavior.

## Key Components

### 1. **ODEsolver_Graaf.m**
This is the core function of the repository. It calculates the time derivatives of the state variables, which include:
- Molar fractions of components (CO2, CO, H2, CH3OH, H2O, N2)
- Pressure
- Temperature
- Adsorbed water concentration

#### Key Steps:
1. **Initialization**: Extracts input parameters and initializes variables.
2. **Momentum Balance**: Calculates the velocity profile across the reactor.
3. **Energy and Mass Balances**: Discretizes the governing equations for pressure, temperature, and molar fractions.
4. **Reaction Rates**: Computes reaction rates using the Graaf et al. model.
5. **Heat Transfer**: Calculates heat transfer coefficients and energy balance.
6. **Adsorption**: Models water adsorption on the solid electrolyte.
7. **Output**: Returns the time derivatives of all state variables.

### 2. **Discretization Functions**
- **dss020.m**: Approximates the first derivative for convection terms.
- **dss044.m**: Approximates the second derivative for diffusion terms.

### 3. **Reaction Kinetics**
The reaction rates are based on the Graaf et al. model (1988) and include:
- CO + 2H2 → CH3OH
- CO2 + H2 → CO + H2O
- CO2 + 3H2 → CH3OH + H2O

### 4. **Heat Transfer**
The `heat_transfer.m` script calculates:
- Wall-to-bed heat transfer coefficient
- Effective thermal conductivity

### 5. **Adsorption**
The adsorption model includes:
- Adsorbent saturation capacity
- Mass transfer coefficients
- Adsorption rate equations

### 6. **Physical Properties**
- **viscosity.m**: Calculates gas mixture viscosity and thermal conductivity.
- **RKS.m**: Computes fugacity coefficients using the Redlich-Kwong-Soave equation of state.

### 7. **Optional Full Cycle Simulation**
- **PSA_script.m**: Simulates the full cycle of the SE reactor model using the Pressure-Swing-Absorption (PSA) logic. This includes adsorption and desorption phases to enhance the reactor's performance.
- **SERP_PDE.m**: Provides a detailed spatial and temporal resolution of the reactor's behavior by solving the system of PDEs for the SE reactor model.
- **equlibria.m**: Calculates equilibrium constnats for the PSA_script.m
- **jacobian_SERP.m**: Auxiliary file for numerical estimation of Jacobian for the ODE solver for PSA_script.m

## How to Use

1. **Set Up Parameters**:
   - Define the reactor dimensions, operating conditions, and initial values in `calculation_script.m`.
   - Adjust global parameters in the `ODEsolver_Graaf.m` file if needed.

2. **Run the Simulation**:
   - Execute `calculation_script.m` in MATLAB.
   - The script will call the ODE solver and other functions to compute the results.

3. **Analyze Results**:
   - The output includes profiles of molar fractions, pressure, temperature, and other variables along the reactor length.

## Dependencies
This code requires MATLAB and the following toolboxes:
- Optimization Toolbox (optional, for parameter fitting)

## Contact
For questions or issues, please contact Pavel Maksimov (maksimov.ps95@gmail.com).
