#pragma once
#include "pmns_matrix.hpp"
#include "electron_density.hpp"

class Hamiltonian {
public:
    Hamiltonian(const PMNSMatrix& pmns, const ElectronDensity& elecDensity);
    
    // Calculate Hamiltonian at given position and energy
    Matrix3c calculate(double r, double energy) const;
    
    // Calculate vacuum oscillation Hamiltonian
    Matrix3c calculateVacuum(double energy) const;
    
    // Calculate matter potential matrix
    Matrix3c calculateMatterPotential(double r) const;
    
private:
    const PMNSMatrix& pmns;
    const ElectronDensity& electronDensity;
};