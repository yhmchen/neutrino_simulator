#include "simulator.hpp"
#include "constants.hpp"
#include <cmath>
#include <iostream>

NeutrinoSimulator::NeutrinoSimulator(const std::string& densityFile)
    : electronDensity(densityFile), pmns(nullptr), hamiltonian(nullptr) {
    
    // Initialize with default parameters
    setParameters(Constants::theta12, Constants::theta13, Constants::theta23,
                 Constants::deltaCP, Constants::deltaM21_2, Constants::deltaM31_2);
}

void NeutrinoSimulator::setParameters(double theta12, double theta13, double theta23, 
                                     double deltaCP, double dm21, double dm31) {
    // Delete old objects if they exist
    if (pmns) delete pmns;
    if (hamiltonian) delete hamiltonian;
    
    // Create new objects with updated parameters
    pmns = new PMNSMatrix(theta12, theta13, theta23, deltaCP, dm21, dm31);
    hamiltonian = new Hamiltonian(*pmns, electronDensity);
}

Matrix3c NeutrinoSimulator::evolveState(double energy, double pathStart, double pathEnd, int steps) {
    using namespace MatrixUtils;
    
    // Initialize evolution operator as identity
    Matrix3c U = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    
    // Calculate step size
    double step_size = (pathEnd - pathStart) / steps;
    
    // Evolve state step by step
    for (int i = 0; i < steps; i++) {
        double r = pathStart + (i + 0.5) * step_size; // Use midpoint
        
        // Calculate Hamiltonian at this position
        Matrix3c H = hamiltonian->calculate(r, energy);
        
        // Evolve for step_size: U_step = exp(-i * H * step_size / ħ)
        double conversion = Constants::hbar * Constants::c / Constants::solarRadius;
        Complex factor(0.0, -step_size / conversion);
        
        Matrix3c H_step = multiply(factor, H);
        Matrix3c U_step = expm(H_step);
        
        // Update evolution operator
        U = multiply(U_step, U);
    }
    
    return U;
}

double NeutrinoSimulator::calculateProbability(int alpha, int beta, double energy, 
                                             double pathStart, double pathEnd) {
    // If calculating vacuum oscillation over long distance (Sun to Earth),
    // use the MSW formula for better accuracy
    if (pathStart >= 1.0 && pathEnd >= Constants::sunEarthDistance / Constants::solarRadius) {
        return mswProbability(alpha, beta, energy, pathEnd - pathStart);
    }
    
    // Calculate evolution operator
    int steps = 1000; // Number of steps for numerical integration
    Matrix3c U = evolveState(energy, pathStart, pathEnd, steps);
    
    // Calculate probability: P(α→β) = |U_βα|²
    Complex amplitude = U[beta][alpha];
    double probability = std::norm(amplitude);
    
    return probability;
}

std::vector<std::vector<double>> NeutrinoSimulator::calculateAllProbabilities(
    double energy, double pathStart, double pathEnd) {
    
    std::vector<std::vector<double>> probs(3, std::vector<double>(3));
    for (int alpha = 0; alpha < 3; alpha++) {
        for (int beta = 0; beta < 3; beta++) {
            probs[alpha][beta] = calculateProbability(alpha, beta, energy, pathStart, pathEnd);
        }
    }
    
    return probs;
}

std::vector<std::vector<std::vector<double>>> NeutrinoSimulator::simulateSunToEarth(
    double energy, const std::vector<double>& radii) {
    
    std::vector<std::vector<std::vector<double>>> results;
    
    // Start from Sun center (r = 0)
    double pathStart = 0.0;
    
    // Calculate probabilities at each specified radius
    for (double r : radii) {
        auto probs = calculateAllProbabilities(energy, pathStart, r);
        results.push_back(probs);
    }
    
    return results;
}

double NeutrinoSimulator::mswProbability(int alpha, int beta, double energy, double L) {
    // This is an implementation based on the formulas in your image
    // Only implemented for electron neutrino survival probability P(νe→νe)
    if (alpha == 0 && beta == 0) {
        // Convert energy to eV
        double E = energy * Constants::MeV_to_eV;
        
        // MSW parameters from PMNS matrix
        double theta12 = Constants::theta12;
        double dm21_2 = Constants::deltaM21_2;
        
        // Calculate oscillation phase
        double delta_m = dm21_2 / (2.0 * E);
        double L_meters = L * Constants::solarRadius;
        double phase = delta_m * L_meters;
        
        // Calculate survival probability using the formula from your image
        double sin2_2theta = std::pow(std::sin(2.0 * theta12), 2);
        double sin2_phase = std::pow(std::sin(phase), 2);
        
        double probability = 1.0 - sin2_2theta * sin2_phase;
        return probability;
    }
    else if (alpha == 0 && beta == 1) {
        // P(νe→νμ) calculation
        double E = energy * Constants::MeV_to_eV;
        double theta12 = Constants::theta12;
        double theta23 = Constants::theta23;
        double dm21_2 = Constants::deltaM21_2;
        
        double delta_m = dm21_2 / (2.0 * E);
        double L_meters = L * Constants::solarRadius;
        double phase = delta_m * L_meters;
        
        double sin2_2theta12 = std::pow(std::sin(2.0 * theta12), 2);
        double cos2_theta23 = std::pow(std::cos(theta23), 2);
        double sin2_phase = std::pow(std::sin(phase), 2);
        
        double probability = sin2_2theta12 * cos2_theta23 * sin2_phase;
        return probability;
    }
    else if (alpha == 0 && beta == 2) {
        // P(νe→ντ) calculation
        double E = energy * Constants::MeV_to_eV;
        double theta12 = Constants::theta12;
        double theta23 = Constants::theta23;
        double dm21_2 = Constants::deltaM21_2;
        
        double delta_m = dm21_2 / (2.0 * E);
        double L_meters = L * Constants::solarRadius;
        double phase = delta_m * L_meters;
        
        double sin2_2theta12 = std::pow(std::sin(2.0 * theta12), 2);
        double sin2_theta23 = std::pow(std::sin(theta23), 2);
        double sin2_phase = std::pow(std::sin(phase), 2);
        
        double probability = sin2_2theta12 * sin2_theta23 * sin2_phase;
        return probability;
    }
    
    // For other transitions, use numerical calculation
    return calculateProbability(alpha, beta, energy, 1.0, L);
}