#pragma once

#include <vector>
#include <complex>
#include <string>
#include <array>

namespace neutrino {

struct SimulationResult {
    std::vector<double> r_vals;              // 距離
    std::vector<std::vector<double>> probs;  // 概率 [距離點][中微子flavor]
};

/**
 * 模擬中微子振盪
 * @param E_nu 中微子能量 (MeV)
 * @param solar_density_profile_path 太陽電子密度檔案路徑
 * @param r_start 起始距離
 * @param r_end 結束距離
 * @param num_points 距離點數量
 * @return 模擬結果
 */

// ✅ 舊版
SimulationResult simulate_oscillation(
    double E_nu,
    const std::string& solar_density_profile_path,
    double r_start,
    double r_end,
    int num_points
);

// ✅ 新版：r_vals 完全由 Python 傳入
SimulationResult simulate_custom_rvals(
    double E_nu,
    const std::vector<double>& input_r_vals,
    const std::vector<double>& density_r_vals,
    const std::vector<double>& Ne_vals);

// 將這段代碼添加到simulator.hpp的namespace neutrino內部

// 電子密度與哈密頓量相關函數
void compute_hamiltonian(std::array<std::complex<double>, 9>& H_standard, double E_nu, double Ve);
void compute_scaled_hamiltonian(std::array<std::complex<double>, 9>& H_scaled, double E_nu, double Ve, double scale_factor);
double estimate_density_gradient(double r, const std::vector<double>& r_vals, const std::vector<double>& Ne_vals);
double interpolate_resonance_position(double r1, double r2, double Ve1, double Ve2, double res_val);

// 共振混合矩陣相關函數
void compute_resonance_mixing_12(std::array<std::complex<double>, 9>& U, std::array<std::complex<double>, 9>& U_dag);
void compute_resonance_mixing_13(std::array<std::complex<double>, 9>& U, std::array<std::complex<double>, 9>& U_dag);

// 共振跳躍相關函數
void transform_to_resonance_eigenbasis_12(const std::array<std::complex<double>, 9>& rho, std::array<std::complex<double>, 9>& rho_eigen);
void apply_jump_probability_12(std::array<std::complex<double>, 9>& rho_eigen, double P_jump);
void transform_back_from_resonance_eigenbasis_12(const std::array<std::complex<double>, 9>& rho_eigen, std::array<std::complex<double>, 9>& rho);

void transform_to_resonance_eigenbasis_13(const std::array<std::complex<double>, 9>& rho, std::array<std::complex<double>, 9>& rho_eigen);
void apply_jump_probability_13(std::array<std::complex<double>, 9>& rho_eigen, double P_jump);
void transform_back_from_resonance_eigenbasis_13(const std::array<std::complex<double>, 9>& rho_eigen, std::array<std::complex<double>, 9>& rho);

// 數值積分相關函數
void normalize_density_matrix(std::array<std::complex<double>, 9>& rho);
bool check_numerical_stability(const std::array<std::complex<double>, 9>& psi);
void ensure_hermiticity(std::array<std::complex<double>, 9>& rho);
void rk4_step(double t, double h, std::array<std::complex<double>, 9>& y, double E_nu, const std::vector<double>& r_vals, const std::vector<double>& Ne_vals);
void enhanced_rk4_step(double& t, double& h, std::array<std::complex<double>, 9>& rho, double E_nu, const std::vector<double>& r_vals, const std::vector<double>& Ne_vals);
void adaptive_rk4_step(double& t, double& h, std::array<std::complex<double>, 9>& y, double E_nu, const std::vector<double>& r_vals, const std::vector<double>& Ne_vals, double tolerance = 1.0e-10);
void scaled_rk4_step(double& t, double h, std::array<std::complex<double>, 9>& rho, double E_nu, const std::vector<double>& r_vals, const std::vector<double>& Ne_vals, double scale_factor);
void ultra_high_energy_rk4_step(double t, double h, double rho, double E_nu, const std::vector<double>& r_vals, const std::vector<double>& Ne_vals);
void integrate_to_position(double& t, double next_r, double& h, std::array<std::complex<double>, 9>& rho, double E_nu, const std::vector<double>& r_vals, const std::vector<double>& Ne_vals);

// 還需要添加這個函數聲明，它在simulator.cpp中被調用，但未聲明
void neutrino_evolution(double t, const std::array<std::complex<double>, 9>& y, double E_nu, const std::vector<double>& r_vals, const std::vector<double>& Ne_vals, std::array<std::complex<double>, 9>& dydt);

// 共振跳躍函數
void apply_resonance_jump_12(std::array<std::complex<double>, 9>& rho, double P_jump);
void apply_resonance_jump_13(std::array<std::complex<double>, 9>& rho, double P_jump);

} // namespace neutrino