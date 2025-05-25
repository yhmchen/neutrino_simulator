#include "../include/simulator.hpp"
#include "../include/constants.hpp"
#include "../include/electron_density.hpp"
#include "../include/hamiltonian.hpp"
#include <vector>
#include <complex>
#include <array>
#include <cmath>
#include <iostream>
#include <fstream> // 加入 fstream 標頭檔

namespace neutrino {

// 簡化的 RK4 求解器
void rk4_step(double t, double h, std::array<std::complex<double>, 9>& y,
              double E_nu, const std::vector<double>& r_vals, 
              const std::vector<double>& Ne_vals) {

    std::cout << "Entering rk4_step() function at t = " << t << std::endl;  // 輸出
    
    std::array<std::complex<double>, 9> k1, k2, k3, k4, y_temp;
    
    // k1 = f(t, y)
    neutrino_evolution(t, y, E_nu, r_vals, Ne_vals, k1);
    std::cout << "k1[0] = " << k1[0] << std::endl;  // 輸出
    
    // k2 = f(t + h/2, y + h*k1/2)
    for (int i = 0; i < 9; i++) {
        y_temp[i] = y[i] + h * k1[i] / 2.0;
    }
    neutrino_evolution(t + h/2, y_temp, E_nu, r_vals, Ne_vals, k2);
    
    // k3 = f(t + h/2, y + h*k2/2)
    for (int i = 0; i < 9; i++) {
        y_temp[i] = y[i] + h * k2[i] / 2.0;
    }
    neutrino_evolution(t + h/2, y_temp, E_nu, r_vals, Ne_vals, k3);
    
    // k4 = f(t + h, y + h*k3)
    for (int i = 0; i < 9; i++) {
        y_temp[i] = y[i] + h * k3[i];
    }
    neutrino_evolution(t + h, y_temp, E_nu, r_vals, Ne_vals, k4);
    
    // y_{n+1} = y_n + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    for (int i = 0; i < 9; i++) {
        y[i] = y[i] + h/6.0 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
        // 檢查波函數是否包含 NaN
        if (std::isnan(y[i].real()) || std::isnan(y[i].imag())) {
            std::cerr << "Wave function psi contains NaN at t = " << t << ", index = " << i << std::endl;
            break;
        }
    }

    std::cout << "Exiting rk4_step() function at t = " << t << std::endl;  // 輸出
}

SimulationResult simulate_oscillation(
    double E_nu,
    const std::string& solar_density_profile_path,
    double r_start,
    double r_end,
    int num_points
) {
    std::ofstream outfile("hamiltonian_output.txt"); // 打開文件
    //outfile << "Entering simulate_oscillation() function with E_nu = " << E_nu << std::endl;
    //std::cout << "Entering simulate_oscillation() function with E_nu = " << E_nu << std::endl;  // 輸出
    SimulationResult result;
    
    // 讀取太陽電子密度數據
    std::vector<double> density_r_vals, Ne_vals;
    if (!solar_electron_density(solar_density_profile_path, density_r_vals, Ne_vals)) {
        std::cerr << "無法讀取電子密度數據，模擬終止。" << std::endl;
        outfile << "無法讀取電子密度數據，模擬終止。" << std::endl;
        std::cout << "Exiting simulate_oscillation() function (failed to read density)" << std::endl;
        outfile << "Exiting simulate_oscillation() function (failed to read density)" << std::endl;
        return result;
    }
    std::cout << "density_r_vals.size() = " << density_r_vals.size() << std::endl;  // 輸出
    std::cout << "Ne_vals.size() = " << Ne_vals.size() << std::endl;  // 輸出
    
    // 生成距離點
    result.r_vals.resize(num_points);
    double h = (r_end - r_start) / (num_points - 1);
    for (int i = 0; i < num_points; i++) {
        result.r_vals[i] = r_start + i * h;
    }
    
    // 初始條件 (單位矩陣)
    std::array<std::complex<double>, 9> psi;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            psi[i * 3 + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // 輸出初始條件
    //std::cout << "Initial psi:" << std::endl;
    //for (const auto& val : psi) {
        // std::cout << val << " ";
    //}
    //std::cout << std::endl;

    // 初始化結果容器
    result.probs.resize(num_points, std::vector<double>(3, 0.0));
    
    // 計算初始點的概率
    for (int j = 0; j < 3; j++) {
        result.probs[0][j] = std::norm(psi[j * 3]);  // |psi[j][0]|^2
    }
    
    // 數值積分
    double t = r_start;
    double step_size = h;
    
    // 如果距離範圍很大，可能需要調整步長
    if (r_end - r_start > 100.0) {
        // 使用較小的步長以保證太陽內部的精度
        step_size = std::min(1.0e-8, h/1000);
    }
    
    for (int i = 1; i < num_points; i++) {
        // 從上一個點到當前點進行多步積分
        double next_t = result.r_vals[i];
        while (t < next_t) {
            double current_step = std::min(step_size, next_t - t);
            rk4_step(t, current_step, psi, E_nu, density_r_vals, Ne_vals);
            t += current_step;
        }
        //outfile.close(); // 關閉文件

        // 計算當前點的概率
        for (int j = 0; j < 3; j++) {
            result.probs[i][j] = std::norm(psi[j * 3]);  // |psi[j][0]|^2
        }
    }
    
    std::cout << "Exiting simulate_oscillation() function successfully with E_nu = " << E_nu << std::endl;  // 輸出
    //outfile << "Exiting simulate_oscillation() function successfully with E_nu = " << E_nu << std::endl;
    //outfile.close(); // 關閉文件
    return result;
}

} // namespace neutrino