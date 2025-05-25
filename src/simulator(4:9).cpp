#include "../include/simulator.hpp"
#include "../include/constants.hpp"
#include "../include/electron_density.hpp"
#include "../include/hamiltonian.hpp"
#include <vector>
#include <complex>
#include <array>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <numeric>

namespace neutrino {

// 波函數規範化函數 - 確保波函數總概率為1
void normalize_wavefunction(std::array<std::complex<double>, 9>& psi) {
    double norm_squared = 0.0;
    
    // 計算總概率
    for (int i = 0; i < 9; i++) {
        norm_squared += std::norm(psi[i]);
    }
    
    // 規範化
    if (norm_squared > 0.0 && std::isfinite(norm_squared)) {
        double normalization_factor = 1.0 / std::sqrt(norm_squared);
        for (int i = 0; i < 9; i++) {
            psi[i] *= normalization_factor;
        }
    } else {
        // 處理NaN或Inf情況
        throw std::runtime_error("Wave function normalization failed due to NaN or Inf values");
    }
}

// 檢查數值穩定性
bool check_numerical_stability(const std::array<std::complex<double>, 9>& psi) {
    for (const auto& val : psi) {
        if (!std::isfinite(val.real()) || !std::isfinite(val.imag())) {
            return false;
        }
    }
    return true;
}

// 改進的RK4步進函數
void rk4_step(double t, double h, std::array<std::complex<double>, 9>& y,
              double E_nu, const std::vector<double>& r_vals, 
              const std::vector<double>& Ne_vals) {

    std::array<std::complex<double>, 9> k1, k2, k3, k4, y_temp;
    
    // k1 = f(t, y)
    neutrino_evolution(t, y, E_nu, r_vals, Ne_vals, k1);
    
    // 檢查k1的穩定性
    for (const auto& val : k1) {
        if (!std::isfinite(val.real()) || !std::isfinite(val.imag())) {
            throw std::runtime_error("RK4 computation failed: k1 contains NaN or Inf");
        }
    }
    
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
    }
    
    // 檢查結果是否穩定
    if (!check_numerical_stability(y)) {
        throw std::runtime_error("RK4 computation failed: result contains NaN or Inf");
    }
    
    // 規範化波函數以保持數值穩定性
    normalize_wavefunction(y);
}

// 適應性步長RK4
void adaptive_rk4_step(double &t, double &h, std::array<std::complex<double>, 9>& y,
                       double E_nu, const std::vector<double>& r_vals, 
                       const std::vector<double>& Ne_vals, double tolerance = 1.0e-6) {
    
    // 儲存原始狀態以便需要回退
    std::array<std::complex<double>, 9> y_original = y;
    double original_h = h;
    
    try {
        // 用一個步長 h 計算
        std::array<std::complex<double>, 9> y_single = y_original;
        rk4_step(t, h, y_single, E_nu, r_vals, Ne_vals);
        
        // 用兩個步長 h/2 計算
        std::array<std::complex<double>, 9> y_double = y_original;
        rk4_step(t, h/2, y_double, E_nu, r_vals, Ne_vals);
        rk4_step(t + h/2, h/2, y_double, E_nu, r_vals, Ne_vals);
        
        // 計算誤差估計
        double error = 0.0;
        for (int i = 0; i < 9; i++) {
            error += std::norm(y_single[i] - y_double[i]);
        }
        error = std::sqrt(error);
        
        // 決定是否接受這一步
        if (error <= tolerance) {
            // 接受這一步，計算下一步的步長
            y = y_double;  // 使用更精確的結果
            t += h;        // 更新時間
            
            // 調整下一步的步長
            double h_new = h * std::min(2.0, std::max(0.5, 0.9 * std::pow(tolerance / std::max(error, 1.0e-15), 0.2)));
            h = std::min(std::max(h_new, 1.0e-12), 1.0e-4);  // 限制步長在合理範圍內
        } 
        
            else {
            // 拒絕這一步，減小步長再試
            double h_new = h * std::max(0.1, 0.9 * std::pow(tolerance / std::max(error, 1.0e-15), 0.25));
            h = std::max(h_new, 1.0e-12);  // 確保步長不會太小
            
            // 不更新t或y，將在下一次迭代中用新步長重試
            throw std::runtime_error("Step rejected, reducing step size to " + std::to_string(h));
        }
        
    } catch (const std::exception& e) {
        // 如果發生錯誤，恢復原始狀態並減小步長
        y = y_original;
        
        // 如果不是因為誤差拒絕，就更激進地減小步長
        if (std::string(e.what()).find("Step rejected") == std::string::npos) {
            h = original_h * 0.1;
        }
        
        // 如果步長太小，放棄
        if (h < 1.0e-12) {
            throw std::runtime_error("Step size too small (< 1.0e-12), integration failed");
        }
    }
}


//一開始模擬4/9

SimulationResult simulate_oscillation(
    double E_nu,
    const std::string& solar_density_profile_path,
    double r_start,
    double r_end,
    int num_points
) {
    //std::ofstream debug_file("debug_simulation.txt");
    //debug_file << "Starting simulation with E_nu = " << E_nu << " MeV" << std::endl;
    //debug_file << "Integration range: " << r_start << " to " << r_end << " with " << num_points << " points" << std::endl;
    
    SimulationResult result;
    
    try {
        // [1]讀取太陽電子密度數據 - 使用現有函數，不修改
        std::vector<double> density_r_vals, Ne_vals;
        if (!solar_electron_density(solar_density_profile_path, density_r_vals, Ne_vals)) {
            throw std::runtime_error("Failed to read electron density data from " + solar_density_profile_path);
        }

        // 確保密度數據是按照半徑遞增排序的
        std::vector<size_t> indices(density_r_vals.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), 
                 [&density_r_vals](size_t i1, size_t i2) { 
                     return density_r_vals[i1] < density_r_vals[i2]; 
                 });
        
        std::vector<double> sorted_r_vals, sorted_Ne_vals;
        for (size_t i = 0; i < indices.size(); i++) {
            sorted_r_vals.push_back(density_r_vals[indices[i]]);
            sorted_Ne_vals.push_back(Ne_vals[indices[i]]);
        }

        
        density_r_vals = sorted_r_vals;
        Ne_vals = sorted_Ne_vals;

        /* 設定距離點的數量 4/7
        num_points = density_r_vals.size();
        result.r_vals.resize(num_points);
        result.probs.resize(num_points, std::vector<double>(3, 0.0));
        
        // 使用從檔案讀取的距離值
        for (int i = 0; i < num_points; i++) {
            result.r_vals[i] = density_r_vals[i];
        }

        // 生成距離點
        result.r_vals.resize(num_points);
        double dr = (r_end - r_start) / (num_points - 1);
        for (int i = 0; i < num_points; i++) {
            result.r_vals[i] = r_start + i * dr;
        }
        

        // 將 num_points 更新為總共的點數
        // [2] 建立 r_vals
        std::vector<double> r_vals;
        num_points = r_vals.size();
        result.r_vals = r_vals;
        result.probs.resize(num_points, std::vector<double>(3, 0.0));

        // 先放入太陽內部（直接用密度檔案裡的距離值）
        for (double r : density_r_vals) {
            if (r >= r_start && r <= 1.0) {  // 只要是 r_start 到 R_sun 的值都收
                r_vals.push_back(r);
            }
        }

        // 再補太陽外部的點（均勻模擬延伸）前先檢查
        if (num_points <= r_vals.size()) {
            std::cerr << "錯誤：num_points 太小，內部點就超過了！" << std::endl;
            exit(1);
        }

        // 再補太陽外部的點（均勻模擬延伸）
        double dr = (r_end - 1.0) / (num_points - r_vals.size());
        for (int i = 1; i < (num_points - r_vals.size() + 1); ++i) {
            r_vals.push_back(1.0 + i * dr);  // 從 R_sun 往外推
        }

        
        
        //debug_file << "Loaded electron density data: " << density_r_vals.size() << " points" << std::endl;
        
        // 輸出一些電子密度數據點進行驗證
        //if (!density_r_vals.empty()) {
            //debug_file << "First few density points:" << std::endl;
            //for (size_t i = 0; i < std::min(size_t(5), density_r_vals.size()); i++) {
                //debug_file << "  r = " << density_r_vals[i] << ", Ne = " << Ne_vals[i] << std::endl;
            //}
        //}
        
        
        
        // 初始條件 (單位矩陣) - 確保正確初始化
        std::array<std::complex<double>, 9> psi;
        for (int i = 0; i < 9; i++) {
            psi[i] = std::complex<double>(0.0, 0.0);
        }
        
        // 設置初始條件：假設系統開始於純電子微中子的狀態
        psi[0] = std::complex<double>(1.0, 0.0); // \( \psi_{ee} = 1 \)
        //psi[4] = std::complex<double>(0.0, 0.0); // \( \psi_{\mu\mu} = 0 \)
        //psi[8] = std::complex<double>(0.0, 0.0); // \( \psi_{\tau\tau} = 0 \)
        
        /* 計算初始點的概率
        result.probs[0][0] = std::norm(psi[0]); // \( P(\nu_e \to \nu_e) \)
        result.probs[0][1] = std::norm(psi[1]); // \( P(\nu_e \to \nu_\mu) \)
        result.probs[0][2] = std::norm(psi[2]); // \( P(\nu_e \to \nu_\tau) \)
            */
        /*正確設置對角元素為1
        for (int i = 0; i < 3; i++) {
            psi[i * 3 + i] = std::complex<double>(1.0, 0.0);
        }
        
        // 驗證初條件
        //debug_file << "Initial psi matrix:" << std::endl;
        //for (int i = 0; i < 3; i++) {
            //for (int j = 0; j < 3; j++) {
                //debug_file << psi[i * 3 + j] << " ";
            //}
            //debug_file << std::endl;
        //}
        
        // 初始化結果容器
        result.probs.resize(num_points, std::vector<double>(3, 0.0));
        
       
        /* 計算初始點的概率
        for (int j = 0; j < 3; j++) {
            result.probs[0][j] = std::norm(psi[j * 3]);  // |psi[j][0]|^2
        }
        */
        //debug_file << "Initial probabilities: ";
        //for (int j = 0; j < 3; j++) {
            //debug_file << result.probs[0][j] << " ";
        //}
        //debug_file << std::endl;
        
        // 數值積分
        // 3/28:double t = r_start;
        double t = result.r_vals[0]; // 從檔案讀取的第一個距離值
        // 根據問題特性選擇適當的初始步長
       // 3/28:double h = std::min(1.0e-7, dr / 1000.0);  // 非常保守的初始步長
        double h = std::min(1.0e-4, (result.r_vals[1] - result.r_vals[0]) / 50.0); // 根據檔案中的距離間隔設定步長
        //debug_file << "Starting integration with initial step size: " << h << std::endl;
        
        // 積分循環
        for (int i = 1; i < num_points; i++) {
            double next_r = result.r_vals[i];
            //debug_file << "Integrating to r = " << next_r << std::endl;
            
            int step_count = 0;
            int max_steps = 10000;  // 防止無限循環
            
            // 儲存上一個成功點的狀態，用於在完全失敗時回退
            std::array<std::complex<double>, 9> last_successful_psi = psi;
            
            while (t < next_r && step_count < max_steps) {
                try {
                    // 確保不會超過目標點
                    double remaining = next_r - t;
                    if (h > remaining) h = remaining;
                    
                    // 使用適應性步長積分
                    adaptive_rk4_step(t, h, psi, E_nu, density_r_vals, Ne_vals);
                    
                    step_count++;
                    
                    // 成功完成一步，更新最後成功狀態
                    last_successful_psi = psi;
                    
                    // 定期輸出調試信息
                    //if (step_count % 1000 == 0) {
                        //debug_file << "  Step " << step_count << ": t = " << t 
                                  //<< ", h = " << h << ", |psi|^2 = ";
                        //double norm_squared = 0.0;
                        //for (const auto& val : psi) {
                            //norm_squared += std::norm(val);
                        //}
                        //debug_file << norm_squared << std::endl;
                    //}
                    
                } catch (const std::exception& e) {
                    if (std::string(e.what()).find("Step rejected") != std::string::npos) {
                        // 這是正常的適應性步長調整，不是錯誤
                        continue;
                    }
                    
                    //debug_file << "Error during integration: " << e.what() << std::endl;
                    
                    // 如果步長已經很小且仍然失敗，可能需要特殊處理
                    if (h < 1.0e-10) {
                        //debug_file << "Fatal error: Integration failed with very small step size." << std::endl;
                        //debug_file << "Using last successful state and continuing..." << std::endl;
                        
                        // 使用上一個成功的狀態
                        psi = last_successful_psi;
                        // 直接跳到下一個輸出點
                        t += (next_r - t) * 0.01;
                        h = std::max(h, 1.0e-11);
                        break;
                    }
                    
                    // 否則減小步長並繼續
                    h *= 0.1;
                    //debug_file << "Reduced step size to " << h << std::endl;
                }
            }
            
            if (step_count >= max_steps) {
                //debug_file << "Warning: Reached maximum step count at r = " << next_r << std::endl;
                //debug_file << "Using last successful state and continuing..." << std::endl;
                t = next_r;  // 強制進入下一個點
            }
            
            // 確保波函數規範化
            try {
                normalize_wavefunction(psi);
            } catch (const std::exception& e) {
                //debug_file << "Normalization failed at r = " << next_r << ": " << e.what() << std::endl;
                //debug_file << "Resetting to unit matrix and continuing..." << std::endl;
                
                if (i > 1) {
                    //debug_file << "Resetting to last successful state..." << std::endl;
                    psi = last_successful_psi;
                    normalize_wavefunction(psi);
                } else {
                    //debug_file << "Resetting to initial condition..." << std::endl;
                    for (int k = 0; k < 9; k++) {
                        psi[k] = std::complex<double>(0.0, 0.0);
                    }
                    psi[0] = std::complex<double>(1.0, 0.0);
                }
            }
            
            // 計算當前點的概率
            //for (int j = 0; j < 3; j++) {
              //  result.probs[i][j] = std::norm(psi[j * 3]);  // |psi[j][0]|^2
            //}
            // 計算當前點的概率 - 修正計算方式
            result.probs[i][0] = std::norm(psi[0]); // P(νe → νe)
            result.probs[i][1] = std::norm(psi[1]); // P(νe → νμ)
            result.probs[i][2] = std::norm(psi[2]); // P(νe → ντ)

            double r = result.r_vals[i];
            /*
            double Pee = norm(psi[0]);
            double Pemu = norm(psi[1]);
            double Petau = norm(psi[2]);
            */
            std::cout << "r = " << r << ", P(νe -> νe) = " << result.probs[i][0]
                      << ", P(νe -> νμ) = " << result.probs[i][1]
                      << ", P(νe -> ντ) = " << result.probs[i][2] << std::endl;
            
        }    
            /* 輸出結果
            for (int i = 0; i < num_points; i++) {
                std::cout << "r = " << result.r_vals[i]
                          << ", P(νe -> νe) = " << result.probs[i][0]
                          << ", P(νe -> νμ) = " << result.probs[i][1]
                          << ", P(νe -> ντ) = " << result.probs[i][2]
                          << std::endl;
        }*/
            //debug_file << "At r = " << next_r << ", probabilities: "
                      //<< result.probs[i][0] << ", " 
                      //<< result.probs[i][1] << ", " 
                     // << result.probs[i][2] << std::endl;
        //}
        
        //debug_file << "Simulation completed successfully" << std::endl;
        std::ofstream output_file("output_data.txt");
        if (output_file.is_open()) {
            for (int i = 0; i < num_points; ++i) {
                output_file << result.r_vals[i] << " "
                            << result.probs[i][0] << " "  // P(νe -> νe)
                            << result.probs[i][1] << " "  // P(νe -> νμ)
                            << result.probs[i][2] << std::endl; // P(νe -> ντ)
            }
            output_file.close();
            std::cout << "Simulation results written to output_data.txt" << std::endl;
        } else {
            std::cerr << "Unable to open file for writing: output_data.txt" << std::endl;}
        
    } catch (const std::exception& e) {
        std::cerr << "Error in simulation: " << e.what() << std::endl;
        //debug_file << "Simulation failed: " << e.what() << std::endl;
        // 返回空結果
        result.r_vals.clear();
        result.probs.clear();
    }
    
    //debug_file.close();
    return result;
}

} // namespace neutrino