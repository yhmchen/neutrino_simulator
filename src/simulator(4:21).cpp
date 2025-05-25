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
#include <iomanip>

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

void normalize_density_matrix(std::array<std::complex<double>, 9>& rho) {
    double trace = std::real(rho[0] + rho[4] + rho[8]);  // rho_ee + rho_mumu + rho_tautau
    if (std::abs(trace) > 1e-12) {
        for (auto& x : rho) {
            x /= trace;
        }
    } else {
        // 處理 trace 非常接近零的情況
        // 可以選擇拋出例外或設定為預設值
        throw std::runtime_error("Trace of density matrix is close to zero!");
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

/* 改進的RK4步進函數
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
}*/

 //改進的RK4步進函數 commutator：-i (Hρ - ρH)
void rk4_step(double t, double h, std::array<std::complex<double>, 9>& y,
              double E_nu, const std::vector<double>& r_vals, 
              const std::vector<double>& Ne_vals) {

    //std::cout << "[rk4] 進入 rk4_step(), t = " << t << ", h = " << h << std::endl;
    std::array<std::complex<double>, 9> k1, k2, k3, k4, y_temp;
    
    // 如果能量非常高，可能需要更小的步長
    if (E_nu > 1e5 && h > 1e-8) {
        std::cout << "⚠️ 高能微中子可能需要更小的步長: 調整 h = " << h << " → ";
        h = std::min(h, 1e-8); // 為高能情況調整步長
        std::cout << h << std::endl;
    }
    
    // k1 = f(t, y)
    neutrino_evolution(t, y, E_nu, r_vals, Ne_vals, k1);
    
    // 檢查k1的穩定性
    double k1_norm = 0.0;
    for (const auto& val : k1) {
        k1_norm += std::norm(val);
        if (!std::isfinite(val.real()) || !std::isfinite(val.imag())) {
            throw std::runtime_error("RK4 computation failed: k1 contains NaN or Inf");
        }
    }
    k1_norm = std::sqrt(k1_norm);
    // 如果演化太微弱，可能需要增加物質效應或調整步長
    if (k1_norm < 1e-12 && E_nu > 1e5) {
        std::cout << "⚠️ 演化率非常小 |k1| = " << k1_norm << " @ r = " << t << std::endl;
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
    
    // 計算當前步的變化
    double step_change = 0.0;
    for (int i = 0; i < 9; i++) {
        step_change += std::norm(h/6.0 * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]));
    }
    step_change = std::sqrt(step_change);
    
    if (step_change < 1e-12 && E_nu > 1e5) {
        std::cout << "⚠️ 步長變化非常小: |Δρ| = " << step_change << " @ r = " << t << std::endl;
    }

    // 檢查結果是否穩定
    if (!check_numerical_stability(y)) {
        throw std::runtime_error("RK4 computation failed: result contains NaN or Inf");
    }
    // ✅ 新增：印出一點點 psi 看變化
    //std::cout << "[RK4] y[0][0] = " << y[0] << std::endl;
    
    // 規範化波函數以保持數值穩定性
    normalize_density_matrix(y);
}

// 適應性步長RK4
void adaptive_rk4_step(
    double &t, 
    double &h, 
    std::array<std::complex<double>, 9>& y,
    double E_nu, const std::vector<double>& r_vals, 
    const std::vector<double>& Ne_vals, double tolerance = 1.0e-8) {
    
    
                        
    //std::cout << "[RK4] → adaptive_rk4_step() 開始，t = " << t << ", h = " << h << ", E = " << E_nu << std::endl;
    // 儲存原始狀態以便需要回退
    std::array<std::complex<double>, 9> y_original = y;
    double original_h = h;
    
    try {
        // << "[RK4] 單步計算前，呼叫 rk4_step()" << std::endl;
        // 用一個步長 h 計算
        std::array<std::complex<double>, 9> y_single = y_original;
        rk4_step(t, h, y_single, E_nu, r_vals, Ne_vals);
        //std::cout << "[RK4] rk4_step() 成功完成！" << std::endl;

        
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
            if (error > 1e-6) {
                std::cout << "[RK4] error = " << error << ", step rejected." << std::endl;
            }
            double h_new = h * std::max(0.1, 0.9 * std::pow(tolerance / std::max(error, 1.0e-15), 0.25));
            h = std::max(h_new, 1.0e-12);  // 確保步長不會太小
            
            // 不更新t或y，將在下一次迭代中用新步長重試
            throw std::runtime_error("Step rejected, reducing step size to " + std::to_string(h));
        }
        
    } catch (const std::exception& e) {
        //std::cout << "[RK4] 發生未預期例外" << std::endl;
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



// ✅ 新版本：模擬積分用 Python 傳入的 r_vals，完全不再從密度檔案決定 r 點
struct SimulationResult;

SimulationResult simulate_custom_rvals(
    double E_nu,
    const std::vector<double>& input_r_vals,               // ✅ Python 傳入的距離點
    const std::vector<double>& density_r_vals,             // ✅ 從密度檔案讀取的 r 值
    const std::vector<double>& Ne_vals                     // ✅ 對應的 Ne(r)
) {
    //std::cout << "[C++] 進入 simulate_custom_rvals(), E_nu = " << E_nu << std::endl;

    std::cout << "==== 開始模擬能量為 " << E_nu << " eV 的微中子 ====" << std::endl;
    
    // 檢查能量是否過高，這可能導致數值問題
    if (E_nu > 1e5) {
        std::cout << "⚠️ 注意：模擬高能微中子 E = " << E_nu << " eV, 可能需要特別注意數值穩定性" << std::endl;
        
        // 計算預期的共振區
        const double delta21 = delta_m21_squared / (2.0 * E_nu);
        double resonance_Ve = delta21 * std::pow(std::cos(theta12), 2.0);
        std::cout << "MSW 共振條件 Ve ≈ " << resonance_Ve << " eV" << std::endl;
        
        // 檢查密度點是否足夠覆蓋共振區
        bool found_resonance = false;
        for (size_t i = 0; i < density_r_vals.size(); i++) {
            double Ve = std::sqrt(2.0) * 1.166e-23 * Ne_vals[i];
            if (std::abs(Ve - resonance_Ve) / resonance_Ve < 0.2) {
                found_resonance = true;
                std::cout << "✓ 找到可能的共振點: r = " << density_r_vals[i] 
                          << ", Ve = " << Ve << " eV" << std::endl;
            }
        }

        if (!found_resonance) {
            std::cout << "⚠️ 未找到符合共振條件的密度點，可能無法正確捕捉到 MSW 效應" << std::endl;
        }
    }

    SimulationResult result;
    int num_points = input_r_vals.size();

    result.r_vals = input_r_vals;
    result.probs.resize(num_points, std::vector<double>(3, 0.0));

    /* ✅ 在這裡加：檢查傳進來的 Ne 是否合理
    std::cout << "===== DEBUG: First 10 Ne values =====" << std::endl;
    for (int i = 0; i < std::min(10, (int)density_r_vals.size()); ++i) {
        std::cout << "[DEBUG] r = " << density_r_vals[i]
                  << ", Ne = " << Ne_vals[i] << " cm^-3" << std::endl;
    }*/
   
    // === 初始化 psi 為密度矩陣：|νe⟩⟨νe| ===
    //std::cout << "[C++] 準備初始化 psi..." << std::endl;
    std::array<std::complex<double>, 3> flavor_state = {
        std::complex<double>(1.0 , 0.0),
        std::complex<double>(0.0 , 0.0),
        std::complex<double>(0.0 , 0.0)
    };
    
    

    std::array<std::complex<double>, 9> psi{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            psi[i * 3 + j] = flavor_state[i] * std::conj(flavor_state[j]);
        }
    }

    // 初始化波函數 psi（3x3 矩陣展開為長度 9）
    //std::array<std::complex<double>, 9> psi{};
    //psi[0] = std::complex<double>(1.0, 0.0);  // 初始為電子微中子

    // 確認初始密度矩陣已正確設置
    std::cout << "初始密度矩陣: ρ_ee = " << psi[0] << ", ρ_μμ = " << psi[4] 
    << ", ρ_ττ = " << psi[8] << std::endl;

    // 打印初始機率
    std::cout << "初始機率: P_ee = " << result.probs[0][0] 
              << ", P_μμ = " << result.probs[0][1]
              << ", P_ττ = " << result.probs[0][2] << std::endl;

    // 初始點的機率
    //std::cout << "[C++] psi 初始化完成，psi[0][0] = " << psi[0] << std::endl;
    //normalize_density_matrix(psi);
    normalize_wavefunction(psi);
    result.probs[0][0] = std::real(psi[0]);
    result.probs[0][1] = std::real(psi[4]);
    result.probs[0][2] = std::real(psi[8]);
    

    double t = input_r_vals[0];
    double h = 1e-7;  // 初始步長，可再優化
    // 高能需要更小的初始步長
    if (E_nu > 1e5) {
        h = 1e-9;
        std::cout << "調整為較小的初始步長: h = " << h << std::endl;
    }
    /*
    std::cout << "[DEBUG] density_r_vals size = " << density_r_vals.size() << std::endl;
    std::cout << "[DEBUG] Ne_vals size         = " << Ne_vals.size() << std::endl;
    std::cout << "[DEBUG] input_r_vals[0]      = " << input_r_vals[0] << std::endl;
    std::cout << "[DEBUG] input_r_vals[1]      = " << input_r_vals[1] << std::endl;
    std::cout << "[DEBUG] E_nu                 = " << E_nu << std::endl;

    // 👇 嘗試提前呼叫 electron_potential 看會不會爆
    double ve_test = electron_potential(input_r_vals[1], density_r_vals, Ne_vals);
    std::cout << "[DEBUG] Ve 試算成功：Ve = " << ve_test << std::endl;


    std::cout << "[C++] 開始進入演化迴圈，第一個 r = " << input_r_vals[1] << std::endl;*/
    double res_val = (delta_m21_squared * std::cos(2 * theta12)) / (2.0 * E_nu);
    double prev_Ve = 0.0;
    bool in_resonance_region = false;
    //這邊都在除錯印各種資料
    for (int i = 1; i < num_points; ++i) {
        double next_r = input_r_vals[i];

        // ✅ 改為插值：確保 Ne 對應 next_r（終點）
        double G_F = 1.166e-23; // 單位：eV^-2
        double Ve = electron_potential(next_r, density_r_vals, Ne_vals); 
        double Ne = Ve / (std::sqrt(2.0) * G_F); // 還原電子密度 (可選擇是否保留這行)
        std::cout << "[DEBUG] E = " << E_nu << " eV, Ne = " << Ne
                  << " cm^-3, Ve = " << Ve << " eV @ r = " << next_r << std::endl;

        double ratio = Ve / E_nu;
        std::cout << "[DEBUG] Ve / E_nu = " << ratio << " @ r = " << next_r << std::endl;
        
        if ((prev_Ve - res_val) * (Ve - res_val) < 0) {
            std::cout << "✅ MSW 共振發生在 r ≈ " << next_r << std::endl;
            std::cout << "    |ψ_e|² = " << std::norm(psi[0]) << std::endl;
        }
        prev_Ve = Ve;
        
        double mass_term = delta_m21_squared / (2.0 * E_nu);
        std::cout << "V_e = " << Ve << ", Delta m^2 / 2E = " << mass_term << std::endl;
        std::cout << "r = " << next_r << ", V_e = " << Ve << ", |ψ_e|^2 = " << std::norm(psi[0]) << std::endl;
        std::cout << "[DEBUG] 共振條件 Ve ≈ " << res_val << " eV" << std::endl;
        std::cout << "[DEBUG] prev_Ve = " << prev_Ve << ", Ve = " << Ve << std::endl;

        // 這裡看Ve和共振條件delta_m^2/2E 的差距
        std::cout << "r = " << next_r << ", Ve = " << Ve
                  << ", Delta = " << mass_term
                  << ", Ve - Delta = " << Ve - mass_term << std::endl;
        
        const double resonance_ve = mass_term;
        if (std::abs(Ve - resonance_ve) / resonance_ve < 0.05) {  // 允許 5% 相對誤差
            std::cout << "⚠️ 在共振區域內，Ve ≈ " << Ve << ", Δ = " << resonance_ve << ", 步長 h = " << h << std::endl;
            h = std::min(h, 1e-10);  // ✅ 強制縮短步長以精準處理共振演化
            std::cout << "✅ 壓縮步長進入 MSW 區域：h = " << h << std::endl;
            std::cout << "🟡 [共振區] r = " << next_r << ", Ve = " << Ve 
            << ", Δ = " << resonance_ve 
            << ", |ψ_e|² = " << std::norm(psi[0]) << std::endl;
            std::cout << "[共振區] psi = ("
                      << psi[0] << "), (" << psi[4] << "), (" << psi[8] << ")" << std::endl;
            
            
        }
        
        if (std::abs(next_r - 0.547) < 1e-4) {
            std::cout << "🔥 進入 MSW resonance 點 r ≈ 0.547" << std::endl;
            std::cout << "ψ = (" << psi[0] << "), (" << psi[4] << "), (" << psi[8] << ")" << std::endl;
            std::cout << "P = " 
                      << std::real(psi[0]) << ", " 
                      << std::real(psi[4]) << ", " 
                      << std::real(psi[8]) << std::endl;
        }
        
        
        int max_steps = 100000;
        int step_count = 0;

        std::array<std::complex<double>, 9> last_successful_psi = psi;

        while (t < next_r && step_count < max_steps) {
            double remaining = next_r - t;
            if (h > remaining) h = remaining;

            try {
                adaptive_rk4_step(t, h, psi, E_nu, density_r_vals, Ne_vals);
                last_successful_psi = psi;
                step_count++;
            } catch (...) {
                h *= 0.5;
                if (h < 1e-10) {
                    psi = last_successful_psi;
                    break;
                }
            }
        }

        // 規範化與機率儲存
        normalize_wavefunction(psi);
        //normalize_density_matrix(psi);
        result.probs[i][0] = std::real(psi[0]);
        result.probs[i][1] = std::real(psi[4]);
        result.probs[i][2] = std::real(psi[8]);

        // === DEBUG: 變化率與最終波函數狀態 ===
        if (next_r != t) {
            double dVe_dr = (Ve - prev_Ve) / (next_r - t);
            std::cout << "[DEBUG] dVe/dr = " << dVe_dr << " eV / R_sun @ r = " << next_r << std::endl;
        }
        
        std::cout << "[DEBUG] Final psi = (" << psi[0] << "), (" << psi[4] << "), (" << psi[8] << ")" << std::endl;


        std::cout << "[DEBUG] Final psi = " << psi[0] << ", " << psi[4] << ", " << psi[8] << std::endl;
        std::cout << "[DEBUG] Final P = " 
                  << std::real(psi[0]) << ", " 
                  << std::real(psi[4]) << ", " 
                  << std::real(psi[8]) << std::endl;


        double P_sum = result.probs[i][0] + result.probs[i][1] + result.probs[i][2];
        if (std::abs(P_sum - 1.0) > 1e-6) {
            std::cout << "⚠️ Probability not normalized at r = " << next_r 
                      << ": P_total = " << P_sum << std::endl;
        }

        // 檢查演化是否正常 - 如果高能但機率幾乎沒變，可能有問題
        if (E_nu > 1e5 && i > 5 && 
            std::abs(result.probs[i][0] - result.probs[0][0]) < 1e-6) {
            std::cout << "⚠️ 警告：高能微中子但機率變化很小，可能有數值問題" << std::endl;
        }

        std::ofstream debug_out("cpp_debug.log", std::ios::app); // append mode
        debug_out << std::setprecision(10)
          << next_r << " " << std::real(psi[0]) << std::endl;
        //debug_out << std::setprecision(10)
        //  << result.r_vals[i] << " " << result.probs[i][0] << std::endl;
        //debug_out << "r = " << next_r << ", P_ee = " << std::norm(psi[0]) << std::endl;

    }
    const double theta_12 = 33.44 * M_PI / 180.0;
    double expected_Pee = std::pow(std::cos(theta_12), 2);  // 你自己用常數算出來
    std::cout << "[DEBUG] 預期 MSW 後 Pee ≈ " << expected_Pee << std::endl;
    std::cout << "模擬結果最終 Pee = " << result.probs.back()[0] << std::endl;
    
    // 檢查結果是否合理
    if (std::abs(result.probs.back()[0] - expected_Pee) < 0.1) {
        std::cout << "✅ 結果與理論預期接近!" << std::endl;
    } else {
        std::cout << "⚠️ 結果與理論預期有顯著差異!" << std::endl;
    }

    

    return result;
}
}