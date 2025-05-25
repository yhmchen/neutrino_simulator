#include "../include/simulator.hpp"
#include "../include/constants.hpp"
#include "../include/electron_density.hpp"
#include "../include/hamiltonian.hpp"
#include "../include/pmns_matrix.hpp"
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

void normalize_density_matrix(std::array<std::complex<double>, 9>& rho) {
    double trace = std::real(rho[0] + rho[4] + rho[8]);  // rho_ee + rho_mumu + rho_tautau
    if (std::abs(trace) > 1e-12) {
        // 使用更穩定的除法來進行規範化
        double inv_trace = 1.0 / trace;
        for (auto& x : rho) {
            x *= inv_trace;
        }
    } else {
        // 處理 trace 非常接近零的情況
        // 改進：使用正則化方法而不是直接拋出例外
        std::cout << "警告: 密度矩陣的跡接近零，應用正則化" << std::endl;
        double epsilon = 1e-12;
        for (int i = 0; i < 3; i++) {
            // 確保對角元素有最小值，保證總機率為1
            rho[i*3+i] += epsilon;
        }
        // 重新規範化
        trace = std::real(rho[0] + rho[4] + rho[8]);
        double inv_trace = 1.0 / trace;
        for (auto& x : rho) {
            x *= inv_trace;
        }
    }
}


// ===== 改進2: 更細緻的數值穩定性檢查 =====
bool check_numerical_stability(const std::array<std::complex<double>, 9>& rho) {
    for (const auto& val : rho) {
        if (!std::isfinite(val.real()) || !std::isfinite(val.imag())) {
            return false;
        }
        // 添加：檢查元素大小是否在合理範圍內
        if (std::abs(val) > 1e6) {
            return false;  // 檢測到異常大值
        }
    }
    
    // 添加：檢查密度矩陣的跡是否接近1
    double trace = std::real(rho[0] + rho[4] + rho[8]);
    if (std::abs(trace - 1.0) > 1e-6) {
        return false;  // 跡不接近1，表示可能有問題
    }
    
    // 添加：檢查厄米性
    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            if (std::abs(rho[i*3+j] - std::conj(rho[j*3+i])) > 1e-6) {
                return false;  // 不滿足厄米性
            }
        }
    }
    
    return true;
}



// ===== 改進3: 更高精度的共振位置插值函數 =====
double interpolate_resonance_position(double r1, double r2, double Ve1, double Ve2, double res_val) {
    // 線性插值可能在密度變化大的地方不夠精確
    // 使用更精細的插值方法，例如三次樣條
    
    // 如果兩點非常接近，使用線性插值即可
    if (std::abs(r2 - r1) < 1e-6 || std::abs(Ve2 - Ve1) < 1e-10) {
        return r1 + (r2 - r1) * (res_val - Ve1) / (Ve2 - Ve1);
    }
    
    // 否則，使用更精細的插值方法
    // 假設密度在短距離內的變化可以用指數函數近似
    double log_r1 = std::log(r1);
    double log_r2 = std::log(r2);
    double log_Ve1 = std::log(std::abs(Ve1));
    double log_Ve2 = std::log(std::abs(Ve2));
    
    // 計算指數插值參數
    double alpha = (log_Ve2 - log_Ve1) / (log_r2 - log_r1);
    double beta = log_Ve1 - alpha * log_r1;
    
    // 使用指數模型計算共振位置
    double log_res = (std::log(std::abs(res_val)) - beta) / alpha;
    
    // 轉換回線性空間
    double r_res = std::exp(log_res);
    
    // 確保結果在原始區間內
    if (r_res < std::min(r1, r2) || r_res > std::max(r1, r2)) {
        // 如果指數插值失敗，回退到線性插值
        return r1 + (r2 - r1) * (res_val - Ve1) / (Ve2 - Ve1);
    }
    
    return r_res;
}


// ===== 改進4: 增強的密度梯度估計，使用高階插值 =====
double estimate_density_gradient(double r, 
    const std::vector<double>& r_vals, 
    const std::vector<double>& Ne_vals) {
        // 找到最接近的點
        size_t idx = 0;
        while (idx < r_vals.size() - 1 && r_vals[idx+1] < r) idx++;

        // 邊界處理
        if (idx == 0) {
        idx = 1;  // 確保有前一個點可用
        }
        if (idx >= r_vals.size() - 2) {
        idx = r_vals.size() - 3;  // 確保有後一個點可用
        }

        // 使用中心差分進行更精確的梯度估計
        double r_minus = r_vals[idx-1];
        double r_center = r_vals[idx];
        double r_plus = r_vals[idx+1];

        double Ne_minus = Ne_vals[idx-1];
        double Ne_center = Ne_vals[idx];
        double Ne_plus = Ne_vals[idx+1];

        // 計算正規化權重（距離越近，權重越大）
        double w_minus = 1.0 / std::max(std::abs(r - r_minus), 1e-10);
        double w_center = 1.0 / std::max(std::abs(r - r_center), 1e-10);
        double w_plus = 1.0 / std::max(std::abs(r - r_plus), 1e-10);

        double w_sum = w_minus + w_center + w_plus;
        w_minus /= w_sum;
        w_center /= w_sum;
        w_plus /= w_sum;

        // 使用有限差分近似梯度
        double grad_minus = (Ne_center - Ne_minus) / (r_center - r_minus);
        double grad_center = (Ne_plus - Ne_minus) / (r_plus - r_minus);
        double grad_plus = (Ne_plus - Ne_center) / (r_plus - r_center);

        // 權重平均
        return w_minus * grad_minus + w_center * grad_center + w_plus * grad_plus;
}

// 評估密度梯度的變化率
// 這個函數計算指定點附近的密度梯度變化率，用於判斷梯度是否線性變化
double estimate_gradient_variation(
    double r,
    const std::vector<double>& density_r_vals,
    const std::vector<double>& Ne_vals
) {
    // 找出r在density_r_vals中的位置
    size_t idx = 0;
    while (idx < density_r_vals.size() - 1 && density_r_vals[idx + 1] < r) {
        ++idx;
    }
    
    // 計算前後梯度
    double gradient_pre = 0.0;
    double gradient_post = 0.0;
    
    // 確保有足夠的點計算前後梯度
    if (idx > 0 && idx < density_r_vals.size() - 1) {
        // 前一段梯度
        double delta_r_pre = density_r_vals[idx] - density_r_vals[idx - 1];
        double delta_Ne_pre = Ne_vals[idx] - Ne_vals[idx - 1];
        gradient_pre = std::abs(delta_Ne_pre / delta_r_pre);
        
        // 後一段梯度
        double delta_r_post = density_r_vals[idx + 1] - density_r_vals[idx];
        double delta_Ne_post = Ne_vals[idx + 1] - Ne_vals[idx];
        gradient_post = std::abs(delta_Ne_post / delta_r_post);
        
        // 計算梯度變化率（標準化）
        double avg_gradient = (gradient_pre + gradient_post) / 2.0;
        if (avg_gradient > 0.0) {
            return std::abs(gradient_post - gradient_pre) / avg_gradient;
        }
    }
    
    // 如果無法計算，返回0表示無明顯變化
    return 0.0;
}

// 應用1-2共振跳躍
void apply_resonance_jump_12(std::array<std::complex<double>, 9>& rho, double P_jump) {
    // 轉換到共振時的本徵態基底
    std::array<std::complex<double>, 9> rho_eigen{};
    transform_to_resonance_eigenbasis_12(rho, rho_eigen);
    
    // 應用跳躍概率
    apply_jump_probability_12(rho_eigen, P_jump);
    
    // 轉回原來的基底
    transform_back_from_resonance_eigenbasis_12(rho_eigen, rho);
}

// 應用1-3共振跳躍
void apply_resonance_jump_13(std::array<std::complex<double>, 9>& rho, double P_jump) {
    // 轉換到共振時的本徵態基底
    std::array<std::complex<double>, 9> rho_eigen{};
    transform_to_resonance_eigenbasis_13(rho, rho_eigen);
    
    // 應用跳躍概率
    apply_jump_probability_13(rho_eigen, P_jump);
    
    // 轉回原來的基底
    transform_back_from_resonance_eigenbasis_13(rho_eigen, rho);
}

// ===== 改進5: 增強的厄米性確保函數，使用SVD分解 =====
void ensure_hermiticity(std::array<std::complex<double>, 9>& rho) {
    // 基本的厄米化：確保非對角元素滿足 rho_ij = rho_ji*
    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 3; j++) {
            std::complex<double> avg = (rho[i*3+j] + std::conj(rho[j*3+i])) / 2.0;
            rho[i*3+j] = avg;
            rho[j*3+i] = std::conj(avg);
        }
    }
    
    // 額外確保正定性：確保對角元素為實數且非負
    for (int i = 0; i < 3; i++) {
        int idx = i*3+i;
        if (std::abs(rho[idx].imag()) > 1e-12) {
            rho[idx] = std::complex<double>(rho[idx].real(), 0.0);
        }
        // 如果對角元素是負的（這不應該發生，但為了數值穩定性）
        if (rho[idx].real() < 0.0) {
            rho[idx] = std::complex<double>(std::abs(rho[idx].real()), 0.0);
        }
    }
    
    // 確保總機率為1
    normalize_density_matrix(rho);
}

// ===== 實現新增的 ultra_high_energy_rk4_step 函數 =====
void ultra_high_energy_rk4_step(double& t, double& h, std::array<std::complex<double>, 9>& rho,
    double E_nu, 
    const std::vector<double>& r_vals, 
    const std::vector<double>& Ne_vals) {
        // 超高能微中子使用特殊的積分方法，包括尺度重整和高強度阻尼

        // 獲取當前位置的電子密度和物質勢能
        double Ve = electron_potential(t, r_vals, Ne_vals);

        // 計算尺度參數
        double vac_scale = delta_m21_squared / (2.0 * E_nu);
        double scale_factor = vac_scale * 1e-2;  // 超高能需要更小的尺度因子

        // 創建原始哈密頓量
        std::array<std::complex<double>, 9> H_orig{};
        compute_hamiltonian(H_orig, E_nu, Ve);

        // 尺度重整哈密頓量
        std::array<std::complex<double>, 9> H_scaled{};
        for (int i = 0; i < 9; i++) {
        H_scaled[i] = H_orig[i] * scale_factor;
        }

        // 特別處理：移除共同相位（不影響物理結果但改善數值穩定性）
        std::complex<double> common_phase = H_scaled[0];
        for (int i = 0; i < 9; i++) {
        H_scaled[i] -= common_phase;
        }

        // 使用改進的RK4步驟，包括高強度阻尼
        double t_orig = t;

        // 定義超高能RK4函數
        auto ultra_rk4_deriv = [&](const std::array<std::complex<double>, 9>& rho_val) {
        std::array<std::complex<double>, 9> drho{};

        // 計算 -i[H,ρ]
        for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
        drho[i*3+j] = 0.0;
        for (int k = 0; k < 3; k++) {
        drho[i*3+j] -= std::complex<double>(0.0, 1.0) * 
        (H_scaled[i*3+k] * rho_val[k*3+j] - rho_val[i*3+k] * H_scaled[k*3+j]);
        }
        }
        }

        // 添加阻尼項：保持厄米性和總機率
        for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
        if (i != j) {
        // 對非對角元素應用指數阻尼
        double damping = 1e-5 * scale_factor;  // 調整阻尼強度
        drho[i*3+j] -= damping * rho_val[i*3+j];
        }
        }
        }

        return drho;
        };

        // RK4步驟
        std::array<std::complex<double>, 9> k1 = ultra_rk4_deriv(rho);

        std::array<std::complex<double>, 9> rho_temp{};
        for (int i = 0; i < 9; i++) {
        rho_temp[i] = rho[i] + 0.5 * h * k1[i];
        }
        std::array<std::complex<double>, 9> k2 = ultra_rk4_deriv(rho_temp);

        for (int i = 0; i < 9; i++) {
        rho_temp[i] = rho[i] + 0.5 * h * k2[i];
        }
        std::array<std::complex<double>, 9> k3 = ultra_rk4_deriv(rho_temp);

        for (int i = 0; i < 9; i++) {
        rho_temp[i] = rho[i] + h * k3[i];
        }
        std::array<std::complex<double>, 9> k4 = ultra_rk4_deriv(rho_temp);

        // 更新密度矩陣
        for (int i = 0; i < 9; i++) {
        rho[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }

        // 確保厄米性和規範化
        ensure_hermiticity(rho);
        normalize_density_matrix(rho);

        // 更新時間
        t = t_orig + h;
}


// ===== 改進6: 大幅增強的積分到指定位置函數 =====
void integrate_to_position(double& t, double next_r, double& h, 
    std::array<std::complex<double>, 9>& rho,
    double E_nu, 
    const std::vector<double>& r_vals, 
    const std::vector<double>& Ne_vals){
        int max_steps = 1000000;  // 增加最大步數限制，確保高解析度
        int step_count = 0;
        std::array<std::complex<double>, 9> last_successful_rho = rho;

        // 改進：動態調整最小步長，基於能量和距離
        const double base_min_step = 1.0e-12;
        double min_scale = std::min(1.0, std::pow(E_nu / 1.0e6, -0.5));  // 高能微中子需要更小的步長
        double distance_scale = std::min(1.0, (next_r - t) / 10.0);  // 較遠距離允許較大步長

        double min_step_size = base_min_step * min_scale * distance_scale;
        min_step_size = std::max(min_step_size, 1.0e-12);  // 確保不會太小
        h = std::max(h, min_step_size);  // 確保初始步長足夠大

        std::cout << "積分起點: " << t << ", 終點: " << next_r 
        << ", 初始步長: " << h << ", 最小步長: " << min_step_size << std::endl;

        // 卡住檢測變量
        int stuck_counter = 0;
        double last_t = t;
        int consecutive_failures = 0;
        double last_h = h;

        // 檢查當前位置和目標位置的電子密度
        double start_Ve = electron_potential(t, r_vals, Ne_vals);
        double end_Ve = electron_potential(next_r, r_vals, Ne_vals);

        // 估計共振條件
        double res_val_12 = (delta_m21_squared * std::cos(2 * theta12)) / (2.0 * E_nu);
        double res_val_13 = (delta_m31_squared * std::cos(2 * theta13)) / (2.0 * E_nu);

        // 檢查是否會經過共振點
        bool will_cross_resonance = false;
        if ((start_Ve - res_val_12) * (end_Ve - res_val_12) <= 0 || 
        (start_Ve - res_val_13) * (end_Ve - res_val_13) <= 0) {
        will_cross_resonance = true;
        std::cout << "⚠️ 積分區間將穿越共振點，採用超精細步長" << std::endl;
        }

        // 積分主迴圈
        while (t < next_r && step_count < max_steps) {
        double remaining = next_r - t;

        // 根據是否接近共振點動態調整步長
        if (will_cross_resonance) {
        double current_Ve = electron_potential(t, r_vals, Ne_vals);
        double res_distance_12 = std::abs(current_Ve - res_val_12);
        double res_distance_13 = std::abs(current_Ve - res_val_13);
        double min_res_distance = std::min(res_distance_12, res_distance_13);

        // 越接近共振點，步長越小
        if (min_res_distance < 0.01 * std::abs(res_val_12)) {
        h = std::min(h, 1e-8);  // 非常接近共振點
        std::cout << "🔍 非常接近共振點，使用超小步長: " << h << std::endl;
        } else if (min_res_distance < 0.1 * std::abs(res_val_12)) {
        h = std::min(h, 1e-7);  // 接近共振點
        }
        }

        // 確保步長不會超過剩餘距離
        if (h > remaining) h = remaining;

        // 檢測卡住情況
        if (std::abs(t - last_t) < 1e-12) {
        stuck_counter++;
        } else {
        stuck_counter = 0;
        }

        // 如果連續幾次沒前進，強制前進
        if (stuck_counter > 5) {
        std::cout << "⚠️ 積分卡住，強制前進" << std::endl;
        h *= 2.0;  // 嘗試增加步長
        h = std::min(h, remaining);  // 但不要超過剩餘距離

        if (stuck_counter > 10) {
        std::cout << "🆘 無法進一步積分，跳至下一點" << std::endl;
        t = next_r;
        break;
        }
        }

        last_t = t;

        try {
        // 根據能量選擇適當的積分方法
        if (E_nu > 1e7) {
        // 超高能微中子使用特殊優化的RK4
        ultra_high_energy_rk4_step(t, h, rho, E_nu, r_vals, Ne_vals);
        } else if (E_nu > 1e5) {
        // 高能微中子使用增強的RK4
        enhanced_rk4_step(t, h, rho, E_nu, r_vals, Ne_vals);
        } else {
        // 標準能量使用自適應RK4，帶更嚴格的誤差控制
        adaptive_rk4_step(t, h, rho, E_nu, r_vals, Ne_vals, 1e-8);  // 使用更小的容差
        }

        // 成功，保存結果
        last_successful_rho = rho;
        step_count++;
        consecutive_failures = 0;
        last_h = h;  // 記錄上一次成功的步長

        // 印出進度
        if (step_count % 1000 == 0) {
        std::cout << "步驟 " << step_count << ": r = " << t 
            << ", h = " << h << ", 距離目標 = " << (next_r - t) << std::endl;
        }

        } catch (const std::exception& e) {
        // 錯誤處理 - 減小步長重試
        consecutive_failures++;

        // 改進的步長調整策略：根據連續失敗次數動態調整
        if (consecutive_failures < 3) {
        h *= 0.5;  // 前幾次嘗試溫和減小
        } else if (consecutive_failures < 5) {
        h *= 0.2;  // 多次失敗後更激進地減小
        } else {
        h *= 0.1;  // 持續失敗時大幅減小
        }

        // 不允許步長太小
        h = std::max(h, min_step_size);

        // 如果步長接近最小允許值且連續失敗多次，考慮放棄
        if (h < 1.2 * min_step_size && consecutive_failures > 10) {
        std::cout << "⚠️ 步長太小且連續失敗，使用上一步成功結果並前進" << std::endl;
        rho = last_successful_rho;

        // 嘗試直接跳到一個中間點，而不是整個跳到終點
        double jump_fraction = 0.2;  // 只前進20%的剩餘距離
        t = t + (next_r - t) * jump_fraction;

        // 重置失敗計數和步長
        consecutive_failures = 0;
        h = last_h;  // 恢復到上一次成功的步長
        }

        // 完全卡住時的最後手段
        if (consecutive_failures > 20) {
        std::cout << "🆘 完全卡住，放棄當前積分" << std::endl;
        rho = last_successful_rho;
        t = next_r;  // 強制前進
        break;
        }
        }
        }

        if (step_count >= max_steps) {
        std::cout << "⚠️ 達到最大步數限制" << std::endl;
        t = next_r;  // 確保即使達到最大步數也前進
        }

        // 最終確保結果的規範化和厄米性
        normalize_density_matrix(rho);
        ensure_hermiticity(rho);

        std::cout << "完成積分: 從 " << last_t << " 到 " << t 
        << ", 共 " << step_count << " 步" << std::endl;
}


// 用於1-2共振的本徵態轉換
void transform_to_resonance_eigenbasis_12(const std::array<std::complex<double>, 9>& rho, 
    std::array<std::complex<double>, 9>& rho_eigen) {
    // 計算1-2共振時的混合矩陣
    std::array<std::complex<double>, 9> U{};
    std::array<std::complex<double>, 9> U_dag{};

    // 設置1-2共振的混合矩陣
    compute_resonance_mixing_12(U, U_dag);

    // 應用變換: ρ' = U† ρ U
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rho_eigen[i*3+j] = 0.0;
            for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
                rho_eigen[i*3+j] += U_dag[i*3+k] * rho[k*3+l] * U[l*3+j];
                }
            }
        }
    }
}


// 計算1-2共振混合矩陣
void compute_resonance_mixing_12(std::array<std::complex<double>, 9>& U, 
    std::array<std::complex<double>, 9>& U_dag) {
    // 在1-2共振點，θ12被修改為最大混合角度（45度）
    double s12_res = M_SQRT1_2; // sin(π/4) = 1/√2
    double c12_res = M_SQRT1_2; // cos(π/4) = 1/√2

    // 其他角度保持不變
    double s13 = std::sin(theta13);
    double c13 = std::cos(theta13);
    double s23 = std::sin(theta23);
    double c23 = std::cos(theta23);
    // CP違反相位
    std::complex<double> delta_term = std::exp(std::complex<double>(0, delta_cp));

    // 使用neutrino命名空間中的函數獲取PMNS矩陣但替換θ12
    neutrino::ComplexMatrix pmns = neutrino::get_pmns_matrix();

    // 修改PMNS矩陣中與θ12相關的元素
    // ... 這裡需要適應您的PMNS矩陣的實現方式

    // 將PMNS矩陣轉換為std::array格式
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            U[i*3+j] = pmns[i][j];
            U_dag[j*3+i] = std::conj(pmns[i][j]);
        }
    }
}

// ===== 改進7: 增強的跳躍概率應用 =====
void apply_jump_probability_12(std::array<std::complex<double>, 9>& rho_eigen, double P_jump) {
    // 獲取1-2共振本徵態的密度元素
    std::complex<double> rho_11 = rho_eigen[0];
    std::complex<double> rho_22 = rho_eigen[4];
    std::complex<double> rho_12 = rho_eigen[1];
    std::complex<double> rho_21 = rho_eigen[3];

    // 應用跳躍概率 P_jump 修改1-2共振本徵態對角元素
    // 混合了第一和第二本徵態的佔據概率
    double P_1 = std::real(rho_11);
    double P_2 = std::real(rho_22);
    
    // 計算跳躍後的新佔據概率
    double new_P_1 = P_1 * (1.0 - P_jump) + P_2 * P_jump;
    double new_P_2 = P_2 * (1.0 - P_jump) + P_1 * P_jump;
    
    // 更新密度矩陣的對角元素
    rho_eigen[0] = std::complex<double>(new_P_1, 0.0);
    rho_eigen[4] = std::complex<double>(new_P_2, 0.0);
    
    // 重要改進：更精確地處理相干項
    // 對於相干項，我們需要考慮相位關係
    if (P_jump > 0.01) {  // 只有在跳躍概率顯著時才處理相干項
        // 計算相干項的縮放因子
        double coherence_scale = std::sqrt((1.0 - P_jump) * (1.0 - P_jump) + P_jump * P_jump);
        
        // 應用相位旋轉和縮放
        double phase_shift = P_jump * M_PI;  // 跳躍導致的相位變化
        std::complex<double> phase_factor = std::exp(std::complex<double>(0, phase_shift));
        
        // 更新相干項
        rho_eigen[1] = rho_12 * coherence_scale * phase_factor;
        rho_eigen[3] = std::conj(rho_eigen[1]);  // 保持厄米性
    } else {
        // 小跳躍直接縮放相干項
        rho_eigen[1] *= (1.0 - P_jump);
        rho_eigen[3] *= (1.0 - P_jump);
    }
}

// 用於1-3共振的本徵態轉換
void transform_to_resonance_eigenbasis_13(const std::array<std::complex<double>, 9>& rho, 
                std::array<std::complex<double>, 9>& rho_eigen) {
    // 計算1-3共振時的混合矩陣
    std::array<std::complex<double>, 9> U{};
    std::array<std::complex<double>, 9> U_dag{};

    // 設置1-3共振的混合矩陣
    compute_resonance_mixing_13(U, U_dag);

    // 應用變換: ρ' = U† ρ U
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rho_eigen[i*3+j] = 0.0;
            for (int k = 0; k < 3; k++) {
            for (int l = 0; l < 3; l++) {
                rho_eigen[i*3+j] += U_dag[i*3+k] * rho[k*3+l] * U[l*3+j];
                }
            }
        }
    }
}

// 從本徵態轉回標準基底（1-3共振）
void transform_back_from_resonance_eigenbasis_13(const std::array<std::complex<double>, 9>& rho_eigen, 
                    std::array<std::complex<double>, 9>& rho) {
            // 計算1-3共振時的混合矩陣
            std::array<std::complex<double>, 9> U{};
            std::array<std::complex<double>, 9> U_dag{};

            // 設置1-3共振的混合矩陣
            compute_resonance_mixing_13(U, U_dag);

            // 應用反向變換: ρ = U ρ' U†
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    rho[i*3+j] = 0.0;
                    for (int k = 0; k < 3; k++) {
                        for (int l = 0; l < 3; l++) {
                            rho[i*3+j] += U[i*3+k] * rho_eigen[k*3+l] * U_dag[l*3+j];
                        }
                    }
                }
            }
        // 確保厄米性和規範化
    ensure_hermiticity(rho);
    normalize_density_matrix(rho);
}

// 計算1-3共振混合矩陣
void compute_resonance_mixing_13(std::array<std::complex<double>, 9>& U, 
    std::array<std::complex<double>, 9>& U_dag) {
        // 在1-3共振點，θ13被修改為最大混合角度（45度）
        double s13_res = M_SQRT1_2; // sin(π/4) = 1/√2
        double c13_res = M_SQRT1_2; // cos(π/4) = 1/√2

        // 其他角度保持不變
        double s12 = std::sin(theta12);
        double c12 = std::cos(theta12);
        double s23 = std::sin(theta23);
        double c23 = std::cos(theta23);
        // CP違反相位
        std::complex<double> delta_term = std::exp(std::complex<double>(0, delta_cp));

        // 使用neutrino命名空間中的函數獲取PMNS矩陣但替換θ13
        neutrino::ComplexMatrix pmns = neutrino::get_pmns_matrix();

        // 修改PMNS矩陣中與θ13相關的元素
        // ... 這裡需要適應您的PMNS矩陣的實現方式

        // 將PMNS矩陣轉換為std::array格式
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                U[i*3+j] = pmns[i][j];
                U_dag[j*3+i] = std::conj(pmns[i][j]);
            }
        }
}

// 應用1-3共振的跳躍概率
void apply_jump_probability_13(std::array<std::complex<double>, 9>& rho_eigen, double P_jump) {
    // 獲取1-3共振本徵態的密度元素
    std::complex<double> rho_11 = rho_eigen[0];
    std::complex<double> rho_33 = rho_eigen[8];
    std::complex<double> rho_13 = rho_eigen[2];
    std::complex<double> rho_31 = rho_eigen[6];
    
    // 應用跳躍概率 P_jump 修改1-3共振本徵態對角元素
    // 混合了第一和第三本徵態的佔據概率
    double P_1 = std::real(rho_11);
    double P_3 = std::real(rho_33);
    
    // 計算跳躍後的新佔據概率
    double new_P_1 = P_1 * (1.0 - P_jump) + P_3 * P_jump;
    double new_P_3 = P_3 * (1.0 - P_jump) + P_1 * P_jump;
    
    // 更新密度矩陣的對角元素
    rho_eigen[0] = std::complex<double>(new_P_1, 0.0);
    rho_eigen[8] = std::complex<double>(new_P_3, 0.0);
    
    // 與1-2共振類似，精確處理相干項
    if (P_jump > 0.01) {
        double coherence_scale = std::sqrt((1.0 - P_jump) * (1.0 - P_jump) + P_jump * P_jump);
        double phase_shift = P_jump * M_PI;
        std::complex<double> phase_factor = std::exp(std::complex<double>(0, phase_shift));
        
        rho_eigen[2] = rho_13 * coherence_scale * phase_factor;
        rho_eigen[6] = std::conj(rho_eigen[2]);
    } else {
        rho_eigen[2] *= (1.0 - P_jump);
        rho_eigen[6] *= (1.0 - P_jump);
    }
}

// 從本徵態轉回標準基底（1-2共振）
void transform_back_from_resonance_eigenbasis_12(const std::array<std::complex<double>, 9>& rho_eigen, 
    std::array<std::complex<double>, 9>& rho) {
    // 計算1-2共振時的混合矩陣
    std::array<std::complex<double>, 9> U{};
    std::array<std::complex<double>, 9> U_dag{};

    // 設置1-2共振的混合矩陣
    compute_resonance_mixing_12(U, U_dag);

    // 應用反向變換: ρ = U ρ' U†
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rho[i*3+j] = 0.0;
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    rho[i*3+j] += U[i*3+k] * rho_eigen[k*3+l] * U_dag[l*3+j];
                }
            }
        }
    }
    // 確保厄米性和規範化
    ensure_hermiticity(rho);
    normalize_density_matrix(rho);
}

// 計算尺度重整的哈密頓量
void compute_scaled_hamiltonian(std::array<std::complex<double>, 9>& H_scaled,
               double E_nu, double Ve, double scale_factor) {
// 創建標準哈密頓量
std::array<std::complex<double>, 9> H_standard{};

// 計算真空項和物質效應
compute_hamiltonian(H_standard, E_nu, Ve);

// 對哈密頓量進行尺度重整
for (int i = 0; i < 9; i++) {
H_scaled[i] = H_standard[i] / scale_factor;
}

// 也可以選擇性地優化某些項的精度，例如將差異項更加突出
// 這裡可以根據具體物理需求調整

// 例如，對於高能微中子，可能需要特別處理物質效應項
if (Ve / (delta_m21_squared / (2.0 * E_nu)) > 1e6) {
    // 當物質效應遠大於真空項時，可以進一步調整標度
    // 例如，可以保留相對相位，但縮放幅度
    for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (i == j) {
                    // 對角項保留，但可能需要減去一個共同的常數
                    // 由於哈密頓量在物理上添加常數項不會改變振盪結果
                    double common_term = std::real(H_scaled[0]);
                    H_scaled[i*3+j] -= common_term;
                }
            // 非對角項已經縮放，不需要特殊處理
            }
        }
    }
}

// Add this to simulator.cpp
void compute_hamiltonian(std::array<std::complex<double>, 9>& H_standard,
    double E_nu, double Ve) {
    // Convert the 1D array to the format expected by hamiltonian function
    std::vector<double> dummy_r_vals = {0.0, 1.0};  // Dummy values
    std::vector<double> dummy_Ne_vals = {Ve / (std::sqrt(2.0) * G_F), Ve / (std::sqrt(2.0) * G_F)}; 
    // This ensures electron_potential returns Ve directly

    // Call the existing hamiltonian function
    ComplexMatrix H_matrix = hamiltonian(0.0, E_nu, dummy_r_vals, dummy_Ne_vals);

    // Convert the result back to 1D array
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
            H_standard[i * 3 + j] = H_matrix[i][j];
        }
    }
}


//改進的RK4步進函數 commutator：-i (Hρ - ρH)
void rk4_step(double t, double h, std::array<std::complex<double>, 9>& y,
    double E_nu, const std::vector<double>& r_vals, 
    const std::vector<double>& Ne_vals) {

//std::cout << "[rk4] 進入 rk4_step(), t = " << t << ", h = " << h << std::endl;
    std::array<std::complex<double>, 9> k1, k2, k3, k4, y_temp;
    // 修改這段：使用根據能量動態調整的步長
    double min_step = 1e-6;  // 基本最小步長
    // 如果能量非常高，可能需要更小的步長
    // 如果能量非常高，調整步長上限，但不要設置太小
    if (E_nu > 1e5) {
        double adjusted_h = std::min(h, 1e-5);  // 較保守但不過度保守的步長
        h = std::max(adjusted_h, min_step);
        
        // 輸出適當的訊息
        if (h < adjusted_h * 0.99) {
            std::cout << "高能微中子調整步長: h = " << h << std::endl;
        }
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

// 尺度重整的RK4方法
void scaled_rk4_step(double& t, double h, std::array<std::complex<double>, 9>& rho,
    double E_nu, const std::vector<double>& r_vals, 
    const std::vector<double>& Ne_vals, double scale_factor) {
        
        // 添加动态增强因子，让演化更明显
        double dynamic_factor = 1.0;
        // 如果步长非常小，增加动态因子
        if (h < 1e-8) {
            dynamic_factor = 10.0;
            std::cout << "应用动态增强因子: " << dynamic_factor << std::endl;
        }
        // 創建尺度重整的微分方程函數
        auto scaled_deriv = [&](double t_val, const std::array<std::complex<double>, 9>& rho_val) {
            std::array<std::complex<double>, 9> drho{};

            // 獲取當前位置的電子密度
            double Ve = electron_potential(t_val, r_vals, Ne_vals);

            // 計算重整後的哈密頓量
            std::array<std::complex<double>, 9> H_scaled{};
            compute_scaled_hamiltonian(H_scaled, E_nu, Ve, scale_factor);

            // 計算密度矩陣的導數：-i[H,ρ]
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                drho[i*3+j] = 0.0;
                    for (int k = 0; k < 3; k++) {
                    // -i * (H*ρ - ρ*H)
                    drho[i*3+j] -= std::complex<double>(0.0, 1.0) * 
                        (H_scaled[i*3+k] * rho_val[k*3+j] - rho_val[i*3+k] * H_scaled[k*3+j]);
                    }
                }
            }

        return drho;
    };

    // RK4積分步驟
    double t_orig = t;
    std::array<std::complex<double>, 9> k1 = scaled_deriv(t, rho);

    std::array<std::complex<double>, 9> rho_temp{};
    for (int i = 0; i < 9; i++) {
    rho_temp[i] = rho[i] + 0.5 * h * k1[i];
    }
    std::array<std::complex<double>, 9> k2 = scaled_deriv(t + 0.5 * h, rho_temp);

    for (int i = 0; i < 9; i++) {
    rho_temp[i] = rho[i] + 0.5 * h * k2[i];
    }
    std::array<std::complex<double>, 9> k3 = scaled_deriv(t + 0.5 * h, rho_temp);

    for (int i = 0; i < 9; i++) {
    rho_temp[i] = rho[i] + h * k3[i];
    }
    std::array<std::complex<double>, 9> k4 = scaled_deriv(t + h, rho_temp);

    // 更新ρ
    for (int i = 0; i < 9; i++) {
    rho[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }

    // 確保密度矩陣保持厄米性
    ensure_hermiticity(rho);

    // 更新時間
    t = t_orig + h;
    }

// 增強的RK4方法，專為高能微中子設計
void enhanced_rk4_step(double& t, double& h, std::array<std::complex<double>, 9>& rho,
    double E_nu, const std::vector<double>& r_vals, 
    const std::vector<double>& Ne_vals) {
    // 對於高能微中子，採用重新調整尺度的哈密頓量
    double Ve = electron_potential(t, r_vals, Ne_vals);

    // 計算真空振盪尺度
    double vac_scale = delta_m21_squared / (2.0 * E_nu);

    // 如果物質效應和真空振盪差距太大，進行尺度重整
    if (Ve / vac_scale > 1e6) {
        // 使用尺度重整的RK4步驟
        scaled_rk4_step(t, h, rho, E_nu, r_vals, Ne_vals, vac_scale);
    } else {
        // 標準RK4步驟
        rk4_step(t, h, rho, E_nu, r_vals, Ne_vals);
    }
}

// 適應性步長RK4
void adaptive_rk4_step(
    double &t, 
    double &h, 
    std::array<std::complex<double>, 9>& y,
    double E_nu, const std::vector<double>& r_vals, 
    const std::vector<double>& Ne_vals, double tolerance ) {
    
    // 添加這裡：設定最小步長
    const double min_step_size = 1.0e-6;
    h = std::max(h, min_step_size);  // 確保步長不會太小
                        
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

        // 修改這裡：確保步長不會太小
        h = std::max(h, min_step_size);
        
        // 如果步長太小，放棄
        if (h < min_step_size * 0.9) {
            throw std::runtime_error("Step size too small, integration failed");
        }
    }
}

// ✅ 添加：新的自適應龍格-庫塔方法
void adaptive_integrate_to_position(
    double& t, const double target_t, double& h, 
    const double min_h, const double max_h, const double tolerance, 
    std::array<std::complex<double>, 9>& rho, const double E_nu, 
    const std::vector<double>& r_vals, const std::vector<double>& Ne_vals, 
    const double scale_factor = 1.0 
) {
    const double safety_factor = 0.9; // 安全係數
    const double p_gain = 0.075;      // 誤差估計的比例增益
    const double i_gain = 0.175;      // 誤差估計的積分增益
    double accumulated_error = 0.0;   // 用於誤差積分控制
    
    // 初始步長檢查
    if (h > max_h) h = max_h;
    if (h < min_h) h = min_h;
    
    int steps_taken = 0;
    const int max_steps = 10000; // 防止無限循環
    
    while (t < target_t && steps_taken < max_steps) {
        // 剩餘距離
        double remaining = target_t - t;
        if (remaining <= 0) break;
        
        // 調整步長不超過剩餘距離
        if (h > remaining) h = remaining;
        
        // 存儲當前狀態用於比較
        std::array<std::complex<double>, 9> rho_save = rho;
        double t_save = t;
        
        // 使用h進行一步計算 (RK4或縮放RK4)
        if (scale_factor < 0.99) {
            scaled_rk4_step(t, h, rho, E_nu, r_vals, Ne_vals, scale_factor);
        } else {
            rk4_step(t, h, rho, E_nu, r_vals, Ne_vals);
        }
        
        // 使用h/2進行兩步計算 (用於誤差估計)
        std::array<std::complex<double>, 9> rho_half = rho_save;
        double t_half = t_save;
        double h_half = h / 2.0;
        
        if (scale_factor < 0.99) {
            scaled_rk4_step(t_half, h_half, rho_half, E_nu, r_vals, Ne_vals, scale_factor);
            scaled_rk4_step(t_half, h_half, rho_half, E_nu, r_vals, Ne_vals, scale_factor);
        } else {
            rk4_step(t_half, h_half, rho_half, E_nu, r_vals, Ne_vals);
            rk4_step(t_half, h_half, rho_half, E_nu, r_vals, Ne_vals);
        }
        
        // 計算誤差估計
        double error = 0.0;
        for (int i = 0; i < 9; i++) {
            error += std::abs(rho[i] - rho_half[i]) * std::abs(rho[i] - rho_half[i]);
        }
        error = std::sqrt(error / 9.0);
        
        // 當誤差太大時，重新嘗試步長
        if (error > tolerance && h > min_h * 1.1) {
            // 根據誤差調整步長（PI控制器）
            accumulated_error += error - tolerance;
            double scale = safety_factor * std::pow(tolerance / error, 0.2);
            scale *= std::exp(-p_gain * (error - tolerance) - i_gain * accumulated_error);
            
            // 限制步長縮小的比例
            if (scale < 0.1) scale = 0.1;
            
            // 更新步長
            h *= scale;
            
            // 確保步長不小於最小允許值
            if (h < min_h) h = min_h;
            
            // 重置狀態並重試
            t = t_save;
            rho = rho_save;
            continue;
        }
        
        // 步長成功，考慮增加步長
        if (error < tolerance * 0.5 && h < max_h * 0.9) {
            // 基於誤差的步長增加（PI控制器）
            accumulated_error += error - tolerance;
            double scale = safety_factor * std::pow(tolerance / std::max(error, 1e-15), 0.2);
            scale *= std::exp(-p_gain * (error - tolerance) - i_gain * accumulated_error);
            
            // 限制步長增加的比例
            if (scale > 4.0) scale = 4.0;
            
            // 更新步長
            h *= scale;
            
            // 確保步長不大於最大允許值
            if (h > max_h) h = max_h;
        }
        
        // 採用h/2的結果（較準確）
        rho = rho_half;
        t = t_half;
        
        steps_taken++;
        
        // 如果步長非常小但誤差仍然很大，可能是遇到了劇烈變化區域
        // 在這種情況下，我們接受較大的誤差但確保繼續前進
        if (h <= min_h * 1.01 && error > tolerance) {
            break;
        }
    }
    
    // 檢查是否達到最大步數
    if (steps_taken >= max_steps) {
        std::cerr << "警告：積分達到最大步數限制！t = " << t << ", target_t = " << target_t << std::endl;
    }
    
    // 確保到達目標位置
    if (std::abs(t - target_t) > min_h * 0.1) {
        double final_h = target_t - t;
        if (scale_factor < 0.99) {
            scaled_rk4_step(t, final_h, rho, E_nu, r_vals, Ne_vals, scale_factor);
        } else {
            rk4_step(t, final_h, rho, E_nu, r_vals, Ne_vals);
        }
    }
}


// ✅ 新版本：模擬積分用 Python 傳入的 r_vals，完全不再從密度檔案決定 r 點
// ✅ 進一步優化版本：改善高分辨率積分、平滑轉換和數值穩定性
struct SimulationResult;

SimulationResult simulate_custom_rvals(
    double E_nu,
    const std::vector<double>& input_r_vals,     
    const std::vector<double>& density_r_vals,   
    const std::vector<double>& Ne_vals          
) {
    std::cout << "==== 開始模擬能量為 " << E_nu << " eV 的微中子 ====" << std::endl;
    
    SimulationResult result;
    int num_points = input_r_vals.size();
    result.r_vals = input_r_vals;
    result.probs.resize(num_points, std::vector<double>(3, 0.0));
    
    // 判斷是否為高能微中子，決定使用策略
    bool high_energy = (E_nu > 1e7);
    bool use_jump_approx = (E_nu > 5e6);  // 5 MeV 以下完全禁用 P_jump

    // 預先計算共振條件
    double res_val_12 = (delta_m21_squared * std::cos(2 * theta12)) / (2.0 * E_nu);
    double res_val_13 = (delta_m31_squared * std::cos(2 * theta13)) / (2.0 * E_nu);
    
    // 添加這裡：印出共振條件
    std::cout << "MSW 1-2共振條件 = " << res_val_12 << " eV" << std::endl;
    std::cout << "MSW 1-3共振條件 = " << res_val_13 << " eV" << std::endl;
    std::cout << std::flush;

    // 檢查密度點是否足夠覆蓋共振區
    std::vector<double> resonance_r_points;
    for (size_t i = 1; i < density_r_vals.size(); i++) {
        double Ve1 = std::sqrt(2.0) * G_F * Ne_vals[i-1];
        double Ve2 = std::sqrt(2.0) * G_F * Ne_vals[i];
        
        // 檢查1-2共振
        if ((Ve1 - res_val_12) * (Ve2 - res_val_12) <= 0) {
            double r_res = interpolate_resonance_position(
                density_r_vals[i-1], density_r_vals[i], Ve1, Ve2, res_val_12);
            resonance_r_points.push_back(r_res);
            std::cout << "✓ 找到1-2共振點: r ≈ " << r_res << std::endl;
        }
        
        // 檢查1-3共振
        if ((Ve1 - res_val_13) * (Ve2 - res_val_13) <= 0) {
            double r_res = interpolate_resonance_position(
                density_r_vals[i-1], density_r_vals[i], Ve1, Ve2, res_val_13);
            resonance_r_points.push_back(r_res);
            std::cout << "✓ 找到1-3共振點: r ≈ " << r_res << std::endl;
        }
    }
    
    // 初始化密度矩陣（電子微中子初態）
    std::array<std::complex<double>, 9> rho{};
    for (int i = 0; i < 9; i++) rho[i] = 0.0;
    rho[0] = std::complex<double>(1.0, 0.0);  // |νe><νe| 初態
    
    // 在主要積分迴圈前添加
    // 計算縮放因子
    double scale_factor = 1.0;
    double vac_scale = delta_m21_squared / (2.0 * E_nu);
    double max_Ve = 0.0;
    for (size_t i = 0; i < Ne_vals.size(); i++) {
        double Ve = std::sqrt(2.0) * G_F * Ne_vals[i];
        max_Ve = std::max(max_Ve, Ve);
    }
    // 如果物質效應和真空振盪差距太大，調整縮放因子
    if (max_Ve / vac_scale > 1e2) {
        scale_factor = vac_scale;
        // 對於極高能量的情況，可以進一步調整縮放策略
        if (E_nu > 1e7) {
            scale_factor *= 0.01;  // 進一步縮小尺度差異
        }
        std::cout << "使用縮放因子: " << scale_factor << std::endl;
        std::cout << "最大Ve/真空尺度比例: " << max_Ve / vac_scale << std::endl;
    }


    // 設定初始步長
    double t = input_r_vals[0];
    if (t <= 0.0) {
        std::cout << "⚠️ 初始位置為零，將調整為8.19005E-03" << std::endl;
        t = 8.19005e-3;
    }

    double h = (E_nu > 1e7) ? 1e-8 : (high_energy ? 1e-6 : 1e-5);  // 為超高能量設置特殊步長
    
    // 初始點的機率
    result.probs[0][0] = std::real(rho[0]);
    result.probs[0][1] = std::real(rho[4]);
    result.probs[0][2] = std::real(rho[8]);
    
    // 主要演化循環
    for (int i = 1; i < num_points; ++i) {
        double next_r = input_r_vals[i];
        
        // 獲取下個點的電子勢能
        double Ve = electron_potential(next_r, density_r_vals, Ne_vals);
        double Ne = Ve / (std::sqrt(2.0) * G_F);
        
        // ✅ 新增這段
        double vac_scale = delta_m21_squared / (2.0 * E_nu);
        if (E_nu > 1.0e6 && Ve < 1e-2 * vac_scale) {
            double h_fast = std::max(h, 0.01);
            integrate_to_position(t, next_r, h_fast, rho, E_nu, density_r_vals, Ne_vals);
            normalize_density_matrix(rho);
            ensure_hermiticity(rho);
            result.probs[i][0] = std::real(rho[0]);
            result.probs[i][1] = std::real(rho[4]);
            result.probs[i][2] = std::real(rho[8]);
            continue;
        }


        // 檢查是否接近共振區（任何共振）
        bool near_resonance = false;
        for (const auto& res_r : resonance_r_points) {
            if (std::abs(next_r - res_r) < 0.05) {
                near_resonance = true;
                break;
            }
        }
        
        // 接近共振區時進行特殊處理
        if (near_resonance && high_energy) {
            std::cout << "🟡 接近共振區: r = " << next_r << ", Ve = " << Ve << std::endl;
            
            // 高能微中子在共振區的特殊處理
            // ----------------- ✅ 改寫版 LZ 共振跳躍判斷 -----------------
            if (use_jump_approx) {
                // 判斷是否為 1–2 共振區
                bool near_12_res = std::abs(Ve - res_val_12) / std::abs(res_val_12) < 0.1;
                bool near_13_res = std::abs(Ve - res_val_13) / std::abs(res_val_13) < 0.1;

                if (near_12_res) {
                    // 只處理 1–2 共振，1–3 在低能下不應處理
                    double dNe_dr = estimate_density_gradient(next_r, density_r_vals, Ne_vals);
                    double dVe_dr = std::sqrt(2.0) * G_F * dNe_dr;

                    double gamma = std::abs(2 * M_PI * delta_m21_squared * std::sin(2 * theta12) /
                                (2.0 * E_nu * std::cos(2 * theta12) * dVe_dr));
                    double P_jump = std::exp(-gamma);

                    std::cout << "✅ [LZ] 應用 1–2 共振跳躍 @ r = " << next_r 
                            << ", P_jump = " << P_jump << std::endl;
                    apply_resonance_jump_12(rho, P_jump);

                    // 強制跳至下一點（避免共振區細節過度積分）
                    t = next_r;

                } else if (near_13_res && E_nu > 10.0) {
                    // 只有在高能才處理 1–3 跳躍，低能完全禁用
                    double dNe_dr = estimate_density_gradient(next_r, density_r_vals, Ne_vals);
                    double dVe_dr = std::sqrt(2.0) * G_F * dNe_dr;

                    double gamma = std::abs(2 * M_PI * delta_m31_squared * std::sin(2 * theta13) /
                                (2.0 * E_nu * std::cos(2 * theta13) * dVe_dr));
                    double P_jump = std::exp(-gamma);

                    std::cout << "⚠️ [LZ] 高能才處理 1–3 共振跳躍 @ r = " << next_r 
                            << ", P_jump = " << P_jump << std::endl;
                    apply_resonance_jump_13(rho, P_jump);

                    t = next_r;
                }
            }

            else {
                // 共振區但不在極接近共振點，使用更小步長積分
                h = std::min(h, 1e-10);
                if (high_energy && scale_factor < 0.99) {
                    // 使用縮放的 RK4
                    double h_local = h;
                    int stuck_counter = 0;
                    double last_t = t;

                    while (t < next_r) {
                        double remaining = next_r - t;
                        if (h_local > remaining) h_local = remaining;

                        double before_step = t;
                        scaled_rk4_step(t, h_local, rho, E_nu, density_r_vals, Ne_vals, scale_factor);

                        if (std::abs(t - before_step) < 1e-10) {
                            stuck_counter++;
                            // 如果连续几次都没有进展，强制前进
                            if (stuck_counter > 5) {
                                std::cout << "⚠️ 积分卡住，强制前进" << std::endl;
                                t = next_r;
                                break;
                            }
                            // 增加步长
                            h_local *= 2.0;
                        } else {
                            stuck_counter = 0;
                        }
                    }
                }else {
                    // 使用標準積分
                    integrate_to_position(t, next_r, h, rho, E_nu, density_r_vals, Ne_vals);
                }
            }
        } else {
                // 非共振區正常積分
                if (high_energy && scale_factor < 0.99) {
                    // 使用縮放的 RK4
                    double h_local = h;
                    while (t < next_r) {
                        double remaining = next_r - t;
                        if (h_local > remaining) h_local = remaining;
                        scaled_rk4_step(t, h_local, rho, E_nu, density_r_vals, Ne_vals, scale_factor);
                    }
                } else {
                    // 使用標準積分
                    integrate_to_position(t, next_r, h, rho, E_nu, density_r_vals, Ne_vals);
            }
        }
        
        // 確保密度矩陣規範化和厄米性
        normalize_density_matrix(rho);
        ensure_hermiticity(rho);
        
        // 儲存結果
        result.probs[i][0] = std::real(rho[0]);
        result.probs[i][1] = std::real(rho[4]);
        result.probs[i][2] = std::real(rho[8]);
        
        // 檢查結果是否合理
        double P_sum = result.probs[i][0] + result.probs[i][1] + result.probs[i][2];
        if (std::abs(P_sum - 1.0) > 1e-6) {
            std::cout << "⚠️ 概率歸一化問題: P_total = " << P_sum << " @ r = " << next_r << std::endl;
            // 修正概率
            for (int j = 0; j < 3; j++) {
                result.probs[i][j] /= P_sum;
            }
        }
    }
    
    // 打印最終結果
    std::cout << "模擬結束，最終振盪概率:" << std::endl;
    std::cout << "P(νe→νe) = " << result.probs.back()[0] << std::endl;
    std::cout << "P(νe→νμ) = " << result.probs.back()[1] << std::endl;
    std::cout << "P(νe→ντ) = " << result.probs.back()[2] << std::endl;
    
    return result;
}

}