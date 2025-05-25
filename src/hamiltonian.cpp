#include "../include/hamiltonian.hpp"
#include "../include/constants.hpp"
#include "../include/pmns_matrix.hpp"
#include "../include/electron_density.hpp"
#include <complex>
#include <iostream>
#include <fstream>

namespace neutrino {

ComplexMatrix hamiltonian(double r, 
                         double E_nu, 
                         const std::vector<double>& r_vals, 
                         const std::vector<double>& Ne_vals) {
    

    // 真空項
    ComplexMatrix H_vac = {{
        {{0.0, 0.0, 0.0}},
        {{0.0, delta_m21_squared , 0.0}},
        {{0.0, 0.0, delta_m31_squared}}
    }};

    // PMNS 矩陣
    ComplexMatrix U = get_pmns_matrix();
    /*
    std::cout << "PMNS Matrix:" << std::endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << U[i][j] << " ";
        }
        std::cout << std::endl;
    }*/
    // 轉換到味道基底 (U * H_vac * U^†)
    ComplexMatrix H_vac_mass = {{}};

    // 先計算 H_vac * U^†
    ComplexMatrix H_vac_Uconj = {{}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            H_vac_Uconj[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                H_vac_Uconj[i][j] += H_vac[i][k] * std::conj(U[j][k]);
            }
        }
    }

    // 計算 U * (H_vac * U^†)
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            H_vac_mass[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                H_vac_mass[i][j] += U[i][k] * H_vac_Uconj[k][j];
            }
            // ⬅️ 真正的質量哈密頓量
            H_vac_mass[i][j] /= (2.0 * E_nu);
        }
    }
    double V_e = electron_potential(r, r_vals, Ne_vals);
    // 物質效應項

    // **新增這行來檢查 V_e 是否為 0**
    //std::cout << "V_e at r = " << r << " is " << V_e << std::endl;

    ComplexMatrix H_matter = {{
        {{V_e, 0.0, 0.0}},
        {{0.0, 0.0, 0.0}},
        {{0.0, 0.0, 0.0}}
    }};
    
    // 转换 H_matter 到质量基底
    ComplexMatrix H_matter_mass = {{}};
    /*
    if (r <= 1.0) {  // 太陽內部
        double V_e = electron_potential(r, r_vals, Ne_vals);
        //std::cout << "V_e = " << V_e << std::endl;  // 輸出
        H_matter[0][0] = V_e;
    }
    */
      // 先计算 H_matter * U^†
    ComplexMatrix H_matter_Uconj = {{}};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                H_matter_Uconj[i][j] = 0.0;
                for (int k = 0; k < 3; k++) {
                    H_matter_Uconj[i][j] += H_matter[i][k] * std::conj(U[j][k]);
                }
            }  
        }          
    // 计算 U * (H_matter * U^†)
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                H_matter_mass[i][j] = 0.0;
                for (int k = 0; k < 3; k++) {
                    H_matter_mass[i][j] += U[i][k] * H_matter_mass[k][j];
                }
            }
        }




    // 總哈密頓量 (H_mass + H_matter)
    ComplexMatrix H_total = {{}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            H_total[i][j] = H_vac_mass[i][j] + H_matter[i][j];
            
        }
    }
    /* 
    std::cout << "[DEBUG] E_nu = " << E_nu 
    << ", H_total[0][0] = " << H_total[0][0]
    << ", H_total[1][1] = " << H_total[1][1]
    << ", H_total[2][2] = " << H_total[2][2]
    << std::endl;
    */
   
    /*
    printMatrix(H_vac, "H_vac");
    printMatrix(H_matter, "H_matter");
    printMatrix(H_total, "H_total");*/
    double mass_term = delta_m21_squared / (2 * E_nu);
    const double resonance_ve = mass_term;
    /*
    if (std::abs(V_e - resonance_ve) / resonance_ve < 0.05) {
        std::cout << "🟡 共振區 r = " << r << std::endl;
        printMatrix(H_total, "H_total");
    }
    if (std::abs(r - 0.547) < 1e-4) {
        std::cout << "🔥 進入 MSW resonance 點 r ≈ 0.547" << std::endl;
        printMatrix(H_total, "H_total");
    }*/
    
    
    
    double res_condition = mass_term * cos(2*theta12);

    const double c12 = std::cos(theta12);
    // 檢查物質效應是否足夠強烈（高能情況下）
    double MSW_resonance = mass_term * c12 * c12;
    //std::cout << "r = " << r << ", Ve = " << V_e << " eV, MSW 共振條件 = " << MSW_resonance << " eV" << std::endl;
    
    // 測試: 如果能量很高，手動提高 Ve 來看到效果 (僅用於測試)
    if (E_nu > 1e5 && V_e < MSW_resonance * 0.1) {
        std::cout << "⚠️ 對於高能微中子 E = " << E_nu << " eV，Ve = " << V_e << " 可能太小" << std::endl;
    }
    //std::cout << "MSW resonance occurs when Ve ≈ " << res_condition << " eV" << std::endl;

    //std::cout << "Mass term = " << mass_term << ", V_e = " << V_e << std::endl;
    // 在這裡輸出哈密頓量
    /*
    std::cout << "H_total at r = " << r << ":" << std::endl;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::cout << H_total[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::flush; // 刷新輸出流
    */
    return H_total;
}




void neutrino_evolution(double t, 
                      const std::array<std::complex<double>, 9>& psi,
                      double E_nu,
                      const std::vector<double>& r_vals,
                      const std::vector<double>& Ne_vals,
                      std::array<std::complex<double>, 9>& dpsi_dt) {

    //std::cout << "Entering neutrino_evolution() function at t = " << t << ", E_nu = " << E_nu << std::endl;

    // 將一維數組重新排列為 3x3 矩陣
    ComplexMatrix psi_matrix = {{}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            psi_matrix[i][j] = psi[i * 3 + j];
        }
    }
    
    // 計算哈密頓量
    ComplexMatrix H = hamiltonian(t, E_nu, r_vals, Ne_vals);
    // 驗證哈密頓量非零
    bool H_is_zero = true;
    for (int i = 0; i < 3 && H_is_zero; i++) {
        for (int j = 0; j < 3 && H_is_zero; j++) {
            if (std::abs(H[i][j]) > 1e-10) {
                H_is_zero = false;
            }
        }
    }
    
    if (H_is_zero) {
        std::cout << "警告: 哈密頓量在 r = " << t << " 處接近零!" << std::endl;
    }
    
    /* DEBUG: 印出哈密頓量以診斷問題
    if (std::fmod(t, 0.5) < 0.01) {  // 每隔0.5單位打印一次
        std::cout << "===== 哈密頓量 @ r = " << t << " =====" << std::endl;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                std::cout << H[i][j] << " ";
            }
            std::cout << std::endl;
        }

        // 檢查密度矩陣
        std::cout << "===== 密度矩陣 @ r = " << t << " =====" << std::endl;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                std::cout << psi_matrix[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }*/
    // 🔹 新增：印出 Hamiltonian
    //printMatrix(H, "H (Hamiltonian at t = " + std::to_string(t) + ")");
    // 輸出 psi_matrix
    /*
    std::cout << "psi_matrix:" << std::endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << psi_matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // 輸出 H
    std::cout << "H:" << std::endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << H[i][j] << " ";
        }
        std::cout << std::endl;
    }
        */

    /*計算演化 (-i * H * psi_matrix)
    ComplexMatrix dpsi_dt_matrix = {{}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            dpsi_dt_matrix[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                dpsi_dt_matrix[i][j] -= std::complex<double>(0, 1) * H[i][k] * psi_matrix[k][j];
            }
            //std::cout << "dpsi_dt_matrix[" << i << "][" << j << "] = " << dpsi_dt_matrix[i][j] << std::endl;
        }
    }*/
    // 🔧 修正：改為使用 commutator：-i (Hρ - ρH)
    ComplexMatrix dpsi_dt_matrix = {{}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            std::complex<double> term1 = 0.0;
            std::complex<double> term2 = 0.0;
            for (int k = 0; k < 3; k++) {
                term1 += H[i][k] * psi_matrix[k][j];  // H * ρ
                term2 += psi_matrix[i][k] * H[k][j];  // ρ * H
            }
            dpsi_dt_matrix[i][j] = -std::complex<double>(0.0, 1.0) * (term1 - term2);  // -i[H, ρ]
        }
    }

    // 確認演化不為零
    bool evolution_is_zero = true;
    for (int i = 0; i < 3 && evolution_is_zero; i++) {
        for (int j = 0; j < 3 && evolution_is_zero; j++) {
            if (std::abs(dpsi_dt_matrix[i][j]) > 1e-10) {
                evolution_is_zero = false;
            }
        }
    }
    
    if (evolution_is_zero) {
        std::cout << "警告: 演化在 r = " << t << " 處接近零!" << std::endl;
    }
     
   
     // 🔹 新增：計算電子中微子的生存機率 P_ee
    //double P_ee = std::norm(psi_matrix[0][0]);  

    // 🔹 印出 P_ee 來檢查
    //std::cout << "P_ee at t = " << t << " is " << P_ee << std::endl;
     
    // 將結果轉回一維數組
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            dpsi_dt[i * 3 + j] = dpsi_dt_matrix[i][j];
        }
    }
    
 }

} // namespace neutrino