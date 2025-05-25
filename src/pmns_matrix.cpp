#include "../include/pmns_matrix.hpp"
#include "../include/constants.hpp"
#include <cmath>

namespace neutrino {

ComplexMatrix get_pmns_matrix() {
    // 計算三角函數值
    double c12 = std::cos(theta12), s12 = std::sin(theta12);
    double c13 = std::cos(theta13), s13 = std::sin(theta13);
    double c23 = std::cos(theta23), s23 = std::sin(theta23);
    
    // CP違反相位
    std::complex<double> delta_term = std::exp(std::complex<double>(0, delta_cp));
    
    // 三個旋轉矩陣
    ComplexMatrix U3 = {{
        {{c12, s12, 0.0}},
        {{-s12, c12, 0.0}},
        {{0.0, 0.0, 1.0}}
    }};
    
    ComplexMatrix U2 = {{
        {{c13, 0.0, s13 * std::conj(delta_term)}},
        {{0.0, 1.0, 0.0}},
        {{-s13 * delta_term, 0.0, c13}}
    }};
    
    ComplexMatrix U1 = {{
        {{1.0, 0.0, 0.0}},
        {{0.0, c23, s23}},
        {{0.0, -s23, c23}}
    }};
    
    // 計算 U3 * U2
    ComplexMatrix U32 = {{}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            U32[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                U32[i][j] += U3[i][k] * U2[k][j];
            }
        }
    }
    
    // 計算 U3 * U2 * U1
    ComplexMatrix U = {{}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            U[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                U[i][j] += U32[i][k] * U1[k][j];
            }
        }
    }
  
    return U;
}

    void printMatrix(const ComplexMatrix& M, const std::string& name){
        std::cout << name << " = [" << std::endl;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                std::cout << M[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "]" << std::endl;
    }




} // namespace neutrino