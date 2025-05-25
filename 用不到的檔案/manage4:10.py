import numpy as np
import matplotlib.pyplot as plt
import time
import multiprocessing as mp
import os
import neutrino_simulator  # 導入 C++ 擴展模組

try:
    import neutrino_simulator
    print("C++ 模組成功導入")
except ImportError as e:
    print(f"C++ 模組導入失敗：{e}")
    exit()  # 如果導入失敗，直接退出程式

import sys
sys.path.append('/Users/hcjhuang/Documents/各種資料/neutrino_simulator')

def plot_neutrino_oscillations(valid_results):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # 使用更多的绘图样式和颜色
    colors = plt.cm.viridis(np.linspace(0, 1, len(valid_results)))
    #line_styles = ['-', '--', '-.', ':']  # 提供一些线型
    #markers = ['o', 's', '^', 'D']      # 提供一些标记
    
    # 全距离图
    for (E_nu, r_vals, probs), color in zip(valid_results, colors):
        # P(νe -> νe)
        survival_probs = [prob[0] for prob in probs]
        ax1.plot(r_vals, survival_probs,label=f'{E_nu:.2f} MeV (νe -> νe)',color='red',markersize=4,alpha=0.7)
        
        # P(νe -> νμ)
        oscillation_probs_mu = [prob[1] for prob in probs]
        ax1.plot(r_vals, oscillation_probs_mu,label=f'{E_nu:.2f} MeV (νe -> νμ)',color='green',markersize=4,alpha=0.7)
        
        # P(νe -> ντ)
        oscillation_probs_tau = [prob[2] for prob in probs]
        ax1.plot(r_vals, oscillation_probs_tau,  label=f'{E_nu:.2f} MeV (νe -> ντ)', color='blue',markersize=4,alpha=0.7)
    
    ax1.set_xlabel('Distance from the center of the sun to the Earth', fontsize=10)
    ax1.set_ylabel('Neutrino Oscillation Probability', fontsize=10)
    ax1.set_title('Neutrino survival probability from the center of the sun to the Earth', fontsize=12)
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend()
        
    # 太阳内部图
    for (E_nu, r_vals, probs), color in zip(valid_results, colors):
        sun_indices = [i for i, r in enumerate(r_vals) if r <= 1.0]
        if sun_indices:
                r_sun_only = [r_vals[i] for i in sun_indices]
                
                # P(νe -> νe)
                probs_sun_only_ee = [probs[i][0] for i in sun_indices]
                ax2.plot(r_sun_only, probs_sun_only_ee,label=f'{E_nu:.2f} MeV (νe -> νe)',color='red',markersize=4,alpha=0.7)

                # P(νe -> νμ)
                probs_sun_only_mu = [probs[i][1] for i in sun_indices]
                ax2.plot(r_sun_only, probs_sun_only_mu,label=f'{E_nu:.2f} MeV (νe -> νμ)', color='green',markersize=4, alpha=0.7)

                # P(νe -> ντ)
                probs_sun_only_tau = [probs[i][2] for i in sun_indices]
                ax2.plot(r_sun_only, probs_sun_only_tau,label=f'{E_nu:.2f} MeV (νe -> ντ)', color='blue',markersize=4, alpha=0.7)
        
        
        ax2.set_xlabel('Distance from the center of the sun (solar radius)', fontsize=10)
        ax2.set_ylabel('Neutrino Oscillation Probability', fontsize=10)
        ax2.set_title('Neutrino Oscillation Probabilities Inside the Sun', fontsize=12)
        ax2.grid(True, linestyle='--', alpha=0.7)
        ax2.legend()
        
        plt.tight_layout()
        plt.savefig("neutrino_oscillation.png", dpi=300, bbox_inches='tight')
        plt.show()
        print("圖形已保存為 neutrino_oscillation.png")

def simulate_single_energy(args):
    #print(f"simulate_single_energy()函數被呼叫，參數為：{args}")
    """為單一能量模擬中微子振盪"""
    E_nu,density_r_vals, r_vals = args
    print(f"simulate_single_energy()函數被呼叫，參數為：{args}")
    # 將 E_nu 轉換為 Python 的 float 類型
    E_nu = float(E_nu)
    '''
    # 用 C++ 實作的模擬器
    r_start = r_vals_space[0]
    r_end = r_vals_space[-1]
    num_points = len(r_vals_space)
    '''
    
    start_time = time.time()
    try:
        result = neutrino_simulator.simulate_custom_rvals(
            E_nu, 
            r_vals.tolist()
            density_r_vals.tolist(),
            Ne_vals.tolist()
        )
        """
        # 新增：印出詳細的模擬結果
        print("Simulation Result:")
        print(f"  Energy: {E_nu:.2f} MeV")
        print(f"  Distance range: {r_start:.2f} to {r_end:.2f}")
        print(f"  Number of points: {num_points}")
        print(f"  R values: {result.r_vals}")
        print(f"  Probabilities: {result.probs}")  # 印出生存機率
        
        duration = time.time() - start_time
        print(f"Energy {E_nu:.2f} MeV simulation completed in {duration:.2f} seconds")
        """
        return E_nu, result.r_vals, result.probs
    except Exception as e:
        print(f"Energy {E_nu:.2f} MeV simulation failed: {str(e)}")
        return E_nu, None, None


def main():
    # 讀取太陽電子密度數據文件路徑
    base_dir = os.path.dirname(os.path.abspath(__file__))
    density_file = input(f"請輸入太陽電子密度數據文件的路徑 (默認: {base_dir}/BP2000 electron density.txt): ")
    if not density_file:
        density_file = os.path.join(base_dir, "BP2000 electron density.txt")
    
    # 確認文件存在
    if not os.path.exists(density_file):
        print(f"找不到文件：{density_file}")
        return
    
    # 設定中微子能量範圍 (MeV)
    energies_input = input("請輸入中微子能量，以逗號分隔 (默認: 1.0,5.0,10.0 MeV): ")
    if energies_input:
        try:
            energies = np.array([float(e.strip()) for e in energies_input.split(",")])
        except ValueError:
            print("能量輸入格式錯誤，使用默認值")
            energies = np.array([1.0, 5.0, 10.0])
    else:
        energies = np.array([1.0, 5.0, 10.0])
    
    '''4/8
    # 設定距離範圍
    r_sun = np.linspace(0, 1.00135, 500)  # 太陽內部
    
    # 太陽表面到地球軌道的距離
    r_earth = neutrino_simulator.R_earth_orbit / neutrino_simulator.R_sun
    r_outer = np.linspace(1.00135, r_earth, 10000)
    
    # 合併距離點，並確保數值唯一
    r_vals_space = np.concatenate([r_sun, r_outer])
    r_vals_space = np.sort(np.unique(r_vals_space))
    '''
    # === 讀取密度檔案，並取得太陽內部的距離點 ===
    density_r_vals, Ne_vals=[],[]
    with open(density_file, 'r') as f:
        lines = f.readlines()
    for line in lines[2:]:  # 跳過前兩行
        if not line.strip():
            continue
        if "BP2000" in line:
            continue
        parts = line.strip().split()
        try:
            r_val = float(parts[0])
            if r_val <= 1.0:
                density_r_vals.append(r_val)
        except ValueError:
            continue  # 忽略無法轉換的行
            
    '''
    #density_r_vals = np.array([r for r in density_r_vals if r <= 1.0])  # 太陽內部定義為 r <= 1.0

    print(f"密度檔案中太陽內部點數：{len(density_r_vals)}")

    # === 太陽外部：從 1.00135 太陽半徑 到 地球軌道之間建立距離點 ===
    r_earth = neutrino_simulator.R_earth_orbit / neutrino_simulator.R_sun
    r_outer_full = np.linspace(1.00135, r_earth, 10000)  # 原本精細的點

    # === Subsample 外部的點：每 step 個取一個，減少數值積分開銷 ===
    outer_step = 5  # 可調整：值越大，點越少 → 越快；越小，點越多 → 越精細
    r_outer_subsampled = r_outer_full[::outer_step]
    print(f"太陽外部 subsampled 點數：{len(r_outer_subsampled)}")

    # === 合併所有距離點 ===
    r_vals_space = np.concatenate([density_r_vals, r_outer_subsampled])
    r_vals_space = np.sort(np.unique(r_vals_space))  # 確保距離是遞增且無重複

    print(f"總距離點數：{len(r_vals_space)}")
    '''
    '''
    # 👉 subsample：每 10 點取一點，加速模擬
    step = 5  # 可以調整成 5、20，看你要多快
    r_vals_space_subsampled = r_vals_space[::step]
    '''
    print(f"模擬 {len(energies)} 個不同能量值的中微子振盪")
    print(f"距離範圍: 0 到 {r_earth:.2f} 太陽半徑")
    
    # === 組合距離點 ===
r_sun = density_r_vals  # 來自檔案，確保太陽內部點精確
r_earth = neutrino_simulator.R_earth_orbit / neutrino_simulator.R_sun
r_outer = np.linspace(1.00135, r_earth, 10000)
r_outer_sub = r_outer[::10]  # 可調整 subsample 精度
r_vals = np.sort(np.unique(np.concatenate([r_sun, r_outer_sub])))

# === 打包模擬參數 ===
simulation_args = [(E, r_vals, np.array(density_r_vals), np.array(Ne_vals)) for E in energies]
    # 準備多處理參數
    num_cores = max(1, mp.cpu_count() - 1)
    print(f"使用 {num_cores} 個 CPU 核心進行並行計算")
    
    # 打包參數
    simulation_args = [(E_nu, density_file, r_vals_space ) for E_nu in energies]
    
    # 使用多處理加速計算
    with mp.Pool(num_cores) as pool:
        results = pool.map(simulate_single_energy, simulation_args)
    print(f"simulation_args 列表：{simulation_args}")

    with mp.Pool(num_cores) as pool:
        results = pool.map(simulate_single_energy, simulation_args)
    print(f"results 列表：{results}")
    # 繪圖
    valid_results = [(E_nu, r_vals, probs) for E_nu, r_vals, probs in results if r_vals is not None]
    
    if valid_results:
        plot_neutrino_oscillations(valid_results) 
    else:
        print("所有模擬都失敗了，請檢查錯誤訊息。")

if __name__ == "__main__":
    main()