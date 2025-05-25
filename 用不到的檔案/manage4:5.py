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
    # 創建用於全距離對數圖的圖表
    plt.figure(figsize=(12, 8))
    ax_full = plt.subplot(111)
    
    # 確保所有能量的結果都有相同的顏色映射
    colors = {'ee': 'red', 'emu': 'green', 'etau': 'blue'}
    
    # 檢查結果數據
    print("有效結果數量:", len(valid_results))
    
    for i, (E_nu, r_vals, probs) in enumerate(valid_results):
        # 先檢查數據
        print(f"能量 {E_nu:.2f} MeV 的數據:")
        print(f"  距離點數量: {len(r_vals)}")
        print(f"  r_vals範圍: {min(r_vals):.6f} 到 {max(r_vals):.6f}")
        print(f"  概率數據點數量: {len(probs)}")
        
        # 確保 r_vals 和 probs 長度相同
        if len(r_vals) != len(probs):
            print(f"警告：距離點數 ({len(r_vals)}) 與概率數據點數 ({len(probs)}) 不匹配！")
            continue
            
        # 確保 r_vals 不含零值（對數尺度需要）
        valid_indices = [i for i, r in enumerate(r_vals) if r > 0]
        if len(valid_indices) < len(r_vals):
            print(f"注意：過濾了 {len(r_vals) - len(valid_indices)} 個零或負值距離點")
        
        r_vals_filtered = [r_vals[i] for i in valid_indices]
        probs_filtered = [probs[i] for i in valid_indices]
        
        if not r_vals_filtered:
            print("警告：過濾後沒有有效的距離點!")
            continue
            
        # 繪製電子中微子存活概率 (νe -> νe)
        survival_probs = [prob[0] for prob in probs_filtered]
        ax_full.plot(r_vals_filtered, survival_probs, color=colors['ee'], linewidth=2, 
                    label=f'{E_nu:.2f} MeV (νe -> νe)')
        
        # 繪製電子中微子到mu中微子的轉換概率 (νe -> νμ)
        oscillation_probs_mu = [prob[1] for prob in probs_filtered]
        ax_full.plot(r_vals_filtered, oscillation_probs_mu, color=colors['emu'], linewidth=2, 
                    label=f'{E_nu:.2f} MeV (νe -> νμ)')
        
        # 繪製電子中微子到tau中微子的轉換概率 (νe -> ντ)
        oscillation_probs_tau = [prob[2] for prob in probs_filtered]
        ax_full.plot(r_vals_filtered, oscillation_probs_tau, color=colors['etau'], linewidth=2, 
                    label=f'{E_nu:.2f} MeV (νe -> ντ)')
    
    # 使用對數尺度以更好地顯示全範圍
    ax_full.set_xscale('log')
    
    # 設定軸範圍和標籤
    r_earth = neutrino_simulator.R_earth_orbit / neutrino_simulator.R_sun
    ax_full.set_xlim(0.01, r_earth)  # 從0.01開始避免對數尺度的0值問題
    ax_full.set_ylim(0, 1.05)        # 確保y軸範圍合適
    
    ax_full.set_xlabel('Distance from the center of the sun (solar radius, log scale)', fontsize=12)
    ax_full.set_ylabel('Neutrino Oscillation Probability', fontsize=12)
    ax_full.set_title('Full Distance Neutrino Oscillation Probabilities (Sun to Earth)', fontsize=14)
    ax_full.grid(True, linestyle='--', alpha=0.7)
    
    # 添加圖例
    ax_full.legend(loc='best')
    
    plt.tight_layout()
    plt.savefig("neutrino_oscillation_full_distance_log.png", dpi=300, bbox_inches='tight')
    plt.show()
    print("全距離對數圖已保存為 neutrino_oscillation_full_distance_log.png")
    
    # 也繪製太陽內部圖 (僅 r_sun，r ≤ 1.0)
    plt.figure(figsize=(12, 6))
    ax_sun = plt.subplot(111)
    
    for i, (E_nu, r_vals, probs) in enumerate(valid_results):
        # 只選取太陽內部的數據點 (r <= 1.0)
        sun_indices = [i for i, r in enumerate(r_vals) if 0 < r <= 1.0]
        
        if sun_indices:
            r_sun_only = [r_vals[i] for i in sun_indices]
            
            # P(νe -> νe)
            probs_sun_only_ee = [probs[i][0] for i in sun_indices]
            ax_sun.plot(r_sun_only, probs_sun_only_ee, color=colors['ee'], linewidth=2, 
                        label=f'{E_nu:.2f} MeV (νe -> νe)')

            # P(νe -> νμ)
            probs_sun_only_mu = [probs[i][1] for i in sun_indices]
            ax_sun.plot(r_sun_only, probs_sun_only_mu, color=colors['emu'], linewidth=2,
                        label=f'{E_nu:.2f} MeV (νe -> νμ)')

            # P(νe -> ντ)
            probs_sun_only_tau = [probs[i][2] for i in sun_indices]
            ax_sun.plot(r_sun_only, probs_sun_only_tau, color=colors['etau'], linewidth=2,
                        label=f'{E_nu:.2f} MeV (νe -> ντ)')
    
    # 明確設定太陽內部圖的 x 軸範圍
    ax_sun.set_xlim(0, 1.0)
    ax_sun.set_ylim(0, 1.05)
    
    ax_sun.set_xlabel('Distance from the center of the sun (solar radius)', fontsize=12)
    ax_sun.set_ylabel('Neutrino Oscillation Probability', fontsize=12)
    ax_sun.set_title('Neutrino Oscillation Probabilities Inside the Sun (r ≤ 1.0)', fontsize=14)
    ax_sun.grid(True, linestyle='--', alpha=0.7)
    
    # 添加圖例
    ax_sun.legend(loc='best')
    
    plt.tight_layout()
    plt.savefig("neutrino_oscillation_sun_interior.png", dpi=300, bbox_inches='tight')
    plt.show()
    print("太陽內部圖已保存為 neutrino_oscillation_sun_interior.png")
    
    # 額外添加：檢查太陽外部的數據
    print("\n檢查太陽外部數據:")
    for i, (E_nu, r_vals, probs) in enumerate(valid_results):
        # 只選取太陽外部的數據點 (r > 1.0)
        out_indices = [i for i, r in enumerate(r_vals) if r > 1.0]
        print(f"能量 {E_nu:.2f} MeV，太陽外部數據點數量: {len(out_indices)}")
        if out_indices:
            max_r = max([r_vals[i] for i in out_indices])
            print(f"  最大距離: {max_r:.2f} 太陽半徑")

def simulate_single_energy(args):
    """為單一能量模擬中微子振盪"""
    E_nu, density_file, r_vals_space = args
    print(f"simulate_single_energy()函數被呼叫，參數為：{args}")
    # 將 E_nu 轉換為 Python 的 float 類型
    E_nu = float(E_nu)

    # 使用 C++ 實作的模擬器
    r_start = r_vals_space[0]
    r_end = r_vals_space[-1]
    num_points = len(r_vals_space)
    
    start_time = time.time()
    try:
        result = neutrino_simulator.simulate_oscillation(
            E_nu, 
            density_file, 
            r_start, 
            r_end, 
            num_points
        )
        # 新增：印出詳細的模擬結果
        print("Simulation Result:")
        print(f"  Energy: {E_nu:.2f} MeV")
        print(f"  Distance range: {r_start:.2f} to {r_end:.2f}")
        print(f"  Number of points: {num_points}")
        print(f"  R values: {result.r_vals}")
        print(f"  Probabilities: {result.probs}")  # 印出生存機率
        
        duration = time.time() - start_time
        print(f"Energy {E_nu:.2f} MeV simulation completed in {duration:.2f} seconds")
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
    
    # 設定距離範圍
    r_sun = np.linspace(0, 1.00135, 500)  # 太陽內部
    
    # 太陽表面到地球軌道的距離
    r_earth = neutrino_simulator.R_earth_orbit / neutrino_simulator.R_sun
    r_outer = np.linspace(1.00135, r_earth, 1000)
    
    # 合併距離點，並確保數值唯一
    r_vals_space = np.concatenate([r_sun, r_outer])
    r_vals_space = np.sort(np.unique(r_vals_space))
    
    print(f"模擬 {len(energies)} 個不同能量值的中微子振盪")
    print(f"距離範圍: 0 到 {r_earth:.2f} 太陽半徑")
    
    # 準備多處理參數
    num_cores = max(1, mp.cpu_count() - 1)
    print(f"使用 {num_cores} 個 CPU 核心進行並行計算")
    
    # 打包參數
    simulation_args = [(E_nu, density_file, r_vals_space) for E_nu in energies]
    
    # 使用多處理加速計算
    # 注意：這裡有一個重複計算，我們刪除了第一次 pool.map 調用，保留第二次
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