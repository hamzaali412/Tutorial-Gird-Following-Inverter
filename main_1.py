from gfi_analysis import SharedRegistry, GFILoader, GFISimulator, GFIAnalyzer, GFIPlotter

def main():
    registry = SharedRegistry()
    loader = GFILoader(registry)
    simulator = GFISimulator(registry)
    analyzer = GFIAnalyzer(registry)
    plotter = GFIPlotter(registry)

    print("--- Q1: PLL Re-synchronization Settling Time ---")
    print(loader.load_parameters("params.json"))
    print(simulator.set_mode("PQ"))
    print(simulator.set_disturbance(phase_jump_angle=0.4))
    print(simulator.run_simulation(method="full"))
    
    print(analyzer.check_stability())
    print(analyzer.get_settling_time())
    print(plotter.plot_results(save_path="main_1_results.png"))

if __name__ == "__main__":
    main()
