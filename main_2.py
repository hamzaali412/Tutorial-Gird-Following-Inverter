from gfi_analysis import SharedRegistry, GFILoader, GFISimulator, GFIAnalyzer, GFIPlotter

def main():
    registry = SharedRegistry()
    loader = GFILoader(registry)
    simulator = GFISimulator(registry)
    analyzer = GFIAnalyzer(registry)
    plotter = GFIPlotter(registry)

    print("--- Q2: Max Phase Jump Limit ---")
    print(loader.load_parameters("params.json"))
    print(simulator.set_mode("PQ"))
    
    # Test a very large phase jump (e.g. 1.0 rad) to check stability
    print(simulator.set_disturbance(phase_jump_angle=1.0))
    print(simulator.run_simulation(method="full"))
    
    print(analyzer.check_stability())
    print(plotter.plot_results(save_path="main_2_results.png"))

if __name__ == "__main__":
    main()
