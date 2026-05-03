from gfi_analysis import SharedRegistry, GFILoader, GFISimulator, GFIAnalyzer, GFIPlotter

def main():
    registry = SharedRegistry()
    loader = GFILoader(registry)
    simulator = GFISimulator(registry)
    analyzer = GFIAnalyzer(registry)
    plotter = GFIPlotter(registry)

    print("--- Q9: Weak-Grid Stability (High Line Impedance) ---")
    print(loader.load_parameters("params.json"))
    print(simulator.set_mode("PQ"))
    
    # Base parameters for R_line_pu and X_line_pu usually loaded from JSON
    # For testing, we just update them manually. Let's multiply them by 5.
    p = registry.get("gfi:params")
    print(simulator.set_parameter("R_line_pu", p["R_line_pu"] * 5))
    print(simulator.set_parameter("X_line_pu", p["X_line_pu"] * 5))
    
    print(simulator.set_disturbance(phase_jump_angle=0.1))
    print(simulator.run_simulation(method="full"))
    
    print(analyzer.check_stability())
    print(analyzer.get_settling_time())
    print(plotter.plot_results(save_path="main_9_results.png"))

if __name__ == "__main__":
    main()
