from gfi_analysis import SharedRegistry, GFILoader, GFISimulator, GFIAnalyzer, GFIPlotter

def main():
    registry = SharedRegistry()
    loader = GFILoader(registry)
    simulator = GFISimulator(registry)
    analyzer = GFIAnalyzer(registry)
    plotter = GFIPlotter(registry)

    print("--- Q10: Combined Phase Jump and Voltage Sag ---")
    print(loader.load_parameters("params.json"))
    print(simulator.set_mode("PQ"))
    print(simulator.set_disturbance(phase_jump_angle=0.3, Vmag_dist=0.92))
    print(simulator.run_simulation(method="full"))
    
    print(analyzer.check_stability())
    print(analyzer.get_settling_time())
    print(analyzer.get_overshoot("P_expr"))
    print(analyzer.get_overshoot("Q_expr"))
    print(plotter.plot_results(save_path="main_10_results.png"))

if __name__ == "__main__":
    main()
