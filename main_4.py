from gfi_analysis import SharedRegistry, GFILoader, GFISimulator, GFIAnalyzer, GFIPlotter

def main():
    registry = SharedRegistry()
    loader = GFILoader(registry)
    simulator = GFISimulator(registry)
    analyzer = GFIAnalyzer(registry)
    plotter = GFIPlotter(registry)

    print("--- Q4: Volt-Var Reactive Power under Voltage Sag ---")
    print(loader.load_parameters("params.json"))
    print(simulator.set_mode("VOLT-VAR", K_droop_v=20.0, Vdb=0.01))
    print(simulator.set_disturbance(Vmag_dist=0.90))
    print(simulator.run_simulation(method="full"))
    
    print(analyzer.check_stability())
    print(analyzer.get_overshoot("Q_expr"))
    print(plotter.plot_results(save_path="main_4_results.png"))

if __name__ == "__main__":
    main()
