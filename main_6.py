from gfi_analysis import SharedRegistry, GFILoader, GFISimulator, GFIAnalyzer, GFIPlotter

def main():
    registry = SharedRegistry()
    loader = GFILoader(registry)
    simulator = GFISimulator(registry)
    analyzer = GFIAnalyzer(registry)
    plotter = GFIPlotter(registry)

    print("--- Q6: Volt-Var Deadband Sensitivity ---")
    print(loader.load_parameters("params.json"))
    print(simulator.set_mode("VOLT-VAR", K_droop_v=20.0, Vdb=0.05))
    print(simulator.set_disturbance(Vmag_dist=1.025))
    print(simulator.run_simulation(method="full"))
    
    print(analyzer.check_stability())
    print(analyzer.get_overshoot("Q_expr"))
    print(plotter.plot_results(save_path="main_6_results.png"))

if __name__ == "__main__":
    main()
