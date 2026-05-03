from gfi_analysis import SharedRegistry, GFILoader, GFISimulator, GFIAnalyzer, GFIPlotter

def main():
    registry = SharedRegistry()
    loader = GFILoader(registry)
    simulator = GFISimulator(registry)
    analyzer = GFIAnalyzer(registry)
    plotter = GFIPlotter(registry)

    print("--- Q3: Frequency Droop Coefficient Effect ---")
    print(loader.load_parameters("params.json"))
    print(simulator.set_mode("FS", K_droop_f=40.0, f_db=0.0005))
    print(simulator.run_simulation(method="full"))
    
    print(analyzer.check_stability())
    print(analyzer.get_overshoot("P_expr"))
    print(plotter.plot_results(save_path="main_3_results.png"))

if __name__ == "__main__":
    main()
