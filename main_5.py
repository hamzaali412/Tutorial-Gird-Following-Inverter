from gfi_analysis import SharedRegistry, GFILoader, GFISimulator, GFIAnalyzer, GFIPlotter

def main():
    registry = SharedRegistry()
    loader = GFILoader(registry)
    simulator = GFISimulator(registry)
    analyzer = GFIAnalyzer(registry)
    plotter = GFIPlotter(registry)

    print("--- Q5: Solver Comparison (Full vs Schur) ---")
    print(loader.load_parameters("params.json"))
    
    print("\nRunning Full Newton-Raphson...")
    print(simulator.run_simulation(method="full"))
    
    print("\nRunning Schur Complement...")
    print(simulator.run_simulation(method="schur"))
    
    print("\n" + analyzer.compare_solvers())
    print(plotter.plot_results(save_path="main_5_results.png"))

if __name__ == "__main__":
    main()
