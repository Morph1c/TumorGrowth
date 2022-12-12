/*
* implement a grid solver for reactoion-diffusion 
* equation through finite difference method
*/

class GridSolver{
	private:
		double delta_x;
		int grid_len;
		vecto<vector<double>> mesh
		int initialized;
	public:
		GridSolver(double size, int len);
		void set_initial_condition();
		void solve();
}
