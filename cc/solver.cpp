/*
* implementing method for gridsolver class
*/

#include "automata.cpp"
#include <iostream>


#define gamma 1


// initialize mesh with 0 everywhere
GridSolver::GridSolver(double size, int len, double time_size){
	delta_x = size;
	delta_t = time_size;
	mesh_size = (int) (len / size); 
	for(int i = 0; i < mesh_size; i++){
		vector<int> temp;
		for(int j = 0; j < mesh_size; j++){
			temp.push_back(0);
		}
		mesh.push_back(temp);
	}
}

// set default initial condition
// now we use dirichlet uniform condition on bottom line
void GridSolver::set_initial_condition(){
	for(int i = 0; i < mesh_size; i++){
		for(int j = 0; j < mesh_size; j++){
			if(i == 0){
				mesh[i][j] = 1
			}
		}
	}

	initialied = 0;
}

// auxiliary function for solver
double k(int i, int j, int type, TumorAutomata cancer_cells, CellularAutomata normal_cells){
	float lambda;
	if(type == 0){
		lambda = lambda_N;
	}
	else{
		lambda = lambda_M;
	}

	return (alpha * alpha) * normal_cells[i][j] + lambda * (alpha * alpha) * cancer_cells[i][j];

}

// this method simulate a PDE reaction-diffusion on a singular delta_t
// in a square grid regular with dirichlet condition
// using finite difference method
void GridSolver::solve(TumorAutomata cancer_cells, CellularAutomata normal_cells){
	for(int i = 1; i < mesh_size - 1; i++){
			for(int j = 1; j < mesh_size - 1; j++){
				//cout << "updating: " << i << " " << j << endl;
				mesh[i][j] = gamma * (mesh[i + 1][j] + mesh[i - 1][j] + mesh[i][j + 1] + mesh[i][j - 1] - 4 * mesh[i][j]) + mesh[i][j] - delta_t * mesh[i][j] * k(i, j, cancer_cells, normal_cells);
			}
		}
}
