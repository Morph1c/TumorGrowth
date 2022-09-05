/*
* porting di cs2.py in C++ per migliorarne 
* le prestazioni di calcolo
* dopo la risoluzione numerica del modello i dati verranno comunque
* passati ad un file python per il post-processing dei dati in una gif
* author: Morph1c
* date: 04/09/2022
*/

#include <iostream>
#include <vector>
#include <math.h>
#include </Users/andrea/desktop/ftensor/src/FTensor.hpp>

// constant definition
#define L_grid 50
#define alpha 3/L_grid
#define MAX_ITER_TIME 10
#define V 1
#define delta_x 1
#define delta_t 0.02
#define delta_s 1
#define lambda_N 150
#define lambda_M 10
#define theta_div 0.3
#define theta_del 0.01

using namespace std;

class CellularAutomata{
	private:
		int grid_size;
	public:
		CellularAutomata(int lenght_square){
			grid_size = lenght_square;
		}
		void start_evolution(int i, int j, vector<vector<vector<double> > > cancer_cells);
		bool is_intern(int k, int i, int j);
		int choose_border_migration_cells(int k, int i, int j);
		void cellular_division(int k, int i, int j, vector<vector<vector<double> > > cancer_cells,  vector<vector<vector<double> > > normal_cells, vector<vector<vector<double> > > death_cells);
		void cellular_mytosis(int k, int i, int j, vector<vector<vector<double> > > cancer_cells, vector<vector<vector<double> > > necrotic_cells);

};

// FUNCTION DEFINITON
double k_N(int k, int i, int j, FTensor::Tensor3<double, MAX_ITER_TIME, L_grid, L_grid> cancer_cells, FTensor::Tensor3<double, MAX_ITER_TIME, L_grid, L_grid> normal_cells){

	return (alpha * alpha) * normal_cells(k, i, j) + lambda_N * (alpha * alpha) * cancer_cells(k, i, j);
}

double k_M(int k, int i, int j, FTensor::Tensor3<double, MAX_ITER_TIME, L_grid, L_grid> cancer_cells, FTensor::Tensor3<double, MAX_ITER_TIME, L_grid, L_grid> normal_cells){

	return (alpha * alpha) * normal_cells(k, i, j) + lambda_M * (alpha * alpha) * cancer_cells(k, i, j);
}



int main(void){

	// variabili per la simulazione
	int cells_number = (int) (L_grid / delta_x);
	cout << "len grid: " << cells_number << endl;
	double gamma = (V * delta_t) / (delta_x * delta_x);

	// definisco tensori utili per la simulazione
	FTensor::Tensor3 <double, MAX_ITER_TIME, L_grid, L_grid> cancer_cells;
	FTensor::Tensor3 <double, MAX_ITER_TIME, L_grid, L_grid> normal_cells;
	FTensor::Tensor3 <double, MAX_ITER_TIME, L_grid, L_grid> death_cells;
	FTensor::Tensor3 <double, MAX_ITER_TIME, L_grid, L_grid> N;
	FTensor::Tensor3 <double, MAX_ITER_TIME, L_grid, L_grid> M;



	for(int k = 0; k < MAX_ITER_TIME - 1; k++){
		cout << "Time: " << k << endl;
		for(int i = 1; i < L_grid - 1; i++){
			for(int j = 1; j < L_grid - 1; j++){
				cout << "updating: " << i << " " << j << endl;
				N(k + 1, i, j) = gamma * (N(k, i + 1, j) + N(k, i-1, j) + N(k, i, j + 1) + N(k, i, j-1) - 4*N(k, i, j)) + N(k, i, j) - delta_t*N(k, i, j)*k_N(k, i, j, cancer_cells, normal_cells);
                M(k + 1, i, j) = gamma * (M(k, i + 1, j) + M(k, i-1, j) + M(k, i, j + 1) + M(k, i, j-1) - 4*M(k, i, j)) + M(k, i, j) - delta_t*M(k, i, j)*k_M(k, i, j, cancer_cells, normal_cells);
			}
		}
	}
	
	return 0;
}