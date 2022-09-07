/*
* porting di cs2.py in C++ per migliorarne 
* le prestazioni di calcolo
* dopo la risoluzione numerica del modello i dati verranno comunque
* passati ad un file python per il post-processing dei dati in una gif
* author: Morph1c
* date: 04/09/2022
*/

#include <iostream>
#include <math.h>
#include <fstream>
#include "automata.cpp"

using namespace std;

// FUNCTION DEFINITON
double k_N(int i, int j, int cancer_cells[L_grid][L_grid], int normal_cells[L_grid][L_grid]){

	return (alpha * alpha) * normal_cells[i][j] + lambda_N * (alpha * alpha) * cancer_cells[i][j];
}

double k_M(int i, int j,  int cancer_cells[L_grid][L_grid], int normal_cells[L_grid][L_grid]){

	return (alpha * alpha) * normal_cells[i][j] + lambda_M * (alpha * alpha) * cancer_cells[i][j];
}

void write_array(double array[L_grid][L_grid]){
	ofstream myfile;
  	myfile.open ("example.txt", fstream::app);
	for(int i = 0; i < L_grid; i++){
		for(int j = 0; j < L_grid; j++){
			myfile << array[i][j] << " ";
		}
		myfile << endl;
	}
	myfile << endl;
	myfile.close();
}

void fill_number(int array[L_grid][L_grid], int value){
	for(int i = 0; i < L_grid; i++){
		for(int j = 0; j < L_grid; j++){
			array[i][j] = value;
		}
	}
}

int main(void){

	// variabili per la simulazione
	int cells_number = (int) (L_grid / delta_x);
	cout << "len grid: " << cells_number << endl;
	double gamma = (V * delta_t) / (delta_x * delta_x);

	// definisco tensori utili per la simulazione
	int cancer_cells[L_grid][L_grid];
	int normal_cells[L_grid][L_grid];
	int death_cells[L_grid][L_grid];
	fill_number(cancer_cells, 0);
	fill_number(normal_cells, 0);
	fill_number(death_cells, 0);

	double N[L_grid][L_grid];
	double M[L_grid][L_grid];

	TumorAutomata<int> A(cancer_cells);
	A.print_state();
	cout << "after" << endl;
	A.start_evolution(5, 5);
	A.print_state();

	// condizioni iniziali nutrimento
	for(int i = 0; i < L_grid; i++){
		for(int j = 0; j < L_grid; j++){
			if(i == 0){
				N[i][j] = 1;
				M[i][j] = 1;
			}
			else{
				N[i][j] = 0;
				M[i][j] = 0;
			}
		}
	}
	

	write_array(N);
	
	for(int k = 1; k < MAX_ITER_TIME - 1; k++){
		//cout << "Time: " << k << endl;
		for(int i = 1; i < L_grid - 1; i++){
			for(int j = 1; j < L_grid - 1; j++){
				//cout << "updating: " << i << " " << j << endl;
				N[i][j] = gamma * (N[i + 1][j] + N[i - 1][j] + N[i][j + 1] + N[i][j - 1] - 4 * N[i][j]) + N[i][j] - delta_t * N[i][j] * k_N(i, j, cancer_cells, normal_cells);
                M[i][j] = gamma * (M[i + 1][j] + M[i - 1][j] + M[i][j + 1] + M[i][j - 1] - 4 * M[i][j]) + M[i][j] - delta_t * M[i][j] * k_M(i, j, cancer_cells, normal_cells);
			}
		}
		write_array(N);
	}
	
	return 0;
}