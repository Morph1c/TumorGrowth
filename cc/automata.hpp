// constant definition
#define L_grid 10
#define alpha 3/L_grid
#define MAX_ITER_TIME 1000
#define V 1
#define delta_x 1
#define delta_t 0.02
#define delta_s 1
#define lambda_N 150
#define lambda_M 10
#define theta_div 0.3
#define theta_del 0.01

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
#include <vector>
//#include "automata.cpp"

using namespace std;


// CLASSES DEFINITION

// Class prototype for normal/death cellular automata
class CellularAutomata{
	private:
		//int L_grid;
		vector<vector<int>> cells_state;	
	public:
		CellularAutomata(); // constructor for no input
		CellularAutomata(int i, int j);
		void print_state();
		void change_value(int i, int j, int k);

};

// Class prototype for tumor cellular automata that is derived
// from the master class CellularAutomata
/*
class TumorAutomata{
	private:
		int L_grid;
		vector<vector<int>> tumor_state		
	public:
	   TumorAutomata(); // same null constructor
	   TumorAutomata(int i, int j, int seed_height, int seed_len);
	   void start_evolution(int i, int j);
	   bool is_intern(int k, int i, int j);
	   int choose_border_migration_cells(int k, int i, int j);
	   void cellular_division(int k, int i, int j, int cancer_cells[L_grid][L_grid],  int normal_cells[L_grid], int death_cells[L_grid][L_grid]);
	   void cellular_mytosis(int k, int i, int j, int cancer_cells[L_grid][L_grid], int necrotic_cells[L_grid][L_grid]);

};
*/





