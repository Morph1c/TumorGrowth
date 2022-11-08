#include <iostream>
#include "automata.hpp"

using namespace std;

// FUNCTIONS CLASSES IMPLEMENTATION
// case with no input, initialize with null arrays
CellularAutomata::CellularAutomata(){
	for(int i = 0; i < L_grid, i++){
		vector<int> temp;
		for(int j = 0, j < L_grid; j++){
			temp.push_back(0);
		}
		cells_state.push_back(temp);
	}
}

// IMPLEMENT AUTOMATIC SEED INITIALIZATION THROUGH RECTANGLE
// case with initial seed:
// (i, j) central position
// (2*seed_height + 1, 2*seed_len + 1) size of tumor rectangle

CellularAutomata:CellularAutomata(int i, int j, int seed_height, int seed_len){
	for(int k = 0; k < L_grid, k++){
		vector<int> temp;
		for(int l = 0, k < L_grid; k++){
			if(k == i && l == j){
				temp.push_back(1);
			}
		}
	}
}

// Change value of a cell
void CellularAutomata::change_value(int i, int j, int k){
	cells_state[i][j] = 1; 
}
