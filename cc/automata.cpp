#include <iostream>
#include "automata.hpp"

using namespace std;

// FUNCTIONS CLASSES IMPLEMENTATION
template<typename T>
int CellularAutomata<T>::print_state(){
	for(int i = 0; i < L_grid; i++){
		for(int j = 0; j < L_grid; j++){
			cout << cells_state[i][j] << " ";
		}
		cout << endl;
	}

	return 0;
}

template<typename T>
void TumorAutomata<T>::start_evolution(int i, int j){
	cells_state[i][j] = 1;
}
