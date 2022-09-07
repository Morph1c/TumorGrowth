// constant definition
#define L_grid 10
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

// CLASSES DEFINITION

// Class prototype for normal/death cellular automata
template<typename T>
class CellularAutomata{
	private:
		T cells_state[L_grid][L_grid];
	public:
		CellularAutomata(T states[L_grid][L_grid]){
			for(int i = 0; i < L_grid; i++){
				for(int j = 0; j < L_grid; j++){
					cells_state[i][j] = states[i][j];
				}
			}
		}
		int print_state();

};

// Class prototype for tumor cellular automata that is derived
// from the master class CellularAutomata

template<typename T>
class TumorAutomata: public CellularAutomata<T>{
	private:
		T cells_state[L_grid][L_grid];
	public:
		TumorAutomata(T states[L_grid][L_grid]): CellularAutomata<T>(states){
			for(int i = 0; i < L_grid; i++){
				for(int j = 0; j < L_grid; j++){
					cells_state[i][j] = states[i][j];
				}
			}
		}
		void start_evolution(int i, int j);
		bool is_intern(int k, int i, int j);
		int choose_border_migration_cells(int k, int i, int j);
		void cellular_division(int k, int i, int j, int cancer_cells[L_grid][L_grid],  int normal_cells[L_grid], int death_cells[L_grid][L_grid]);
		void cellular_mytosis(int k, int i, int j, int cancer_cells[L_grid][L_grid], int necrotic_cells[L_grid][L_grid]);

};


