##
# SIMULAZIONE NUMERICA DI UN AUTOMA CELLULARE
# CHE MIMA LA CRESCITA TUMORALE 
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from math import exp
#from numba import jit
# Lunghezza rettangolo e numero di iterazioni
L = 50
max_iter_time = 15000# delta_t * max_iter_time = T

v = 1 # velocità propagazione

# parametri di discretizzazione spaziale e temporale
delta_x = 1
delta_t = 0.02
gamma = (v * delta_t) / (delta_x ** 2) # < 0.5 per CFL o 0.25

# fisso la scala evolutiva dell'automa, determino fattore delta_s rispetto a delta_t
delta_s = 1
s = 10

# fisso parametri evolutivi dell'automa
# per valori piccoli di alpha e lambda_N corrispondono ad una richiesta bassa
# di nutrimento e quindi il tumore evolve più lentamente ossia
# esso è meno aggressivo
alpha = 3/L # contenuto nell'intervallo 1/L e 4/L
lambda_N = 150 # contenuto nell'intervallo [0, 200]
lambda_M = 10
theta_div = 0.3
theta_del = 0.01
theta_mov = 3

##
# Inizializzo i valori a 0 di tutte le componenti
# M, N nutrienti principali e secondari
# _cells vari tipi di cellule presenti nella simulazione
M = np.zeros((max_iter_time, L, L), dtype = np.longdouble)
N = np.zeros((max_iter_time, L, L), dtype = np.longdouble)
cancer_cells = np.zeros((max_iter_time, L, L), dtype = np.uint)
normal_cells = np.ones((max_iter_time, L, L), dtype = np.uint) # c'è una cellula normale in ogni sito
necrotic_cells = np.zeros((max_iter_time, L, L), dtype = np.uint)
num_cancer_cells = np.zeros((1, int(max_iter_time / s)))

# condizioni al bordo
top = 1.0 # fisicamente è un apporto di nutrimentp continuo
left = 0.0
bottom = 0.0
right = 0.0

# imposto le condizioni iniziali di Dirichlet ossia un apporto costante di nutrientri 
# nella barra inferiore sia per M che N 
N[:, 0, :] = top
M[:, 0, :] = top
#N[:, L-1, :] = top
#M[:, L-1, :] = top


# parte lineare nelle PDE che viene aggiornato a seconda del tempo e della posizione
def k_N(k, i, j):
    return (alpha**2) * normal_cells[k, i, j] + lambda_N * (alpha**2) * cancer_cells[k, i, j]
def k_M(k, i, j):
    return (alpha**2) * normal_cells[k, i, j] + lambda_M * (alpha**2) * cancer_cells[k, i, j]

# funzioni ausiliari per ricerca di vicini e scelta di un vicino

# Restituisce True se la cellula tumorale è dentro al tumore ossia
# se tutti i suoi 4 vicini(a croce) sono cellule tumorali, se è vicino al bordo
# controlla solo le celle vicini disbonibili
def is_intern(k, i, j):
    for c in range(i - 1, i + 2):
        for l in range(j - 1, j + 2):
            if cancer_cells[k, c, l] == 0:
                return False

    return True

# Restituisce un vicino ammissibile dove non vi è una cellula tumorale
def choose_border_migration_cells(k, i, j):
    possible_position = []
    for c in range(i - 1, i + 2):
        for l in range(j - 1, j + 2):
            if cancer_cells[k, c, l] == 0: # devo aggiungere il caso in cui la cellula è al bordo del dominio e devo escludere alcune celle perchè non definite
                                           # potrebbe dare errore index out of range
                possible_position.append((c, l))
    lun = len(possible_position)
    if lun  > 1:
        choose = np.random.randint(0, lun - 1)
    else:
        choose = 0            
    return possible_position[choose]

# Effettua una divisione cellulare
def cellular_division(k, i, j, cancer_cells, normal_cells, necrotic_cells):
    if is_intern(k, i, j): # la cellula si trova all'interno del tumore
        cancer_cells[k:, i, j] = cancer_cells[k - 1, i, j] + 1
    else: # la cellula si trova al bordo del tumore
        new_i, new_j = choose_border_migration_cells(k, i, j) # ritorna una lista di due elementi
        print("sito scelto per la divisione cellulare scelto (%d, %d) per la cellula tumorale alla pos (%d, %d)" %(new_i, new_j, i, j))
        cancer_cells[k:, new_i, new_j] = 1 #poichè in tale posizione prima non c'era nessuna cellula tumorale
        if normal_cells[k, new_i, new_j] == 1:
            normal_cells[k:, new_i, new_j] = 0
        if necrotic_cells[k, new_i, new_j] == 1:
            necrotic_cells[k:, new_i, new_j] = 0

# Effettua una morte cellulare
def cellular_mytosis(k, i, j, cancer_cells, necrotic_cells):
    cancer_cells[k:, i, j] = cancer_cells[k - 1, i, j] - 1
    if cancer_cells[k, i, j] == 0:
        necrotic_cells[k:, i, j] = 1

# conta tutte le cellule tumorali presenti
def count_cancer_cells(k, cancer_cells):
    num = 0
    for i in range(0, L):
        for j in range(0, L):
            if cancer_cells[k, i, j] >= 1:
                num += 1

    return num


def evolution(N, M, cancer_cells, necrotic_cells):
    p = 0.5 # probabilità che la cellula selezionata compia la divisione o la morte
    #p2 = 0.3
    p_div = 0
    p_del = 0
    p_mov = 0

    for k in range(0, max_iter_time-1, 1): # EVOLUZIONE TEMPORALE, ogni iterazione passa un tempo delta_t
        for i in range(1, L - 1, delta_x):
            for j in range(1, L - 1, delta_x):
                # AGGIORNO I NUTRIENTI AL TEMPO K+1
                #print(k_N(i, j))
                N[k + 1, i, j] = gamma * (N[k, i + 1, j] + N[k, i-1, j] + N[k, i, j + 1] + N[k, i, j-1] - 4*N[k, i, j]) + N[k, i, j] - delta_t*N[k][i][j]*k_N(k, i, j)
                M[k + 1, i, j] = gamma * (M[k, i + 1, j] + M[k, i-1, j] + M[k, i, j + 1] + M[k, i, j-1] - 4*M[k, i, j]) + M[k, i, j] - delta_t*M[k][i][j]*k_M(k, i, j)
        # impongo le condizioni di bordo
        for l in range(0, L, delta_x):
            # condizione periodica sui due assi x di bordo
            N[k + 1, l, 0] = (N[k + 1, l, 1] + N[k + 1, l, L - 2])/2
            N[k + 1, l, L - 1] = N[k + 1, l, 0]

            M[k + 1, l, 0] = (M[k + 1, l, 1] + M[k + 1, l, L - 2])/2
            M[k + 1, l, L - 1] = M[k + 1, l, 0]
            # condizione di neuman derivata in x=L == 0
            N[k + 1, L - 1, l] = N[k + 1, L - 2, l]
            M[k + 1, L - 1, l] = M[k + 1, L - 2, l]
        
        print("iteration %d" %(k))
        if k % s == 0: # ossia sono nella scala temporale evolutiva del CA, che ha scala temporale delta_s / delta_t 
            for i in range(1, L - 1, delta_x): # i bordi del rettangolo non vengono modificati per problemi di aggiornamento regole
                for j in range(1, L - 1, delta_x):
                    if np.random.binomial(1, p, 1) == 1 and cancer_cells[k, i, j] >= 1: # setaccio tutte le cellule tumorali e se ve ne sono in un sito allora provo ad effettuare un evoluzione
                        # la cellula tumorale è stata selezionata per evolvere
                        c = np.random.binomial(1, p, 1)[0] # con probabilità 1/2 scelgo un'azione o l'altra
                        if c == 1: #l'azione selezionata è la divisione cellulare
                            p_div = 1 - exp(-1 * ((N[k, i, j] / (cancer_cells[k, i, j] * theta_div))**2))
                            print("p_div at time %d is %.3f nutriment at (%d, %d) = %.3f" % (k, p_div, i, j, N[k, i, j]))
                            is_divide = np.random.binomial(1, p_div, 1)[0] # non so se è giusto perchè la tratto come una variabile di bernoulli ma scelgo di effettura la divisione con prob p_div
                            if is_divide == 1: # la cellula si divide
                                print("l'automa si divide")
                                cellular_division(k, i, j, cancer_cells, normal_cells, necrotic_cells)
                                print("dopo la divisione cellulare le cellule tumorali sono %d" % count_cancer_cells(k, cancer_cells))

                        elif c == 0: #azione selezionata è morte cellulare
                            p_del = exp(-1 * ((M[k, i, j] / (cancer_cells[k, i, j] * theta_del)) ** 2))
                            print("p_del at time %d is %.3f nutriment at (%d, %d) = %.3f and in this site cancer cells are %d" % (k, p_del, i, j, M[k, i, j], cancer_cells[k, i, j]))
                            is_death = np.random.binomial(1, p_del, 1)[0] # stesso dubbio di p_div
                            if is_death == 1 and count_cancer_cells(k, cancer_cells) > 1: # la seconda condizione è necessaria poichè altrimenti quasi sempre la simulazione inizia
                                                                                          # con la morte dell'unica cellula tumorale
                                print("muore una cellula tumorale") 
                                cellular_mytosis(k, i, j, cancer_cells, necrotic_cells)
                        ## DA IMPLEMENTARE ULTERIORE AZIONE
                        #else: # c == 2 e quindi compio una migrazione cellulare
                        #    p_mov = 1 -  exp(-1 * ((M[k, i, j] / (cancer_cells[k, i, j] * theta_mov)) ** 2))
                        #    is_mov = np.random.binomial(1, p_del, 1)[0]

            # Tengo conto del numero ci cellule tumorali che crescono per vedere se verifica la legge di gompertz
            num_cancer_cells[0, int(k/s) - 1] = count_cancer_cells(k, cancer_cells)

## 
# Routine per il plot della heatmap ad ogni iterazione k
def plotheatmap(u_k, k, m):
    # elimina eventuali altre figure dal plot
    plt.clf()

    plt.title(f"Cancer cells at t = {k*delta_t:.3f} unit time")
    plt.xlabel("x")
    plt.ylabel("y")

    # This is to plot u_k (u at time-step k)
    plt.pcolormesh(u_k, cmap=plt.cm.jet, vmin=0, vmax=m)
    plt.colorbar()

    return plt

# Inizializza l'evoluzione dell'automa cellulare con un seme
# che è la cellula tumorale più o meno al centro
if __name__ == "__main__":
    print("In")
    cancer_cells[0:, int(L/2), int(L/2)] = 1
    #evo_jit = jit(nopython = True)(evolution)
    evolution(N, M, cancer_cells, necrotic_cells)

##
# Plot di controllo su tutte le cellule
#print(cancer_cells[max_iter_time-1, :, :])
#plt.imshow(normal_cells[max_iter_time-1, :, :], interpolation='none')
#plt.title("normal cells")
#plt.show()
#plt.imshow(cancer_cells[max_iter_time-1, :, :], interpolation='none')
#plt.title("Cancer cells")
#plt.show()
#plt.imshow(necrotic_cells[max_iter_time-1, :, :], interpolation='none')
#plt.title("necrotic cells")
#plt.show()

plt.plot(num_cancer_cells[0])
#plt.title("necrotic cells")
plt.show()


# iterazione necessaria per il plot dell'animazione
def animate_nutr(k):
    plotheatmap(M[k, :, :], k, 2)

def animate_cancer(k):
    plotheatmap(cancer_cells[k, :, :], k, 10)

#anim = animation.FuncAnimation(plt.figure(), animate_nutr, interval=20, frames=max_iter_time, repeat=False)
#anim.save("nutr.gif")

anim = animation.FuncAnimation(plt.figure(), animate_nutr, interval=20, frames=max_iter_time, repeat=False)
anim.save("reactdiffTest9.gif")
#anim = animation.FuncAnimation(plt.figure(), animate_cancer, interval=20, frames=max_iter_time, repeat=False)
#anim.save("cancerTest8.gif")
#fig, axes = plt.subplots(ncols=2)
#ax1, ax2 = axes

#im1 = ax1.matshow(cancer_cells[int(max_iter_time / 2), :, :], vmin = 0, vmax = 10, cmap = "viridis")
#im2 = ax2.matshow(cancer_cells[max_iter_time - 1, :, :], vmin = 0, vmax = 10, cmap = "viridis")

#fig.colorbar(im1, ax=ax1)
#fig.colorbar(im2, ax=ax2)
#plt.tight_layout()

#plt.show()

print("Done!")