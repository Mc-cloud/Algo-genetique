from genetic_algo.core.algogenetique import Individu
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
from matplotlib.widgets import Slider
from genetic_algo.dna.Traj3D import Traj3D
import numpy as np
from .resultsmanager import load_simulation_data

def get_indicators(coords):
    dist = np.linalg.norm(coords[0]-coords[-2])
    v1 = coords[-1]-coords[-2]
    v2 = coords[1] -coords[0]
    v1,v2 = v1/np.linalg.norm(v1),v2/np.linalg.norm(v2)
    return dist,np.linalg.norm(v1-v2),np.dot(v2,v1)

def get_trajectories(indiv_list:list,dna_seq:str):
    trajectories=[]
    traj = Traj3D()
    for indiv in indiv_list:
        traj.compute(dna_seq,indiv.Rot_table)
        xyz = np.array(traj.getTraj())
        trajectories.append(xyz)
    return np.array(trajectories)

def plot_with_slider(trajectories, block=True):
    """Mets de la documentation dans tes fonctions, Clément, je t'en supplie…
    Au moins les plus importantes ! 
    Là j'ai aucune idée de ce qu'est censé être une trajectoire à priori, une instance de Traj3D ?"""
    num_steps = len(trajectories)
    
    # Nombre de points dans une trajectoire (pour aller chercher le dernier)
    num_points = trajectories.shape[1] 

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    plt.subplots_adjust(bottom=0.25) 

    # 1. LA LIGNE PRINCIPALE (Bleue)
    line, = ax.plot(
        trajectories[0, :, 0], 
        trajectories[0, :, 1], 
        trajectories[0, :, 2], 
        color='b', lw=1, alpha=0.6 # Un peu transparent pour voir les points
    )

    # 2. LE POINT DE DÉPART (Vert) - Index 0
    # Note: On utilise le slicing 0:1 pour garder une forme de liste/array [valeur]
    start_pt, = ax.plot(
        trajectories[0, 0:1, 0], 
        trajectories[0, 0:1, 1], 
        trajectories[0, 0:1, 2], 
        color='green', marker='o', markersize=4, linestyle='' 
    )

    # 3. LE POINT D'ARRIVÉE (Rouge) - Index -1 (le dernier)
    end_pt, = ax.plot(
        trajectories[0, -1:, 0], 
        trajectories[0, -1:, 1], 
        trajectories[0, -1:, 2], 
        color='red', marker='o', markersize=4, linestyle='' 
    )

    # Titre initial
    title = ax.set_title(f"Étape : 0")

    # --- FIXER LES AXES ---
    all_x = trajectories[:, :, 0].flatten()
    all_y = trajectories[:, :, 1].flatten()
    all_z = trajectories[:, :, 2].flatten()

    ax.set_xlim(all_x.min(), all_x.max())
    ax.set_ylim(all_y.min(), all_y.max())
    ax.set_zlim(all_z.min(), all_z.max())
    
    # --- Slider ---
    ax_slider = plt.axes([0.2, 0.1, 0.65, 0.03], facecolor='lightgoldenrodyellow')
    slider = Slider(ax_slider, 'Étape', 0, num_steps - 1, valinit=0, valstep=1)

    # --- Fonction de mise à jour ---
    def update(val):
        step = int(slider.val)
        
        line.set_data(trajectories[step, :, 0], trajectories[step, :, 1])
        line.set_3d_properties(trajectories[step, :, 2])
        
        start_pt.set_data(trajectories[step, 0:1, 0], trajectories[step, 0:1, 1])
        start_pt.set_3d_properties(trajectories[step, 0:1, 2])
        
        end_pt.set_data(trajectories[step, -1:, 0], trajectories[step, -1:, 1])
        end_pt.set_3d_properties(trajectories[step, -1:, 2])
        
        title.set_text(f"Étape : {step}")
        fig.canvas.draw_idle()

    slider.on_changed(update)
    plt.show(block=block)

def save_trajectory_gif(trajectories, filename="gifs/etapes.gif", fps=10):
    folder = os.path.dirname(filename)
    if folder and not os.path.exists(folder):
        os.makedirs(folder)
        print(f"Dossier '{folder}' créé.")

    num_steps = len(trajectories)
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # La ligne ADN
    line, = ax.plot([], [], [], color='b', lw=1, alpha=0.6)
    # Point de départ (Vert)
    start_pt, = ax.plot([], [], [], color='green', marker='o', markersize=8, linestyle='')
    # Point d'arrivée (Rouge)
    end_pt, = ax.plot([], [], [], color='red', marker='o', markersize=8, linestyle='')

    title = ax.set_title("Iter 0")

    all_x = trajectories[:, :, 0].flatten()
    all_y = trajectories[:, :, 1].flatten()
    all_z = trajectories[:, :, 2].flatten()

    ax.set_xlim(all_x.min(), all_x.max())
    ax.set_ylim(all_y.min(), all_y.max())
    ax.set_zlim(all_z.min(), all_z.max())


    def update(frame):
        line.set_data(trajectories[frame, :, 0], trajectories[frame, :, 1])
        line.set_3d_properties(trajectories[frame, :, 2])
        
        start_pt.set_data(trajectories[frame, 0:1, 0], trajectories[frame, 0:1, 1])
        start_pt.set_3d_properties(trajectories[frame, 0:1, 2])
        
        end_pt.set_data(trajectories[frame, -1:, 0], trajectories[frame, -1:, 1])
        end_pt.set_3d_properties(trajectories[frame, -1:, 2])
        
        title.set_text(f"Iter {frame}")
        
        return line, start_pt, end_pt, title

    ani = animation.FuncAnimation(
        fig, 
        update, 
        frames=num_steps, 
        interval=1000/fps, 
        blit=False 
    )

    ani.save(filename, writer='pillow', fps=fps)
    
    plt.close(fig) 
    print(f"GIF sauvegardé avec succès : {filename}")

def plot_best_worst(filename,dna_seq,title=""):
    _,best_list,worst_list,_ = load_simulation_data(filename,dna_seq)
    N= [i for i in range(len(best_list))]
    plt.plot(N,best_list, label = "best table")
    plt.plot(N,worst_list, label = "worst table")
    plt.yscale("log")
    plt.xlabel("générations")
    plt.ylabel("score ")
    plt.title(title)
    plt.legend()
    plt.show()

def plot_best_multiple(filenames,dna_seq,labels,title=""):
    if len(labels) != len(filenames):
        labels = [os.path.basename(f) for f in filenames]
    for i in range(len(filenames)):
        filename = filenames[i]
        _,best_list,worst_list,_ = load_simulation_data(filename,dna_seq)
        N= [i for i in range(len(best_list))]
        plt.plot(N,best_list, label = f"best {labels[i]}")
        #plt.plot(N,worst_list, label = f"worst {labels[i]}")
    plt.yscale("log")
    plt.xlabel("générations")
    plt.ylabel("score ")
    plt.title(title)
    plt.legend()
    plt.show()
    
def plot_three_indicators(dist,norm,ps,title=""):
    fig, axs = plt.subplots(3,1, figsize = (12,14), sharex = True)
    generations = range(len(dist))

    axs[0].plot(generations, dist)
    axs[0].set_ylabel('Distance')
    axs[0].set_title('Évolution de la distance de fermeture')
    axs[0].grid(True, linestyle='--')
    # axs[0].legend()

    axs[1].plot(generations, norm)
    axs[1].set_ylabel('$|\\vec{v}{start} - \\vec{v}{end}|$')
    axs[1].set_title('Continuité : Norme de la différence des directions')
    axs[1].grid(True, linestyle='--')

    axs[2].plot(generations,ps)
    axs[2].set_ylabel('$\\vec{v}{start} \\cdot \\vec{v}{end}$')
    axs[2].set_xlabel('Génération')
    axs[2].set_title('Alignement : Produit scalaire des directions')
    axs[2].grid(True, linestyle='--')

    fig.suptitle(title)
    #plt.tight_layout()
    plt.show()
    #plt.savefig('Evolution_metrique_genetique.png')

def plot_multiple_three_indicators(filenames,dna_seq,labels,title=""):
    if len(labels) != len(filenames):
        labels = [os.path.basename(f) for f in filenames]
    fig, axs = plt.subplots(3,1, figsize = (12,14), sharex = True)

    axs[0].set_ylabel('Distance')
    axs[0].set_title('Évolution de la distance de fermeture')
    axs[0].grid(True, linestyle='--')
    # axs[0].legend()

    axs[1].set_ylabel('$|\\vec{v}{start} - \\vec{v}{end}|$')
    axs[1].set_title('Continuité : Norme de la différence des directions')
    axs[1].grid(True, linestyle='--')

    axs[2].set_ylabel('$\\vec{v}{start} \\cdot \\vec{v}{end}$')
    axs[2].set_xlabel('Génération')
    axs[2].set_title('Alignement : Produit scalaire des directions')
    axs[2].grid(True, linestyle='--')

    fig.suptitle(title)
    for i in range(len(filenames)):
        filename = filenames[i]
        indiv_list,_,_,_ = load_simulation_data(filename,dna_seq)
        dist,norm,ps = [],[],[]
        for indiv in indiv_list:
            traj_res = Traj3D()
            traj_res.compute(dna_seq,indiv.Rot_table)
            d,n,p = get_indicators(traj_res.getTraj())
            dist.append(d)
            norm.append(n)
            ps.append(p)
        axs[0].plot(range(len(dist)), dist, label = f"{labels[i]}")
        axs[1].plot(range(len(norm)), norm, label = f"{labels[i]}")
        axs[2].plot(range(len(ps)), ps, label = f"{labels[i]}")
    axs[0].legend()
    axs[1].legend()
    axs[2].legend()
    plt.show()
        
def plot_average_multiple(filename_groups, dna_seq, labels, title=""):
    """
    filename_groups : Liste de listes de noms de fichiers. 
                      Ex: [['simulA_1.pkl', 'simulA_2.pkl'], ['simulB_1.pkl', ...]]
    labels          : Liste des noms pour la légende (un par groupe)
    """
    
    if len(labels) != len(filename_groups):
        print("Attention : nombre de labels différent du nombre de groupes")
    
    # Boucle sur chaque GROUPE de simulations (ex: Modèle A, Modèle B...)
    for i in range(len(filename_groups)):
        group = filename_groups[i]
        all_runs_best = [] 
        
        for filename in group:
            _, best_list, _, _ = load_simulation_data(filename, dna_seq)
            all_runs_best.append(best_list)
        
        min_len = min(len(run) for run in all_runs_best)
        all_runs_best = [run[:min_len] for run in all_runs_best]
        
        arr = np.array(all_runs_best)
        
        # 4. Calcul de la moyenne (axis=0 signifie "colonne par colonne", donc génération par génération)
        mean_curve = np.mean(arr, axis=0)
        
        # 5. Plot
        N = range(len(mean_curve))
        plt.plot(N, mean_curve, label=labels[i])
        

    plt.yscale("log")
    plt.xlabel("Générations")
    plt.ylabel("Score Moyen")
    plt.title(title)
    plt.legend()
    plt.show()