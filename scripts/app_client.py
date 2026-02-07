import streamlit as st
import os
import sys
import matplotlib.pyplot as plt
import json
import io
import contextlib

# --- Setup Paths ---
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
src_path = os.path.join(project_root, 'src')
data_raw_path = os.path.join(project_root, 'data', 'raw')

if src_path not in sys.path:
    sys.path.insert(0, src_path)

# --- Imports ---
from genetic_algo.dna.RotTable import RotTable
from genetic_algo.core.algogenetique import AlgoGenetique
from genetic_algo.dna.Traj3D import Traj3D

# --- Page Config ---
st.set_page_config(page_title="Plasmid Optimizer", layout="wide")
st.title('🧬 Plasmid Rotation Optimizer')

# --- Helper: Terminal Capture Class ---
class StreamlitOutput:
    """Captures print() output and updates a Streamlit container."""
    def __init__(self, parent_container):
        self.parent = parent_container
        self.buffer = io.StringIO()
        
    def write(self, text):
        self.buffer.write(text)
        # Update the container with the full content so far
        self.parent.code(self.buffer.getvalue(), language="text")

    def flush(self):
        pass

# --- Initialize Session State ---
# This keeps data alive so moving the slider doesn't restart the algorithm
if 'results' not in st.session_state:
    st.session_state.results = None
if 'logs' not in st.session_state:
    st.session_state.logs = ""

# --- SIDEBAR ---
with st.sidebar:
    st.header("1. Input Data")

    data_source = st.radio("Data Source", ["📁 Select from Database", "⬆️ Upload File"])
    
    seq = None
    filename_display = ""

    if data_source == "📁 Select from Database":
        if os.path.exists(data_raw_path):
            fasta_files = [f for f in os.listdir(data_raw_path) if f.endswith('.fasta')]
            if fasta_files:
                selected_file = st.selectbox("Choose a Plasmid", fasta_files)
                file_path = os.path.join(data_raw_path, selected_file)
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    seq = "".join([line.strip() for line in lines if not line.startswith('>')])
                filename_display = selected_file
            else:
                st.warning("No .fasta files found in 'data/raw'.")
        else:
            st.error(f"Directory not found: {data_raw_path}")

    else: # Upload Mode
        uploaded_file = st.file_uploader("Upload .fasta file", type=["fasta"])
        if uploaded_file is not None:
            stringio = uploaded_file.getvalue().decode("utf-8")
            seq = "".join([line.strip() for line in stringio.splitlines() if not line.startswith('>')])
            filename_display = uploaded_file.name

    st.divider()

    st.header("2. Rotation Table")
    rot_table_dir = os.path.join(src_path, "genetic_algo", "dna")
    try:
        available_tables = [f for f in os.listdir(rot_table_dir) if f.endswith('.json')]
    except FileNotFoundError:
        available_tables = ["table.json"]
    
    selected_table = st.selectbox("Initial Table", available_tables)

    st.divider()
    
    st.header("3. Algorithm Parameters")
    nb_indiv = st.slider("Population Size", 10, 1000, 100, 10)
    nb_gen = st.slider("Generations", 1, 200, 20, 1)
    
    selection_method = st.selectbox("Selection Method", 
        ["elitiste", "tournament", "roulette", "rang_reel", "rang_geo", "roulette_exp"]
    )
    taux_selec = st.slider("Selection Rate", 0.01, 0.99, 0.3)
    
    st.header("4. Fitness Parameters")
    nb_append = st.number_input("Overlap (bases)", min_value=1, value=2)
    nb_cuts = st.number_input("Additional Cuts", min_value=0, value=0)

    run_btn = st.button("🚀 Run Optimization", type="primary")

# --- MAIN LOGIC ---

# 1. RUNNING THE ALGORITHM
if run_btn:
    if seq is None:
        st.error("❌ Please select or upload a FASTA file first.")
    else:
        table_path = os.path.join(rot_table_dir, selected_table)
        
        # Create a placeholder for the "Terminal Output"
        st.subheader("🖥️ Execution Log")
        log_placeholder = st.empty()
        
        # Capture stdout
        capture = StreamlitOutput(log_placeholder)
        
        with st.spinner(f"Running Genetic Algorithm on {filename_display}..."):
            # Redirect stdout to our capturer
            with contextlib.redirect_stdout(capture):
                Best_indiv_list, Best_scores, Worst_scores = AlgoGenetique(
                    table_path, 
                    seq, 
                    nb_individus=nb_indiv, 
                    nb_generations=nb_gen, 
                    taux_selec=taux_selec, 
                    selection_type=selection_method, 
                    nb_cuts=nb_cuts, 
                    nb_append=nb_append
                )
        
        # Save logs and results to Session State (so they persist)
        st.session_state.logs = capture.buffer.getvalue()
        st.session_state.results = {
            "best_list": Best_indiv_list,
            "best_scores": Best_scores,
            "worst_scores": Worst_scores,
            "seq": seq,
            "filename": filename_display,
            "method": selection_method
        }
        st.success("Optimization Complete!")

# 2. DISPLAYING RESULTS (If they exist in memory)
if st.session_state.results is not None:
    res = st.session_state.results
    
    # Show the saved logs in an expander
    with st.expander("Show Execution Logs (Terminal Output)", expanded=False):
        st.code(st.session_state.logs, language="text")

    # --- TABS ---
    tab1, tab2, tab3 = st.tabs(["📈 Fitness Evolution", "🧬 3D Trajectory", "💾 Data"])

    with tab1:
        fig, ax = plt.subplots()
        ax.plot(res["best_scores"], label="Best Score", color='green')
        ax.plot(res["worst_scores"], label="Worst Score", linestyle="--", color='red')
        ax.set_xlabel("Generation")
        ax.set_ylabel("Fitness Score")
        ax.set_yscale("log")
        ax.legend()
        ax.grid(True, alpha=0.3)
        st.pyplot(fig)

    with tab2:
        st.markdown("### Evolution of the Structure")
        st.caption("Use the slider to see the best shape found at each generation.")
        
        # --- NEW LOGIC: Slider selects GENERATION, not STEP ---
        num_generations = len(res["best_list"]) - 1
        
        # Slider for Generation
        gen_idx = st.slider("Generation", 0, num_generations, num_generations)
        
        # Retrieve the best individual for THIS specific generation
        current_best_indiv = res["best_list"][gen_idx]
        
        # Compute trajectory for this individual
        traj = Traj3D()
        traj.compute(res["seq"], current_best_indiv.Rot_table)
        xyz = traj.getTraj()
        
        # Get Score
        score = res["best_scores"][gen_idx]
        st.info(f"🧬 Generation: **{gen_idx}** | Score: **{score:.4f}**")
        
        # Plot
        fig3d = plt.figure(figsize=(10, 8))
        ax3d = fig3d.add_subplot(111, projection='3d')
        
        ax3d.plot(xyz[:,0], xyz[:,1], xyz[:,2], linewidth=2, color='blue')
        ax3d.scatter(xyz[0,0], xyz[0,1], xyz[0,2], color='green', s=50, label="Start")
        ax3d.scatter(xyz[-1,0], xyz[-1,1], xyz[-1,2], color='red', s=50, label="End")

        # To keep the camera stable while moving the slider, we set fixed limits
        # based on the FINAL best individual (or you could calculate global max/min)
        final_indiv = res["best_list"][-1]
        traj_final = Traj3D()
        traj_final.compute(res["seq"], final_indiv.Rot_table)
        xyz_final = traj_final.getTraj()
        
        all_x, all_y, all_z = xyz_final[:,0], xyz_final[:,1], xyz_final[:,2]
        pad = 50 # padding
        ax3d.set_xlim(all_x.min()-pad, all_x.max()+pad)
        ax3d.set_ylim(all_y.min()-pad, all_y.max()+pad)
        ax3d.set_zlim(all_z.min()-pad, all_z.max()+pad)
        
        ax3d.legend()
        st.pyplot(fig3d)

    with tab3:
        st.write(f"**Final Best Score:** {res['best_scores'][-1]}")
        
        best_candidate = res["best_list"][-1]
        json_str = json.dumps(best_candidate.Rot_table.rot_table, indent=4)
        
        st.download_button(
            label="⬇️ Download Optimal Rotation Table (.json)",
            data=json_str,
            file_name=f"optimized_{res['filename']}_{res['method']}.json",
            mime="application/json"
        )