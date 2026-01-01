# ==========================================
# requirements
# ==========================================
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.stats import norm
from scipy.stats.qmc import LatinHypercube
from dataclasses import dataclass
import time

# ==========================================
# 0. pagina iniziale
# ==========================================
st.set_page_config(
    page_title="Buffer Calculator",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# SETUP GRAFICO MATPLOTLIB
plt.style.use('dark_background')
plt.rcParams.update({
    'font.size': 8,
    'figure.facecolor': '#000000',
    'axes.facecolor': '#000000',
    'savefig.facecolor': '#000000',
    'text.color': 'white',
    'axes.labelcolor': 'white',
    'xtick.color': 'white',
    'ytick.color': 'white',
    'grid.color': '#333333',
    'grid.alpha': 0.4,
    'axes.spines.top': False,
    'axes.spines.right': False
})

# CSS PERSONALIZZATO
st.markdown("""
<style>
    .stApp { background-color: #000000 !important; color: #E0E0E0 !important; }
    [data-testid="stSidebar"] { background-color: #111111 !important; border-right: 1px solid #333; }
    .stNumberInput input, .stTextInput input, .stSelectbox div, div[data-baseweb="select"] { 
        background-color: #222 !important; color: white !important; border-color: #444 !important;
    }
    h1, h2, h3, h4, h5, h6, p, label, span, li, div[data-testid="stMarkdownContainer"] > p { color: #FFFFFF !important; }
    [data-testid="stExpander"] { background-color: #1A1A1A !important; border: 1px solid #333; }
    div[data-testid="stStatusWidget"] { background-color: #222 !important; }
    div[data-testid="stDataFrame"] { background-color: #111 !important; }
    [data-testid="stMetricValue"] { color: #00E5FF !important; }
    [data-testid="stMetricLabel"] { color: #AAAAAA !important; }
    button[data-baseweb="tab"] { background-color: #111 !important; color: #aaa !important; border-bottom: 2px solid #333 !important; }
    button[data-baseweb="tab"][aria-selected="true"] { color: #00E5FF !important; border-bottom: 2px solid #00E5FF !important;}
</style>
""", unsafe_allow_html=True)

# ==========================================
# 1. Calcoli classici e stocastici
# ==========================================

@dataclass
class VariableConfig:
    mean: float
    std: float
    name: str

@dataclass
class SimulationParams:
    mass_NaOAc: VariableConfig
    V_tot: VariableConfig
    Conc_HCl: VariableConfig
    Vol_HCl: VariableConfig
    pH_measured: VariableConfig
    pKa_ref: VariableConfig
    Temp: VariableConfig

class ChemEngine:
    R_GAS = 8.314
    T_REF_K = 298.15
    DELTA_H_KW = 55836.0
    DELTA_H_AC = -410.0
    MW_NaOAc = 82.034

    @staticmethod
    def davies_gamma(I, z):
        I = np.maximum(I, 1e-9)
        A = 0.509
        sqrt_I = np.sqrt(I)
        return 10.0**(-A * (z**2) * ((sqrt_I / (1.0 + sqrt_I)) - 0.3 * I))

    @staticmethod
    def vant_hoff(K_ref, delta_H, T_k):
        return K_ref * np.exp((-delta_H / ChemEngine.R_GAS) * ((1/T_k) - (1/ChemEngine.T_REF_K)))

    @staticmethod
    def calc_beta(a_H, Kw, Ka, C_tot):
        H_conc = a_H 
        OH_conc = Kw / H_conc
        denom = (Ka + H_conc)**2
        num = C_tot * Ka * H_conc
        buffer_term = num / np.maximum(denom, 1e-12)
        beta = 2.303 * (H_conc + OH_conc + buffer_term)
        return beta

    @staticmethod
    def solve_system(df):
        T_K = df['Temp'] + 273.15
        Ka_T = ChemEngine.vant_hoff(10**(-df['pKa_ref']), ChemEngine.DELTA_H_AC, T_K)
        Kw_T = ChemEngine.vant_hoff(1e-14, ChemEngine.DELTA_H_KW, T_K)

        C_Na = (df['mass_NaOAc'] / ChemEngine.MW_NaOAc) / df['V_tot']
        C_Cl = (df['Conc_HCl'] * df['Vol_HCl']) / df['V_tot']
        
        I = C_Na.copy()
        a_H_final = None
        
        for _ in range(15):
            a_H = 10**(-df['pH_measured'])
            g1 = ChemEngine.davies_gamma(I, 1)
            c_H3O = a_H / g1
            c_OH = (Kw_T / a_H) / g1
            r = Ka_T / (a_H * g1)
            c_Ac = C_Na * (r / (1 + r))
            I_new = 0.5 * (C_Na + C_Cl + c_Ac + c_H3O + c_OH)
            if np.allclose(I, I_new, atol=1e-8): break
            I = I_new
            a_H_final = a_H

        beta_final = ChemEngine.calc_beta(a_H_final, Kw_T, Ka_T, C_Na)
        return I, beta_final

    @staticmethod
    def calc_theoretical_ph(df):
        T_K = df['Temp'] + 273.15
        Ka_T = ChemEngine.vant_hoff(10**(-df['pKa_ref']), ChemEngine.DELTA_H_AC, T_K)
        n_salt = df['mass_NaOAc'] / ChemEngine.MW_NaOAc
        n_acid_added = df['Conc_HCl'] * df['Vol_HCl']
        n_A_minus = n_salt - n_acid_added 
        n_HA = n_acid_added
        C_A = n_A_minus / df['V_tot']
        C_HA = n_HA / df['V_tot']
        I_calc = 0.5 * ((n_salt/df['V_tot']) + (n_acid_added/df['V_tot']) + C_A) 
        g1 = ChemEngine.davies_gamma(I_calc, 1)
        log_ratio = np.log10(np.maximum(C_A, 1e-10) / np.maximum(C_HA, 1e-10))
        log_gamma = np.log10(g1)
        pKa_T = -np.log10(Ka_T)
        ph_theoretical = pKa_T + log_ratio + log_gamma
        return ph_theoretical

class SimulationManager:
    def __init__(self, params):
        self.params = params
        self.vars = {k: v for k, v in params.__dict__.items()}
        self.keys = list(self.vars.keys())

    def get_deterministic(self):
        data = {k: [v.mean] for k, v in self.vars.items()}
        return pd.DataFrame(data)

    def generate_lhs(self, n):
        sampler = LatinHypercube(d=len(self.keys), seed=42)
        sample = sampler.random(n=n)
        data = {}
        for i, k in enumerate(self.keys):
            cfg = self.vars[k]
            data[k] = norm.ppf(sample[:, i], loc=cfg.mean, scale=cfg.std)
        return pd.DataFrame(data)

    def generate_mc(self, n):
        rng = np.random.default_rng(int(time.time()))
        data = {k: rng.normal(v.mean, v.std, n) for k, v in self.vars.items()}
        return pd.DataFrame(data)

    def check_convergence(self, data, tolerance_percent=0.01):
        try:
            data_values = np.array(data).flatten()
            cumsum = np.cumsum(data_values)
            steps = np.arange(1, len(data_values) + 1)
            cum_means = cumsum / steps
            final_mean = cum_means[-1]
            rel_error = np.abs((cum_means - final_mean) / final_mean) * 100
            
            not_converged_indices = np.where(rel_error > tolerance_percent)[0]
            if len(not_converged_indices) == 0: step_conv = 1 
            elif not_converged_indices[-1] == len(data_values) - 1: step_conv = "Non raggiunta"
            else: step_conv = not_converged_indices[-1] + 1
            return cum_means, rel_error, step_conv
        except:
            return np.zeros(len(data)), np.zeros(len(data)), "Errore"

    def create_animation(self, data_values, total_sim_steps, n_frames=30, title="Simulazione"):
        try:
            bar_color = '#00E5FF' 
            line_color = '#FF00FF'
            plt.close('all')
            fig, ax = plt.subplots(figsize=(5, 4))
            fig.patch.set_facecolor('black')
            ax.set_facecolor('black')
            data_min, data_max = min(data_values), max(data_values)
            bins = np.linspace(data_min, data_max, 40)
            
            def update(frame):
                ax.clear()
                ax.set_facecolor('black')
                # FIX TRONCAMENTO: Assicuriamo che l'ultimo frame prenda tutto
                chunk_size = max(1, total_sim_steps // n_frames)
                end_idx = min((frame + 1) * chunk_size, len(data_values))
                if frame == n_frames - 1:
                    end_idx = len(data_values) # Forza l'ultimo punto
                
                curr_data = data_values[:end_idx]
                ax.hist(curr_data, bins=bins, color=bar_color, density=True, alpha=0.8, edgecolor='black')
                if len(curr_data) > 0:
                    mu = np.mean(curr_data)
                    ax.axvline(mu, color=line_color, linestyle='--', linewidth=2)
                
                ax.set_title(f"{title}: {len(curr_data)} / {total_sim_steps}", color='white', fontweight='bold')
                ax.tick_params(colors='white')
                ax.set_xlim(data_min, data_max)
                ax.spines['bottom'].set_color('#555')
                ax.spines['left'].set_color('#555')

            ani = animation.FuncAnimation(fig, update, frames=n_frames, repeat=False)
            gif_path = f"temp_sim_{title.replace(' ', '')}.gif"
            ani.save(gif_path, writer='pillow', fps=10)
            plt.close('all')
            return gif_path
        except Exception as e:
            return None

    def plot_static_histogram(self, data, title, color_bar):
        fig, ax = plt.subplots(figsize=(5, 3.5))
        ax.hist(data, bins=40, color=color_bar, alpha=0.7, density=True, edgecolor='black')
        mu = np.mean(data)
        ax.axvline(mu, color='white', linestyle='--', linewidth=1, label=f'Media: {mu:.4f}')
        ax.set_title(title)
        ax.set_xlabel("Forza Ionica [M]")
        ax.legend()
        return fig

# ==========================================
# 2. Interfaccia utente
# ==========================================

def main():
    st.title("Analisi Stocastica Avanzata per Tamponi Acetato 🧪")
    st.markdown("")

    # --- SIDEBAR ---
    st.sidebar.header("📝 Setup Sperimentale")
    
    with st.sidebar.expander("⚖️ Masse e Molarità", expanded=True):
        mean_mass = st.number_input("Massa NaOAc (g)", value=0.8200, format="%.4f")
        std_mass = st.number_input("Err. Massa (g)", value=0.0010, format="%.4f")
        mean_conc_hcl = st.number_input("[HCl] (M)", value=0.100, format="%.3f")
        std_conc_hcl = st.number_input("Err. [HCl]", value=0.002, format="%.3f")

    with st.sidebar.expander("💧 Volumi", expanded=False):
        mean_vtot = st.number_input("V Totale Soluzione (L)", value=0.1000, format="%.4f")
        std_vtot = st.number_input("Err. V Tot (L)", value=0.0005, format="%.4f")
        mean_vol_hcl = st.number_input("Vol HCl Aggiunto (L)", value=0.0500, format="%.4f")
        std_vol_hcl = st.number_input("Err. Vol HCl (L)", value=0.0002, format="%.4f")

    with st.sidebar.expander("🌡️ Condizioni", expanded=False):
        mean_ph = st.number_input("pH Misurato", min_value=0.0, max_value=14.0, value=4.70, format="%.2f")
        std_ph = st.slider("Incertezza pH", 0.01, 0.10, 0.02, 0.01)
        mean_temp = st.number_input("Temp (°C)", value=25.0, format="%.1f")
        std_temp = st.slider("Incertezza T (°C)", 0.1, 5.0, 2.0, 0.1)

    st.sidebar.divider()
    st.sidebar.subheader("⚙️ Impostazioni Simulazione")
    # SLIDER CON MAX 200.000 per evitare limiti
    n_sim_mc = st.sidebar.slider("Step Monte Carlo (MC)", 5000, 200000, 10000, 5000)
    n_sim_lhs = st.sidebar.slider("Step Latin Hypercube (LHS)", 5000, 200000, 20000, 5000)
    tol_conv = st.sidebar.select_slider("Soglia Convergenza (%)", options=[0.1, 0.05, 0.01, 0.005], value=0.01)
    
    if st.sidebar.button("🚀 Avvia Simulazione", type="primary", use_container_width=True):
        
        config = SimulationParams(
            mass_NaOAc=VariableConfig(mean_mass, std_mass, "Massa NaOAc"),
            V_tot=VariableConfig(mean_vtot, std_vtot, "Volume Tot"),
            Conc_HCl=VariableConfig(mean_conc_hcl, std_conc_hcl, "[HCl]"),
            Vol_HCl=VariableConfig(mean_vol_hcl, std_vol_hcl, "Vol HCl"),
            pH_measured=VariableConfig(mean_ph, std_ph, "pH Misurato"),
            pKa_ref=VariableConfig(4.76, 0.01, "pKa (Ref)"),
            Temp=VariableConfig(mean_temp, std_temp, "Temperatura")
        )
        sim = SimulationManager(config)

        # CALCOLI
        with st.status("Calcolo in corso...", expanded=True) as status:
            df_det = sim.get_deterministic()
            res_I_det, res_Beta_det = ChemEngine.solve_system(df_det)
            res_I_det, res_Beta_det = res_I_det[0], res_Beta_det[0]
            
            # Monte Carlo
            st.write(f"🔹 Monte Carlo ({n_sim_mc} step)...")
            df_mc = sim.generate_mc(n_sim_mc)
            res_I_mc, res_Beta_mc = ChemEngine.solve_system(df_mc)
            
            # LHS
            st.write(f"🔹 Latin Hypercube ({n_sim_lhs} step)...")
            df_lhs = sim.generate_lhs(n_sim_lhs)
            res_I_lhs, res_Beta_lhs = ChemEngine.solve_system(df_lhs)
            
            st.write("🔹 Calcolo pH Teorico...")
            ph_theory = ChemEngine.calc_theoretical_ph(df_lhs)
            df_lhs['pH_calc'] = ph_theory
            df_lhs['I_final'] = res_I_lhs
            df_lhs['Beta_final'] = res_Beta_lhs
            
            st.write("🔹 Verifica Convergenza (LHS)...")
            cum_means, rel_error, step_conv = sim.check_convergence(res_I_lhs, tol_conv)
            
            status.update(label="✅ Completato!", state="complete", expanded=False)

        # 1. KPI PRINCIPALI
        st.divider()
        k1, k2, k3, k4 = st.columns(4)
        k1.metric("Forza Ionica", f"{np.mean(res_I_lhs):.4f} M")
        k2.metric("Capacità Tamp.", f"{np.mean(res_Beta_lhs):.4f} M")
        
        if isinstance(step_conv, (int, float, np.integer)):
            k3.metric("Convergenza $I$", f"Step {int(step_conv)}", delta="Ottimale")
        else:
            k3.metric("Convergenza $I$", "Non Raggiunta", delta="-Instabile", delta_color="inverse")
            
        k4.empty()

        # Tabella Riepilogo Statistico
        summary_data = {
            "Metodo": ["Classico", "Monte Carlo", "Latin Hypercube (LHS)"],
            "Step": [1, n_sim_mc, n_sim_lhs],
            "Media $I$ (M)": [res_I_det, np.mean(res_I_mc), np.mean(res_I_lhs)],
            "Dev. Std $I$ (M)": [0.0, np.std(res_I_mc), np.std(res_I_lhs)],
            "Media $\beta$ (M)": [res_Beta_det, np.mean(res_Beta_mc), np.mean(res_Beta_lhs)]
        }
        df_sum = pd.DataFrame(summary_data)
        format_dict = {"Media $I$ (M)": "{:.6f}", "Dev. Std $I$ (M)": "{:.6f}", "Media $\beta$ (M)": "{:.5f}"}
        st.dataframe(df_sum.style.format(format_dict), use_container_width=True, hide_index=True)

        st.divider()

        # =========================================================
        # SEZIONE 1: VALIDAZIONE pH
        # =========================================================
        with st.expander("⚖️ 1. Validazione Modello (pH)", expanded=True):
            col_ph_kpi, col_ph_plot = st.columns([4, 6])
            
            with col_ph_kpi:
                st.markdown("##### 🔍 Confronto Medie")
                kp1, kp2 = st.columns(2)
                mean_meas = np.mean(df_lhs['pH_measured'])
                mean_calc = np.mean(df_lhs['pH_calc'])
                delta_mean = mean_meas - mean_calc
                kp1.metric("pH Misurato (Input)", f"{mean_meas:.3f}")
                kp2.metric("pH Teorico (Calc)", f"{mean_calc:.3f}", delta=f"Diff: {delta_mean:.3f}", delta_color="off")
                st.info("Confronto tra pH misurato e calcolato tramite la stechiometria classica")

            with col_ph_plot:
                st.markdown("##### 📊 Distribuzione Errore")
                delta_vals = df_lhs['pH_measured'] - df_lhs['pH_calc']
                fig_delta, ax_delta = plt.subplots(figsize=(5, 2.5))
                ax_delta.hist(delta_vals, bins=40, color='#FFA500', alpha=0.7, density=True, edgecolor='black')
                ax_delta.axvline(0, color='white', linestyle='--', linewidth=1)
                ax_delta.set_xlabel("Distribuzione pH stocastico")
                st.pyplot(fig_delta)
                plt.close(fig_delta)

        # =========================================================
        # SEZIONE 2: ANALISI DI CONVERGENZA (FIX TRONCAMENTO)
        # =========================================================
        with st.expander("🎯 2. Analisi di Convergenza mediante LHS", expanded=False):
            c_conv1, c_conv2 = st.columns(2)
            with c_conv1:
                st.caption(f"Stabilizzazione Media (Step totali: {n_sim_lhs})")
                fig_conv, ax_conv = plt.subplots(figsize=(5, 3.5))
                # FIX: Assicuriamo che l'asse X copra tutto
                step_plot = max(1, len(cum_means) // 2000)
                x_ax = np.arange(1, len(cum_means)+1)[::step_plot]
                y_ax = cum_means[::step_plot]
                
                # Se il slicing ha tagliato l'ultimo punto, lo aggiungiamo manualmente
                if x_ax[-1] != len(cum_means):
                    x_ax = np.append(x_ax, len(cum_means))
                    y_ax = np.append(y_ax, cum_means[-1])
                
                ax_conv.plot(x_ax, y_ax, color="#00FF00", linewidth=1.5)
                ax_conv.axhline(res_I_det, color="#FF00FF", linestyle="--", alpha=0.7, label="Valore Classico")
                
                if isinstance(step_conv, (int, float, np.integer)):
                    ax_conv.axvline(step_conv, color="white", linestyle=":", alpha=0.8)
                
                ax_conv.set_xlim(0, n_sim_lhs) # FORZA LIMITE ASSE X
                ax_conv.set_xlabel("Iterazioni")
                ax_conv.set_ylabel("Forza Ionica [M]")
                st.pyplot(fig_conv); plt.close(fig_conv)

            with c_conv2:
                st.caption(f"Decadimento Errore % (Soglia: {tol_conv}%)")
                fig_err, ax_err = plt.subplots(figsize=(5, 3.5))
                start_cut = 100
                if n_sim_lhs > start_cut:
                    x_ax_err = np.arange(start_cut, len(rel_error))[::step_plot]
                    y_ax_err = rel_error[start_cut:][::step_plot]
                    
                    if x_ax_err[-1] < len(rel_error) + start_cut - 1:
                        x_ax_err = np.append(x_ax_err, len(rel_error) + start_cut)
                        y_ax_err = np.append(y_ax_err, rel_error[-1])
                        
                    ax_err.plot(x_ax_err, y_ax_err, color="#00E5FF", linewidth=1)
                    ax_err.axhline(tol_conv, color="red", linestyle="--", alpha=0.8)
                    ax_err.set_yscale('log')
                    ax_err.set_xlim(0, n_sim_lhs) # FORZA LIMITE ASSE X
                    ax_err.set_xlabel("Iterazioni")
                    st.pyplot(fig_err); plt.close(fig_err)

        # =========================================================
        # SEZIONE 3: DISTRIBUZIONE E ANIMAZIONI (SEPARATE MC / LHS)
        # =========================================================
        with st.expander("📊 3. Distribuzione Statistica (MC vs LHS)", expanded=False):
            
            tab_mc, tab_lhs, tab_torn = st.tabs(["🎲 Monte Carlo (MC)", "📐 Latin Hypercube (LHS)", "🌪️ Sensibilità"])
            
            with tab_mc:
                c1, c2 = st.columns([4, 6])
                with c1:
                    st.subheader("🎥 Animazione MC")
                    gif_mc = sim.create_animation(res_I_mc, total_sim_steps=n_sim_mc, n_frames=30, title="MC")
                    if gif_mc: st.image(gif_mc, use_container_width=True)
                with c2:
                    st.subheader("📊 Istogramma Finale MC")
                    fig_mc = sim.plot_static_histogram(res_I_mc, f"Distribuzione Monte Carlo ({n_sim_mc} step)", "#00E5FF")
                    st.pyplot(fig_mc); plt.close(fig_mc)
            
            with tab_lhs:
                c1, c2 = st.columns([4, 6])
                with c1:
                    st.subheader("🎥 Animazione LHS")
                    gif_lhs = sim.create_animation(res_I_lhs, total_sim_steps=n_sim_lhs, n_frames=30, title="LHS")
                    if gif_lhs: st.image(gif_lhs, use_container_width=True)
                with c2:
                    st.subheader("📊 Istogramma Finale LHS")
                    fig_lhs = sim.plot_static_histogram(res_I_lhs, f"Distribuzione LHS ({n_sim_lhs} step)", "#00FF00")
                    st.pyplot(fig_lhs); plt.close(fig_lhs)
            
            with tab_torn:
                st.subheader("🌪️ Tornado Plot (LHS)")
                try:
                    inputs_only = df_lhs.drop(['I_final', 'Beta_final', 'pH_calc'], axis=1)
                    corrs = inputs_only.corrwith(df_lhs['I_final']).sort_values()
                    fig_torn, ax_torn = plt.subplots(figsize=(8, 4))
                    colors = ['#FF2B2B' if x < 0 else '#00FF00' for x in corrs]
                    corrs.plot(kind='barh', color=colors, ax=ax_torn, width=0.7)
                    ax_torn.set_xlabel("Correlazione Pearson con Forza Ionica")
                    st.pyplot(fig_torn); plt.close(fig_torn)
                except:
                    st.warning("Errore calcolo sensibilità")

        # =========================================================
        # SEZIONE 4: SCATTER PLOTS (A TENDINA)
        # =========================================================
        with st.expander("🔍 4. Dettaglio Correlazioni (Scatter)", expanded=False):
            def plot_scatter(data, col_x, col_y, color_point):
                fig, ax = plt.subplots(figsize=(5, 3.5))
                ax.scatter(data[col_x], data[col_y], alpha=0.05, s=2, color=color_point, rasterized=True)
                z = np.polyfit(data[col_x], data[col_y], 1)
                p = np.poly1d(z)
                ax.plot(data[col_x], p(data[col_x]), color='white', linestyle='--', linewidth=1, alpha=0.5)
                ax.set_xlabel(col_x)
                return fig

            tab_I, tab_Beta = st.tabs(["📈 Forza Ionica ($I$)", "🛡️ Capacità Tamp. ($\beta$)"])
            with tab_I:
                c1, c2, c3 = st.columns(3)
                c1.pyplot(plot_scatter(df_lhs, 'mass_NaOAc', 'I_final', '#00E5FF'))
                c2.pyplot(plot_scatter(df_lhs, 'V_tot', 'I_final', '#00E5FF'))
                c3.pyplot(plot_scatter(df_lhs, 'pH_measured', 'I_final', '#00E5FF'))
                plt.close('all') 
            with tab_Beta:
                c1, c2, c3 = st.columns(3)
                c1.pyplot(plot_scatter(df_lhs, 'mass_NaOAc', 'Beta_final', '#FFD700'))
                c2.pyplot(plot_scatter(df_lhs, 'V_tot', 'Beta_final', '#FFD700'))
                c3.pyplot(plot_scatter(df_lhs, 'pH_measured', 'Beta_final', '#FFD700'))
                plt.close('all')

        # EXPORT
        st.divider()
        st.subheader("💾 Export Dati")
        c1, c2 = st.columns(2)
        csv_sum = df_sum.to_csv(index=False).encode('utf-8')
        c1.download_button("📄 Scarica Report", csv_sum, "report.csv", "text/csv", type="primary")
        csv_raw = df_lhs.to_csv(index=False).encode('utf-8')
        c2.download_button("📈 Scarica Dati Grezzi (LHS)", csv_raw, "raw_data.csv", "text/csv")

    else:
        st.markdown("<h2 style='text-align: center;'>👋 Benvenuto nel Simulatore</h2>", unsafe_allow_html=True)

if __name__ == "__main__":
    main()

##########################
#streamlit run nomefile.py
##########################