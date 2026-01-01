# 🧪 Simulatore Buffer Ultimate (Acetate Edition)

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![Streamlit](https://img.shields.io/badge/Streamlit-1.28%2B-red)
![Status](https://img.shields.io/badge/Status-Stable-success)
![License](https://img.shields.io/badge/License-MIT-green)

Un simulatore stocastico avanzato per l'analisi della propagazione dell'incertezza nella preparazione di tamponi acetato. Sviluppato per fornire un'analisi rigorosa della **Forza Ionica ($I$)**, della **Capacità Tamponante ($\beta$)** e della validazione del **pH**.

## ✨ Caratteristiche Principali

Il software combina la chimica analitica con metodi statistici avanzati:

* **🔬 Motore Chimico Rigoroso**:
    * Correzione termodinamica delle costanti ($K_a$, $K_w$) tramite equazione di Van't Hoff.
    * Calcolo dei coefficienti di attività ($\gamma$) con l'equazione di Davies.
    * Solver iterativo per il calcolo esatto della Forza Ionica e del pH all'equilibrio.
* **🎲 Analisi Stocastica**:
    * **Monte Carlo (MC)**: Per animazioni fluide e visualizzazione della distribuzione.
    * **Latin Hypercube Sampling (LHS)**: Per una copertura ottimale dello spazio delle probabilità e alta precisione statistica.
* **🎯 Analisi di Convergenza**: Algoritmi automatici per determinare quando la simulazione raggiunge la stabilità statistica (errore relativo < 0.01%).
* **⚖️ Validazione Sperimentale**: Confronto automatico tra pH misurato (input sperimentale) e pH teorico (stechiometrico) per identificare bias sistematici.
* **📊 Visualizzazione Avanzata**:
    * Interfaccia "Dark Mode" ottimizzata.
    * Animazioni GIF della distribuzione in tempo reale.
    * Tornado Plots per l'analisi di sensibilità (quali variabili impattano di più?).
    * Scatter plots per l'analisi delle correlazioni.

## 🚀 Installazione e Avvio

### Prerequisiti
Assicurati di avere Python installato (versione 3.8 o superiore).

1.  **Clona il repository (o scarica i file):**
    ```bash
    git clone [https://github.com/tuo-username/simulatore-buffer.git](https://github.com/tuo-username/simulatore-buffer.git)
    cd simulatore-buffer
    ```

2.  **Installa le dipendenze:**
    È consigliato usare un ambiente virtuale (`venv`).
    ```bash
    pip install streamlit numpy pandas matplotlib scipy seaborn
    ```
    *Nota: `scipy` deve essere aggiornato (>= 1.7.0) per supportare `LatinHypercube`.*

3.  **Avvia l'applicazione:**
    ```bash
    streamlit run app.py
    ```

## 🖥️ Interfaccia Utente

### Sidebar (Setup)
Inserisci i dati sperimentali con le relative incertezze (deviazione standard):
* **Reagenti**: Massa NaOAc, Concentrazione HCl.
* **Volumi**: Volume totale soluzione, Volume HCl aggiunto.
* **Condizioni**: pH misurato, Temperatura.
* **Impostazioni Simulazione**: Numero di step (MC/LHS) e soglia di convergenza.

### Dashboard Analitica
1.  **KPI Cards**: Risultati immediati per $I$, $\beta$, pH Teorico e stato della Convergenza.
2.  **Validazione**: Confronto tra pH misurato e teorico.
3.  **Convergenza**: Grafici di stabilizzazione della media e decadimento dell'errore.
4.  **Distribuzione**: Istogrammi e GIF animate.
5.  **Export**: Download dei report in formato CSV.

## 🧪 Background Scientifico

Il simulatore utilizza le seguenti equazioni fondamentali:

**1. Equazione di Davies (Attività):**
$$ \log \gamma = -0.509 z^2 \left( \frac{\sqrt{I}}{1 + \sqrt{I}} - 0.3I \right) $$

**2. Capacità Tamponante ($\beta$):**
$$ \beta = 2.303 \left( [H^+] + [OH^-] + \frac{C_{tot} \cdot K_a \cdot [H^+]}{(K_a + [H^+])^2} \right) $$

**3. Equazione di Van't Hoff (Temperatura):**
$$ \ln \left( \frac{K_2}{K_1} \right) = \frac{-\Delta H}{R} \left( \frac{1}{T_2} - \frac{1}{T_1} \right) $$

## 📂 Struttura del Codice

* `app.py`: File principale contenente l'intera logica.
    * `ChemEngine`: Classe statica per i calcoli chimico-fisici.
    * `SimulationManager`: Gestisce la generazione di numeri casuali (LHS/MC) e l'analisi statistica.
    * `main()`: Gestisce l'interfaccia Streamlit e il flusso di lavoro.

## 🤝 Contribuire

I contributi sono benvenuti! Sentiti libero di aprire una *issue* o inviare una *pull request*.

1.  Forka il progetto
2.  Crea il tuo feature branch (`git checkout -b feature/NuovaFeature`)
3.  Committa le modifiche (`git commit -m 'Aggiunta nuova feature'`)
4.  Pusha sul branch (`git push origin feature/NuovaFeature`)
5.  Apri una Pull Request

## 📜 Licenza

Distribuito sotto licenza MIT. Vedi `LICENSE` per maggiori informazioni.

---
*Autore: [Il Tuo Nome] - Tesi Triennale/Progetto di Ricerca - 2026*
