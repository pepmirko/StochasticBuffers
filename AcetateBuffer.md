## Monte Carlo Estimation of the Ionic Strength of an Acetate Buffer

This script estimates the ionic strength of an acetate buffer prepared from
sodium acetate and hydrochloric acid using a Monte Carlo approach.

**Experimental conditions (nominal):**

- Buffer components: NaOAc (sodium acetate) + HCl (hydrochloric acid)  
- pH: 4.00 ± 0.02  
- Total acetate concentration: \( C_T \approx 100\ \mu\text{M} \) (to be verified)

**Model assumptions:**

1. Initially, all acetate is present as the deprotonated species Ac⁻ (from NaOAc).
2. Upon addition of HCl, a fraction of Ac⁻ is protonated to HAc.
3. The Henderson–Hasselbalch equation is used to obtain \([ \mathrm{Ac}^- ]\) and
   \([ \mathrm{HAc} ]\) from \( C_T \), pH, and pKₐ.
4. All relevant ionic species in solution are included: Na⁺, Cl⁻, Ac⁻, H₃O⁺, and OH⁻.
5. The ionic strength is evaluated via a standard Monte Carlo approach by sampling
   experimental parameters (mass, volumes, pH, pKₐ, etc.) from normal distributions.
   A typical simulation uses \( N = 2 \times 10^5 \) samples.

The script outputs:

- A text file with summary statistics of the ionic strength distribution.
- A CSV file with all sampled parameters and the corresponding ionic strength.
