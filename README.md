# relSig

`relSig` is a small **MATLAB** package that helps you decide *which* out of several simultaneously tested regression effects remain reliable once you control the **family‑wise error rate (FWE)**. It combines ordinary‑least‑squares (OLS) modelling with a bootstrap‑calibrated version of the **Holm–Šídák** multiple‑comparison procedure.

---

## Key features

| What | Where in the code |
|------|-------------------|
| OLS estimation that returns *p*-values, F‑statistics and residuals | `OLS.m` |
| Berens–Holm–Šídák step‑down correction | `berens_holm_sidak.m` |
| Automatic search for the effective number of tests *m′* that keeps the bootstrap FWE at α&nbsp;=&nbsp;0.05 | `fnc2min_BHS.m` & `relSig_FWE.m` |
| End‑to‑end helper that wraps everything into a single call | `relSig_FWE.m` |
| Fully‑worked reproducible example | `relSig_Example.m` |

---

## Quick start

```bash
git clone https://github.com/Sam-Berens/relSig.git
cd relSig
matlab -batch "relSig_Example"
```

The example script will

1. simulate multivariate data with a known covariance matrix,  
2. fit an OLS model,  
3. obtain raw *p*-values for each response variable,  
4. call `relSig_FWE`, and  
5. print a logical vector `sig` whose `true` entries mark effects that survive FWE control.

---

## Usage in your own code
```matlab
[p,fStat,Bhat,~,Err] = OLS(Y,X,H);   % raw tests  (m tests)
ErrSigmaHat         = cov(Err);      % estimated Σ̂
sig = relSig_FWE(ErrSigmaHat,X,H,p); % 1×m logical result
```
| Argument | Meaning |
|----------|---------|
| `Y` | *n × m* matrix of responses |
| `X` | *n × p* design matrix (include a column of ones for an intercept) |
| `H` | Contrast matrix that defines the hypothesis **H β&nbsp;=&nbsp;0** |
| `p` | 1 × *m* vector of raw *p*-values from `OLS` |
| output `sig` | 1 × *m* logical mask: **`true` = significant at FWE α = 0.05** |

### Optional tweaks

Open **`relSig_FWE.m`** and modify
```matlab
nIter = 1e5;   % number of bootstrap samples
alpha = 0.05;  % target family‑wise error rate
```
to trade speed for precision or change the error‑rate threshold.

---

## How it works (method summary)

1. **OLS step** – obtain test statistics and residual covariance Σ̂.  
2. **Null bootstrap** – repeatedly draw *n*‑sample data from **N**(0, Σ̂), refit the same model and store the vector of null *p*-values.  
3. **Find effective m′** – choose *m′* that makes the empirical FWE equal to the nominal α (root‑finding via `fminbnd`).  
4. **Correct** – apply Berens–Holm–Šídák with that *m′* to the original *p*-values.
The result is a step‑down procedure whose FWE is calibrated in finite samples instead of relying on large‑sample or independence assumptions.

---

## Requirements

* MATLAB **R2020a** or later (earlier versions may work but are untested)  
* Statistics & Machine Learning Toolbox – for `mvnrnd`, `fcdf`, etc.

---

## Contributing

Contributions, suggestions, and improvements are welcome! Please fork the repository and submit a pull request with your changes.

---

## License

This project is licensed under the [CC BY 4.0 License](https://creativecommons.org/licenses/by/4.0/legalcode.en).
You are free to share and adapt the material under the following terms:
- **Attribution:** You must give appropriate credit, provide a link to the license, and indicate if changes were made.
- **No additional restrictions:** You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.

---

## Contact
For questions or feedback, please contact Sam Berens at [s.berens@sussex.ac.uk](mailto:s.berens@sussex.ac.uk).
