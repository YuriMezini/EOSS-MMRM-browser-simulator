# EOSS MMRM browser simulator

Browser version of the EOSS_MMRM power simulation workflow, with scripts for comparing the JavaScript fit against R `mmrm`.

Jing Dai and Jurgen Mezinaj

This project follows the estimand-oriented MMRM setup from the original EOSS_MMRM Shiny app (PHUSE US Connect 2024, PAP_AS03). Source for the original R app: https://github.com/Kong-WayneState/EOSS_MMRM-V1.0

## Files

- `eoss_mmrm_simulator_final.html` — open this in a browser
- `codes/` — simulation engine and paired JS–R evaluation scripts

The browser analysis is an iterative GLS fit with EM-type covariance updates. It is not REML. For exact inference, use R `mmrm` with Kenward–Roger.

## Using the HTML app

1. Open `eoss_mmrm_simulator_final.html` in Chrome, Firefox, Safari, or Edge.
2. Choose an example (treatment policy, hypothetical, placebo-like EFU, or null).
3. Set sample size and number of simulations, then run.

No installation is required for the HTML app.

## Using the paired scripts

The scripts in `codes/` regenerate shared datasets and fit both engines on the same observed data. You will need Node.js and R with the `mmrm` package.

See `codes/README.md` for a short file list.
