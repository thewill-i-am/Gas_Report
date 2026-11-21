# GAS Report

[Open the English report](https://thewill-i-am.github.io/Gas_Report/)

Interactive benchmarks for uniform and guided Grover Adaptive Search (GAS):

- Plasma, S1 and Shotgun, each at panel sizes **k = 3, 5 and 7**.
- Best objective per dataset, panel size and preparation, selected by mean LDA/SVM balanced accuracy over 20 seeds.
- Frozen random, QAOA and published Bias Field comparators.
- Collapsible results, biomarker frequencies, exact panels, QUBO construction and preparation guides.
- Offline LaTeX rendering with embedded KaTeX and accessible MathML.

The main report is in English. Its archive link opens the original Spanish round-2 report.

## Interpretation

Results are exploratory: feature selection precedes LOOCV and objective selection is retrospective. Published BF-DCQO results use k=5 accuracy and are identified separately from balanced accuracy. These results do not establish clinical validity or quantum advantage.

## Files and publishing

- `index.html`: entry point to the English report.
- `report/top1.html`: self-contained English report.
- `report/top1.json`: selected configurations and aggregate biomarker/panel frequencies.
- `report/index.html`: original round-2 archive.

GitHub Pages publishes the root of the `main` branch. To update the site, replace the report files and push to `main`.

Source scripts and experimental provenance are maintained in [Quantum-Feature-Selection-](https://github.com/thewill-i-am/Quantum-Feature-Selection-/tree/de92848/presentations/gas_objective_round2).
KaTeX's MIT license is embedded in the report.
