# Supporting repository for Uncovering the flexibility of CDR loops in antibodies and TCRs through large-scale molecular dynamics

## 📝 About
Supporting repository for the manuscript "Uncovering the flexibility of CDR loops in antibodies and TCRs through large-scale molecular dynamics" ([Cagiada M., Spoendlin F.C., et al., bioRxiv 2025](https://www.biorxiv.com)).

This repository provides notebooks, data, and analysis scripts to reproduce the figures presented in the manuscript. Additionally, the repository includes CV3-Fv pipeline scripts for implementing our setup in CALVADOS 3 (COMING SOON...).

---

## 🗂️	Repository layout

The contents of this repository allow you to reproduce the results and figures that appear in the manuscript.
- `notebooks`: Contains Jupyter notebooks used to perform the analyses and generate the figures included in the manuscript
- `src`: Contains all datasets and supplementary files necessary for running the analysis notebooks and reproducing the results.
- `figures`: Includes figures generated as outputs from the analysis notebooks.
- `cv3-pipeline`: (coming soon): Scripts and instructions for reproducing our CV3-Fv molecular dynamics simulation setup.

**N.B.**: due to their large size, the all-atom simulations used in this analysis are not included in the current repository. Simulations can be found at: [https://doi.org/10.5281/zenodo.17525665](https://doi.org/10.5281/zenodo.17525665).

---

## ⚙️ CV3-Fv pipeline

Coming soon...

----

## 📚 Data

- CV3-Fv simulations for experimentally resolved structures from the FlAbDab database are accessible through the dedicated webserver [https://opig.stats.ox.ac.uk/webapps/flabdab](https://opig.stats.ox.ac.uk/webapps/flabdab).
- FlAbDab and FTCRDab simulations for experimentally resolved structures can also be downloaded as from Zenodo [https://doi.org/10.5281/zenodo.17552276](https://doi.org/10.5281/zenodo.17552276).
- Ensembles of representative structures for all the antibodies in FlAbDab, which were predicted using ABodyBuilder2, are available on Zenodo [https://doi.org/10.5281/zenodo.17555649](https://doi.org/10.5281/zenodo.17555649).
- All-atom MD simulations for the two antibodies and one TCR, along with associated metadata, are available at: [https://doi.org/10.5281/zenodo.17525665](https://doi.org/10.5281/zenodo.17525665).

---- 
## 📝 Reference this work

If you use our model please cite:

Cagiada, M., Spoendlin, F.C., Ifashe K. & Deane C.M (2025). Uncovering the flexibility of CDR loops in antibodies and TCRs through large-scale molecular dynamics. In bioRxiv (p.). https://doi.org/

```text
@ARTICLE{Cagiada2025-cv,
  title    = "Uncovering the flexibility of CDR loops in antibodies and TCRs through large-scale molecular dynamics",
  author   = "Cagiada, Matteo and Spoendlin, Fabian. C. and Ifashe, King and Deane, Charlotte M.",
  journal  = "bioRxiv",
  pages    = "",
  month    =  ,
  year     =  2025,
  language = "en"
```

Also if you use the CALVADOS 3 suite for our CV3-Fv setup, please cite:

- von Bülow, S., Yasuda, I., Cao, F., Schulze, T. K., Trolle, A. I., Rauh, A. S., Crehuet, R., Lindorff-Larsen, K., & Tesei, G. (2025). Software package for simulations using the coarse-grained CALVADOS model. In arXiv [q-bio.BM]. arXiv. http://arxiv.org/abs/2504.10408

## 🙌 Acknowledgements  

The research was supported by a Novo Nordisk Foundation Postdoctoral Fellowship (NNF23OC0082912) awarded to MC, by  research funding from the UK Engineering and Physical Sciences Research Council (EPSCR) [Grant number EP/S024093/1], Roche and the Royal Commission for the Exhibition of 1851 awarded to FS and KI.

----

## 📜 License  
This project is licensed under the **MIT License**. See [LICENSE](LICENSE) for details.

## 📬 Contact  
For questions or support with this repository, please use the GitHub issue tab or reach out to us via email:  

📧 Matteo Cagiada:  [matteo.cagiada@bio.ku.dk](mailto:matteo.cagiada@bio.ku.dk) / [matteo.cagiada@stats.ox.ac.uk](mailto:matteo.cagiada@stats.ox.ac.uk)

📧 Charlotte M. Deane: [deane@stats.ox.ac.uk](mailto:deane@stats.ox.ac.uk)

