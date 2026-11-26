# calibration-lmm

Code of the article *Finite-sample calibration of regional mixed-model association tests with an
estimated variance ratio* by Kossi N. Doulabe and Lajmi Lakhal-Chaieb (2026), in
[`calibration/`](calibration/README.md). The scripts compute the exact finite-sample distributions
of the Wald, score, likelihood-ratio and Bayes-factor statistics of a candidate variant when the
variance ratio of a regional polygenic background is estimated, the candidate-specific index
$\nu_j$ and the plug-in and protected $p$-values. They run the published mixed-model and
fine-mapping software on the same data and regenerate every table and figure of the article. The
[README of the folder](calibration/README.md) gives the requirements, the data and the scripts
behind each table and figure.

## Citation

Please cite the article (see [`CITATION.cff`](CITATION.cff)):

```
@article{Doulabe2026Calibration,
  author  = {Doulabe, Kossi N. and Lakhal-Chaieb, Lajmi},
  title   = {Finite-sample calibration of regional mixed-model association
             tests with an estimated variance ratio},
  note    = {Preprint},
  year    = {2026}
}
```

## License

MIT, see [`LICENSE`](LICENSE). The patch `calibration/comparators/gemma_arm64_debug.patch`
modifies GEMMA and follows its GPL-3 licence.
