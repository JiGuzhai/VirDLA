# Scope and provenance audit

The package code was authored for this release after inspecting the APSCS production v11 scorer and the CBC manuscript's assay specification. No production implementation file was copied.

## Framework and algorithms

VirDLA is the **assay framework**: it defines the reporter construct, reference configuration and biological context, assembles the model input and returns relative reporter activity. **Readout algorithms plug into this framework.** Prediction-based readout (PB) and Trunk Activation Readout (TAR) are two such algorithms, and model-specific adapters connect the framework to a given backend.

This release ships the **framework** (a generic cassette assembler, an assay specification, a backend protocol and simple correlation utilities) together with **one reference readout algorithm, PB** (a model-output reporter-region mean and a test/reference ratio). It does not ship the production readout algorithm (TAR) or the backend layer.

## Deliberately excluded

- **Readout algorithm:** the v11 trunk extraction and local RNA-seq head path, TAR algorithms and parameters.
- **Backend layer:** model-specific adapters, pretrained weights, ontology/track routing and routing calibration.
- **Assay and data assets:** production T1/H5 cassette sequences, Atlas matrices, benchmark datasets and third-party source.
- **Release hygiene:** private service files and credentials.

The included mock backend is synthetic.

## Third-party provenance of the excluded cassette sequences

The two production cassettes referenced by the manuscript are not distributed here because they derive from third-party material with terms independent of this package's license:

- **T1** (`pGL4.10_DLAmonitor`) derives from the Promega **pGL4.10** dual-luciferase vector. Its backbone and sequence are Promega product material; callers who need it should obtain it under Promega's own terms (see the pGL4.10 vector product documentation).
- **H5** derives from the reporter plasmids of the SHARPR-MPRA benchmark study: Ernst, J. et al. Genome-scale high-resolution mapping of activating and repressive nucleotides in regulatory regions. *Nature Biotechnology* **34**, 1180–1190 (2016). doi:10.1038/nbt.3678. Experimental data: GEO accession GSE71279. Attribution and any reuse conditions follow that publication and dataset.

The insulator shared by T1 and H5 (cHS4) is itself a third-party element; users publishing or redistributing it should attribute its source.

Callers must supply their own legally usable reporter and reference sequences. `AssaySpec` deliberately ships without any T1/H5 cassette sequence, and no third-party sequence is relicensed by this package.

The PB implementation demonstrates the assay contract and is not a numerical reproduction of the production v11 TAR route. Exact reproduction of the manuscript's reported results requires its execution environment, assay sequences, model access and data.

Only files inside this package are released under GPL-3.0-or-later. External model/API terms remain independent.
