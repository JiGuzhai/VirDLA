# VirDLA

VirDLA is a virtual dual-luciferase assay interface for sequence-to-function
model backends. It fixes a reporter construct and biological context, executes
a user-supplied backend, and converts aligned model outputs into relative
test-to-reference reporter activity.

This repository contains the VirDLA assay framework and one reference readout
algorithm, the prediction-based readout (PB). Readouts are interchangeable
within the framework. The package does not contain model-specific adapters,
pretrained models, context-to-track maps, production T1 or H5 cassette
sequences, or APromoter Atlas values, so it does not reproduce the reported
Atlas snapshot by itself.

The accompanying manuscript is titled **“VirDLA: A Zero-Shot Virtual Assay
Harness for Multi-Omics Foundation Models.”** Supporting data are archived at
[Zenodo](https://doi.org/10.5281/zenodo.22763539). The APromoter web service and
precomputed Atlas are available at [apromoter.bio](https://apromoter.bio).

## Installation

```bash
pip install .
```

## Quick Start

```python
from virdla import AssaySpec, Prediction, VirtualAssay


class ExampleBackend:
    name = "synthetic"

    def predict(self, sequence, context, tracks):
        values = [2.0 if base == "G" else 1.0 for base in sequence]
        return Prediction(1, {track: values for track in tracks})

spec = AssaySpec(test_reporter="GGGG", reference_reporter="AAAA", reference_promoter="ACGT")
assay = VirtualAssay(spec, ExampleBackend(), {"example_cell": ("example_track",)})
result = assay.run("ACGTACGT", "example_cell")

print(result.activity)
```

The example backend generates synthetic values and is not a biological model.
See `examples/mock_assay.py` for the complete offline example. Real model
backends and API services must be obtained and used under their own terms.

## Output alignment

`Prediction.values` contains one position-aligned array per requested track.
`Prediction.bin_size` gives the number of construct bases covered by each
value. Reporter means are weighted by the number of reporter bases overlapping
each bin, so partially overlapping boundary bins do not receive full weight.

## Third-party material

The production cassettes referenced by the manuscript are not distributed here.
**T1** (`pGL4.10_DLAmonitor`) derives from the Promega pGL4.10 dual-luciferase
vector; **H5** derives from the reporter plasmids of Ernst et al.,
*Nature Biotechnology* **34**, 1180–1190 (2016), doi:10.1038/nbt.3678
(data: GEO [GSE71279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71279)).
Both remain subject to their own terms and are not relicensed by this package.
Callers must supply their own legally usable reporter and reference sequences.

## License

The package is released under the GNU General Public License, version 3 or
later (GPL-3.0-or-later). That license applies only to the files in this
package; it does not license external models, model outputs, data or services.
See `PROVENANCE.md` for the release boundary.

---

*An earlier version of this repository has been superseded and is no longer maintained.*

