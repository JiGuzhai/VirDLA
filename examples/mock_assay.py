"""Offline example with synthetic predictions, not a biological model."""

from virdla import AssaySpec, Prediction, VirtualAssay


class MockBackend:
    name = "synthetic"

    def predict(self, sequence: str, context: str, tracks: tuple[str, ...]) -> Prediction:
        # Position-aligned dummy values only demonstrate the backend contract.
        values = [2.0 if base == "G" else 1.0 for base in sequence]
        return Prediction(bin_size=1, values={track: values for track in tracks})


spec = AssaySpec(test_reporter="GGGG", reference_reporter="AAAA",
                 reference_promoter="ACGT", insulator="NN")

# 1. Assemble construct directly from spec (OOP method)
construct = spec.assemble("ACGT")
print("Assembled construct:", construct.sequence)

# 2. Run single assay
assay = VirtualAssay(spec, MockBackend(), {"example_cell": ("example_track",)})
result = assay.run("ACGT", "example_cell")
print("Single result:", result.to_dict())

# 3. Batch assay across multiple candidates
batch_results = assay.run_batch(["ACGT", "GGCC"], ["example_cell"])
print(f"Batch evaluated {len(batch_results)} results:")
for r in batch_results:
    print(f"  Candidate: {r.candidate} | Context: {r.context} | Activity: {r.activity:.2f}")

