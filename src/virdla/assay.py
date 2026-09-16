"""Reporter-construct assembly and public prediction-based readout."""

from dataclasses import dataclass
from math import isfinite
from typing import Mapping

from .backend import Prediction, PredictionBackend

DNA = frozenset("ACGTN")


def _dna(value: str, field: str, *, empty: bool = False) -> str:
    sequence = value.strip().upper()
    if (not sequence and not empty) or set(sequence) - DNA:
        raise ValueError(f"{field} must contain DNA bases A, C, G, T or N")
    return sequence


@dataclass(frozen=True)
class AssaySpec:
    """Fixed cassette components supplied by the assay designer.

    No production T1/H5 sequences are shipped with this package. A caller may
    provide its own legally usable reporter and reference sequences.
    """

    test_reporter: str
    reference_reporter: str
    reference_promoter: str
    insulator: str = ""
    window: int | None = None

    def assemble(self, candidate: str) -> "Construct":
        """Assemble a candidate promoter with this assay specification."""
        return build_construct(candidate, self)


@dataclass(frozen=True)
class Construct:
    sequence: str
    test_span: tuple[int, int]
    reference_span: tuple[int, int]

    @classmethod
    def from_spec(cls, candidate: str, spec: AssaySpec) -> "Construct":
        """Assemble a construct given a candidate sequence and assay specification."""
        return build_construct(candidate, spec)


def build_construct(candidate: str, spec: AssaySpec) -> Construct:
    candidate = _dna(candidate, "candidate")
    test = _dna(spec.test_reporter, "test_reporter")
    reference = _dna(spec.reference_reporter, "reference_reporter")
    promoter = _dna(spec.reference_promoter, "reference_promoter")
    insulator = _dna(spec.insulator, "insulator", empty=True)
    cassette = candidate + test + insulator + promoter + reference
    if spec.window is not None:
        if not isinstance(spec.window, int) or spec.window < len(cassette):
            raise ValueError("window must be an integer at least as long as the cassette")
        left = (spec.window - len(cassette)) // 2
        cassette = "N" * left + cassette + "N" * (spec.window - len(cassette) - left)
    else:
        left = 0
    test_start = left + len(candidate)
    reference_start = test_start + len(test) + len(insulator) + len(promoter)
    return Construct(cassette, (test_start, test_start + len(test)),
                     (reference_start, reference_start + len(reference)))


def _region_mean(prediction: Prediction, span: tuple[int, int], tracks: tuple[str, ...]) -> float:
    if not isinstance(prediction.bin_size, int) or prediction.bin_size <= 0:
        raise ValueError("prediction bin_size must be a positive integer")
    start = span[0] // prediction.bin_size
    stop = (span[1] + prediction.bin_size - 1) // prediction.bin_size
    weighted_sum = 0.0
    total_weight = 0
    for track in tracks:
        if track not in prediction.values:
            raise KeyError(f"backend did not return requested track {track!r}")
        values = prediction.values[track]
        if stop > len(values):
            raise ValueError(f"track {track!r} does not cover reporter span")
        for index in range(start, stop):
            value = float(values[index])
            if not isfinite(value):
                raise ValueError("reporter signal is non-finite")
            bin_start = index * prediction.bin_size
            bin_stop = bin_start + prediction.bin_size
            overlap = max(0, min(span[1], bin_stop) - max(span[0], bin_start))
            weighted_sum += value * overlap
            total_weight += overlap
    if total_weight == 0:
        raise ValueError("reporter signal is empty")
    return weighted_sum / total_weight


@dataclass(frozen=True)
class Result:
    candidate: str
    context: str
    backend: str
    readout: str
    activity: float
    test_signal: float
    reference_signal: float

    def to_dict(self) -> dict[str, float | str]:
        """Convert result to a plain dictionary for tabular serialization."""
        return {
            "candidate": self.candidate,
            "context": self.context,
            "backend": self.backend,
            "readout": self.readout,
            "activity": self.activity,
            "test_signal": self.test_signal,
            "reference_signal": self.reference_signal,
        }


class VirtualAssay:
    """Run a fixed assay specification through a replaceable PB backend."""

    def __init__(self, spec: AssaySpec, backend: PredictionBackend,
                 context_tracks: Mapping[str, tuple[str, ...]]):
        self.spec = spec
        self.backend = backend
        self.context_tracks = dict(context_tracks)

    def run(self, candidate: str, context: str) -> Result:
        tracks = self.context_tracks.get(context)
        if not tracks:
            raise KeyError(f"no tracks configured for context {context!r}")
        construct = self.spec.assemble(candidate)
        prediction = self.backend.predict(construct.sequence, context, tracks)
        test = _region_mean(prediction, construct.test_span, tracks)
        reference = _region_mean(prediction, construct.reference_span, tracks)
        if reference == 0:
            raise ZeroDivisionError("reference reporter signal is zero")
        activity = test / reference
        if not isfinite(activity):
            raise ValueError("activity is non-finite")
        return Result(_dna(candidate, "candidate"), context, self.backend.name,
                      "PB", activity, test, reference)

    def run_batch(self, candidates: list[str], contexts: list[str]) -> list[Result]:
        """Evaluate multiple candidate sequences across multiple biological contexts.

        Parameters
        ----------
        candidates : list[str]
            List of candidate promoter DNA sequences.
        contexts : list[str]
            List of biological contexts (e.g., cell types or tissues).

        Returns
        -------
        list[Result]
            List of virtual assay results for each (candidate, context) pair.
        """
        results: list[Result] = []
        for cand in candidates:
            for ctx in contexts:
                results.append(self.run(cand, ctx))
        return results
