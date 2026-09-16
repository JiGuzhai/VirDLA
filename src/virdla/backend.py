"""Interface for legally obtained, user-supplied model predictions."""

from dataclasses import dataclass
from typing import Mapping, Protocol, Sequence


@dataclass(frozen=True)
class Prediction:
    """Position-aligned model outputs for each requested track.

    ``values`` maps track names to one scalar per genomic bin. The first bin
    begins at position zero of ``Construct.sequence``. Backend adapters must
    provide outputs in model space and document their units.
    """

    bin_size: int
    values: Mapping[str, Sequence[float]]


class PredictionBackend(Protocol):
    name: str

    def predict(self, sequence: str, context: str, tracks: tuple[str, ...]) -> Prediction:
        """Return aligned prediction tracks for one assembled construct."""
