import unittest

from virdla import AssaySpec, Prediction, VirtualAssay, build_construct, pearson, spearman


class FixedBackend:
    name = "fixed"

    def predict(self, sequence, context, tracks):
        values = [3.0 if base == "G" else 1.0 for base in sequence]
        return Prediction(1, {track: values for track in tracks})


class AlternativeBackend(FixedBackend):
    name = "alternative"

    def predict(self, sequence, context, tracks):
        values = [5.0 if base == "G" else 1.0 for base in sequence]
        return Prediction(1, {track: values for track in tracks})


class AssayTests(unittest.TestCase):
    def setUp(self):
        self.spec = AssaySpec("GG", "AA", "AC", "NN", window=20)

    def test_construct_and_spans(self):
        construct = build_construct("ACGT", self.spec)
        self.assertEqual(len(construct.sequence), 20)
        self.assertEqual(construct.sequence[construct.test_span[0]:construct.test_span[1]], "GG")
        self.assertEqual(construct.sequence[construct.reference_span[0]:construct.reference_span[1]], "AA")

    def test_reference_ratio(self):
        result = VirtualAssay(self.spec, FixedBackend(), {"cell": ("track",)}).run("ACGT", "cell")
        self.assertEqual(result.test_signal, 3.0)
        self.assertEqual(result.reference_signal, 1.0)
        self.assertEqual(result.activity, 3.0)
        self.assertEqual(result.readout, "PB")

    def test_partial_bins_are_overlap_weighted(self):
        class CoarseBackend:
            name = "coarse"

            def predict(self, sequence, context, tracks):
                # The 2-bp test reporter spans one base in each of two bins.
                return Prediction(5, {"track": [2.0, 4.0, 1.0]})

        spec = AssaySpec("GG", "AA", "AC", "NNN")
        result = VirtualAssay(spec, CoarseBackend(), {"cell": ("track",)}).run("ACGT", "cell")
        self.assertEqual(result.test_signal, 3.0)
        self.assertEqual(result.reference_signal, 1.0)
        self.assertEqual(result.activity, 3.0)

    def test_backend_is_replaceable(self):
        a = VirtualAssay(self.spec, FixedBackend(), {"cell": ("track",)}).run("ACGT", "cell")
        b = VirtualAssay(self.spec, AlternativeBackend(), {"cell": ("track",)}).run("ACGT", "cell")
        self.assertEqual((a.candidate, a.context), (b.candidate, b.context))
        self.assertEqual((a.activity, b.activity), (3.0, 5.0))

    def test_missing_track_fails(self):
        assay = VirtualAssay(self.spec, FixedBackend(), {"cell": ("track",)})
        with self.assertRaises(KeyError):
            assay.run("ACGT", "other")

    def test_invalid_construct_fails(self):
        with self.assertRaises(ValueError):
            build_construct("AX", self.spec)
        with self.assertRaises(ValueError):
            build_construct("ACGT", AssaySpec("GG", "AA", "AC", window=4))

    def test_correlations(self):
        self.assertAlmostEqual(pearson([1, 2, 3], [2, 4, 6]), 1.0)
        self.assertAlmostEqual(spearman([1, 3, 2], [10, 30, 20]), 1.0)
        self.assertAlmostEqual(spearman([1, 2, 2, 4], [4, 2, 2, 1]), -1.0)


    def test_oop_assemble_and_batch(self):
        c1 = self.spec.assemble("ACGT")
        c2 = build_construct("ACGT", self.spec)
        self.assertEqual(c1.sequence, c2.sequence)
        self.assertEqual(c1.test_span, c2.test_span)

        assay = VirtualAssay(self.spec, FixedBackend(), {"cell": ("track",)})
        res = assay.run("ACGT", "cell")
        d = res.to_dict()
        self.assertEqual(d["candidate"], "ACGT")
        self.assertEqual(d["activity"], 3.0)

        batch = assay.run_batch(["ACGT", "GGCC"], ["cell"])
        self.assertEqual(len(batch), 2)
        self.assertEqual(batch[0].candidate, "ACGT")
        self.assertEqual(batch[1].candidate, "GGCC")


if __name__ == "__main__":
    unittest.main()
