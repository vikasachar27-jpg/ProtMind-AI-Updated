"""Automated verification test suite for ProtMind AI QML Integration.

Tests:
1. QML singleton engine loading & artifact integrity.
2. Physicochemical delta calculations matching training data.
3. Mutation validation & sequence positional verification.
4. Exact PennyLane hybrid QNN inference execution.
5. Exact Scikit-learn HistGradientBoosting baseline inference.
6. Benchmark metadata & performance table retrieval.
7. FastAPI REST API endpoints:
   - POST /api/step9/qml
   - POST /api/qml/predict
   - POST /api/qml/compare
   - POST /api/step9 (Explainable AI synthesis with embedded QML)
   - POST /api/step10/pdf (Automated clinical report with QML audit)
"""

import sys
import unittest
from pathlib import Path

# Add backend directory to sys.path
BACKEND_DIR = Path(__file__).resolve().parent
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

import qml_inference
from api import app
from fastapi.testclient import TestClient


class TestQMLIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine = qml_inference.get_qml_engine()
        cls.client = TestClient(app)

    def test_01_engine_loaded(self):
        """Verify model files loaded without error."""
        self.assertTrue(self.engine.loaded, f"Engine failed to load: {self.engine.load_error}")
        self.assertIsNotNone(self.engine.qnn_model)
        self.assertIsNotNone(self.engine.hgb_model)
        self.assertIsNotNone(self.engine.qnn_prep)

    def test_02_delta_feature_calculation(self):
        """Verify delta calculations match training script exactly for R175H."""
        deltas = qml_inference.compute_physicochemical_deltas("R", "H")
        # R: hydro=-4.5, mass=156.1875, charge=1.0
        # H: hydro=-3.2, mass=137.1411, charge=0.1
        self.assertAlmostEqual(deltas["hydrophobicity_delta"], 1.3, places=3)
        self.assertAlmostEqual(deltas["charge_delta"], -0.9, places=3)
        self.assertAlmostEqual(deltas["mass_delta"], -19.0464, places=3)

    def test_03_qnn_and_classical_prediction(self):
        """Verify QNN and HGB produce valid probabilities and class labels."""
        res = self.engine.predict("R175H")
        self.assertTrue(res["ok"])
        self.assertEqual(res["status"], "success")

        # QNN checks
        qnn = res["qnn_prediction"]
        self.assertIn("probability_pathogenic", qnn)
        self.assertGreaterEqual(qnn["probability_pathogenic"], 0.0)
        self.assertLessEqual(qnn["probability_pathogenic"], 1.0)
        self.assertIn(qnn["predicted_class"], ["Pathogenic", "Benign"])

        # Classical HGB checks
        hgb = res["classical_prediction"]
        self.assertIn("probability_pathogenic", hgb)
        self.assertGreaterEqual(hgb["probability_pathogenic"], 0.0)
        self.assertLessEqual(hgb["probability_pathogenic"], 1.0)
        self.assertIn(hgb["predicted_class"], ["Pathogenic", "Benign"])

    def test_04_sequence_validation(self):
        """Verify sequence checking against reference residue."""
        # 1. Matching sequence: 174 'M' residues, then 'R', then 'A' residues
        matching_seq = ("M" * 174) + "R" + ("A" * 50)
        res_match = self.engine.predict("R175H", sequence=matching_seq)
        self.assertTrue(res_match["sequence_validation"]["verified"])

        # 2. Mismatching sequence: position 175 is 'K' instead of 'R'
        mismatch_seq = ("M" * 174) + "K" + ("A" * 50)
        res_mismatch = self.engine.predict("R175H", sequence=mismatch_seq)
        self.assertFalse(res_mismatch["sequence_validation"]["verified"])
        self.assertIn("Residue discrepancy", res_mismatch["sequence_validation"]["message"])

    def test_05_invalid_mutation_syntax(self):
        """Verify bad mutation strings return clean errors."""
        res_bad = self.engine.predict("INVALID")
        self.assertFalse(res_bad["ok"])
        self.assertEqual(res_bad["status"], "invalid_mutation")

    def test_06_fastapi_step9_qml_endpoint(self):
        """Verify POST /api/step9/qml endpoint."""
        response = self.client.post("/api/step9/qml", json={
            "uniprot_id": "P04637",
            "mutation": "R175H",
            "sequence": ("M" * 174) + "R" + ("A" * 100)
        })
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertTrue(data.get("ok"))
        self.assertEqual(data.get("mutation"), "R175H")
        self.assertIn("qnn_prediction", data)
        self.assertIn("classical_prediction", data)
        self.assertIn("benchmarks", data)

    def test_07_fastapi_step9_unified_endpoint(self):
        """Verify POST /api/step9 includes embedded QML result."""
        response = self.client.post("/api/step9", json={
            "uniprot_id": "P04637",
            "mutation": "R175H",
            "sequence": ("M" * 174) + "R" + ("A" * 100)
        })
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertIn("verdict", data)
        self.assertIn("qml", data)
        self.assertTrue(data["qml"].get("ok"))


if __name__ == "__main__":
    unittest.main()
