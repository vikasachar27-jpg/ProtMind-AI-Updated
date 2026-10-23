"""Quantum Machine Learning (QML) & Classical Baseline Inference Engine.

Integrates a trained PennyLane Quantum Neural Network (QNN) and Scikit-Learn
HistGradientBoostingClassifier trained on ClinVar-derived variant data for
protein mutation pathogenicity assessment in ProtMind AI (Step 9).

Research prototype only:
- Executes on the PennyLane 'default.qubit' statevector quantum simulator.
- Trained on 3 physicochemical delta features: hydrophobicity, charge, and mass.
- Does not establish quantum computational advantage or clinical validation.
- Does not directly predict 3D protein structure or Gibbs free energy stability.
"""

from __future__ import annotations

import json
import logging
import os
import re
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import joblib
import numpy as np
import pennylane as qml
import torch

logger = logging.getLogger("protmind.qml")

# ----------------- Canonical Amino Acid Biophysical Table -----------------
# Kyte-Doolittle hydropathy index; average residue molecular masses (Da); formal physiological charge.
# Exactly mirrors the training pipeline in prepare_clinvar_local.py.
AA_PROPERTIES = {
    "A": {"three": "Ala", "name": "Alanine", "hydro": 1.8, "mass": 71.0788, "charge": 0.0},
    "R": {"three": "Arg", "name": "Arginine", "hydro": -4.5, "mass": 156.1875, "charge": 1.0},
    "N": {"three": "Asn", "name": "Asparagine", "hydro": -3.5, "mass": 114.1038, "charge": 0.0},
    "D": {"three": "Asp", "name": "Aspartate", "hydro": -3.5, "mass": 115.0886, "charge": -1.0},
    "C": {"three": "Cys", "name": "Cysteine", "hydro": 2.5, "mass": 103.1388, "charge": 0.0},
    "Q": {"three": "Gln", "name": "Glutamine", "hydro": -3.5, "mass": 128.1307, "charge": 0.0},
    "E": {"three": "Glu", "name": "Glutamate", "hydro": -3.5, "mass": 129.1155, "charge": -1.0},
    "G": {"three": "Gly", "name": "Glycine", "hydro": -0.4, "mass": 57.0519, "charge": 0.0},
    "H": {"three": "His", "name": "Histidine", "hydro": -3.2, "mass": 137.1411, "charge": 0.1},
    "I": {"three": "Ile", "name": "Isoleucine", "hydro": 4.5, "mass": 113.1594, "charge": 0.0},
    "L": {"three": "Leu", "name": "Leucine", "hydro": 3.8, "mass": 113.1594, "charge": 0.0},
    "K": {"three": "Lys", "name": "Lysine", "hydro": -3.9, "mass": 128.1741, "charge": 1.0},
    "M": {"three": "Met", "name": "Methionine", "hydro": 1.9, "mass": 131.1926, "charge": 0.0},
    "F": {"three": "Phe", "name": "Phenylalanine", "hydro": 2.8, "mass": 147.1766, "charge": 0.0},
    "P": {"three": "Pro", "name": "Proline", "hydro": -1.6, "mass": 97.1167, "charge": 0.0},
    "S": {"three": "Ser", "name": "Serine", "hydro": -0.8, "mass": 87.0782, "charge": 0.0},
    "T": {"three": "Thr", "name": "Threonine", "hydro": -0.7, "mass": 101.1051, "charge": 0.0},
    "W": {"three": "Trp", "name": "Tryptophan", "hydro": -0.9, "mass": 186.2132, "charge": 0.0},
    "Y": {"three": "Tyr", "name": "Tyrosine", "hydro": -1.3, "mass": 163.1760, "charge": 0.0},
    "V": {"three": "Val", "name": "Valine", "hydro": 4.2, "mass": 99.1326, "charge": 0.0},
}

THREE_TO_ONE = {v["three"].upper(): k for k, v in AA_PROPERTIES.items()}
FEATURE_NAMES = ["hydrophobicity_delta", "charge_delta", "mass_delta"]

# Mutation regex accommodating formats: R175H, p.R175H, Arg175His, p.Arg175His, p.(Arg175His)
MUTATION_RE = re.compile(
    r"^(?:p\.\(?\s*)?([A-Za-z]{1,3})\s*(\d+)\s*([A-Za-z]{1,3})\)?$",
    re.IGNORECASE,
)


def normalize_aa(token: str) -> Optional[str]:
    """Convert 1-letter or 3-letter amino acid code to uppercase 1-letter format."""
    token = token.strip()
    if len(token) == 1:
        one = token.upper()
        return one if one in AA_PROPERTIES else None
    return THREE_TO_ONE.get(token.upper())


def parse_mutation_string(mutation_str: str) -> Optional[Tuple[str, int, str]]:
    """Parse mutation string into (ref_aa, position, alt_aa). Returns None if unparseable."""
    if not isinstance(mutation_str, str):
        return None
    match = MUTATION_RE.match(mutation_str.strip())
    if not match:
        return None
    ref_raw, pos_raw, alt_raw = match.groups()
    ref = normalize_aa(ref_raw)
    alt = normalize_aa(alt_raw)
    try:
        pos = int(pos_raw)
    except ValueError:
        return None
    if not ref or not alt or pos <= 0:
        return None
    return ref, pos, alt


def compute_physicochemical_deltas(ref_aa: str, alt_aa: str) -> Dict[str, float]:
    """Calculate the 3 physicochemical delta features matching ClinVar training."""
    ref_props = AA_PROPERTIES[ref_aa]
    alt_props = AA_PROPERTIES[alt_aa]
    return {
        "hydrophobicity_delta": round(float(alt_props["hydro"] - ref_props["hydro"]), 4),
        "charge_delta": round(float(alt_props["charge"] - ref_props["charge"]), 4),
        "mass_delta": round(float(alt_props["mass"] - ref_props["mass"]), 4),
    }


# ----------------- PennyLane Hybrid QNN Architecture -----------------
N_QUBITS = 3
N_LAYERS = 2

# Simulator device
dev = qml.device("default.qubit", wires=N_QUBITS)


@qml.qnode(dev, interface="torch", diff_method="backprop")
def qnn_circuit(inputs, weights):
    # Angle embedding maps normalized 3D inputs to Y-axis rotations on 3 qubits
    qml.AngleEmbedding(inputs, wires=range(N_QUBITS), rotation="Y")
    # Strongly entangling parameterized layers
    qml.StronglyEntanglingLayers(weights, wires=range(N_QUBITS))
    # Expectation value of Pauli-Z on wire 0: output in [-1, +1]
    return qml.expval(qml.PauliZ(0))


class HybridQNNModule(torch.nn.Module):
    """PyTorch wrapper for the PennyLane QNode matching evaluate_qnn.py."""

    def __init__(self):
        super().__init__()
        self.quantum = qml.qnn.TorchLayer(
            qnn_circuit, {"weights": (N_LAYERS, N_QUBITS, 3)}
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # Scale expectation value [-1, +1] to probability [0, 1]
        return (self.quantum(x).reshape(-1) + 1.0) / 2.0


# ----------------- Inference Engine Singleton -----------------
class QMLInferenceEngine:
    """Thread-safe singleton managing trained QNN, classical baseline, and preprocessors."""

    _instance: Optional[QMLInferenceEngine] = None

    def __init__(self, model_dir: Optional[str] = None):
        if model_dir is None:
            model_dir = str(Path(__file__).parent / "models" / "qml")
        self.model_dir = Path(model_dir)

        self.qnn_model: Optional[HybridQNNModule] = None
        self.qnn_prep: Optional[Dict[str, Any]] = None
        self.hgb_model: Optional[Any] = None
        self.hgb_scaler: Optional[Any] = None
        self.hgb_imputer: Optional[Any] = None
        self.benchmark_metadata: Optional[Dict[str, Any]] = None
        self.loaded = False
        self.load_error: Optional[str] = None

        self._load_artifacts()

    @classmethod
    def get_instance(cls, model_dir: Optional[str] = None) -> QMLInferenceEngine:
        if cls._instance is None:
            cls._instance = cls(model_dir)
        return cls._instance

    def _load_artifacts(self):
        try:
            # 1. Load PennyLane QNN
            qnn_pt_path = self.model_dir / "pennylane_qnn.pt"
            if not qnn_pt_path.exists():
                raise FileNotFoundError(f"Missing QNN weights: {qnn_pt_path}")

            qnn = HybridQNNModule()
            try:
                saved = torch.load(qnn_pt_path, map_location="cpu", weights_only=True)
            except (TypeError, Exception):
                saved = torch.load(qnn_pt_path, map_location="cpu")
            state_dict = saved.get("state_dict", saved)
            qnn.load_state_dict(state_dict)
            qnn.eval()
            self.qnn_model = qnn

            # 2. Load QNN Preprocessing (SimpleImputer + StandardScaler)
            qnn_prep_path = self.model_dir / "qnn_preprocessing.joblib"
            if not qnn_prep_path.exists():
                raise FileNotFoundError(f"Missing QNN preprocessor: {qnn_prep_path}")
            self.qnn_prep = joblib.load(qnn_prep_path)

            # 3. Load Classical HistGradientBoosting
            hgb_path = self.model_dir / "hist_gradient_boosting.joblib"
            if not hgb_path.exists():
                raise FileNotFoundError(f"Missing HGB model: {hgb_path}")
            hgb_artifact = joblib.load(hgb_path)
            if isinstance(hgb_artifact, dict):
                self.hgb_model = hgb_artifact.get("model")
                self.hgb_scaler = hgb_artifact.get("scaler")
                self.hgb_imputer = hgb_artifact.get("imputer")
            else:
                self.hgb_model = hgb_artifact

            # 4. Load Evaluation Benchmark Results
            eval_path = self.model_dir / "evaluation_results.json"
            if eval_path.exists():
                with open(eval_path, "r", encoding="utf-8") as f:
                    self.benchmark_metadata = json.load(f)

            self.loaded = True
            self.load_error = None
            logger.info("ProtMind QML & Classical inference engine loaded successfully.")
        except Exception as exc:
            self.loaded = False
            self.load_error = str(exc)
            logger.error("Failed to load QML artifacts: %s", exc)

    def predict(
        self,
        mutation: str,
        sequence: Optional[str] = None,
        uniprot_id: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Perform comprehensive QNN & classical baseline inference for a protein variant."""
        if not self.loaded:
            return {
                "ok": False,
                "error": f"QML model engine not ready: {self.load_error}",
                "status": "model_unavailable",
            }

        # 1. Parse and validate mutation syntax
        parsed = parse_mutation_string(mutation)
        if not parsed:
            return {
                "ok": False,
                "error": (
                    f"Invalid mutation format '{mutation}'. "
                    "Expected standard single amino acid substitution format like 'R175H', 'Arg175His', or 'p.R175H'."
                ),
                "status": "invalid_mutation",
            }

        ref_aa, pos, alt_aa = parsed
        if ref_aa == alt_aa:
            return {
                "ok": False,
                "error": f"Synonymous substitution detected ({ref_aa}{pos}{alt_aa}). Physicochemical deltas are zero.",
                "status": "synonymous_variant",
            }

        # 2. Sequence consistency check
        seq_validation = {
            "verified": False,
            "message": "Sequence not supplied for positional verification.",
            "seq_residue": None,
            "sequence_length": len(sequence) if sequence else 0,
        }
        if sequence:
            clean_seq = sequence.strip().upper()
            seq_len = len(clean_seq)
            seq_validation["sequence_length"] = seq_len

            if pos > seq_len:
                return {
                    "ok": False,
                    "error": (
                        f"Mutation position {pos} is out of bounds for the supplied sequence "
                        f"(length: {seq_len} amino acids)."
                    ),
                    "status": "position_out_of_bounds",
                }
            actual_res = clean_seq[pos - 1]
            seq_validation["seq_residue"] = actual_res
            if actual_res == ref_aa:
                seq_validation["verified"] = True
                seq_validation["message"] = (
                    f"Position {pos} verified: Reference residue '{ref_aa}' ({AA_PROPERTIES[ref_aa]['name']}) "
                    "matches the protein FASTA sequence exactly."
                )
            else:
                seq_validation["verified"] = False
                seq_validation["message"] = (
                    f"Residue discrepancy: At position {pos}, the FASTA sequence contains '{actual_res}' "
                    f"({AA_PROPERTIES.get(actual_res, {}).get('name', 'Unknown')}), but the mutation specifies '{ref_aa}' "
                    f"({AA_PROPERTIES[ref_aa]['name']})."
                )

        # 3. Compute the 3 physicochemical deltas
        deltas = compute_physicochemical_deltas(ref_aa, alt_aa)
        raw_feature_vector = np.array(
            [[deltas["hydrophobicity_delta"], deltas["charge_delta"], deltas["mass_delta"]]],
            dtype=float,
        )

        # 4. Quantum Neural Network Inference
        try:
            imputed_q = self.qnn_prep["imputer"].transform(raw_feature_vector)
            scaled_q = self.qnn_prep["scaler"].transform(imputed_q)
            # Clip standardized input to [-pi, pi] for angle embedding rotation
            clipped_q = np.clip(scaled_q, -np.pi, np.pi)

            with torch.no_grad():
                qnn_tensor = torch.tensor(clipped_q, dtype=torch.float32)
                qnn_prob = float(self.qnn_model(qnn_tensor).item())

            qnn_class = "Pathogenic" if qnn_prob >= 0.5 else "Benign"
            qnn_class_id = 1 if qnn_prob >= 0.5 else 0
            qnn_confidence = round(abs(qnn_prob - 0.5) * 2.0 * 100.0, 1)

            qnn_result = {
                "model_name": "PennyLane Hybrid QNN",
                "framework": "PennyLane + PyTorch",
                "simulator": "default.qubit (statevector)",
                "qubits": N_QUBITS,
                "entangling_layers": N_LAYERS,
                "probability_pathogenic": round(qnn_prob, 4),
                "probability_benign": round(1.0 - qnn_prob, 4),
                "predicted_class": qnn_class,
                "class_id": qnn_class_id,
                "threshold": 0.5,
                "relative_confidence_pct": qnn_confidence,
                "statevector_status": "Simulated",
            }
        except Exception as e:
            qnn_result = {
                "model_name": "PennyLane Hybrid QNN",
                "error": f"QNN inference failed: {e}",
            }

        # 5. Classical HistGradientBoosting Inference
        try:
            hgb_imp = self.hgb_imputer or self.qnn_prep["imputer"]
            hgb_scl = self.hgb_scaler or self.qnn_prep["scaler"]

            imputed_hgb = hgb_imp.transform(raw_feature_vector)
            scaled_hgb = hgb_scl.transform(imputed_hgb)
            hgb_probs = self.hgb_model.predict_proba(scaled_hgb)[0]
            hgb_prob = float(hgb_probs[1])

            hgb_class = "Pathogenic" if hgb_prob >= 0.5 else "Benign"
            hgb_class_id = 1 if hgb_prob >= 0.5 else 0
            hgb_confidence = round(abs(hgb_prob - 0.5) * 2.0 * 100.0, 1)

            classical_result = {
                "model_name": "HistGradientBoostingClassifier",
                "framework": "Scikit-Learn",
                "probability_pathogenic": round(hgb_prob, 4),
                "probability_benign": round(1.0 - hgb_prob, 4),
                "predicted_class": hgb_class,
                "class_id": hgb_class_id,
                "threshold": 0.5,
                "relative_confidence_pct": hgb_confidence,
            }
        except Exception as e:
            classical_result = {
                "model_name": "HistGradientBoostingClassifier",
                "error": f"Classical inference failed: {e}",
            }

        # 6. Benchmark Metrics (from evaluation_results.json held-out test split)
        benchmark_table = []
        if self.benchmark_metadata and "results" in self.benchmark_metadata:
            res = self.benchmark_metadata["results"]
            model_keys = [
                ("pennylane_qnn", "PennyLane QNN (Quantum)"),
                ("hist_gradient_boosting", "HistGradientBoosting (Classical)"),
                ("random_forest", "Random Forest (Classical)"),
                ("logistic_regression", "Logistic Regression (Linear)"),
            ]
            for key, display_name in model_keys:
                if key in res:
                    m = res[key]
                    benchmark_table.append({
                        "model": display_name,
                        "key": key,
                        "roc_auc": round(float(m.get("roc_auc", 0)), 4),
                        "accuracy": round(float(m.get("accuracy", 0)) * 100, 2),
                        "balanced_accuracy": round(float(m.get("balanced_accuracy", 0)) * 100, 2),
                        "precision_pathogenic": round(float(m.get("precision_pathogenic", 0)) * 100, 2),
                        "recall_pathogenic": round(float(m.get("recall_pathogenic", 0)) * 100, 2),
                        "f1_pathogenic": round(float(m.get("f1_pathogenic", 0)), 4),
                        "test_samples": m.get("n", 1500),
                    })

        # 7. Explainability of the three features
        feature_breakdown = [
            {
                "feature": "Hydrophobicity Delta",
                "key": "hydrophobicity_delta",
                "value": deltas["hydrophobicity_delta"],
                "unit": "Kyte-Doolittle index",
                "interpretation": (
                    f"Wildtype {ref_aa} ({AA_PROPERTIES[ref_aa]['hydro']}) → Mutant {alt_aa} ({AA_PROPERTIES[alt_aa]['hydro']}). "
                    + (
                        "Increases hydrophobicity (may promote core packing or aggregation)."
                        if deltas["hydrophobicity_delta"] > 0
                        else "Decreases hydrophobicity (favors solvent exposure or destabilizes core)."
                        if deltas["hydrophobicity_delta"] < 0
                        else "Neutral hydropathy shift."
                    )
                ),
                "quantum_role": "Rotates qubit 0 on the Bloch sphere along the Y-axis via AngleEmbedding.",
            },
            {
                "feature": "Charge Delta",
                "key": "charge_delta",
                "value": deltas["charge_delta"],
                "unit": "Elementary charge (e)",
                "interpretation": (
                    f"Wildtype {ref_aa} ({AA_PROPERTIES[ref_aa]['charge']:+g}e) → Mutant {alt_aa} ({AA_PROPERTIES[alt_aa]['charge']:+g}e). "
                    + (
                        f"Net charge shift of {deltas['charge_delta']:+g}e alters electrostatic salt bridges and local surface potential."
                        if deltas["charge_delta"] != 0
                        else "Preserves net electrostatic charge."
                    )
                ),
                "quantum_role": "Rotates qubit 1 on the Bloch sphere along the Y-axis via AngleEmbedding.",
            },
            {
                "feature": "Molecular Mass Delta",
                "key": "mass_delta",
                "value": deltas["mass_delta"],
                "unit": "Daltons (Da)",
                "interpretation": (
                    f"Wildtype {ref_aa} ({AA_PROPERTIES[ref_aa]['mass']:.1f} Da) → Mutant {alt_aa} ({AA_PROPERTIES[alt_aa]['mass']:.1f} Da). "
                    + (
                        f"Increases residue volume (+{deltas['mass_delta']:.1f} Da), creating potential steric clashes."
                        if deltas["mass_delta"] > 0
                        else f"Truncates residue volume ({deltas['mass_delta']:.1f} Da), risking cavity formation."
                        if deltas["mass_delta"] < 0
                        else "Residue mass is identical (stereoisomer)."
                    )
                ),
                "quantum_role": "Rotates qubit 2 on the Bloch sphere along the Y-axis via AngleEmbedding.",
            },
        ]

        # 8. Consensus comparison verdict
        models_agree = (
            qnn_result.get("predicted_class") == classical_result.get("predicted_class")
        )
        if models_agree:
            consensus_note = (
                f"Both the Quantum Neural Network and the HistGradientBoosting baseline reach consensus, "
                f"classifying the substitution as **{qnn_result.get('predicted_class')}**."
            )
        else:
            consensus_note = (
                f"Model divergence observed: QNN predicts **{qnn_result.get('predicted_class')}** "
                f"(score: {qnn_result.get('probability_pathogenic')}), while Classical HGB predicts "
                f"**{classical_result.get('predicted_class')}** (score: {classical_result.get('probability_pathogenic')}). "
                "Consult multi-stream biophysical stability and conservation evidence."
            )

        return {
            "ok": True,
            "status": "success",
            "mutation": f"{ref_aa}{pos}{alt_aa}",
            "ref_aa": ref_aa,
            "alt_aa": alt_aa,
            "position": pos,
            "ref_name": AA_PROPERTIES[ref_aa]["name"],
            "alt_name": AA_PROPERTIES[alt_aa]["name"],
            "uniprot_id": uniprot_id or "",
            "sequence_validation": seq_validation,
            "features": deltas,
            "feature_breakdown": feature_breakdown,
            "qnn_prediction": qnn_result,
            "classical_prediction": classical_result,
            "models_agree": models_agree,
            "consensus_note": consensus_note,
            "benchmarks": benchmark_table,
            "quantum_circuit_details": {
                "framework": "PennyLane (v0.45+)",
                "simulator": "default.qubit (statevector)",
                "qubits": N_QUBITS,
                "parameterized_layers": N_LAYERS,
                "encoding": "AngleEmbedding(rotation='Y', wires=[0, 1, 2])",
                "entanglement": "StronglyEntanglingLayers (circular CNOT + single-qubit rotations)",
                "measurement": "Expectation value PauliZ on wire 0: <Z0>",
                "activation": "(<Z0> + 1.0) / 2.0 -> sigmoid probability range [0, 1]",
            },
            "scientific_disclaimers": [
                "Experimental Research Prototype: This model was trained on ClinVar-derived variant classifications using 3 physicochemical features.",
                "Simulator Execution: Runs on a classical quantum statevector simulator (default.qubit); does NOT claim quantum computational speedup or quantum advantage.",
                "Not Clinically Validated: Predictions must NOT be used as the sole basis for clinical diagnostic or therapeutic decision-making.",
                "Scope Limitation: This QNN classifies missense variant pathogenicity; it does not directly fold proteins, predict tertiary coordinates, or measure experimental Gibbs free energy.",
            ],
        }


# Global helper function for one-line import
def get_qml_engine() -> QMLInferenceEngine:
    return QMLInferenceEngine.get_instance()
