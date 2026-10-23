"""ProtMind AI - FastAPI layer.

Provides REST API endpoints connecting the React frontend to biophysical calculations,
UniProtKB domain extraction, STRING database pathways, CHARMM potential energy models,
Kyte-Doolittle / Zamyatnin stability thermodynamics, and Ramachandran validation.
"""
import json
from typing import Optional

from fastapi import FastAPI, File, Form, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, Response
from pydantic import BaseModel

import backend as bk
from qml_inference import get_qml_engine

# Pre-load QML and Classical models once at server startup
_qml_engine = get_qml_engine()

app = FastAPI(title="ProtMind AI API")
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
    expose_headers=["Content-Disposition"],
)


class Payload(BaseModel):
    sequence: str
    uniprot_id: str
    mutation: str
    source: Optional[str] = ""
    custom_pattern: Optional[str] = ""


class MotifQuery(BaseModel):
    sequence: str
    uniprot_id: Optional[str] = ""
    mutation: Optional[str] = ""
    custom_pattern: Optional[str] = ""


def _mut_pos(mutation: str):
    try:
        return int(mutation[1:-1])
    except (ValueError, IndexError):
        return None


# ---------------- Root status endpoint ----------------
@app.get("/")
def root():
    return {
        "status": "online",
        "service": "ProtMind AI Backend API",
        "version": "2.0",
        "health": "/api/health",
        "documentation": "/docs"
    }


# ---------------- STEP 1: User Input ----------------
@app.post("/api/step1/validate")
async def step1(
    uniprot_id: str = Form(""),
    mutation: str = Form(""),
    sequence_input: str = Form(""),
    file: Optional[UploadFile] = File(None),
):
    errors, final_sequence, source_type = [], "", ""

    if file is not None and file.filename:                      # Priority 1: file upload
        try:
            raw = (await file.read()).decode("utf-8")
            final_sequence = bk.parse_fasta(raw)
            source_type = f"Uploaded file ({file.filename})"
        except Exception:
            errors.append("Could not read the uploaded FASTA file. Ensure it is UTF-8 encoded text.")
    elif sequence_input.strip():                                # Priority 2: text area
        final_sequence = bk.parse_fasta(sequence_input)
        source_type = "Manual input"
    else:
        errors.append("Please upload a FASTA file or paste a protein sequence.")

    if final_sequence and not bk.validate_sequence(final_sequence):
        errors.append("Invalid amino acid characters detected. Only standard 20 amino acids allowed.")
    if not uniprot_id or not bk.validate_uniprot_id(uniprot_id):
        errors.append("Invalid UniProt ID format. Example: P04637")
    if not mutation or not bk.validate_mutation(mutation):
        errors.append("Invalid Mutation format. Example: R175H")

    if errors:
        return {"ok": False, "errors": errors}
    return {
        "ok": True,
        "payload": {
            "sequence": final_sequence.strip().upper(),
            "uniprot_id": uniprot_id.strip().upper(),
            "mutation": mutation.strip().upper(),
            "source": source_type,
        },
        "length": len(final_sequence),
    }


# ---------------- STEP 2: Data Retrieval ----------------
@app.post("/api/step2")
def step2(p: Payload):
    try:
        uniprot = bk.fetch_uniprot_data(p.uniprot_id)
        af_url = bk.fetch_alphafold_url(p.uniprot_id)
    except Exception as e:
        uniprot, af_url = {"Error": f"Could not reach UniProt/AlphaFold ({type(e).__name__})"}, None
    return {
        "uniprot": uniprot,
        "alphafold_url": af_url,
        "mutation": p.mutation,
        "gnomad": {"European": 0.005, "African": 0.001, "Asian": 0.003, "Global": 0.004},
    }


# ---------------- STEP 3: Preprocessing & Motif/Domain Search ----------------
@app.post("/api/step3")
def step3(p: Payload):
    features = bk.extract_features(p.sequence, p.mutation)
    mutant = list(p.sequence)
    if features["Positional Match"]:
        mutant[features["Position"] - 1] = features["Mutant AA"]
    score = bk.calculate_alignment_score(p.sequence, "".join(mutant))
    
    # Run Motif and Domain searching through the protein sequence
    motifs_data = bk.search_motifs_and_domains(
        sequence=p.sequence,
        uniprot_id=p.uniprot_id,
        mutation=p.mutation,
        custom_pattern=p.custom_pattern or ""
    )

    return {
        "features": features,
        "alignment_score": score,
        "motif_domain_analysis": motifs_data
    }


@app.post("/api/step3/motifs")
def search_motifs(q: MotifQuery):
    """Real-time custom regex / motif pattern scanning."""
    return bk.search_motifs_and_domains(
        sequence=q.sequence,
        uniprot_id=q.uniprot_id or "",
        mutation=q.mutation or "",
        custom_pattern=q.custom_pattern or ""
    )


# ---------------- STEP 4: Mutation Analysis & STRING Pathways ----------------
@app.post("/api/step4")
def step4(p: Payload):
    string_data = bk.fetch_string_pathways_and_ppi(p.uniprot_id, limit=8)
    return {
        "ppi": string_data["ppi"],
        "pathways": string_data["pathways"],
        "network_graph": string_data["network_graph"],
        "diagram_urls": string_data["diagram_urls"],
        "pathogenicity": bk.predict_pathogenicity(p.mutation),
    }


# ---------------- STEP 5: Structure Prediction ----------------
@app.post("/api/step5")
def step5(p: Payload):
    try:
        pdb = bk.fetch_pdb_content(bk.fetch_alphafold_url(p.uniprot_id))
    except Exception:
        pdb = ""
    return {
        "pdb": pdb,
        "mut_pos": _mut_pos(p.mutation),
        "mutation": p.mutation,
        "uniprot_id": p.uniprot_id,
    }


# ---------------- STEP 6: Thermodynamic Stability & Potential Energy ----------------
@app.post("/api/step6")
def step6(p: Payload):
    f = bk.extract_features(p.sequence, p.mutation)
    pos = _mut_pos(p.mutation) or 1

    # Fetch Phi/Psi if structural coordinates are accessible
    phi, psi = None, None
    try:
        pdb = bk.fetch_pdb_content(bk.fetch_alphafold_url(p.uniprot_id))
        angles = bk.calculate_phi_psi(pdb, pos)
        if angles:
            phi, psi = angles.get("phi"), angles.get("psi")
    except Exception:
        pass

    stab = bk.predict_thermodynamic_stability(
        wt_aa=f["Wildtype AA"],
        mut_aa=f["Mutant AA"],
        wt_property=f["WT Property"],
        mut_property=f["Mutant Property"],
        phi=phi,
        psi=psi
    )
    return {"features": f, "stability": stab}


# ---------------- STEP 7: Stereochemical Validation (Ramachandran) ----------------
@app.post("/api/step7")
def step7(p: Payload):
    pos = _mut_pos(p.mutation) or 1
    try:
        pdb = bk.fetch_pdb_content(bk.fetch_alphafold_url(p.uniprot_id))
    except Exception:
        pdb = ""
    angles = bk.calculate_phi_psi(pdb, pos)
    if not angles:
        return {"angles": None}

    # Extract all residues for authentic PROCHECK plot & stats
    all_residues = bk.extract_all_phi_psi(pdb)

    # Classify angles against Ramachandran regions
    sec_struct, steric_feas, status_code = bk.classify_ramachandran(angles["phi"], angles["psi"])
    
    # 1. Primary: PROCHECK / SAVES Plot (matching user's uploaded reference)
    fig_procheck = bk.generate_procheck_ramachandran_plot(
        phi=angles["phi"],
        psi=angles["psi"],
        mutation=p.mutation,
        all_residues=all_residues
    )
    
    # 2. Reference: Grayscale density plot (Image 5)
    fig_density = bk.generate_ramachandran_plot(
        phi=angles["phi"],
        psi=angles["psi"],
        mutation=p.mutation
    )

    # 3. PROCHECK statistics distribution
    stats = bk.calculate_procheck_statistics(all_residues)
    target_zone = bk.lookup_procheck_zone(angles["phi"], angles["psi"])
    zone_names = {
        3: "Most Favoured [A,B,L]",
        2: "Additional Allowed [a,b,l,p]",
        1: "Generously Allowed [~a,~b,~l,~p]",
        0: "Disallowed Region"
    }

    return {
        "angles": angles,
        "classification": {
            "secondary_structure": sec_struct,
            "steric_feasibility": steric_feas,
            "status_code": status_code,
            "procheck_zone": zone_names.get(target_zone, "Allowed")
        },
        "figure": json.loads(fig_procheck.to_json()),
        "figure_procheck": json.loads(fig_procheck.to_json()),
        "figure_density": json.loads(fig_density.to_json()),
        "statistics": stats
    }


# ---------------- STEP 8: Therapeutic Targeting ----------------
@app.post("/api/step8")
def step8(p: Payload):
    try:
        gene = bk.fetch_uniprot_data(p.uniprot_id).get("Gene", "Unknown")
    except Exception:
        gene = "Unknown"
    return {"gene": gene, "drugs": bk.fetch_chembl_drugs(p.uniprot_id)}


@app.post("/api/step8/drugs")
def step8_drugs(p: Payload):
    return step8(p)


# ---------------- STEP 9: Explainable AI ----------------
@app.post("/api/step9")
def step9(p: Payload):
    return bk.generate_explainable_ai(p.model_dump())


@app.post("/api/step9/explain")
def step9_explain(p: Payload):
    return step9(p)


# Backward-compatibility alias for legacy clients expecting drug endpoint on step 9
@app.post("/api/step9/drugs")
def step9_drugs_compat(p: Payload):
    return step8(p)


# ---------------- STEP 9: Quantum Machine Learning (QML) and Classical Comparison ----------------
@app.post("/api/step9/qml")
def step9_qml(p: Payload):
    """
    Inference endpoint for Quantum Neural Network (PennyLane default.qubit)
    and HistGradientBoosting classical baseline on ClinVar-derived physicochemical deltas.
    """
    engine = get_qml_engine()
    return engine.predict(mutation=p.mutation, sequence=p.sequence, uniprot_id=p.uniprot_id)


@app.post("/api/qml/predict")
def qml_predict(p: Payload):
    """Alias for direct QNN and classical baseline inference."""
    return step9_qml(p)


@app.post("/api/qml/compare")
def qml_compare(p: Payload):
    """Comparative benchmarking metrics and predictions."""
    return step9_qml(p)


# ---------------- STEP 10: Automated Clinical Report (PDF) ----------------
@app.post("/api/step10/pdf")
def step10(p: Payload):
    pdf = bk.generate_pdf_report(p.model_dump())
    name = f"{p.uniprot_id}_{p.mutation}_Automated_Clinical_Report.pdf"
    return Response(pdf, media_type="application/pdf",
                    headers={"Content-Disposition": f'attachment; filename="{name}"'})


@app.get("/api/health")
def health():
    return JSONResponse({"status": "ok"})


# Mount built React frontend for production full-stack deployment
import os
from fastapi.staticfiles import StaticFiles

_dist_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "frontend", "dist"))
if os.path.exists(_dist_dir):
    app.mount("/", StaticFiles(directory=_dist_dir, html=True), name="frontend")

