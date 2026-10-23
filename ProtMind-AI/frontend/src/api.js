// Thin client for the FastAPI backend (see backend/api.py).

function resolveApiBase() {
  if (import.meta.env.VITE_API_URL) {
    return import.meta.env.VITE_API_URL.replace(/\/+$/, "");
  }
  // When running locally in browser (e.g. localhost:5173), prioritize local backend on port 8000
  if (typeof window !== "undefined") {
    const host = window.location.hostname;
    if (host === "localhost" || host === "127.0.0.1" || host === "0.0.0.0" || host === "::1") {
      return "http://localhost:8000";
    }
  }
  return "https://protmind-ai-react.onrender.com";
}

let BASE = resolveApiBase();
const FALLBACK_BASE = (BASE.includes("localhost") || BASE.includes("127.0.0.1"))
  ? "https://protmind-ai-react.onrender.com"
  : "http://localhost:8000";

async function request(path, options) {
  let url = `${BASE}${path}`;
  try {
    const res = await fetch(url, options);
    // If primary server returns 404 on an upgraded endpoint, attempt fallback server
    if (!res.ok && res.status === 404 && FALLBACK_BASE && BASE !== FALLBACK_BASE) {
      try {
        const altRes = await fetch(`${FALLBACK_BASE}${path}`, options);
        if (altRes.ok) {
          BASE = FALLBACK_BASE;
          return altRes;
        }
      } catch {
        // Fallback network failure, proceed to standard error
      }
    }
    if (!res.ok) throw new Error(`Server error (${res.status}) on ${path}`);
    return res;
  } catch (err) {
    // If primary server failed to connect (e.g. local backend wasn't running), attempt fallback
    if (FALLBACK_BASE && BASE !== FALLBACK_BASE) {
      try {
        const altRes = await fetch(`${FALLBACK_BASE}${path}`, options);
        if (altRes.ok) {
          BASE = FALLBACK_BASE;
          return altRes;
        }
      } catch {
        // Ignore fallback failure
      }
    }
    throw new Error(`Failed to reach backend (${err.message}) on ${path}`);
  }
}

export async function postJson(path, payload) {
  const res = await request(path, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  return res.json();
}

export async function validateInput({ uniprotId, mutation, sequence, file }) {
  const fd = new FormData();
  fd.append("uniprot_id", uniprotId);
  fd.append("mutation", mutation);
  fd.append("sequence_input", sequence);
  if (file) fd.append("file", file);
  const res = await request("/api/step1/validate", { method: "POST", body: fd });
  return res.json();
}

export async function fetchReportPdf(payload) {
  const res = await request("/api/step10/pdf", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  return res.blob();
}

export async function fetchQmlPrediction(payload) {
  return postJson("/api/step9/qml", payload);
}
