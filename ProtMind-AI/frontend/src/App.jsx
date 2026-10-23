import { useState, Component } from "react";
import "./App.css";
import { Gate, Step1, Step2, Step3, Step4, Step5, Step6, Step7, Step8, Step9, Step10 } from "./Steps";
const heroProtein = "/hero.png";

const workflow = [
  {
    n: "01", title: "User Input & Authentication", short: "Secure research workspace", icon: "◉",
    desc: "Create an account, verify the user, and provide a UniProt accession or protein sequence together with the mutation and analysis options."
  },
  {
    n: "02", title: "Data Retrieval", short: "External databases & APIs", icon: "◌",
    desc: "Retrieve protein sequences, structures, variants, population evidence, interactions, literature and drug-related data from trusted sources.",
    chips: ["UniProt", "PDB", "AlphaFold DB", "ClinVar", "dbSNP", "gnomAD", "STRING", "ChEMBL", "PubMed"]
  },
  {
    n: "03", title: "Data Preprocessing & Feature Extraction", short: "Clean, standardize & extract", icon: "✦",
    desc: "Validate sequences, map mutations, perform multiple sequence alignment and extract evolutionary and physicochemical features."
  },
  {
    n: "04", title: "Protein Mutation Analysis", short: "Pathogenicity & functional impact", icon: "◇",
    desc: "Integrate mutation predictors and evidence to assess pathogenicity, functional impact, conservation, population frequency and clinical relevance.",
    chips: ["AlphaMissense", "EVE", "SIFT", "PolyPhen-2"]
  },
  {
    n: "05", title: "Protein Structure Prediction & Visualization", short: "3D structure & mutation context", icon: "⌬",
    desc: "Retrieve or predict a high-confidence 3D protein structure and visualize the mutation in its structural context.",
    chips: ["AlphaFold2", "PDB", "Interactive 3D"]
  },
  {
    n: "06", title: "Thermodynamic & Structural Stability", short: "Potential energy & ΔG", icon: "∿",
    desc: "Evaluate energetic components and Gibbs free-energy changes to quantify structural stability changes caused by mutation.",
    chips: ["Bond energy", "Angle energy", "Dihedral", "Van der Waals", "Electrostatic", "ΔG"]
  },
  {
    n: "07", title: "Stereochemical Validation", short: "Ramachandran & geometry", icon: "⌁",
    desc: "Validate structural quality using stereochemical analysis, including favored, allowed and outlier regions in the Ramachandran plot."
  },
  {
    n: "08", title: "Therapeutic Targeting", short: "Drug & target exploration", icon: "💊",
    desc: "Connect mutation findings with drug targets, drug–target interactions, pathways, clinical evidence and therapeutic candidates.",
    chips: ["ChEMBL", "Pharmacological mechanisms", "Target druggability"]
  },
  {
    n: "09", title: "Explainable Intelligence Summary", short: "Multi-evidence mechanistic reasoning", icon: "🧠",
    desc: "Translate computational predictions and biophysical evidence into transparent mechanistic insights explaining how and why the verdict was issued.",
    chips: ["Consensus verdict", "Quantum ML (QNN)", "Biophysical drivers", "Evidence matrix", "Limitations"]
  },
  {
    n: "10", title: "Automated Clinical Report", short: "Downloadable self-contained report", icon: "📄",
    desc: "Bring the analysis together into a self-contained, traceable clinical report covering Steps 1–9 outputs, completeness check, structural plots, and references.",
    chips: ["PDF", "Steps 1–9 outputs", "Completeness audit", "References"]
  },
];

const stepPages = [
  { title: "User Input", icon: "⚡", sub: "Provide the protein, mutation and FASTA sequence, then validate and save.", Comp: Step1 },
  { title: "Data Retrieval Results", icon: "🔍", sub: "Live UniProt, NIH, gnomAD and AlphaFold evidence for your target.", Comp: Step2 },
  { title: "Preprocessing & Motif/Domain Search", icon: "⚙️", sub: "Physicochemical features, evolutionary alignment, and sequence motif/domain search.", Comp: Step3 },
  { title: "Mutation Analysis & STRING Pathways", icon: "🧬", sub: "Pathogenicity predictors and interactive STRING database pathway network.", Comp: Step4 },
  { title: "Protein Structure Prediction", icon: "🧩", sub: "Interactive AlphaFold 3D structure with the mutation highlighted.", Comp: Step5 },
  { title: "Thermodynamic Stability (ΔG & Potential Energy)", icon: "🔥", sub: "Gibbs free energy (ΔG & ΔΔG), CHARMM potential energy and literature calculations.", Comp: Step6 },
  { title: "Stereochemical Validation (Ramachandran)", icon: "📐", sub: "Backbone Phi/Psi angles on the Ramachandran plot with steric validation.", Comp: Step7 },
  { title: "Therapeutic Targeting", icon: "💊", sub: "ChEMBL drugs, mechanisms, and pharmacological candidates for the target.", Comp: Step8 },
  { title: "Explainable Intelligence Summary", icon: "🧠", sub: "Automated consensus verdict, Quantum ML (QNN), and biophysical drivers.", Comp: Step9 },
  { title: "Automated Clinical Report", icon: "📄", sub: "Self-contained clinical report with complete pipeline outputs and completeness check.", Comp: Step10 },
];

function App() {
  const [page, setPage] = useState("welcome");
  const [authMode, setAuthMode] = useState("login");
  const [expanded, setExpanded] = useState(1);
  const [user, setUser] = useState(null);
  const [payload, setPayload] = useState(null); // was st.session_state['payload']

  const openAuth = (mode) => {
    setAuthMode(mode);
    setPage(mode === "login" ? "login" : "register");
    window.scrollTo({ top: 0, behavior: "smooth" });
  };

  const handleAuth = (payload) => {
    setUser({ name: payload.name || payload.email.split("@")[0], email: payload.email });
    setPage("dashboard");
    window.scrollTo({ top: 0, behavior: "smooth" });
  };

  const startAnalysis = () => {
    if (!user) {
      openAuth("login");
      return;
    }
    setPage("analysis");
    window.scrollTo({ top: 0, behavior: "smooth" });
  };

  return (
    <div className="app-shell">
      <div className="bg-orb orb-one" />
      <div className="bg-orb orb-two" />
      <div className="bg-orb orb-three" />
      <Header page={page} user={user} onAuth={openAuth} onHome={() => setPage("welcome")} onDashboard={() => setPage("dashboard")} onLogout={() => { setUser(null); setPayload(null); setPage("welcome"); }} />

      {page === "welcome" && (
        <WelcomePage expanded={expanded} setExpanded={setExpanded} onAuth={openAuth} onStart={startAnalysis} />
      )}
      {page === "login" && <AuthPage mode="login" onAuth={handleAuth} onSwitch={() => openAuth("register")} onHome={() => setPage("welcome")} />}
      {page === "register" && <AuthPage mode="register" onAuth={handleAuth} onSwitch={() => openAuth("login")} onHome={() => setPage("welcome")} />}
      {page === "dashboard" && <Dashboard user={user} onStart={startAnalysis} onWorkflow={() => setPage("welcome")} />}
      {page === "analysis" && <AnalysisWorkspace payload={payload} setPayload={setPayload} onBack={() => setPage("dashboard")} />}

      <footer className="footer">
        <div><strong>ProtMind <em>AI</em></strong><span>Explainable protein mutation intelligence</span></div>
        <span>From mutation to mechanism to medicine.</span>
      </footer>
    </div>
  );
}

function Header({ page, user, onAuth, onHome, onDashboard, onLogout }) {
  return (
    <header className="site-header">
      <button className="brand" onClick={onHome} aria-label="ProtMind AI home">
        <span className="brand-mark brand-logo"><img src="/protmind-logo.jpeg" alt="ProtMind AI logo" /></span>
        <span><strong>ProtMind <em>AI</em></strong><small>EXPLAINABLE PROTEIN MUTATION INTELLIGENCE</small></span>
      </button>
      <nav>
        <button className={page === "welcome" ? "active" : ""} onClick={onHome}>Home</button>
        <button onClick={onHome}>Workflow</button>
        {user ? <button className={page === "dashboard" ? "active" : ""} onClick={onDashboard}>Dashboard</button> : null}
      </nav>
      <div className="header-actions">
        {user ? (
          <div className="user-menu"><span>{user.name.charAt(0).toUpperCase()}</span><div><b>{user.name}</b><small>Researcher</small></div><button onClick={onLogout}>Log out</button></div>
        ) : (
          <><button className="header-login" onClick={() => onAuth("login")}>Login</button><button className="header-cta" onClick={() => onAuth("register")}>Create Account</button></>
        )}
      </div>
    </header>
  );
}

function WelcomePage({ expanded, setExpanded, onAuth, onStart }) {
  return (
    <main className="welcome-page">
      <section className="hero-section">
        <div className="hero-copy">
          <span className="eyebrow">EXPLAINABLE PROTEIN MUTATION INTELLIGENCE</span>
          <h1>From <span>mutation</span><br />to <strong>mechanism.</strong></h1>
          <p className="hero-lead">An evidence-driven platform for understanding protein mutations, structural consequences, energetic stability and therapeutic possibilities.</p>
          <div className="hero-actions"><button className="primary-btn" onClick={() => onAuth("register")}>Explore ProtMind AI <span>→</span></button></div>
          <div className="hero-points"><span>✓ Explainable analysis</span><span>✓ Evidence-driven workflow</span><span>✓ Research-ready reports</span></div>
        </div>
        <div className="hero-art">
          <div className="art-glow" />
          <img src={heroProtein} alt="Protein molecular visualization" />
          <div className="floating-tag tag-one"><b>10</b><span>analysis stages</span></div>
          <div className="floating-tag tag-two"><b>ΔG</b><span>stability insight</span></div>
          <div className="floating-tag tag-three"><b>AI</b><span>explainable evidence</span></div>
        </div>
      </section>

      <section className="capability-section">
        <div className="section-heading"><div><span className="eyebrow">WHY PROTMIND AI</span><h2>From biological evidence to <em>understandable insight.</em></h2></div><p>ProtMind AI connects sequence, structure, mutation, energetic, population and therapeutic evidence in one research workflow.</p></div>
        <div className="capability-grid">
          {[
            ["🧬", "Understand Mutations", "Assess mutation effects, pathogenicity and functional consequences."],
            ["🔬", "Study Protein Structure", "Visualize structures and locate mutations in their molecular context."],
            ["⚡", "Analyze Stability", "Explore potential-energy components and Gibbs free-energy changes."],
            ["💊", "Map Therapeutics", "Connect findings with drug targets, interactions and clinical evidence."],
            ["🧠", "Explain Mechanisms", "Turn computational predictions and biophysical data into transparent mechanistic insights."],
            ["📄", "Generate Reports", "Compile traceable results, evidence and references into a self-contained report."],
          ].map(([icon,title,desc]) => <div className="capability-card" key={title}><span>{icon}</span><h3>{title}</h3><p>{desc}</p></div>)}
        </div>
      </section>

      <section className="workflow-section" id="workflow">
        <div className="section-heading workflow-heading"><div><span className="eyebrow">PROTMIND AI · ANALYSIS PIPELINE</span><h2>The complete <em>10-step workflow.</em></h2></div><p>Each stage adds another layer of biological evidence, leading from a mutation input to a comprehensive, explainable research report.</p></div>
        <div className="workflow-grid">
          {workflow.map((item, index) => {
            const open = expanded === index;
            return <article className={`workflow-card ${open ? "expanded" : ""}`} key={item.n}>
              <button className="workflow-card-head" onClick={() => setExpanded(open ? -1 : index)}>
                <span className="workflow-number">{item.n}</span><span className="workflow-icon">{item.icon}</span><span className="workflow-title"><strong>{item.title}</strong><small>{item.short}</small></span><span className="workflow-toggle">{open ? "−" : "+"}</span>
              </button>
              {open && <div className="workflow-card-body"><p>{item.desc}</p>{item.chips && <div className="chips">{item.chips.map(c => <span key={c}>{c}</span>)}</div>}<button className="step-link" onClick={onStart}>Open analysis workspace →</button></div>}
            </article>;
          })}
        </div>
      </section>

      <section className="cta-section"><div><span className="eyebrow">READY TO BEGIN?</span><h2>Start with a protein. End with an explanation.</h2><p>Use Login or Create Account at the top right to open your research workspace and begin the analysis pipeline.</p></div></section>
    </main>
  );
}

function AuthPage({ mode, onAuth, onSwitch, onHome }) {
  const [form, setForm] = useState({ name: "", email: "", password: "", confirm: "", human: false });
  const [error, setError] = useState("");
  const isLogin = mode === "login";
  const update = e => setForm(p => ({ ...p, [e.target.name]: e.target.type === "checkbox" ? e.target.checked : e.target.value }));
  const submit = e => {
    e.preventDefault();
    if (!form.email || !form.password) return setError("Please enter your email and password.");
    if (!form.human) return setError("Please confirm that you are human.");
    if (!isLogin && (!form.name || !form.confirm)) return setError("Please complete all account fields.");
    if (!isLogin && form.password !== form.confirm) return setError("Passwords do not match.");
    setError(""); onAuth(form);
  };
  return (
    <main className="auth-page">
      <button className="back-link" onClick={onHome}>← Back to ProtMind AI</button>
      <section className="auth-shell">
        <div className="auth-showcase"><span className="eyebrow">PROTMIND AI</span><h1>From mutation<br />to <em>mechanism.</em></h1><p>A secure workspace for protein mutation intelligence, structural analysis and therapeutic discovery.</p><div className="showcase-points"><span><b>01</b> Secure research workspace</span><span><b>02</b> Evidence-driven analysis</span><span><b>03</b> Explainable computational insights</span></div><div className="auth-orbit" /></div>
        <form className="auth-panel" onSubmit={submit}>
          <div className="auth-icon">↗</div><span className="eyebrow">RESEARCH WORKSPACE</span><h2>{isLogin ? "Welcome back" : "Create your account"}</h2><p>{isLogin ? "Sign in to continue to your protein analysis workspace." : "Create a secure account for your research workspace."}</p>
          {!isLogin && <AuthField name="name" label="Full name" value={form.name} onChange={update} placeholder="e.g. Vikas K R" />}
          <AuthField name="email" label="Email" value={form.email} onChange={update} placeholder="researcher@example.com" type="email" />
          <AuthField name="password" label="Password" value={form.password} onChange={update} placeholder={isLogin ? "Enter your password" : "Create a password"} type="password" />
          {!isLogin && <AuthField name="confirm" label="Confirm password" value={form.confirm} onChange={update} placeholder="Repeat your password" type="password" />}
          <label className="human-check"><input name="human" type="checkbox" checked={form.human} onChange={update}/><span><b>I'm human</b><small>Human verification for research workspace access</small></span><strong>CAPTCHA</strong></label>
          {error && <div className="form-error">{error}</div>}
          <button className="auth-submit" type="submit">{isLogin ? "Login to ProtMind AI" : "Create ProtMind AI Account"} <span>→</span></button>
          <small className="secure-note">🔒 Secure authentication</small>
          <div className="switch-auth">{isLogin ? "Don't have an account?" : "Already have an account?"}<button type="button" onClick={onSwitch}>{isLogin ? "Create Account" : "Login"}</button></div>
        </form>
      </section>
    </main>
  );
}

function AuthField({ name, label, value, onChange, placeholder, type = "text" }) { return <label className="auth-field"><span>{label}</span><input name={name} value={value} onChange={onChange} placeholder={placeholder} type={type}/></label>; }

function Dashboard({ user, onStart, onWorkflow }) {
  return <main className="dashboard-page">
    <section className="dashboard-hero dashboard-hero-pro">
      <div>
        <span className="eyebrow">PROTMIND AI · RESEARCH WORKSPACE</span>
        <h1>Welcome to ProtMind AI.</h1>
        <p>Explore protein mutations through an evidence-driven workflow from sequence to mechanism to medicine.</p>
      </div>
      <div className="dashboard-status"><span className="status-dot"/> <b>Workspace ready</b><small>10 analysis stages available</small></div>
    </section>

    <section className="dashboard-stat-grid">
      <div><span>01</span><b>Protein input</b><small>UniProt or FASTA</small></div>
      <div><span>06</span><b>Analysis modules</b><small>Selectable insights</small></div>
      <div><span>10</span><b>Workflow stages</b><small>From input to report</small></div>
      <div><span>∞</span><b>Research sessions</b><small>Ready for future history</small></div>
    </section>

    <section className="welcome-panel dashboard-welcome-panel">
      <div className="welcome-panel-copy">
        <span className="eyebrow">YOUR PROTEIN INTELLIGENCE WORKSPACE</span>
        <h2>What can you do here?</h2>
        <p>Start with a protein and mutation, then move through structural, energetic, stereochemical, explainable-AI and therapeutic evidence in one connected workspace.</p>
        <div className="mini-grid">
          <span>🧬 <b>Mutation impact</b><small>Pathogenicity & function</small></span>
          <span>🔬 <b>Structure context</b><small>3D mutation mapping</small></span>
          <span>⚡ <b>Stability & ΔG</b><small>Energetic changes</small></span>
          <span>💊 <b>Therapeutic mapping</b><small>Drug & target evidence</small></span>
          <span>🧠 <b>Explainable AI</b><small>Mechanistic insights</small></span>
          <span>📄 <b>Final report</b><small>Self-contained research report</small></span>
        </div>
        <button className="primary-btn dashboard-start" onClick={onStart}>Start New Protein Analysis <span>→</span></button>
      </div>
      <div className="dashboard-visual dashboard-visual-pro">
        <div className="visual-ring ring-a"/><div className="visual-ring ring-b"/>
        <img src={heroProtein} alt="Protein molecular visualization"/>
        <div className="visual-badge"><b>10</b><span>connected analysis stages</span></div>
      </div>
    </section>

    <section className="dashboard-section dashboard-workflow">
      <div className="section-heading"><div><span className="eyebrow">ANALYSIS PIPELINE</span><h2>Your research workflow</h2></div><button className="text-btn" onClick={onWorkflow}>View complete workflow →</button></div>
      <div className="dashboard-steps">{workflow.map(w => <div className="dash-step" key={w.n}><span>{w.n}</span><div><b>{w.title}</b><small>{w.short}</small></div><i>→</i></div>)}</div>
    </section>
  </main>;
}

class ErrorBoundary extends Component {
  constructor(props) {
    super(props);
    this.state = { hasError: false, error: null };
  }
  static getDerivedStateFromError(error) {
    return { hasError: true, error };
  }
  componentDidCatch(error, errorInfo) {
    console.error("Pipeline Step caught in ErrorBoundary:", error, errorInfo);
  }
  render() {
    if (this.state.hasError) {
      return (
        <div className="step-panel" style={{ borderLeft: "4px solid #ef4444", padding: "1.5rem", background: "rgba(239, 68, 68, 0.08)" }}>
          <h3 style={{ color: "#ef4444", marginTop: 0 }}>⚠️ Analysis Step Render Notice</h3>
          <p style={{ color: "#cbd5e1" }}>
            An unexpected error occurred while displaying this step: <code>{this.state.error?.message || "Unknown error"}</code>
          </p>
          <button className="primary-btn" onClick={() => this.setState({ hasError: false, error: null })}>
            🔄 Reload This Step
          </button>
        </div>
      );
    }
    return this.props.children;
  }
}

function AnalysisWorkspace({ payload, setPayload, onBack }) {
  const [step, setStep] = useState(1);
  const page = stepPages[step - 1];
  const Comp = page.Comp;
  const pad = n => String(n).padStart(2, "0");
  const go = n => { setStep(n); window.scrollTo({ top: 0, behavior: "smooth" }); };
  return <main className="analysis-page">
    <div className="analysis-top"><button className="back-link" onClick={onBack}>← Dashboard</button><span className="step-counter">STEP {pad(step)} OF 10</span></div>
    <nav className="stepper" aria-label="Pipeline steps">
      {stepPages.map((p, i) => {
        const n = i + 1, locked = n > 1 && !payload;
        return <button key={n} className={`stepper-dot ${n === step ? "current" : ""} ${n < step ? "done" : ""}`} disabled={locked} title={p.title} onClick={() => go(n)}>{pad(n)}</button>;
      })}
    </nav>
    <section className="analysis-heading"><div className="step-badge">{pad(step)}</div><div><span className="eyebrow">PROTMIND AI · ANALYSIS PIPELINE</span><h1>{page.icon} {page.title}</h1><p>{page.sub}</p></div>{payload && <div className="ready-pill"><span/> {payload.uniprot_id} · {payload.mutation}</div>}</section>
    <div className="step-body">
      <ErrorBoundary key={step}>
        {step === 1
          ? <Step1 payload={payload} onSaved={setPayload} />
          : <Gate payload={payload}><Comp key={step} payload={payload} /></Gate>}
      </ErrorBoundary>
    </div>
    <div className="analysis-actions">
      <button className="secondary-btn" onClick={() => step === 1 ? onBack() : go(step - 1)}>{step === 1 ? "Back" : `← Step ${pad(step - 1)}`}</button>
      {step < 10 && <button className="primary-btn" disabled={!payload} onClick={() => go(step + 1)}>Continue to {stepPages[step].title} <span>→</span></button>}
    </div>
  </main>;
}

export default App;
