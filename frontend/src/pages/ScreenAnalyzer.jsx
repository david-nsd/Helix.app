import { useState, useRef, useCallback } from 'react'
import api from '../api'

const T = {
  bg:      '#020a06',
  surface: '#0a1a0f',
  border:  'rgba(0,255,136,0.12)',
  green:   '#00ff88',
  amber:   '#ffaa00',
  red:     '#ff2244',
  dim:     '#1a4a2a',
  text:    '#c8f5d8',
  muted:   '#4a8a5a',
  font:    'monospace',
}

const card = {
  background: T.surface,
  border: `1px solid ${T.border}`,
  borderRadius: 8,
  padding: '20px 24px',
}

// ─── example data ─────────────────────────────────────────────────────────────

function negBinom(mean, dispersion) {
  const r     = mean / dispersion
  const p     = dispersion / (1 + dispersion)
  const gamma = -Math.log(Math.random() || 1e-10) / r
  return Math.max(0, Math.floor(-Math.log(Math.random() || 1e-10) / Math.max(gamma, 1e-10) * p / (1 - p)))
}

function gaussNoise(std) {
  const u1 = Math.random() || 1e-10
  const u2 = Math.random()
  return std * Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2)
}

function buildExampleCSV() {
  let csv = 'gene,guide,control_count,treatment_count\n'

  // 50 essential genes — strongly depleted
  for (let g = 1; g <= 50; g++) {
    for (let i = 1; i <= 6; i++) {
      const ctrl       = Math.floor(negBinom(500, 0.3))
      const depletion  = 0.03 + Math.random() * 0.08
      const trt        = Math.max(1, Math.floor(ctrl * depletion + gaussNoise(5)))
      csv += `ESS${String(g).padStart(3,'0')},ESS${String(g).padStart(3,'0')}_g${i},${ctrl},${trt}\n`
    }
  }

  // 20 enriched genes — tumor suppressors
  for (let g = 1; g <= 20; g++) {
    for (let i = 1; i <= 6; i++) {
      const ctrl        = Math.floor(negBinom(300, 0.3))
      const enrichment  = 3 + Math.random() * 4
      const trt         = Math.max(1, Math.floor(ctrl * enrichment + gaussNoise(20)))
      csv += `ENR${String(g).padStart(3,'0')},ENR${String(g).padStart(3,'0')}_g${i},${ctrl},${trt}\n`
    }
  }

  // 430 neutral genes
  for (let g = 1; g <= 430; g++) {
    for (let i = 1; i <= 6; i++) {
      const ctrl  = Math.floor(negBinom(400, 0.3))
      const fold  = 0.85 + Math.random() * 0.3
      const trt   = Math.max(1, Math.floor(ctrl * fold + gaussNoise(15)))
      csv += `NEU${String(g).padStart(3,'0')},NEU${String(g).padStart(3,'0')}_g${i},${ctrl},${trt}\n`
    }
  }

  return csv
}

// ─── volcano plot ─────────────────────────────────────────────────────────────

const SVG_W = 800, SVG_H = 320
const PAD_L = 50, PAD_R = 20, PAD_T = 20, PAD_B = 40
const PLOT_W = SVG_W - PAD_L - PAD_R
const PLOT_H = SVG_H - PAD_T - PAD_B

function geneX(lfc)  { return PAD_L + ((lfc + 6) / 12) * PLOT_W }
function geneY(fdr)  {
  const nl = Math.min(-Math.log10(Math.max(fdr, 0.001)), 6)
  return PAD_T + PLOT_H - (nl / 6) * PLOT_H
}

const X_ZERO   = geneX(0)
const X_MINUS1 = geneX(-1)
const X_PLUS1  = geneX(1)
const Y_FDR    = geneY(0.1)

function VolcanoPlot({ genes }) {
  const svgRef = useRef(null)
  const [tooltip, setTooltip] = useState(null)

  if (!genes || genes.length === 0) return null

  const dotColor = g => {
    if (g.classification === 'depleted') return T.red
    if (g.classification === 'enriched') return T.green
    return '#2a5a3a'
  }

  const labeledGenes = [...genes]
    .sort((a, b) => a.fdr - b.fdr)
    .slice(0, 8)

  return (
    <div style={{ position: 'relative', width: '100%', height: '320px', background: '#020a06' }}>
      <svg
        ref={svgRef}
        viewBox="0 0 800 320"
        width="100%"
        height="320"
        style={{ display: 'block' }}
        onMouseLeave={() => setTooltip(null)}
      >
        {/* axes */}
        <line x1={PAD_L} y1={PAD_T} x2={PAD_L} y2={PAD_T + PLOT_H}
          stroke={T.muted} strokeWidth={1} />
        <line x1={PAD_L} y1={PAD_T + PLOT_H} x2={PAD_L + PLOT_W} y2={PAD_T + PLOT_H}
          stroke={T.muted} strokeWidth={1} />

        {/* threshold lines */}
        <line x1={PAD_L} y1={Y_FDR} x2={PAD_L + PLOT_W} y2={Y_FDR}
          stroke="rgba(255,255,255,0.1)" strokeWidth={1} strokeDasharray="4" />
        <line x1={X_ZERO} y1={PAD_T} x2={X_ZERO} y2={PAD_T + PLOT_H}
          stroke="rgba(255,255,255,0.1)" strokeWidth={1} strokeDasharray="4" />
        <line x1={X_MINUS1} y1={PAD_T} x2={X_MINUS1} y2={PAD_T + PLOT_H}
          stroke="rgba(255,255,255,0.07)" strokeWidth={1} strokeDasharray="4" />
        <line x1={X_PLUS1} y1={PAD_T} x2={X_PLUS1} y2={PAD_T + PLOT_H}
          stroke="rgba(255,255,255,0.07)" strokeWidth={1} strokeDasharray="4" />

        {/* dots */}
        {genes.map((g, i) => {
          const cx = geneX(g.mean_log2fc)
          const cy = geneY(g.fdr)
          return (
            <circle
              key={i} cx={cx} cy={cy}
              r={g.classification !== 'neutral' ? 5 : 3}
              fill={dotColor(g)}
              opacity={g.classification === 'neutral' ? 0.4 : 0.85}
              style={{ cursor: 'pointer' }}
              onMouseEnter={() => setTooltip({ gene: g, cx, cy })}
            />
          )
        })}

        {/* gene labels */}
        {labeledGenes.map((g, i) => {
          const cx = geneX(g.mean_log2fc)
          const cy = geneY(g.fdr)
          const anchor = cx > SVG_W / 2 ? 'end' : 'start'
          const dx = cx > SVG_W / 2 ? -7 : 7
          return (
            <text key={i} x={cx + dx} y={cy - 5}
              fill={dotColor(g)} fontSize={9} fontFamily={T.font}
              textAnchor={anchor} opacity={0.9}
            >
              {g.gene}
            </text>
          )
        })}

        {/* x-axis labels */}
        {[-6, -4, -2, 0, 2, 4, 6].map(v => (
          <text key={v} x={geneX(v)} y={SVG_H - 5}
            fill="#1a4a2a" fontSize={9} fontFamily={T.font} textAnchor="middle">
            {v}
          </text>
        ))}

        {/* y-axis labels */}
        {[0, 1, 2, 3, 4, 5, 6].map(v => (
          <text key={v} x={PAD_L - 6} y={geneY(Math.pow(10, -v)) + 3}
            fill="#1a4a2a" fontSize={9} fontFamily={T.font} textAnchor="end">
            {v}
          </text>
        ))}

        {/* axis titles */}
        <text x={PAD_L + PLOT_W / 2} y={SVG_H - 2}
          fill="#1a4a2a" fontSize={10} fontFamily={T.font} textAnchor="middle">
          mean log₂ FC
        </text>
        <text x={12} y={PAD_T + PLOT_H / 2}
          fill="#1a4a2a" fontSize={10} fontFamily={T.font} textAnchor="middle"
          transform={`rotate(-90, 12, ${PAD_T + PLOT_H / 2})`}>
          −log₁₀(FDR)
        </text>

        {/* FDR label */}
        <text x={PAD_L + PLOT_W - 2} y={Y_FDR - 3}
          fill="rgba(255,255,255,0.2)" fontSize={8} fontFamily={T.font} textAnchor="end">
          FDR=0.1
        </text>

        {/* tooltip */}
        {tooltip && (() => {
          const { gene: g, cx, cy } = tooltip
          const bx = cx > SVG_W / 2 ? cx - 138 : cx + 8
          const by = Math.max(cy - 30, PAD_T)
          return (
            <g>
              <rect x={bx} y={by} width={130} height={60} rx={4}
                fill="#0a1a0f" stroke="rgba(0,255,136,0.3)" />
              <text x={bx + 8} y={by + 15} fill={T.green} fontSize={10} fontFamily={T.font} fontWeight={700}>
                {g.gene}
              </text>
              <text x={bx + 8} y={by + 28} fill={T.text} fontSize={9} fontFamily={T.font}>
                LFC: {g.mean_log2fc.toFixed(3)}
              </text>
              <text x={bx + 8} y={by + 40} fill={T.text} fontSize={9} fontFamily={T.font}>
                FDR: {g.fdr < 0.001 ? g.fdr.toExponential(1) : g.fdr.toFixed(3)}
              </text>
              <text x={bx + 8} y={by + 52} fill={dotColor(g)} fontSize={9} fontFamily={T.font}>
                {g.classification} · {g.n_guides} guides
              </text>
            </g>
          )
        })()}
      </svg>

      {/* legend */}
      <div style={{
        position: 'absolute', bottom: 44, right: 24,
        display: 'flex', gap: 12, fontSize: 10, fontFamily: T.font, color: T.muted,
      }}>
        {[['Depleted', T.red], ['Enriched', T.green], ['Neutral', '#2a5a3a']].map(([label, color]) => (
          <span key={label} style={{ display: 'flex', alignItems: 'center', gap: 4 }}>
            <span style={{ width: 8, height: 8, borderRadius: '50%', background: color, display: 'inline-block' }} />
            {label}
          </span>
        ))}
      </div>
    </div>
  )
}

// ─── gene table ───────────────────────────────────────────────────────────────

function GeneTable({ genes, type }) {
  const filtered = genes.filter(g => g.classification === type).slice(0, 20)
  if (filtered.length === 0) return (
    <div style={{ color: T.dim, fontSize: 12, fontFamily: T.font, padding: '12px 0' }}>
      No {type} genes found
    </div>
  )

  const accent = type === 'depleted' ? T.red : T.green

  return (
    <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: 12, fontFamily: T.font }}>
      <thead>
        <tr style={{ borderBottom: `1px solid ${T.border}` }}>
          {['Rank', 'Gene', 'log₂FC', 'FDR', 'Guides'].map(h => (
            <th key={h} style={{ padding: '6px 10px', textAlign: 'left', color: T.muted, fontWeight: 600, fontSize: 10 }}>
              {h}
            </th>
          ))}
        </tr>
      </thead>
      <tbody>
        {filtered.map((g, i) => (
          <tr key={g.gene} style={{ borderBottom: `1px solid rgba(0,255,136,0.05)` }}>
            <td style={{ padding: '5px 10px', color: T.dim }}>{g.rank}</td>
            <td style={{ padding: '5px 10px', color: accent, fontWeight: 600 }}>{g.gene}</td>
            <td style={{ padding: '5px 10px', color: g.mean_log2fc < 0 ? T.red : T.green }}>
              {g.mean_log2fc > 0 ? '+' : ''}{g.mean_log2fc.toFixed(3)}
            </td>
            <td style={{ padding: '5px 10px', color: g.fdr < 0.05 ? T.amber : T.text }}>
              {g.fdr < 0.001 ? g.fdr.toExponential(1) : g.fdr.toFixed(3)}
            </td>
            <td style={{ padding: '5px 10px', color: T.muted }}>{g.n_guides}</td>
          </tr>
        ))}
      </tbody>
    </table>
  )
}

// ─── guide table ──────────────────────────────────────────────────────────────

function GuideTable({ guides }) {
  const [filter, setFilter]   = useState('')
  const [sortKey, setSortKey] = useState('log2fc')
  const [sortAsc, setSortAsc] = useState(true)

  const cols = [
    { key: 'guide',           label: 'Guide'     },
    { key: 'gene',            label: 'Gene'      },
    { key: 'control_count',   label: 'Ctrl'      },
    { key: 'treatment_count', label: 'Trt'       },
    { key: 'log2fc',          label: 'log₂FC'    },
    { key: 'z_score',         label: 'Z-score'   },
  ]

  const toggle = key => {
    if (sortKey === key) setSortAsc(a => !a)
    else { setSortKey(key); setSortAsc(true) }
  }

  const visible = guides
    .filter(g => !filter || g.gene.toLowerCase().includes(filter.toLowerCase())
                         || g.guide.toLowerCase().includes(filter.toLowerCase()))
    .sort((a, b) => {
      const va = a[sortKey], vb = b[sortKey]
      if (typeof va === 'number') return sortAsc ? va - vb : vb - va
      return sortAsc ? String(va).localeCompare(String(vb)) : String(vb).localeCompare(String(va))
    })
    .slice(0, 500)

  return (
    <div>
      <input
        value={filter}
        onChange={e => setFilter(e.target.value)}
        placeholder="Filter by gene or guide…"
        style={{
          background: '#0d1f12', border: `1px solid ${T.border}`, borderRadius: 6,
          color: T.text, fontFamily: T.font, fontSize: 12,
          padding: '6px 12px', width: 260, marginBottom: 12, outline: 'none',
        }}
      />
      <div style={{ maxHeight: 360, overflowY: 'auto', border: `1px solid ${T.border}`, borderRadius: 6 }}>
        <table style={{ width: '100%', borderCollapse: 'collapse', fontSize: 11, fontFamily: T.font }}>
          <thead style={{ position: 'sticky', top: 0, background: '#0d1f12', zIndex: 1 }}>
            <tr>
              {cols.map(({ key, label }) => (
                <th key={key} onClick={() => toggle(key)}
                  style={{ padding: '7px 10px', textAlign: 'left', color: sortKey === key ? T.green : T.muted,
                    fontWeight: 600, fontSize: 10, cursor: 'pointer', userSelect: 'none',
                    borderBottom: `1px solid ${T.border}` }}>
                  {label} {sortKey === key ? (sortAsc ? '↑' : '↓') : ''}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {visible.map((g, i) => (
              <tr key={i} style={{ borderBottom: `1px solid rgba(0,255,136,0.04)` }}>
                <td style={{ padding: '4px 10px', color: T.muted, fontSize: 10, fontFamily: T.font }}>
                  {g.guide.slice(0, 12)}…
                </td>
                <td style={{ padding: '4px 10px', color: T.text }}>{g.gene}</td>
                <td style={{ padding: '4px 10px', color: T.muted }}>{g.control_count.toLocaleString()}</td>
                <td style={{ padding: '4px 10px', color: T.muted }}>{g.treatment_count.toLocaleString()}</td>
                <td style={{ padding: '4px 10px', color: g.log2fc < -1 ? T.red : g.log2fc > 1 ? T.green : T.text }}>
                  {g.log2fc > 0 ? '+' : ''}{g.log2fc.toFixed(3)}
                </td>
                <td style={{ padding: '4px 10px', color: T.muted }}>{g.z_score.toFixed(2)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      {visible.length === 500 && (
        <div style={{ fontSize: 10, color: T.dim, fontFamily: T.font, marginTop: 6 }}>
          Showing first 500 guides. Export CSV for full data.
        </div>
      )}
    </div>
  )
}

// ─── stat card ────────────────────────────────────────────────────────────────

function StatCard({ label, value, accent }) {
  return (
    <div style={{ ...card, flex: '1 1 140px', textAlign: 'center', padding: '14px 16px' }}>
      <div style={{ fontSize: 22, fontFamily: T.font, fontWeight: 700, color: accent || T.green }}>
        {value}
      </div>
      <div style={{ fontSize: 10, color: T.muted, fontFamily: T.font, marginTop: 4, letterSpacing: '0.5px' }}>
        {label}
      </div>
    </div>
  )
}

// ─── qc bar ───────────────────────────────────────────────────────────────────

function QCBar({ label, value, max, accent }) {
  const pct = Math.min((value / max) * 100, 100)
  return (
    <div style={{ marginBottom: 10 }}>
      <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: 4 }}>
        <span style={{ fontSize: 11, fontFamily: T.font, color: T.muted }}>{label}</span>
        <span style={{ fontSize: 11, fontFamily: T.font, color: accent || T.text }}>{value}</span>
      </div>
      <div style={{ height: 4, background: '#0d1f12', borderRadius: 2, overflow: 'hidden' }}>
        <div style={{ width: `${pct}%`, height: '100%', background: accent || T.green,
          borderRadius: 2, transition: 'width 0.4s ease' }} />
      </div>
    </div>
  )
}

// ─── main component ───────────────────────────────────────────────────────────

export default function ScreenAnalyzer() {
  const [file, setFile]           = useState(null)
  const [headers, setHeaders]     = useState([])
  const [colMap, setColMap]       = useState({ guide: '', gene: '', control: '', treatment: '' })
  const [minReads, setMinReads]   = useState(30)
  const [loading, setLoading]     = useState(false)
  const [error, setError]         = useState(null)
  const [results, setResults]     = useState(null)
  const [showGuides, setShowGuides] = useState(false)
  const [isDragging, setIsDragging] = useState(false)
  const fileInputRef              = useRef(null)
  const svgContainerRef           = useRef(null)

  const KEYWORD_MAP = {
    guide:     ['guide', 'sgrna', 'grna', 'spacer', 'sequence'],
    gene:      ['gene', 'target', 'symbol'],
    control:   ['control', 'ctrl', 'plasmid', 'input', 't0'],
    treatment: ['treatment', 'treat', 'screen', 'output', 't1', 'selected'],
  }

  const parseHeaders = useCallback((text) => {
    const firstLine = text.split('\n')[0]
    const sep = firstLine.includes('\t') ? '\t' : ','
    const cols = firstLine.split(sep).map(h => h.trim().replace(/^"|"$/g, ''))
    setHeaders(cols)

    const auto = { guide: '', gene: '', control: '', treatment: '' }
    for (const [role, keywords] of Object.entries(KEYWORD_MAP)) {
      for (const col of cols) {
        if (keywords.some(kw => col.toLowerCase().includes(kw))) {
          auto[role] = col
          break
        }
      }
    }
    setColMap(auto)
  }, [])

  const handleFile = useCallback((f) => {
    if (!f) return
    setFile(f)
    setResults(null)
    setError(null)
    const reader = new FileReader()
    reader.onload = e => parseHeaders(e.target.result)
    reader.readAsText(f)
  }, [parseHeaders])

  const handleDrop = useCallback((e) => {
    e.preventDefault()
    setIsDragging(false)
    const f = e.dataTransfer.files[0]
    if (f && f.name.endsWith('.csv')) handleFile(f)
  }, [handleFile])

  const handleExampleData = async () => {
    const csv = buildExampleCSV()
    const blob = new Blob([csv], { type: 'text/csv' })
    const file = new File([blob], 'example_screen.csv')
    setFile(file)
    setHeaders(['gene', 'guide', 'control_count', 'treatment_count'])
    setColMap({ guide: 'guide', gene: 'gene', control: 'control_count', treatment: 'treatment_count' })
    const formData = new FormData()
    formData.append('file', file)
    formData.append('guide_col', 'guide')
    formData.append('gene_col', 'gene')
    formData.append('control_col', 'control_count')
    formData.append('treatment_col', 'treatment_count')
    formData.append('min_reads', '30')
    setLoading(true)
    setError(null)
    try {
      const resp = await api.post('/api/screen/analyze', formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
      })
      setResults(resp.data)
    } catch (e) {
      setError(e.response?.data?.detail || 'Analysis failed')
    } finally {
      setLoading(false)
    }
  }

  const handleAnalyze = async () => {
    if (!file || !colMap.guide || !colMap.gene || !colMap.control || !colMap.treatment) return
    setLoading(true)
    setError(null)
    try {
      const form = new FormData()
      form.append('file', file)
      form.append('guide_col', colMap.guide)
      form.append('gene_col', colMap.gene)
      form.append('control_col', colMap.control)
      form.append('treatment_col', colMap.treatment)
      form.append('min_reads', String(minReads))
      const { data } = await api.post('/api/screen/analyze', form, {
        headers: { 'Content-Type': 'multipart/form-data' },
      })
      setResults(data)
    } catch (e) {
      setError(e.response?.data?.detail || e.message)
    } finally {
      setLoading(false)
    }
  }

  const exportGeneCSV = () => {
    if (!results) return
    const header = 'rank,gene,n_guides,mean_log2fc,median_log2fc,p_value,fdr,classification\n'
    const rows = results.genes.map(g =>
      `${g.rank},${g.gene},${g.n_guides},${g.mean_log2fc},${g.median_log2fc},${g.p_value},${g.fdr},${g.classification}`
    ).join('\n')
    const blob = new Blob([header + rows], { type: 'text/csv' })
    const a = document.createElement('a')
    a.href = URL.createObjectURL(blob)
    a.download = 'screen_genes.csv'
    a.click()
  }

  const exportVolcanoSVG = () => {
    const svg = svgContainerRef.current?.querySelector('svg')
    if (!svg) return
    const data = new XMLSerializer().serializeToString(svg)
    const blob = new Blob([data], { type: 'image/svg+xml' })
    const a = document.createElement('a')
    a.href = URL.createObjectURL(blob)
    a.download = 'volcano_plot.svg'
    a.click()
  }

  const allMapped = colMap.guide && colMap.gene && colMap.control && colMap.treatment

  const inputStyle = {
    background: '#0d1f12', border: `1px solid ${T.border}`, borderRadius: 6,
    color: T.text, fontFamily: T.font, fontSize: 12,
    padding: '6px 10px', outline: 'none', width: '100%',
  }

  const selectStyle = { ...inputStyle, cursor: 'pointer' }

  const btnPrimary = {
    background: T.green, color: '#020a06', border: 'none', borderRadius: 6,
    fontFamily: T.font, fontSize: 12, fontWeight: 700, cursor: 'pointer',
    padding: '8px 20px', letterSpacing: '0.5px',
  }

  const btnSecondary = {
    background: 'transparent', color: T.muted, border: `1px solid ${T.dim}`,
    borderRadius: 6, fontFamily: T.font, fontSize: 11, cursor: 'pointer',
    padding: '6px 14px',
  }

  const sectionLabel = { fontSize: 10, color: T.muted, fontFamily: T.font,
    letterSpacing: '1px', textTransform: 'uppercase', marginBottom: 6 }

  return (
    <div style={{ maxWidth: 1100, margin: '0 auto' }}>
      {/* header */}
      <div style={{ marginBottom: 28 }}>
        <h1 style={{ color: T.green, fontFamily: T.font, fontSize: 22, fontWeight: 700,
          letterSpacing: '2px', margin: 0 }}>
          CRISPR SCREEN ANALYZER
        </h1>
        <p style={{ color: T.muted, fontFamily: T.font, fontSize: 12, marginTop: 6 }}>
          MAGeCK-inspired pooled screen analysis · log₂FC · RRA · BH FDR correction
        </p>
      </div>

      {/* upload + config */}
      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 16, marginBottom: 16 }}>

        {/* dropzone */}
        <div style={{ ...card }}>
          <div style={sectionLabel}>Input Data</div>
          <div
            onDragOver={e => { e.preventDefault(); setIsDragging(true) }}
            onDragLeave={() => setIsDragging(false)}
            onDrop={handleDrop}
            onClick={() => fileInputRef.current?.click()}
            style={{
              border: `2px dashed ${isDragging ? T.green : T.border}`,
              borderRadius: 8, padding: '28px 20px', textAlign: 'center',
              cursor: 'pointer', transition: 'border-color 0.2s',
              background: isDragging ? 'rgba(0,255,136,0.04)' : 'transparent',
              marginBottom: 12,
            }}
          >
            <input ref={fileInputRef} type="file" accept=".csv"
              style={{ display: 'none' }}
              onChange={e => handleFile(e.target.files[0])} />
            <div style={{ fontSize: 28, marginBottom: 8 }}>📂</div>
            <div style={{ color: file ? T.green : T.muted, fontFamily: T.font, fontSize: 12 }}>
              {file ? file.name : 'Drop CSV or click to browse'}
            </div>
            <div style={{ color: T.dim, fontFamily: T.font, fontSize: 10, marginTop: 4 }}>
              guide_id · gene · control_reads · treatment_reads
            </div>
          </div>
          <button onClick={handleExampleData} disabled={loading} style={{ ...btnSecondary, width: '100%' }}>
            {loading ? 'Analyzing…' : 'Load Example Dataset (500 genes · 3000 guides)'}
          </button>
        </div>

        {/* column mapper */}
        <div style={{ ...card }}>
          <div style={sectionLabel}>Column Mapping</div>
          {(['guide', 'gene', 'control', 'treatment']).map(role => (
            <div key={role} style={{ marginBottom: 10 }}>
              <label style={{ fontSize: 10, color: T.muted, fontFamily: T.font,
                display: 'block', marginBottom: 3, textTransform: 'capitalize' }}>
                {role === 'control' ? 'Control Reads' : role === 'treatment' ? 'Treatment Reads'
                  : role === 'guide' ? 'Guide Column' : 'Gene Column'}
              </label>
              <select
                value={colMap[role]}
                onChange={e => setColMap(m => ({ ...m, [role]: e.target.value }))}
                style={selectStyle}
                disabled={headers.length === 0}
              >
                <option value="">— select column —</option>
                {headers.map(h => <option key={h} value={h}>{h}</option>)}
              </select>
            </div>
          ))}
          <div style={{ display: 'flex', alignItems: 'center', gap: 10, marginTop: 4 }}>
            <label style={{ fontSize: 10, color: T.muted, fontFamily: T.font, whiteSpace: 'nowrap' }}>
              Min reads (ctrl)
            </label>
            <input
              type="number" value={minReads} min={1} max={1000}
              onChange={e => setMinReads(Number(e.target.value))}
              style={{ ...inputStyle, width: 70 }}
            />
          </div>
        </div>
      </div>

      {/* run button */}
      <div style={{ display: 'flex', gap: 12, marginBottom: 24 }}>
        <button
          onClick={handleAnalyze}
          disabled={!file || !allMapped || loading}
          style={{
            ...btnPrimary,
            opacity: (!file || !allMapped || loading) ? 0.4 : 1,
            cursor: (!file || !allMapped || loading) ? 'not-allowed' : 'pointer',
          }}
        >
          {loading ? 'Analyzing…' : 'Run Analysis'}
        </button>
        {results && (
          <>
            <button onClick={exportGeneCSV} style={btnSecondary}>Export Gene CSV</button>
            <button onClick={exportVolcanoSVG} style={btnSecondary}>Export Volcano SVG</button>
          </>
        )}
      </div>

      {error && (
        <div style={{ ...card, border: `1px solid ${T.red}`, color: T.red,
          fontFamily: T.font, fontSize: 12, marginBottom: 16 }}>
          {error}
        </div>
      )}

      {/* results */}
      {results && (() => {
        const { summary, genes, guides, qc_metrics } = results
        return (
          <>
            {/* summary cards */}
            <div style={{ display: 'flex', gap: 12, flexWrap: 'wrap', marginBottom: 16 }}>
              <StatCard label="Total Guides"     value={summary.total_guides.toLocaleString()} />
              <StatCard label="Total Genes"      value={summary.total_genes.toLocaleString()} />
              <StatCard label="Filtered Out"     value={summary.guides_filtered.toLocaleString()} accent={T.amber} />
              <StatCard label="Depleted Genes"   value={summary.depleted_genes} accent={T.red} />
              <StatCard label="Enriched Genes"   value={summary.enriched_genes} accent={T.green} />
              <StatCard label="Neutral Genes"    value={summary.neutral_genes} accent={T.muted} />
              <StatCard label="Mean log₂FC"      value={summary.mean_log2fc.toFixed(3)} accent={T.text} />
              <StatCard label="Std log₂FC"       value={summary.std_log2fc.toFixed(3)} accent={T.text} />
            </div>

            {/* qc metrics */}
            <div style={{ ...card, marginBottom: 16 }}>
              <div style={sectionLabel}>QC Metrics</div>
              <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '0 40px' }}>
                <QCBar label="Gini Coefficient (library evenness)" value={qc_metrics.gini_coefficient}
                  max={1} accent={qc_metrics.gini_coefficient > 0.4 ? T.red : T.green} />
                <QCBar label="Ctrl–Trt Pearson Correlation" value={qc_metrics.correlation}
                  max={1} accent={qc_metrics.correlation > 0.8 ? T.green : T.amber} />
                <div style={{ marginBottom: 8 }}>
                  <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: 4 }}>
                    <span style={{ fontSize: 11, fontFamily: T.font, color: T.muted }}>Zero-count guides</span>
                    <span style={{ fontSize: 11, fontFamily: T.font, color: T.amber }}>{qc_metrics.zero_count_guides}</span>
                  </div>
                </div>
                <div>
                  <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: 4 }}>
                    <span style={{ fontSize: 11, fontFamily: T.font, color: T.muted }}>High-count guides (top 5%)</span>
                    <span style={{ fontSize: 11, fontFamily: T.font, color: T.text }}>{qc_metrics.high_count_guides}</span>
                  </div>
                </div>
              </div>
            </div>

            {/* volcano plot */}
            <div style={{ ...card, marginBottom: 16 }} ref={svgContainerRef}>
              <div style={sectionLabel}>Volcano Plot</div>
              <VolcanoPlot genes={genes} />
            </div>

            {/* gene tables */}
            <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 16, marginBottom: 16 }}>
              <div style={card}>
                <div style={{ ...sectionLabel, color: T.red }}>Top Depleted Genes</div>
                <GeneTable genes={genes} type="depleted" />
              </div>
              <div style={card}>
                <div style={{ ...sectionLabel, color: T.green }}>Top Enriched Genes</div>
                <GeneTable genes={genes} type="enriched" />
              </div>
            </div>

            {/* guide table toggle */}
            <div style={card}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: showGuides ? 16 : 0 }}>
                <div style={sectionLabel}>Guide-level Data</div>
                <button onClick={() => setShowGuides(s => !s)} style={btnSecondary}>
                  {showGuides ? 'Hide' : `Show ${guides.length.toLocaleString()} guides`}
                </button>
              </div>
              {showGuides && <GuideTable guides={guides} />}
            </div>
          </>
        )
      })()}
    </div>
  )
}
