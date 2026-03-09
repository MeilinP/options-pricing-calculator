# Options Pricing Engine & Live Implied Volatility Surface

[![React](https://img.shields.io/badge/React-18-blue)](https://reactjs.org)
[![Python](https://img.shields.io/badge/Python-3.9+-blue)](https://python.org)
[![Streamlit](https://img.shields.io/badge/Streamlit-Live-red)](https://volatility-surface-mp.streamlit.app)
[![Vercel](https://img.shields.io/badge/Vercel-Live-black)](https://options-pricing-calculator.vercel.app)

A full-stack options analytics suite combining a multi-model pricing engine with a live implied volatility surface visualization tool.

| Component | Stack | Link |
|-----------|-------|------|
| Options Pricing Calculator | React, Recharts | [options-pricing-calculator.vercel.app](https://options-pricing-calculator.vercel.app) |
| Live IV Surface | Python, Streamlit, yfinance | [volatility-surface-mp.streamlit.app](https://volatility-surface-mp.streamlit.app) |

---

## Features

**Pricing Models**
- Black-Scholes closed-form solution for European options
- Monte Carlo simulation (Box-Muller GBM sampling, 10,000 paths)
- Binomial tree with American early-exercise support (CRR parameterization)

**Greeks & Risk**
- Full Greeks: Δ, Γ, ν, Θ, ρ
- Newton-Raphson implied volatility solver (converges to 10⁻⁶ tolerance)
- Interactive spot price and volatility sensitivity charts

**Live IV Surface**
- Real-time options chain via yfinance; Polygon.io integration available for institutional data feeds
- Moneyness filter: 92%–108% of spot; IV bounds: 5%–80%
- 3D volatility surface, front-month skew, ATM term structure

---

## File Structure
```
├── src/                          # Pricing engine (Black-Scholes, Monte Carlo, Binomial Tree)
├── App.js                        # React frontend
├── streamlit_app.py              # Live IV surface (Streamlit + yfinance)
├── polygon_volatility_surface.py # Polygon.io integration variant
├── alpaca_volatility_surface.py  # Alpaca API variant
├── live_volatility_surface.py    # Real-time streaming surface
└── requirements.txt
```

## Usage

**React Pricing Calculator**
```bash
npm install
npm start
```

**Streamlit IV Surface**
```bash
pip install -r requirements.txt
streamlit run streamlit_app.py
```
